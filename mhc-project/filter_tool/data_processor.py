import glob
import gzip
import os
import re
import shutil
import sys

import numpy as np
import pandas as pd
import requests
from Bio import SeqIO


def get_base_dir() -> str:
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))


BASE_DIR = get_base_dir()
INPUT_DIR = BASE_DIR
FILE_00_NAME = 'mhc_ligand_full_00.csv'
FASTA_GZ = os.path.join(BASE_DIR, 'human_proteome.fasta.gz')
FASTA_UNZIPPED = os.path.join(BASE_DIR, 'human_proteome.fasta')
UNIPROT_URL = 'https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&includeIsoform=true&query=%28proteome%3AUP000005640%29'


def iedb_strict_gold_mode_v3_with_gene():
    print('启动筛选：当前目录 =', INPUT_DIR)

    if not os.path.exists(FASTA_UNZIPPED):
        print('[0/6] 下载人类蛋白库...')
        try:
            r = requests.get(UNIPROT_URL, stream=True, timeout=120)
            r.raise_for_status()
            with open(FASTA_GZ, 'wb') as f:
                for chunk in r.iter_content(8192):
                    f.write(chunk)
            with gzip.open(FASTA_GZ, 'rb') as f_in, open(FASTA_UNZIPPED, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        except Exception as e:
            print(f'下载失败: {e}')
            input('按回车退出...')
            return

    print('[构建 ID 与 Gene 索引]...')
    id_to_gene_map = {}
    re_gene = re.compile(r'GN=([^\s]+)')

    for record in SeqIO.parse(FASTA_UNZIPPED, 'fasta'):
        parts = record.id.split('|')
        acc = parts[1] if len(parts) >= 2 else record.id
        match = re_gene.search(record.description)
        gene_name = match.group(1) if match else 'Unknown'
        id_to_gene_map[acc.upper()] = gene_name
        base_acc = acc.split('-')[0].upper()
        if base_acc not in id_to_gene_map:
            id_to_gene_map[base_acc] = gene_name

    valid_human_ids = set(id_to_gene_map.keys())

    try:
        print('[1/6] 加载数据...')
        all_csv_files = sorted(glob.glob(os.path.join(INPUT_DIR, 'mhc_ligand_full_*.csv')))
        path_00 = os.path.join(INPUT_DIR, FILE_00_NAME)

        if not all_csv_files or not os.path.exists(path_00):
            print('错误：当前目录未找到 mhc_ligand_full_*.csv，且必须包含 mhc_ligand_full_00.csv')
            input('按回车退出...')
            return

        template_df = pd.read_csv(path_00, header=[0, 1], nrows=0)
        original_multi_columns = template_df.columns
        target_col_count = len(original_multi_columns)

        combined_dfs = []
        for fpath in all_csv_files:
            fname = os.path.basename(fpath)
            try:
                if fname == FILE_00_NAME:
                    df_tmp = pd.read_csv(fpath, header=[0, 1], low_memory=False)
                else:
                    df_tmp = pd.read_csv(fpath, header=None, low_memory=False)
                    if len(df_tmp.columns) > target_col_count:
                        df_tmp = df_tmp.iloc[:, :target_col_count]
                    elif len(df_tmp.columns) < target_col_count:
                        for _ in range(target_col_count - len(df_tmp.columns)):
                            df_tmp[len(df_tmp.columns)] = np.nan
                    df_tmp.columns = original_multi_columns
                combined_dfs.append(df_tmp)
            except Exception as e:
                print(f'跳过 {fname}: {e}')

        df = pd.concat(combined_dfs, ignore_index=True)
        df.columns = [f'{str(col[0]).strip()}_{str(col[1]).strip()}' for col in df.columns]

        print('[2/6] 锁定关键列...')
        seq_col = [c for c in df.columns if 'Epitope_Name' in c][0]
        mhc_col = [c for c in df.columns if 'MHC Restriction_Name' in c][0]
        qual_col = [c for c in df.columns if 'Assay_Qualitative Measurement' in c][0]
        quant_col = [c for c in df.columns if 'Assay_Quantitative measurement' in c][0]
        method_col = [c for c in df.columns if 'Assay_Method' in c][0]
        iri_col = 'Epitope_Source Molecule IRI'
        host_col = [c for c in df.columns if 'Host_Name' in c][0]
        start_col = next((c for c in df.columns if 'Epitope_Start' in c), None)
        end_col = next((c for c in df.columns if 'Epitope_End' in c), None)

        print('[3/6] 基础筛选...')
        df = df[df[qual_col].str.contains('Positive', case=False, na=False)].copy()
        df = df[df[host_col].str.contains('Human|Homo sapiens', case=False, na=False)].copy()

        print('[4/6] 提取 Protein ID 并映射 Gene Name...')

        def extract_protein_id(row):
            iri = str(row[iri_col]).strip()
            if pd.isna(row[iri_col]) or iri == '' or iri.lower() == 'nan' or '/' not in iri:
                return None
            raw_id = iri.split('/')[-1].upper()
            base_id = raw_id.split('-')[0].split('.')[0]
            if raw_id in valid_human_ids:
                return raw_id
            if base_id in valid_human_ids:
                return base_id
            return None

        df['Mapped_Protein_ID'] = df.apply(extract_protein_id, axis=1)
        df_gold = df.dropna(subset=['Mapped_Protein_ID']).copy()
        df_gold['Gene_Name'] = df_gold['Mapped_Protein_ID'].map(id_to_gene_map)
        print('Gold 数据量:', len(df_gold))

        print('[5/6] 智能去重...')
        df_gold['__Method_Score'] = df_gold[method_col].apply(
            lambda x: 1 if any(s in str(x).lower() for s in ['purified', 'fp', 'elisa', 'mass spec']) else 2
        )
        df_gold[quant_col] = pd.to_numeric(df_gold[quant_col], errors='coerce')
        merge_subset = ['Mapped_Protein_ID', seq_col, mhc_col]
        df_gold['Record_Count'] = df_gold.groupby(merge_subset)[seq_col].transform('count')
        df_gold.sort_values(
            by=['Mapped_Protein_ID', seq_col, mhc_col, '__Method_Score', quant_col],
            ascending=True,
            inplace=True,
        )
        df_best = df_gold.drop_duplicates(subset=merge_subset, keep='first')

        print('[6/6] 保存结果...')
        target_cols = [seq_col, start_col, end_col, 'Gene_Name', 'Mapped_Protein_ID', mhc_col, 'Record_Count']
        final_cols = [c for c in target_cols if c and c in df_best.columns]

        mhc_names = df_best[mhc_col].astype(str).str.upper()
        df_i = df_best[mhc_names.str.contains(r'HLA-[ABCEFG]', regex=True)][final_cols]
        df_ii = df_best[mhc_names.str.contains(r'HLA-D', regex=True)][final_cols]

        path_i = os.path.join(INPUT_DIR, 'Class_I_HomoSapiens_Gold_Gene.csv')
        path_ii = os.path.join(INPUT_DIR, 'Class_II_HomoSapiens_Gold_Gene.csv')
        df_i.to_csv(path_i, index=False, encoding='utf-8-sig')
        df_ii.to_csv(path_ii, index=False, encoding='utf-8-sig')

        print('完成！')
        print(f'Class I: {len(df_i)}')
        print(f'Class II: {len(df_ii)}')
        print('输出文件：')
        print(' - Class_I_HomoSapiens_Gold_Gene.csv')
        print(' - Class_II_HomoSapiens_Gold_Gene.csv')
        input('按回车退出...')

    except Exception as e:
        print(f'错误: {e}')
        import traceback
        traceback.print_exc()
        input('按回车退出...')


if __name__ == '__main__':
    iedb_strict_gold_mode_v3_with_gene()
