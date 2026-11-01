# Filter tool

Builds the class-stratified gold tables from the IEDB MHC ligand export.

## Files

| File | Purpose |
| --- | --- |
| `data_processor.py` | Pipeline: load shards → filter → resolve proteins → deduplicate → export |
| `requirements.txt` | Python dependencies |
| `build_exe.bat` | Builds a standalone `MHC_Data_Processor.exe` |

## Usage

Place the shards and `human_proteome.fasta` in one directory, then run:

```bash
pip install -r requirements.txt
python data_processor.py --input-dir .
```

If the FASTA is absent it is downloaded once from UniProt (UP000005640,
isoforms included). Useful options:

```text
--input-dir       directory holding mhc_ligand_full_*.csv (default: script folder)
--output-dir      where to write the gold tables (default: same as --input-dir)
--reference-shard shard defining the column layout (default: mhc_ligand_full_00.csv)
--fasta           use a local proteome FASTA instead of downloading
--no-download     fail if the proteome is missing
--modality        all (default) | eluted | binding
-v, --verbose     verbose logging
```

`--modality` is orthogonal to the evidence ranking used to pick a
representative record: it selects *which kind of evidence* the tables are built
from. The default labels every record in the `Assay_Modality` column without
dropping any; `--modality eluted` keeps only mass-spectrometry elution data,
which is the right input when the question is about presentation rather than
in vitro affinity.

On Windows, `build_exe.bat` builds the standalone executable; copy the input
files next to `MHC_Data_Processor.exe` before running it.

## Output

- `Class_I_HomoSapiens_Gold_Gene.csv` — restrictions matching `HLA-[ABCEFG]`
- `Class_II_HomoSapiens_Gold_Gene.csv` — restrictions matching `HLA-D`

Columns: peptide, start, end, gene name, UniProt accession, HLA restriction,
`Record_Count` (number of independent records supporting the entry), and
`Assay_Modality` (`eluted` / `binding` / `other`).
