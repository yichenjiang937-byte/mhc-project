import os
import sys
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import numpy as np


def get_base_dir() -> str:
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))


BASE_DIR = get_base_dir()
DATA_FILE_PATH_I = os.path.join(BASE_DIR, 'Class_I_HomoSapiens_Gold_Gene.csv')
DATA_FILE_PATH_II = os.path.join(BASE_DIR, 'Class_II_HomoSapiens_Gold_Gene.csv')

st.set_page_config(page_title='MHC Critical Analysis', layout='wide', page_icon='⚖️')

st.markdown("""
<style>
    .stTabs [data-baseweb="tab-list"] { gap: 24px; }
    .stTabs [data-baseweb="tab"] { height: 50px; white-space: pre-wrap; background-color: #f0f2f6; border-radius: 4px 4px 0 0; gap: 1px; }
    .stTabs [aria-selected="true"] { background-color: #ffffff; border-top: 3px solid #2563eb; }
</style>
""", unsafe_allow_html=True)


def file_status(label: str, path: str) -> str:
    return '✅' if os.path.exists(path) else '❌'


with st.sidebar:
    st.header('文件检查')
    st.write(f'Class I: {file_status("Class I", DATA_FILE_PATH_I)}')
    st.write(f'Class II: {file_status("Class II", DATA_FILE_PATH_II)}')
    st.caption('两个 CSV 需要和 exe 放在同一目录。')


@st.cache_data
def load_and_clean_data(path):
    if not os.path.exists(path):
        msg = (
            f'File not found: {os.path.basename(path)} '\
            '请把 Class_I_HomoSapiens_Gold_Gene.csv 和 '\
            'Class_II_HomoSapiens_Gold_Gene.csv 放到 exe 同目录。'
        )
        return None, msg

    try:
        df = pd.read_csv(path, low_memory=False, encoding='utf-8-sig')
    except Exception as e:
        return None, f'读取失败: {e}'

    df.columns = [c.strip() for c in df.columns]

    rename_map = {
        'Gene_Name': 'Gene',
        'Mapped_Protein_ID': 'Protein',
        'Record_Count': 'Count',
        'Epitope_Starting Position': 'Start',
        'Epitope_Ending Position': 'End',
        'MHC Restriction_Name': 'HLA',
        'Epitope_Name': 'Peptide',
    }
    df = df.rename(columns=rename_map)

    if 'Gene' not in df.columns and 'Protein' in df.columns:
        df['Gene'] = df['Protein']

    df['Gene'] = df['Gene'].fillna('Unknown')
    df['Start'] = pd.to_numeric(df['Start'], errors='coerce').fillna(0).astype(int)
    df['End'] = pd.to_numeric(df['End'], errors='coerce').fillna(0).astype(int)
    df['Count'] = pd.to_numeric(df['Count'], errors='coerce').fillna(1).astype(int)

    return df, None


def build_leaderboard(df):
    return (
        df[df['Gene'] != 'Unknown']
        .groupby(['Gene', 'Protein'])
        .agg(Unique_HLA_Count=('HLA', 'nunique'), Total_Reports=('Count', 'sum'))
        .reset_index()
        .sort_values(by=['Unique_HLA_Count', 'Total_Reports'], ascending=[False, False])
    )


def plot_speculative_landscape(df_sub, hla_label, depth_line_color, depth_fill_rgba):
    if df_sub.empty:
        return go.Figure()

    max_len = int(df_sub['End'].max()) + 20
    coverage_depth = np.zeros(max_len)
    hla_sets = [set() for _ in range(max_len)]

    for _, row in df_sub.iterrows():
        s, e = row['Start'], row['End']
        if 0 <= s < max_len and 0 <= e < max_len:
            coverage_depth[s:e + 1] += row['Count']
            for i in range(s, e + 1):
                hla_sets[i].add(row['HLA'])

    hla_counts = [len(s) for s in hla_sets]
    x_axis = list(range(max_len))

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=x_axis, y=coverage_depth,
        name='Total Repeat Count (Depth)',
        fill='tozeroy',
        line=dict(color=depth_line_color, width=0),
        fillcolor=depth_fill_rgba,
        yaxis='y',
        hovertemplate='Position: %{x}<br>Repeat Count: %{y}<extra></extra>',
    ))
    fig.add_trace(go.Scatter(
        x=x_axis, y=hla_counts,
        name=f'{hla_label} Type Coverage Count (Breadth)',
        mode='lines',
        line=dict(color='#ef4444', width=2.8),
        yaxis='y2',
        hovertemplate=f'Position: %{{x}}<br>{hla_label} Type Count: %{{y}}<extra></extra>',
    ))

    fig.update_layout(
        title='<b>📊 Immunogenic Speculative Graph (Depth vs Breadth)</b>'
              f'<br><sup>Fill = Depth | Red line = {hla_label} Breadth</sup>',
        xaxis=dict(title='Amino Acid Position'),
        yaxis=dict(title='Total Detection/Repeat Count (Depth)', showgrid=False),
        yaxis2=dict(
            title=f'{hla_label} Types Recognizing This Position (Count)',
            overlaying='y', side='right', showgrid=False,
            range=[0, max(hla_counts) * 1.2 if hla_counts else 5],
        ),
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1),
        template='plotly_white',
        height=520,
        hovermode='x unified',
    )
    return fig


mhc_class = st.selectbox('Select MHC Class', ['Class I', 'Class II'], index=0)

if mhc_class == 'Class I':
    data_file_path = DATA_FILE_PATH_I
    hla_label = 'HLA I'
    depth_line_color = '#2563eb'
    depth_fill_rgba = 'rgba(37, 99, 235, 0.28)'
    bar_scale = 'Blues'
    pep_cmap = 'Blues'
else:
    data_file_path = DATA_FILE_PATH_II
    hla_label = 'HLA II'
    depth_line_color = '#9333ea'
    depth_fill_rgba = 'rgba(147, 51, 234, 0.28)'
    bar_scale = 'RdPu'
    pep_cmap = 'PuRd'


df, err = load_and_clean_data(data_file_path)
if err:
    st.error(err)
    st.stop()

st.title(f'🧬 {mhc_class} Target Discovery: Depth vs Breadth')

tab1, tab2 = st.tabs(['🌍 Global View (Global Ranking)', '🔬 Deep Dive Analysis'])

with tab1:
    st.subheader('🏆 Global Gene Breadth Ranking')
    st.info(f'The ranking is based on how many types of {hla_label} recognize it, not just the data volume.')
    leaderboard = build_leaderboard(df)

    c1, c2 = st.columns([2, 1])
    with c1:
        fig = px.bar(
            leaderboard.head(20),
            x='Unique_HLA_Count', y='Gene', orientation='h', color='Total_Reports',
            color_continuous_scale=bar_scale,
            labels={'Unique_HLA_Count': f'{hla_label} Type Coverage Count', 'Gene': 'Gene Name'},
        )
        fig.update_layout(yaxis={'categoryorder': 'total ascending'})
        st.plotly_chart(fig, use_container_width=True)
    with c2:
        st.dataframe(leaderboard.head(20)[['Gene', 'Unique_HLA_Count', 'Total_Reports']], use_container_width=True, hide_index=True)

with tab2:
    st.subheader('🔎 Single Protein Deep Dive')
    leaderboard = build_leaderboard(df)
    all_genes = leaderboard['Gene'].tolist()
    col_ctrl1, col_ctrl2 = st.columns([1, 3])
    with col_ctrl1:
        sel_gene = st.selectbox('Step 1: Select Target Gene/Protein', all_genes)

    df_gene_raw = df[df['Gene'] == sel_gene]
    available_hlas = sorted(df_gene_raw['HLA'].unique().tolist())
    with col_ctrl2:
        sel_hlas = st.multiselect(
            f'Step 2: Filter {hla_label} Subtypes',
            available_hlas,
            default=available_hlas,
        )

    if sel_hlas:
        df_final = df_gene_raw[df_gene_raw['HLA'].isin(sel_hlas)]
        st.divider()
        st.markdown('#### 📈 Core Speculative Graph')
        st.plotly_chart(
            plot_speculative_landscape(df_final, hla_label, depth_line_color, depth_fill_rgba),
            use_container_width=True,
        )

        st.divider()
        st.markdown('#### 🧪 Peptide-level inspection')
        peptide_df = (
            df_final.groupby('Peptide')
            .agg(
                Start=('Start', 'min'),
                End=('End', 'max'),
                HLA_Count=('HLA', 'nunique'),
                Total_Reports=('Count', 'sum'),
            )
            .reset_index()
            .sort_values(by=['HLA_Count', 'Total_Reports'], ascending=[False, False])
        )
        st.dataframe(peptide_df, use_container_width=True, hide_index=True)
    else:
        st.warning('请至少选择一个 HLA subtype。')
