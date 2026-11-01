"""Interactive explorer for the class I / class II gold ligand tables.

Two views are offered:

* **Global View** ranks genes by the number of distinct HLA restrictions that
  present peptides from them (breadth), with total report count shown as a
  secondary encoding.
* **Deep Dive** plots, for one protein, peptide depth and HLA breadth along the
  amino-acid sequence, so that allele-restricted peaks can be told apart from
  regions presented across many alleles.

The application expects :file:`Class_I_HomoSapiens_Gold_Gene.csv` and
:file:`Class_II_HomoSapiens_Gold_Gene.csv` next to the script or executable.
"""

from __future__ import annotations

import os
import sys
from typing import Optional, Tuple

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

DATA_FILE_I = "Class_I_HomoSapiens_Gold_Gene.csv"
DATA_FILE_II = "Class_II_HomoSapiens_Gold_Gene.csv"

#: Padding added beyond the last coordinate when building the position axis.
POSITION_MARGIN = 20

#: Rows whose coordinates exceed this are treated as malformed. No known
#: protein is anywhere near this length; the guard only prevents an accidental
#: huge coordinate from allocating a gigantic array.
MAX_POSITION = 1_000_000

DEFAULT_COUNT = 1

TAB_STYLE = """
<style>
    .stTabs [data-baseweb="tab-list"] { gap: 24px; }
    .stTabs [data-baseweb="tab"] { height: 50px; white-space: pre-wrap; background-color: #f0f2f6; border-radius: 4px 4px 0 0; gap: 1px; }
    .stTabs [aria-selected="true"] { background-color: #ffffff; border-top: 3px solid #2563eb; }
</style>
"""

COLUMN_ALIASES = {
    "Gene_Name": "Gene",
    "Mapped_Protein_ID": "Protein",
    "Record_Count": "Count",
    "Epitope_Starting Position": "Start",
    "Epitope_Ending Position": "End",
    "MHC Restriction_Name": "HLA",
    "Epitope_Name": "Peptide",
}


def get_base_dir() -> str:
    """Directory of the script, or of the executable when frozen."""
    if getattr(sys, "frozen", False):
        return getattr(sys, "_MEIPASS", os.path.dirname(sys.executable))
    return os.path.dirname(os.path.abspath(__file__))


BASE_DIR = get_base_dir()
DATA_PATH_I = os.path.join(BASE_DIR, DATA_FILE_I)
DATA_PATH_II = os.path.join(BASE_DIR, DATA_FILE_II)


def file_status(path: str) -> str:
    return "✅" if os.path.exists(path) else "❌"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def load_and_clean_data(path: str) -> Tuple[Optional[pd.DataFrame], Optional[str], int]:
    """Load one gold table and normalise its columns.

    Returns ``(frame, error, dropped)`` where ``dropped`` counts rows whose
    coordinates could not be used for the position axis.
    """
    if not os.path.exists(path):
        message = (
            f"File not found: {os.path.basename(path)}\n\n"
            "Place Class_I_HomoSapiens_Gold_Gene.csv and "
            "Class_II_HomoSapiens_Gold_Gene.csv next to the executable."
        )
        return None, message, 0

    try:
        frame = pd.read_csv(path, low_memory=False, encoding="utf-8-sig")
    except Exception as exc:  # surfaced to the user instead of crashing
        return None, f"Could not read {os.path.basename(path)}: {exc}", 0

    frame.columns = [c.strip() for c in frame.columns]
    frame = frame.rename(columns=COLUMN_ALIASES)

    missing = [c for c in ("Protein", "HLA", "Start", "End") if c not in frame.columns]
    if missing:
        return None, f"Missing required columns: {', '.join(missing)}", 0

    if "Gene" not in frame.columns:
        frame["Gene"] = frame["Protein"]
    frame["Gene"] = frame["Gene"].fillna("Unknown")
    frame["HLA"] = frame["HLA"].astype(str)

    if "Count" not in frame.columns:
        frame["Count"] = DEFAULT_COUNT
    frame["Count"] = (
        pd.to_numeric(frame["Count"], errors="coerce")
        .fillna(DEFAULT_COUNT)
        .astype(int)
    )

    frame["Start"] = pd.to_numeric(frame["Start"], errors="coerce")
    frame["End"] = pd.to_numeric(frame["End"], errors="coerce")

    usable = (
        frame["Start"].notna()
        & frame["End"].notna()
        & (frame["Start"] >= 0)
        & (frame["End"] >= frame["Start"])
        & (frame["End"] <= MAX_POSITION)
    )

    dropped = int((~usable).sum())
    frame = frame.loc[usable].copy()
    frame["Start"] = frame["Start"].astype(int)
    frame["End"] = frame["End"].astype(int)

    return frame, None, dropped


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner=False)
def build_leaderboard(frame: pd.DataFrame) -> pd.DataFrame:
    """Rank genes by HLA breadth, then by total report count."""
    return (
        frame[frame["Gene"] != "Unknown"]
        .groupby(["Gene", "Protein"])
        .agg(Unique_HLA_Count=("HLA", "nunique"), Total_Reports=("Count", "sum"))
        .reset_index()
        .sort_values(by=["Unique_HLA_Count", "Total_Reports"], ascending=[False, False])
    )


def build_profile(
    frame: pd.DataFrame,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return per-position peptide depth and HLA breadth.

    Depth is the summed report count of every peptide covering a position;
    breadth is the number of distinct HLA restrictions covering it. Both are
    accumulated with interval difference arrays, which is equivalent to
    iterating over every covered position but scales with the number of
    peptides instead of the peptide length.
    """
    if frame.empty:
        return np.zeros(0), np.zeros(0)

    starts = frame["Start"].to_numpy(dtype=np.int64, copy=False)
    ends = frame["End"].to_numpy(dtype=np.int64, copy=False)
    counts = frame["Count"].to_numpy(dtype=float, copy=False)
    alleles = frame["HLA"].to_numpy()

    max_len = int(ends.max()) + POSITION_MARGIN
    valid = starts <= ends

    depth_delta = np.zeros(max_len + 1, dtype=float)
    np.add.at(depth_delta, starts[valid], counts[valid])
    np.add.at(depth_delta, ends[valid] + 1, -counts[valid])
    depth = np.cumsum(depth_delta)[:max_len]

    breadth = np.zeros(max_len, dtype=np.int64)
    for allele in np.unique(alleles[valid]):
        selected = valid & (alleles == allele)
        delta = np.zeros(max_len + 1, dtype=np.int64)
        np.add.at(delta, starts[selected], 1)
        np.add.at(delta, ends[selected] + 1, -1)
        breadth += (np.cumsum(delta)[:max_len] > 0)

    return depth, breadth


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def plot_speculative_landscape(
    frame: pd.DataFrame,
    hla_label: str,
    depth_line_color: str,
    depth_fill_rgba: str,
) -> go.Figure:
    """Dual-axis depth/breadth profile for one protein."""
    if frame.empty:
        return go.Figure()

    depth, breadth = build_profile(frame)
    positions = np.arange(depth.size)

    figure = go.Figure()
    figure.add_trace(
        go.Scatter(
            x=positions, y=depth,
            name="Total Repeat Count (Depth)",
            fill="tozeroy",
            line=dict(color=depth_line_color, width=0),
            fillcolor=depth_fill_rgba,
            yaxis="y",
            hovertemplate="Position: %{x}<br>Repeat Count: %{y}<extra></extra>",
        )
    )
    figure.add_trace(
        go.Scatter(
            x=positions, y=breadth,
            name=f"{hla_label} Type Coverage Count (Breadth)",
            mode="lines",
            line=dict(color="#ef4444", width=2.8),
            yaxis="y2",
            hovertemplate=f"Position: %{{x}}<br>{hla_label} Type Count: %{{y}}<extra></extra>",
        )
    )

    breadth_max = int(breadth.max()) if breadth.size else 0
    figure.update_layout(
        title="<b>📊 Immunogenic Speculative Graph (Depth vs Breadth)</b>"
              f"<br><sup>Fill = Depth | Red line = {hla_label} Breadth</sup>",
        xaxis=dict(title="Amino Acid Position"),
        yaxis=dict(title="Total Detection/Repeat Count (Depth)", showgrid=False),
        yaxis2=dict(
            title=f"{hla_label} Types Recognizing This Position (Count)",
            overlaying="y", side="right", showgrid=False,
            range=[0, breadth_max * 1.2 if breadth_max else 5],
        ),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        template="plotly_white",
        height=520,
        hovermode="x unified",
    )
    return figure


# ---------------------------------------------------------------------------
# Application
# ---------------------------------------------------------------------------
def main() -> None:
    st.set_page_config(page_title="MHC Critical Analysis", layout="wide", page_icon="⚖️")
    st.markdown(TAB_STYLE, unsafe_allow_html=True)

    with st.sidebar:
        st.header("文件检查")
        st.write(f"Class I: {file_status(DATA_PATH_I)}")
        st.write(f"Class II: {file_status(DATA_PATH_II)}")
        st.caption("两个 CSV 需要和 exe 放在同一目录。")

    mhc_class = st.selectbox("Select MHC Class", ["Class I", "Class II"], index=0)

    if mhc_class == "Class I":
        data_path = DATA_PATH_I
        hla_label = "HLA I"
        depth_line_color = "#2563eb"
        depth_fill_rgba = "rgba(37, 99, 235, 0.28)"
        bar_scale = "Blues"
    else:
        data_path = DATA_PATH_II
        hla_label = "HLA II"
        depth_line_color = "#9333ea"
        depth_fill_rgba = "rgba(147, 51, 234, 0.28)"
        bar_scale = "RdPu"

    frame, error, dropped = load_and_clean_data(data_path)
    if error:
        st.error(error)
        st.stop()

    if dropped:
        st.warning(f"{dropped} row(s) without usable coordinates were excluded.")

    st.title(f"🧬 {mhc_class} Target Discovery: Depth vs Breadth")

    tab_global, tab_deep = st.tabs(
        ["🌍 Global View (Global Ranking)", "🔬 Deep Dive Analysis"]
    )

    leaderboard = build_leaderboard(frame)

    with tab_global:
        st.subheader("🏆 Global Gene Breadth Ranking")
        st.info(
            f"The ranking is based on how many types of {hla_label} recognize it, "
            "not just the data volume."
        )

        left, right = st.columns([2, 1])
        with left:
            figure = px.bar(
                leaderboard.head(20),
                x="Unique_HLA_Count", y="Gene", orientation="h", color="Total_Reports",
                color_continuous_scale=bar_scale,
                labels={"Unique_HLA_Count": f"{hla_label} Type Coverage Count",
                        "Gene": "Gene Name"},
            )
            figure.update_layout(yaxis={"categoryorder": "total ascending"})
            st.plotly_chart(figure, use_container_width=True)
        with right:
            st.dataframe(
                leaderboard.head(20)[["Gene", "Unique_HLA_Count", "Total_Reports"]],
                use_container_width=True, hide_index=True,
            )

    with tab_deep:
        st.subheader("🔎 Single Protein Deep Dive")
        all_genes = leaderboard["Gene"].tolist()
        if not all_genes:
            st.info("No annotated genes available in this table.")
            st.stop()

        control, selector = st.columns([1, 3])
        with control:
            selected_gene = st.selectbox("Step 1: Select Target Gene/Protein", all_genes)

        gene_rows = frame[frame["Gene"] == selected_gene]
        available_hlas = sorted(gene_rows["HLA"].unique().tolist())
        with selector:
            selected_hlas = st.multiselect(
                f"Step 2: Filter {hla_label} Subtypes",
                available_hlas,
                default=available_hlas,
            )

        if not selected_hlas:
            st.warning("请至少选择一个 HLA subtype。")
            return

        selected_rows = gene_rows[gene_rows["HLA"].isin(selected_hlas)]

        st.divider()
        st.markdown("#### 📈 Core Speculative Graph")
        st.plotly_chart(
            plot_speculative_landscape(
                selected_rows, hla_label, depth_line_color, depth_fill_rgba
            ),
            use_container_width=True,
        )

        st.divider()
        st.markdown("#### 🧪 Peptide-level inspection")
        if "Peptide" not in selected_rows.columns:
            st.info("Peptide sequences are not present in this table.")
            return

        peptide_table = (
            selected_rows.groupby("Peptide")
            .agg(
                Start=("Start", "min"),
                End=("End", "max"),
                HLA_Count=("HLA", "nunique"),
                Total_Reports=("Count", "sum"),
            )
            .reset_index()
            .sort_values(by=["HLA_Count", "Total_Reports"], ascending=[False, False])
        )
        st.dataframe(peptide_table, use_container_width=True, hide_index=True)


if __name__ == "__main__":
    main()
