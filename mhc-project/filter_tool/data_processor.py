"""Build class-stratified gold tables of naturally eluted HLA ligands.

The pipeline reads the IEDB ``mhc_ligand_full`` multi-file export together with
the UniProt human reference proteome (UP000005640) and writes two tables:

    Class_I_HomoSapiens_Gold_Gene.csv
    Class_II_HomoSapiens_Gold_Gene.csv

Each surviving row is one *presentation unit*: a unique
(protein, peptide, HLA restriction) combination, annotated with its gene name
and with the number of independent database records that support it.

Usage
-----
    python data_processor.py --input-dir .
    python data_processor.py --input-dir D:/data --output-dir out --no-download

When no ``--input-dir`` is given the directory containing this script is used,
so the module also works unchanged inside a frozen executable.
"""

from __future__ import annotations

import argparse
import glob
import gzip
import logging
import re
import shutil
import sys
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
import pandas as pd
import requests
from Bio import SeqIO

LOGGER = logging.getLogger("mhc_filter")

PROTEOME_ACCESSION = "UP000005640"
UNIPROT_STREAM_URL = (
    "https://rest.uniprot.org/uniprotkb/stream"
    "?compressed=true&format=fasta&includeIsoform=true"
    f"&query=%28proteome%3A{PROTEOME_ACCESSION}%29"
)

DEFAULT_REFERENCE_SHARD = "mhc_ligand_full_00.csv"
SHARD_GLOB = "mhc_ligand_full_*.csv"
PROTEOME_FASTA = "human_proteome.fasta"
PROTEOME_FASTA_GZ = "human_proteome.fasta.gz"

OUTPUT_CLASS_I = "Class_I_HomoSapiens_Gold_Gene.csv"
OUTPUT_CLASS_II = "Class_II_HomoSapiens_Gold_Gene.csv"

# IEDB restriction labels are parsed by pattern rather than by table lookup so
# that allele families added to the database are handled without code changes.
CLASS_I_PATTERN = r"HLA-[ABCEFG]"
CLASS_II_PATTERN = r"HLA-D"

# Assay methods that count as a direct experimental measurement (ligand
# elution or binding). This is a deliberately coarse two-tier ranking: tier 1
# wins over tier 2 when several records describe the same presentation unit,
# and the reported quantitative value only breaks ties inside a tier.
PRIMARY_METHOD_TERMS: Tuple[str, ...] = ("mass spec", "purified", "fp", "elisa")
METHOD_TIER_PRIMARY = 1
METHOD_TIER_FALLBACK = 2

# Evidence modality. This describes *what kind* of evidence a record represents
# and is kept orthogonal to the representative-selection ranking above: it is
# used for optional stratification, never mixed into the ranking.
MODALITY_ALL = "all"
MODALITY_ELUTED = "eluted"
MODALITY_BINDING = "binding"
MODALITY_OTHER = "other"

#: Assay methods that indicate the ligand was eluted from an MHC molecule and
#: detected by mass spectrometry -- direct evidence of presentation.
ELUTED_TERMS: Tuple[str, ...] = ("mass spectrometry", "coelution")

#: Assay methods that measure structure or T-cell function rather than
#: presentation or binding.
NON_PRESENTATION_TERMS: Tuple[str, ...] = (
    "x-ray crystallography",
    "edman degradation",
    "t cell",
)

# Composite key that defines a presentation unit.
PRESENTATION_UNIT_KEY = ("Mapped_Protein_ID", "Epitope_Name", "MHC Restriction_Name")

# Columns written to the gold tables, in order. Missing ones are skipped.
OUTPUT_COLUMNS = (
    "Epitope_Name",
    "Epitope_Starting Position",
    "Epitope_Ending Position",
    "Gene_Name",
    "Mapped_Protein_ID",
    "MHC Restriction_Name",
    "Record_Count",
    "Assay_Modality",
)

_GENE_NAME_PATTERN = re.compile(r"GN=([^\s]+)")


# ---------------------------------------------------------------------------
# Reference proteome
# ---------------------------------------------------------------------------
def download_proteome(destination: Path, url: str = UNIPROT_STREAM_URL) -> Path:
    """Download and decompress the human reference proteome.

    Returns the path to the uncompressed FASTA. Downloading is skipped when
    the decompressed file is already present.
    """
    fasta_path = destination / PROTEOME_FASTA
    if fasta_path.exists():
        return fasta_path

    gz_path = destination / PROTEOME_FASTA_GZ
    LOGGER.info("Downloading UniProt proteome %s", PROTEOME_ACCESSION)
    with requests.get(url, stream=True, timeout=120) as response:
        response.raise_for_status()
        with gz_path.open("wb") as handle:
            for chunk in response.iter_content(chunk_size=8192):
                handle.write(chunk)

    with gzip.open(gz_path, "rb") as source, fasta_path.open("wb") as target:
        shutil.copyfileobj(source, target)
    return fasta_path


def build_accession_index(fasta_path: Path) -> Dict[str, str]:
    """Map every UniProt accession in the FASTA to its gene name.

    Both the full accession (``P12345-2``) and the canonical accession
    (``P12345``) are indexed so that isoform identifiers reported by IEDB
    resolve to the same protein.
    """
    id_to_gene: Dict[str, str] = {}
    for record in SeqIO.parse(str(fasta_path), "fasta"):
        parts = record.id.split("|")
        accession = parts[1] if len(parts) >= 2 else record.id
        match = _GENE_NAME_PATTERN.search(record.description)
        gene_name = match.group(1) if match else "Unknown"

        id_to_gene[accession.upper()] = gene_name
        canonical = accession.split("-")[0].upper()
        id_to_gene.setdefault(canonical, gene_name)
    return id_to_gene


# ---------------------------------------------------------------------------
# Shard loading
# ---------------------------------------------------------------------------
def _read_reference_columns(shard_path: Path) -> pd.MultiIndex:
    """Read only the two-level header of the reference shard."""
    header = pd.read_csv(shard_path, header=[0, 1], nrows=0)
    return header.columns


def _load_shard(shard_path: Path, template_columns: pd.MultiIndex) -> Optional[pd.DataFrame]:
    """Load one shard, aligning it to the reference column layout.

    Only the reference shard carries a header row; the remaining shards are
    positional and are padded or truncated to the reference width.
    """
    expected = len(template_columns)
    try:
        if shard_path.name == DEFAULT_REFERENCE_SHARD:
            frame = pd.read_csv(shard_path, header=[0, 1], low_memory=False)
        else:
            frame = pd.read_csv(shard_path, header=None, low_memory=False)
            if frame.shape[1] > expected:
                frame = frame.iloc[:, :expected]
            elif frame.shape[1] < expected:
                for _ in range(expected - frame.shape[1]):
                    frame[len(frame.columns)] = np.nan
            frame.columns = template_columns
    except Exception as exc:  # a broken shard must not abort the whole run
        LOGGER.warning("Skipping %s: %s", shard_path.name, exc)
        return None

    if frame.shape[1] != expected:
        LOGGER.warning(
            "%s has %d columns, expected %d", shard_path.name, frame.shape[1], expected
        )
    return frame


def load_ligand_data(input_dir: Path, reference_shard: str) -> pd.DataFrame:
    """Concatenate every shard in ``input_dir`` into a single frame."""
    shard_paths = sorted(Path(p) for p in glob.glob(str(input_dir / SHARD_GLOB)))
    reference_path = input_dir / reference_shard

    if not shard_paths:
        raise FileNotFoundError(f"No {SHARD_GLOB} files in {input_dir}")
    if not reference_path.exists():
        raise FileNotFoundError(f"Reference shard not found: {reference_path}")

    template_columns = _read_reference_columns(reference_path)

    frames = []
    for path in shard_paths:
        LOGGER.info("Reading %s", path.name)
        frame = _load_shard(path, template_columns)
        if frame is not None:
            frames.append(frame)

    if not frames:
        raise RuntimeError("No shard could be read")

    data = pd.concat(frames, ignore_index=True)
    data.columns = [
        f"{str(top).strip()}_{str(bottom).strip()}" for top, bottom in data.columns
    ]
    LOGGER.info("Loaded %s records from %d shards", f"{len(data):,}", len(frames))
    return data


def locate_columns(data: pd.DataFrame) -> Dict[str, Optional[str]]:
    """Resolve every column the pipeline depends on, by substring match."""
    def first(*needles: str) -> Optional[str]:
        return next((c for c in data.columns if all(n in c for n in needles)), None)

    columns = {
        "sequence": first("Epitope_Name"),
        "mhc": first("MHC Restriction_Name"),
        "outcome": first("Assay_Qualitative Measurement"),
        "quantity": first("Assay_Quantitative measurement"),
        "method": first("Assay_Method"),
        "host": first("Host_Name"),
        "iri": "Epitope_Source Molecule IRI",
        "start": first("Epitope_Start"),
        "end": first("Epitope_End"),
    }
    required = ("sequence", "mhc", "outcome", "method", "host")
    missing = [name for name in required if columns[name] is None]
    if missing:
        raise KeyError(f"Required columns not found: {', '.join(missing)}")
    if columns["iri"] not in data.columns:
        raise KeyError(f"Required column not found: {columns['iri']}")
    return columns


# ---------------------------------------------------------------------------
# Selection, resolution, deduplication
# ---------------------------------------------------------------------------
def filter_primary_records(
    data: pd.DataFrame, columns: Dict[str, Optional[str]]
) -> pd.DataFrame:
    """Keep positive assays performed under a human host context."""
    outcome = data[columns["outcome"]].astype(str)
    host = data[columns["host"]].astype(str)

    mask = outcome.str.contains("Positive", case=False, na=False)
    mask &= host.str.contains("Human|Homo sapiens", case=False, na=False)

    filtered = data.loc[mask].copy()
    LOGGER.info("Positive human-host records: %s", f"{len(filtered):,}")
    return filtered


def resolve_protein_ids(
    data: pd.DataFrame,
    iri_column: str,
    valid_accessions: Iterable[str],
) -> pd.Series:
    """Map each IEDB source-molecule IRI onto a validated UniProt accession.

    Records whose accession is not part of the human proteome come back as
    ``NA``. Both the full (``P12345-2``) and canonical (``P12345``) accession
    are accepted.
    """
    valid = set(valid_accessions)

    iri = data[iri_column].astype("string").str.strip()
    is_iri = (
        iri.notna()
        & (iri != "")
        & (iri.str.lower() != "nan")
        & iri.str.contains("/", regex=False, na=False)
    ).fillna(False)

    tail = iri.str.rsplit("/", n=1).str[-1].str.upper().where(is_iri)
    canonical = tail.str.split("-").str[0].str.split(".").str[0]

    resolved = tail.where(tail.isin(valid))
    return resolved.where(resolved.notna(), canonical.where(canonical.isin(valid)))


def score_methods(methods: pd.Series) -> np.ndarray:
    """Two-tier priority: direct experimental measurement vs. everything else."""
    text = methods.astype("string").str.lower()
    pattern = "|".join(re.escape(term) for term in PRIMARY_METHOD_TERMS)
    is_primary = text.str.contains(pattern, regex=True, na=False).to_numpy()
    return np.where(is_primary, METHOD_TIER_PRIMARY, METHOD_TIER_FALLBACK)


def classify_modality(methods: pd.Series) -> np.ndarray:
    """Label each record by the kind of evidence its assay provides.

    ``eluted``
        Mass spectrometry of ligands eluted from MHC molecules -- direct
        evidence that the peptide was processed and presented.
    ``binding``
        In vitro binding assays -- evidence of affinity, not of presentation.
    ``other``
        Structural or functional assays that measure neither.
    """
    text = methods.astype("string").str.lower()

    def matches(terms: Tuple[str, ...]) -> np.ndarray:
        pattern = "|".join(re.escape(term) for term in terms)
        return text.str.contains(pattern, regex=True, na=False).to_numpy()

    labels = np.full(len(text), MODALITY_BINDING, dtype=object)
    labels[matches(NON_PRESENTATION_TERMS)] = MODALITY_OTHER
    labels[matches(ELUTED_TERMS)] = MODALITY_ELUTED
    return labels


def deduplicate(
    data: pd.DataFrame,
    sequence_column: str,
    quantity_column: Optional[str],
) -> pd.DataFrame:
    """Collapse each presentation unit onto its best-supported record.

    ``Record_Count`` is computed *before* the collapse so that the retained row
    still carries the number of independent records behind it.
    """
    key = list(PRESENTATION_UNIT_KEY)
    data = data.assign(
        Record_Count=data.groupby(key)[sequence_column].transform("count")
    )

    order = key + ["__method_tier"]
    if quantity_column and quantity_column in data.columns:
        order.append(quantity_column)
    data = data.sort_values(by=order, ascending=True)

    deduplicated = data.drop_duplicates(subset=key, keep="first")
    LOGGER.info("Presentation units after deduplication: %s", f"{len(deduplicated):,}")
    return deduplicated


def export_class_tables(
    data: pd.DataFrame,
    columns: Dict[str, Optional[str]],
    output_dir: Path,
) -> Tuple[int, int]:
    """Write the class I and class II gold tables."""
    present = [c for c in OUTPUT_COLUMNS if c in data.columns]
    restrictions = data[columns["mhc"]].astype(str).str.upper()

    class_i = data.loc[restrictions.str.contains(CLASS_I_PATTERN, regex=True), present]
    class_ii = data.loc[restrictions.str.contains(CLASS_II_PATTERN, regex=True), present]

    class_i_path = output_dir / OUTPUT_CLASS_I
    class_ii_path = output_dir / OUTPUT_CLASS_II
    class_i.to_csv(class_i_path, index=False, encoding="utf-8-sig")
    class_ii.to_csv(class_ii_path, index=False, encoding="utf-8-sig")

    LOGGER.info("Wrote %s (%s rows)", class_i_path.name, f"{len(class_i):,}")
    LOGGER.info("Wrote %s (%s rows)", class_ii_path.name, f"{len(class_ii):,}")
    return len(class_i), len(class_ii)


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------
def run_pipeline(
    input_dir: Path,
    output_dir: Path,
    reference_shard: str = DEFAULT_REFERENCE_SHARD,
    fasta_path: Optional[Path] = None,
    allow_download: bool = True,
    modality: str = MODALITY_ALL,
) -> Tuple[int, int]:
    """Run the full filtering pipeline and return the two output row counts."""
    output_dir.mkdir(parents=True, exist_ok=True)

    if fasta_path is None:
        fasta_path = input_dir / PROTEOME_FASTA
    if not fasta_path.exists():
        if not allow_download:
            raise FileNotFoundError(
                f"{fasta_path} not found and downloading is disabled"
            )
        fasta_path = download_proteome(input_dir)

    LOGGER.info("Indexing reference proteome")
    id_to_gene = build_accession_index(fasta_path)
    LOGGER.info("Indexed %s accessions", f"{len(id_to_gene):,}")

    data = load_ligand_data(input_dir, reference_shard)
    columns = locate_columns(data)

    data = filter_primary_records(data, columns)
    data["Assay_Modality"] = classify_modality(data[columns["method"]])

    if modality != MODALITY_ALL:
        before = len(data)
        data = data[data["Assay_Modality"] == modality].copy()
        LOGGER.info(
            "Restricted to %s evidence: %s -> %s records",
            modality, f"{before:,}", f"{len(data):,}",
        )

    data["Mapped_Protein_ID"] = resolve_protein_ids(
        data, columns["iri"], id_to_gene.keys()
    )
    data = data.dropna(subset=["Mapped_Protein_ID"]).copy()
    data["Gene_Name"] = data["Mapped_Protein_ID"].map(id_to_gene)
    LOGGER.info("Records with a validated human accession: %s", f"{len(data):,}")

    if columns["quantity"] in data.columns:
        data[columns["quantity"]] = pd.to_numeric(
            data[columns["quantity"]], errors="coerce"
        )
    data["__method_tier"] = score_methods(data[columns["method"]])

    data = deduplicate(data, columns["sequence"], columns["quantity"])
    return export_class_tables(data, columns, output_dir)


def _default_input_dir() -> Path:
    """Directory of the script, or of the executable when frozen."""
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parent


def _pause_if_frozen() -> None:
    """Keep the console window open when the tool is launched by double-click."""
    if getattr(sys, "frozen", False):
        try:
            input("Press Enter to exit...")
        except EOFError:
            pass


def parse_args(argv: Optional[list] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build class-stratified gold tables of eluted HLA ligands.",
    )
    parser.add_argument(
        "--input-dir", type=Path, default=None,
        help="Directory holding the mhc_ligand_full_*.csv shards (default: script folder)",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=None,
        help="Directory for the gold tables (default: same as --input-dir)",
    )
    parser.add_argument(
        "--reference-shard", default=DEFAULT_REFERENCE_SHARD,
        help="Shard that defines the canonical column layout",
    )
    parser.add_argument(
        "--fasta", type=Path, default=None,
        help="Pre-downloaded human proteome FASTA",
    )
    parser.add_argument(
        "--no-download", action="store_true",
        help="Fail instead of downloading the proteome when it is missing",
    )
    parser.add_argument(
        "--modality", choices=[MODALITY_ALL, MODALITY_ELUTED, MODALITY_BINDING],
        default=MODALITY_ALL,
        help=(
            "Restrict the gold tables to one kind of evidence. 'eluted' keeps "
            "mass-spectrometry elution data (direct evidence of presentation); "
            "'binding' keeps in vitro binding assays; 'all' (default) keeps every "
            "record and only labels it in the Assay_Modality column"
        ),
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    return parser.parse_args(argv)


def main(argv: Optional[list] = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(message)s",
    )

    input_dir = (args.input_dir or _default_input_dir()).resolve()
    output_dir = (args.output_dir or input_dir).resolve()
    LOGGER.info("Input directory: %s", input_dir)

    try:
        run_pipeline(
            input_dir=input_dir,
            output_dir=output_dir,
            reference_shard=args.reference_shard,
            fasta_path=args.fasta,
            allow_download=not args.no_download,
            modality=args.modality,
        )
    except Exception:
        LOGGER.exception("Processing failed")
        _pause_if_frozen()
        return 1

    _pause_if_frozen()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
