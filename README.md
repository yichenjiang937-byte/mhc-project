# MHC Ligand Filtering and Visualization

A two-step toolkit for turning the raw IEDB MHC ligand export into an
interactive map of antigen-presentation hotspots.

1. **Filter tool** — reads the raw `mhc_ligand_full_*.csv` shards together with
   the UniProt human reference proteome, keeps high-confidence human ligand
   records, resolves their source proteins, and writes two class-stratified
   "gold" tables.
2. **Visualizer tool** — reads those gold tables and renders, for any protein,
   peptide depth and HLA allele breadth along the amino-acid sequence.

The companion paper is published in the *Proceedings of the 2026 18th
International Conference on Bioinformatics and Biomedical Technology*
(ICBBT 2026, Suzhou, China), ACM, pages 78–82:
[doi:10.1145/3822521.3822541](https://doi.org/10.1145/3822521.3822541).

---

## Repository layout

```text
mhc-project/
├─ README.md
├─ CHANGELOG.md
├─ .gitignore
├─ mhc-project/
│  ├─ readme.md
│  ├─ filter_tool/
│  │  ├─ data_processor.py
│  │  ├─ requirements.txt
│  │  ├─ build_exe.bat
│  │  └─ README.md
│  └─ visualizer_tool/
│     ├─ app_visualization.py
│     ├─ main_visualizer.py
│     ├─ requirements.txt
│     ├─ start_visualizer.bat
│     ├─ build_exe.bat
│     └─ README.md
```

---

## Requirements

Python 3.9 or newer. Install the dependencies of each tool separately:

```bash
pip install -r mhc-project/filter_tool/requirements.txt
pip install -r mhc-project/visualizer_tool/requirements.txt
```

---

## Quick start

### 1. Build the gold tables

Put the shards and the reference FASTA in one directory, then run:

```bash
cd mhc-project/filter_tool
python data_processor.py --input-dir /path/to/data --output-dir /path/to/output
```

If `human_proteome.fasta` is missing it is downloaded once from UniProt
(UP000005640, including isoforms). Use `--no-download` to fail instead, or
`--fasta` to point at a local copy.

Command-line options:

| Option | Description |
| --- | --- |
| `--input-dir` | Directory holding `mhc_ligand_full_*.csv` (default: script folder) |
| `--output-dir` | Directory for the gold tables (default: same as `--input-dir`) |
| `--reference-shard` | Shard that defines the canonical column layout (default `mhc_ligand_full_00.csv`) |
| `--fasta` | Pre-downloaded human proteome FASTA |
| `--no-download` | Do not download the proteome |
| `--modality` | `all` (default), `eluted`, or `binding` — see [Evidence modality](#evidence-modality) |
| `-v`, `--verbose` | Verbose logging |

On Windows, `build_exe.bat` produces a standalone `MHC_Data_Processor.exe`.

### 2. Explore the results

Put both gold tables next to the application, then either double-click
`start_visualizer.bat` or run:

```bash
cd mhc-project/visualizer_tool
python main_visualizer.py
```

The browser opens at `http://127.0.0.1:8501`. Set `MHC_VISUALIZER_PORT` to use
a different port. `build_exe.bat` produces a standalone `MHC_Visualizer.exe`.

---

## Input data

The raw shards come from the IEDB **MHC Ligand Export (v3)**, dataset
`mhc_ligand_full`, distributed as nine CSV files (`mhc_ligand_full_00.csv` …
`mhc_ligand_full_08.csv`). The export uses a two-level header, and only the
first shard carries one; the remaining shards are positional and are aligned to
the reference column layout at load time.

The reference proteome is the UniProt human proteome `UP000005640`, downloaded
with isoforms included.

Raw inputs and the generated gold tables are archived on
[Figshare](https://doi.org/10.6084/m9.figshare.31566175).

---

## How a record is selected

1. **Primary filter** — keep assays with a positive qualitative outcome that
   were performed in a human host context.
2. **Source-protein resolution** — parse the UniProt accession out of the
   `Epitope_Source Molecule IRI`, accepting both full (`P12345-2`) and canonical
   (`P12345`) accessions. Records whose accession is absent from the human
   proteome are dropped.
3. **Evidence accounting** — a presentation unit is the composite key
   *(protein, peptide, HLA restriction)*. Its `Record_Count` is the number of
   independent database records behind it.
4. **Representative selection** — of all records for one presentation unit, the
   one with the strongest experimental evidence is kept, using reported
   quantitative values to break ties within an evidence tier.
5. **Class split** — restrictions matching `HLA-[ABCEFG]` go to the class I
   table, those matching `HLA-D` to the class II table.

### Output schema

| Column | Meaning |
| --- | --- |
| `Epitope_Name` | Peptide sequence |
| `Epitope_Starting Position` | Start coordinate on the source protein |
| `Epitope_Ending Position` | End coordinate on the source protein |
| `Gene_Name` | Gene name parsed from the UniProt `GN=` field |
| `Mapped_Protein_ID` | Validated UniProt accession |
| `MHC Restriction_Name` | HLA restriction as reported by IEDB |
| `Record_Count` | Number of independent records for this presentation unit |
| `Assay_Modality` | `eluted`, `binding` or `other` — the kind of evidence behind the retained record |

### Evidence modality

Records come from two scientifically different kinds of experiment, and mixing
them makes `Record_Count` mean two things at once:

| Modality | IEDB assay examples | What it establishes |
| --- | --- | --- |
| `eluted` | `cellular mhc/mass spectrometry`, `secreted mhc/mass spectrometry`, `coelution` | The peptide was processed and presented by a cell |
| `binding` | `purified mhc/competitive/radioactivity`, `purified mhc/direct/fluorescence`, `binding assay` | The peptide binds the allele in vitro |
| `other` | `x-ray crystallography`, `edman degradation`, T-cell assays | Neither |

The default, `--modality all`, keeps every record and only labels it, which is
the behaviour used for the published tables. `--modality eluted` restricts the
gold tables to elution evidence, which is the appropriate input when the
analysis is about presentation rather than binding affinity. On a single-shard
test this retained 96.6% of class I records, and the top-20 gene ranking
overlapped 19 of 20 entries with the unrestricted run.

---

## Reading the visualization

The deep-dive view plots two quantities per amino-acid position:

- **Depth** (filled area) — summed report count of peptides covering the
  position, i.e. how strongly it is represented in the database.
- **Breadth** (red line) — number of distinct HLA restrictions covering it.

Regions where both rise are pan-allelic candidates; regions where only depth
rises are more consistent with allele-restricted presentation. Deselecting
individual HLA subtypes in the sidebar is a quick way to test whether an
apparent hotspot survives removal of a few alleles.

---

## Citation

```bibtex
@inproceedings{jiang2026mhc,
  title     = {A Lightweight Visual Analytics Tool for Regional Enrichment
               Analysis of Immunopeptidomic Data},
  author    = {Jiang, Yichen and Xie, Menghan and Yuan, Xuye and
               Meng, Jia and Cai, Yuheng},
  booktitle = {Proceedings of the 2026 18th International Conference on
               Bioinformatics and Biomedical Technology},
  pages     = {78--82},
  year      = {2026},
  month     = may,
  address   = {Suzhou, China},
  publisher = {Association for Computing Machinery},
  doi       = {10.1145/3822521.3822541}
}
```

---

## Data availability

- Raw IEDB export and generated gold tables:
  <https://doi.org/10.6084/m9.figshare.31566175>
