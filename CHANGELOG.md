# Changelog

## 2026-11 — Evidence modality

Adds an explicit description of *what kind* of evidence each record represents,
which the previous evidence ranking conflated.

- New `Assay_Modality` column: `eluted` (mass-spectrometry elution),
  `binding` (in vitro binding assays) or `other` (structural / functional
  assays).
- New `--modality {all,eluted,binding}` option. The default, `all`, keeps every
  record and only labels it, so the published tables are reproduced unchanged.
- The modality is intentionally **not** part of the representative-selection
  ranking. The two questions are separate: "which record represents this
  presentation unit" and "what kind of evidence supports it".

Backwards compatibility: the default run produces the same rows in the same
order with identical values in every existing column; the only difference is
the appended `Assay_Modality` column. Verified on `mhc_ligand_full_00.csv`
(59,083 class I and 22,227 class II rows, original columns identical).
`--modality eluted` yields 57,055 class I and 19,567 class II rows on the same
shard.

## 2026-11 — Refactor and hardening

A quality pass over both tools. **Filtering behaviour is unchanged**: the
updated pipeline reproduces the published gold tables exactly (verification
below).

### Filter tool

- Vectorised protein-accession resolution, gene mapping and evidence tiering.
  These now operate on whole columns rather than per-row callables, removing
  the main scaling bottleneck on the full nine-shard export.
- Split the procedural block into named, individually callable stages:
  `load_ligand_data` → `locate_columns` → `filter_primary_records` →
  `resolve_protein_ids` → `deduplicate` → `export_class_tables`.
- Added a command-line interface (`--input-dir`, `--output-dir`, `--fasta`,
  `--no-download`, `--verbose`) so the tool can run unattended.
- The interactive pause is now applied only to the frozen executable, letting
  scripted runs terminate on their own while double-clicked runs keep their
  console open.
- Shard-level read failures are reported; a missing required column raises with
  the column name rather than failing downstream.
- Added `requirements.txt` and `.gitignore`.

### Visualizer

- Depth and breadth are accumulated with interval difference arrays instead of
  a per-peptide Python loop, tying cost to the number of peptides rather than
  to their length.
- Coordinates are validated explicitly; rows that cannot be placed are excluded
  and counted, and the count is surfaced in the UI.
- Missing or renamed columns produce a readable message instead of a `KeyError`.
- The Streamlit page is wrapped in `main()` and the leaderboard is cached, so
  widget re-renders no longer re-aggregate the full table.
- `MHC_VISUALIZER_PORT` overrides the default port.

### Verification

Both versions of the filter were run against `mhc_ligand_full_00.csv`
(599,998 records) with the same reference proteome:

| Output | Rows | Checksum |
| --- | --- | --- |
| `Class_I_HomoSapiens_Gold_Gene.csv` | 59,083 | identical before and after |
| `Class_II_HomoSapiens_Gold_Gene.csv` | 22,227 | identical before and after |

The visualizer's position profile was checked against the previous
implementation on EEF2, KRT1, MYH9, CD74 and PRKDC; depth and breadth arrays
match exactly in every case.
