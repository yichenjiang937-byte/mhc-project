# MHC Ligand Filtering and Visualization

Two connected tools:

- **filter_tool** — turns the raw IEDB `mhc_ligand_full_*.csv` shards plus the
  UniProt human proteome into `Class_I_HomoSapiens_Gold_Gene.csv` and
  `Class_II_HomoSapiens_Gold_Gene.csv`.
- **visualizer_tool** — reads those gold tables and plots peptide depth against
  HLA allele breadth along a protein sequence.

Full documentation, usage and citation: see the [repository README](../README.md).

```text
Raw mhc_ligand_full_*.csv + human_proteome.fasta
                ↓
           filter_tool
                ↓
Class_I_HomoSapiens_Gold_Gene.csv
Class_II_HomoSapiens_Gold_Gene.csv
                ↓
         visualizer_tool
```
