# Visualizer tool

Interactive Streamlit application for the gold tables produced by `filter_tool`.

## Files

| File | Purpose |
| --- | --- |
| `app_visualization.py` | Streamlit application |
| `main_visualizer.py` | Launcher: starts the local server and opens the browser |
| `requirements.txt` | Python dependencies |
| `start_visualizer.bat` | Runs the app with `streamlit run` |
| `build_exe.bat` | Builds a standalone `MHC_Visualizer.exe` |

## Usage

Put both gold tables next to the application:

```text
Class_I_HomoSapiens_Gold_Gene.csv
Class_II_HomoSapiens_Gold_Gene.csv
```

Then either double-click `start_visualizer.bat`, or run:

```bash
pip install -r requirements.txt
python main_visualizer.py
```

The browser opens at `http://127.0.0.1:8501`. Set `MHC_VISUALIZER_PORT` to use
a different port.

On Windows, `build_exe.bat` builds the standalone executable; copy the two CSV
files next to `MHC_Visualizer.exe` before running it.

## Views

- **Global View** — genes ranked by the number of distinct HLA restrictions,
  with total report count as a secondary encoding.
- **Deep Dive** — for one protein, peptide depth (filled area) and HLA breadth
  (red line) along the amino-acid sequence. Filtering HLA subtypes on the left
  is a quick way to test whether a hotspot depends on a few alleles.
