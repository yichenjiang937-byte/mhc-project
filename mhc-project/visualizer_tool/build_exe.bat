@echo off
cd /d "%~dp0"

echo [1/3] Installing packages...
py -m pip install --upgrade pip
py -m pip install streamlit pandas plotly numpy pyinstaller

echo [2/3] Cleaning old build...
if exist build rmdir /s /q build
if exist dist rmdir /s /q dist
if exist MHC_Visualizer.spec del /q MHC_Visualizer.spec

echo [3/3] Building exe...
py -m PyInstaller --noconfirm --onefile --name MHC_Visualizer --add-data "app_visualization.py;." main_visualizer.py

echo.
echo Build done. Put these files into dist folder:
echo   Class_I_HomoSapiens_Gold_Gene.csv
echo   Class_II_HomoSapiens_Gold_Gene.csv
pause
