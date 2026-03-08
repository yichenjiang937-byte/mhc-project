@echo off
cd /d "%~dp0"

echo [1/3] Installing packages...
py -m pip install --upgrade pip
py -m pip install pandas requests biopython numpy pyinstaller

echo [2/3] Cleaning old build...
if exist build rmdir /s /q build
if exist dist rmdir /s /q dist
if exist MHC_Data_Processor.spec del /q MHC_Data_Processor.spec

echo [3/3] Building exe...
py -m PyInstaller --noconfirm --onefile --name MHC_Data_Processor data_processor.py

echo.
echo Build done. Put your mhc_ligand_full_*.csv files into dist folder.
pause
