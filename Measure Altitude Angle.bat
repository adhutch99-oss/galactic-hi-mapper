@echo off
rem ============================================================
rem  Double-click to work out the dish altitude angle from a
rem  couple of tape-measure numbers (measure along the feed boom).
rem ============================================================
title Dish Altitude Angle Helper
cd /d "%~dp0"

rem  Find Python 3.13 (all dependencies are installed there).  Override by
rem  setting HI_PYTHON to a full python.exe path before running this file.
set "PY=%HI_PYTHON%"
if not defined PY for /f "delims=" %%i in ('py -3.13 -c "import sys;print(sys.executable)" 2^>nul') do set "PY=%%i"
if not defined PY for /f "delims=" %%i in ('py -3 -c "import sys;print(sys.executable)" 2^>nul') do set "PY=%%i"
if not defined PY set "PY=python"

"%PY%" "%~dp0angle_calc.py"

echo.
pause
