@echo off
rem ============================================================
rem  Double-click to make the SAME plots as "Make Plots
rem  (Analyze Log).bat", but using ONLY real-hardware rows:
rem  every row tagged Source=Simulation is ignored.
rem
rem  Use THIS one once the real telescope is running, so old
rem  simulation practice rows never mix into your science plots.
rem ============================================================
title Analyze HI Log (Hardware Only)
cd /d "%~dp0"

rem  Find Python 3.13 (all dependencies are installed there).  Override by
rem  setting HI_PYTHON to a full python.exe path before running this file.
set "PY=%HI_PYTHON%"
if not defined PY for /f "delims=" %%i in ('py -3.13 -c "import sys;print(sys.executable)" 2^>nul') do set "PY=%%i"
if not defined PY for /f "delims=" %%i in ('py -3 -c "import sys;print(sys.executable)" 2^>nul') do set "PY=%%i"
if not defined PY set "PY=python"

"%PY%" "%~dp0analyze_log.py" --show --hardware-only

echo.
pause
