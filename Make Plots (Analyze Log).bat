@echo off
rem ============================================================
rem  Double-click to turn galactic_plane_log.csv into plots:
rem    lv_diagram.png      (spiral-arm / longitude-velocity)
rem    rotation_curve.png  (galactic rotation curve)
rem  The plots also pop up on screen.
rem ============================================================
title Analyze HI Log
cd /d "%~dp0"

rem  Find Python 3.13 (all dependencies are installed there).  Override by
rem  setting HI_PYTHON to a full python.exe path before running this file.
set "PY=%HI_PYTHON%"
if not defined PY for /f "delims=" %%i in ('py -3.13 -c "import sys;print(sys.executable)" 2^>nul') do set "PY=%%i"
if not defined PY for /f "delims=" %%i in ('py -3 -c "import sys;print(sys.executable)" 2^>nul') do set "PY=%%i"
if not defined PY set "PY=python"

"%PY%" "%~dp0analyze_log.py" --show

echo.
pause
