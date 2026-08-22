@echo off
rem ============================================================
rem  Double-click this file to open the Galactic Plane HI Mapper.
rem  It runs the app with the correct Python (3.13) no matter
rem  where this folder lives.
rem ============================================================
title Galactic Plane 21-cm HI Mapper
cd /d "%~dp0"

rem  Find Python 3.13 (all dependencies are installed there).  Override by
rem  setting HI_PYTHON to a full python.exe path before running this file.
set "PY=%HI_PYTHON%"
if not defined PY for /f "delims=" %%i in ('py -3.13 -c "import sys;print(sys.executable)" 2^>nul') do set "PY=%%i"
if not defined PY for /f "delims=" %%i in ('py -3 -c "import sys;print(sys.executable)" 2^>nul') do set "PY=%%i"
if not defined PY set "PY=python"

"%PY%" -c "import sys" >nul 2>&1
if errorlevel 1 (
  echo Could not find a usable Python ^(tried: %PY%^).
  echo Install Python 3.13, or set HI_PYTHON to your python.exe path.
  pause
  exit /b 1
)

"%PY%" "%~dp0hydrogen_mapper.py"

if errorlevel 1 (
  echo.
  echo -----------------------------------------------------------
  echo The app closed with an error. Screenshot the lines above
  echo if you need help figuring out what went wrong.
  echo -----------------------------------------------------------
  pause
)
