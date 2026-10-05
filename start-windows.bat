@echo off
rem EZplot launcher for Windows: double-click this file.
rem First run installs uv (a small Python manager, in your user folder, no admin rights)
rem and the app's packages (~2 minutes). Later runs start in a few seconds.
title EZplot
cd /d "%~dp0"
set "PATH=%USERPROFILE%\.local\bin;%PATH%"
where uv >nul 2>nul
if errorlevel 1 (
  echo First run: installing uv ^(one time only^)...
  powershell -NoProfile -ExecutionPolicy ByPass -Command "$env:UV_NO_MODIFY_PATH=1; irm https://astral.sh/uv/install.ps1 | iex"
)
echo Starting EZplot ^(the first start downloads Python and packages, please wait^)...
uv run --frozen --no-dev ezplot %*
if errorlevel 1 (
  echo.
  echo EZplot stopped with an error ^(see above^).
  pause
)
