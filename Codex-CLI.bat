@echo off
setlocal

REM Run in the folder where this .bat lives
cd /d "%~dp0"

codex --yolo -m gpt-5.2-codex ^
  -c model_reasoning_effort="high" ^
  -c model_reasoning_summary_format=experimental %*

pause
