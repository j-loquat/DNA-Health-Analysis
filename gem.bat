@echo off
setlocal

REM Run in the folder where this .bat lives
cd /d "%~dp0"

gemini --yolo %*

pause
