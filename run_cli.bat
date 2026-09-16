@echo off
title Bindora Dock v2.0 - Research Terminal CLI
cd /d "%~dp0"
cls
python bindora_cli.py
if errorlevel 1 (
    echo.
    echo An error occurred. Press any key to exit...
    pause >nul
)
