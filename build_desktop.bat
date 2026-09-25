@echo off
echo ========================================================
echo   Bindora Dock - Desktop Application Build Script
echo ========================================================
echo.

echo [1/3] Generating installer wizard images...
python scripts/generate_installer_assets.py
if errorlevel 1 goto error

echo.
echo [2/3] Building standalone executable with PyInstaller...
python -m PyInstaller -y --clean bindora.spec
if errorlevel 1 goto error

echo.
echo [3/3] Checking for Inno Setup compiler (ISCC.exe)...
where iscc >nul 2>nul
if %errorlevel% equ 0 (
    echo Compiling Inno Setup installer...
    iscc installer.iss
    if errorlevel 1 goto error
    echo.
    echo [SUCCESS] Installer generated in installer\Output\BindoraDock-Setup.exe
) else (
    echo [NOTE] Inno Setup compiler (iscc) not found in PATH.
    echo Standalone application folder is ready at: dist\bindora_launcher\
    echo You can test the app now: dist\bindora_launcher\bindora_launcher.exe
    echo To compile the Setup.exe wizard, install Inno Setup 6 from https://jrsoftware.org/isdl.php
)

echo.
echo ========================================================
echo   Build Finished Successfully!
echo ========================================================
exit /b 0

:error
echo.
echo [ERROR] Build failed! Check the output logs above.
exit /b 1
