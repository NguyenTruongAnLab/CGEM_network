@echo off
setlocal enabledelayedexpansion

set "srcDir=%~dp0..\src"
set "binDir=%~dp0..\bin\Debug"
set "exePath=%~dp0..\bin\Debug\CGEM_Network.exe"

if not exist "%binDir%" mkdir "%binDir%"

echo Building CGEM Network...

set "objs="
for %%f in ("%srcDir%\*.c") do (
    set "obj=%binDir%\%%~nf.o"
    echo Compiling %%f...
    "C:\msys64\mingw64\bin\gcc.exe" -std=c11 -O2 -Wall -Wextra -c "%%f" -o "!obj!"
    if errorlevel 1 goto :error
    set "objs=!objs! !obj!"
)

echo Linking...
"C:\msys64\mingw64\bin\gcc.exe" -o "%exePath%" %objs%
if errorlevel 1 goto :error

echo Build completed successfully.
goto :eof

:error
echo Build failed.
exit /b 1