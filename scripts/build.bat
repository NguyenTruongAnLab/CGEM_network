@echo off
setlocal enabledelayedexpansion

REM ===========================================================================
REM CGEM Network Build Script
REM ===========================================================================
REM 
REM Structure:
REM   src/              - Core files (main, network, init)
REM   src/physics/      - Hydrodynamics + transport (Saint-Venant, advection-dispersion)
REM   src/rive/         - RIVE biogeochemistry + sediment
REM   src/io/           - Input/output files
REM 
REM ===========================================================================

set "rootDir=%~dp0.."
set "srcDir=%rootDir%\src"
set "binDir=%rootDir%\bin\Debug"
set "exePath=%binDir%\CGEM_Network.exe"
set "nloptDir=%rootDir%\external\nlopt"

REM Compiler settings
set "CC=C:\msys64\mingw64\bin\gcc.exe"
set "CFLAGS=-std=c11 -O2 -Wall -Wextra"
set "INCLUDES=-I%srcDir% -I%srcDir%\rive -I%srcDir%\physics -I%srcDir%\io -I%nloptDir%\src\api"
set "LDFLAGS=-L%nloptDir%\build -lnlopt -lm"

if not exist "%binDir%" mkdir "%binDir%"

echo ===========================================================================
echo Building CGEM Network
echo ===========================================================================

set "objs="

REM ---------------------------------------------------------------------------
REM Compile core source files (src/*.c)
REM ---------------------------------------------------------------------------
echo.
echo [1/4] Compiling core modules...
for %%f in ("%srcDir%\*.c") do (
    set "obj=%binDir%\%%~nf.o"
    echo   Compiling %%~nf.c...
    "%CC%" %CFLAGS% %INCLUDES% -c "%%f" -o "!obj!"
    if errorlevel 1 goto :error
    set "objs=!objs! !obj!"
)

REM ---------------------------------------------------------------------------
REM Compile physics module (src/physics/*.c) - Hydrodynamics + Transport
REM ---------------------------------------------------------------------------
echo.
echo [2/4] Compiling physics module (hydrodynamics + transport)...
if exist "%srcDir%\physics\*.c" (
    for %%f in ("%srcDir%\physics\*.c") do (
        set "obj=%binDir%\phys_%%~nf.o"
        echo   Compiling physics/%%~nf.c...
        "%CC%" %CFLAGS% %INCLUDES% -c "%%f" -o "!obj!"
        if errorlevel 1 goto :error
        set "objs=!objs! !obj!"
    )
)

REM ---------------------------------------------------------------------------
REM Compile RIVE biogeochemistry + sediment module (src/rive/*.c)
REM ---------------------------------------------------------------------------
echo.
echo [3/4] Compiling RIVE module (biogeochemistry + sediment)...
if exist "%srcDir%\rive\*.c" (
    for %%f in ("%srcDir%\rive\*.c") do (
        set "obj=%binDir%\rive_%%~nf.o"
        echo   Compiling rive/%%~nf.c...
        "%CC%" %CFLAGS% %INCLUDES% -c "%%f" -o "!obj!"
        if errorlevel 1 goto :error
        set "objs=!objs! !obj!"
    )
)

REM ---------------------------------------------------------------------------
REM Compile I/O module (src/io/*.c)
REM ---------------------------------------------------------------------------
echo.
echo [4/5] Compiling I/O module...
if exist "%srcDir%\io\*.c" (
    for %%f in ("%srcDir%\io\*.c") do (
        set "obj=%binDir%\io_%%~nf.o"
        echo   Compiling io/%%~nf.c...
        "%CC%" %CFLAGS% %INCLUDES% -c "%%f" -o "!obj!"
        if errorlevel 1 goto :error
        set "objs=!objs! !obj!"
    )
)

REM ---------------------------------------------------------------------------
REM Compile optimization/calibration module (src/optimization/*.c)
REM ---------------------------------------------------------------------------
echo.
echo [5/5] Compiling optimization module (calibration)...
if exist "%srcDir%\optimization\*.c" (
    for %%f in ("%srcDir%\optimization\*.c") do (
        set "obj=%binDir%\opt_%%~nf.o"
        echo   Compiling optimization/%%~nf.c...
        "%CC%" %CFLAGS% %INCLUDES% -c "%%f" -o "!obj!"
        if errorlevel 1 goto :error
        set "objs=!objs! !obj!"
    )
)

REM ---------------------------------------------------------------------------
REM Link all object files (with NLopt library)
REM ---------------------------------------------------------------------------
echo.
echo Linking with NLopt library...
"%CC%" -o "%exePath%" %objs% %LDFLAGS%
if errorlevel 1 goto :error

REM Copy NLopt DLL to bin directory if not already there
if not exist "%binDir%\libnlopt.dll" (
    echo Copying libnlopt.dll to bin directory...
    copy "%nloptDir%\build\libnlopt.dll" "%binDir%\" >nul
)
if errorlevel 1 goto :error

echo.
echo ===========================================================================
echo Build completed successfully.
echo Executable: %exePath%
echo ===========================================================================
goto :eof

:error
echo.
echo ===========================================================================
echo BUILD FAILED
echo ===========================================================================
exit /b 1