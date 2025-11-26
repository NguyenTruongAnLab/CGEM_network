param(
    [string]$CaseConfig,
    [Alias("r")]
    [string]$RunCase,
    [string]$MakePath = "C:\msys64\mingw64\bin\mingw32-make.exe",
    [ValidateSet("Debug")]
    [string]$Configuration = "Debug"
)

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$projectRoot = Resolve-Path (Join-Path $scriptDir "..")
Set-Location $projectRoot

if ($RunCase) {
    $CaseConfig = "INPUT/Cases/$RunCase/case_config.txt"
} elseif (!$CaseConfig) {
    $CaseConfig = "INPUT/Cases/SaigonDongNai/case_config.txt"
}

if (!(Test-Path $MakePath)) {
    throw "mingw32-make.exe not found at '$MakePath'. Update -MakePath to match your MSYS2 installation."
}

Write-Host "==> Building CGEM Network ($Configuration)" -ForegroundColor Cyan

$srcDir = Join-Path $projectRoot "src"
$binDir = Join-Path $projectRoot "bin/$Configuration"
$exePath = Join-Path $binDir "CGEM_Network.exe"

# Create bin dir
if (!(Test-Path $binDir)) {
    New-Item -ItemType Directory -Path $binDir | Out-Null
}

# Compile all .c files
$srcFiles = Get-ChildItem -Path $srcDir -Filter "*.c" | Select-Object -ExpandProperty FullName
$objs = @()
foreach ($src in $srcFiles) {
    $obj = [IO.Path]::ChangeExtension($src, ".o").Replace($srcDir, $binDir)
    $objs += $obj
    # Suppress format-truncation warnings that are false positives for our use case
    $compileOutput = & $MakePath.Replace("mingw32-make.exe", "gcc.exe") -std=c11 -O2 -Wall -Wextra -Wno-format-truncation -c $src -o $obj 2>&1
    $compileOutput | Out-Host
    if ($LASTEXITCODE -ne 0) {
        throw "Compilation failed for $src"
    }
}

# Link
$linkOutput = & $MakePath.Replace("mingw32-make.exe", "gcc.exe") -o $exePath $objs 2>&1
$linkOutput | Out-Host
if ($LASTEXITCODE -ne 0) {
    throw "Linking failed."
}

Write-Host "==> Running $exePath $CaseConfig" -ForegroundColor Cyan
& $exePath $CaseConfig
