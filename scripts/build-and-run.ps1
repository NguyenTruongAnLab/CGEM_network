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
$nloptDir = Join-Path $projectRoot "external/nlopt"
$nloptLib = Join-Path $nloptDir "build"
$nloptInclude = Join-Path $nloptDir "src/api"

# Create bin dir
if (!(Test-Path $binDir)) {
    New-Item -ItemType Directory -Path $binDir | Out-Null
}

# Compile all .c files RECURSIVELY (includes physics/, rive/, io/, optimization/)
$srcFiles = Get-ChildItem -Path $srcDir -Filter "*.c" -Recurse | Select-Object -ExpandProperty FullName
$objs = @()
$gcc = $MakePath.Replace("mingw32-make.exe", "gcc.exe")

Write-Host "Compiling $($srcFiles.Count) source files..." -ForegroundColor Yellow

foreach ($src in $srcFiles) {
    # Create output path preserving directory structure
    $relativePath = $src.Substring($srcDir.Length + 1)
    $objName = [IO.Path]::ChangeExtension($relativePath, ".o").Replace("\", "_").Replace("/", "_")
    $obj = Join-Path $binDir $objName
    $objs += $obj
    
    # Compile with include paths for all subdirectories
    $includes = "-I`"$srcDir`" -I`"$srcDir\rive`" -I`"$srcDir\physics`" -I`"$srcDir\io`" -I`"$nloptInclude`""
    $compileCmd = "$gcc -std=c11 -O2 -Wall -Wextra -Wno-format-truncation $includes -c `"$src`" -o `"$obj`""
    $compileOutput = Invoke-Expression $compileCmd 2>&1
    $compileOutput | Out-Host
    if ($LASTEXITCODE -ne 0) {
        throw "Compilation failed for $src"
    }
}

# Link with NLopt library
Write-Host "Linking with NLopt library..." -ForegroundColor Yellow
$objList = $objs -join '" "'
$linkCmd = "$gcc -o `"$exePath`" `"$objList`" -L`"$nloptLib`" -lnlopt -lm"
$linkOutput = Invoke-Expression $linkCmd 2>&1
$linkOutput | Out-Host
if ($LASTEXITCODE -ne 0) {
    throw "Linking failed."
}

# Copy NLopt DLL if not present
$dllSrc = Join-Path $nloptLib "libnlopt.dll"
$dllDst = Join-Path $binDir "libnlopt.dll"
if ((Test-Path $dllSrc) -and !(Test-Path $dllDst)) {
    Write-Host "Copying libnlopt.dll to bin directory..." -ForegroundColor Yellow
    Copy-Item $dllSrc $dllDst
}

Write-Host "==> Running $exePath $CaseConfig" -ForegroundColor Cyan
& $exePath $CaseConfig
