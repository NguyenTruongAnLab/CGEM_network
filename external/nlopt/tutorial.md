# Method 1: Compile using codeblock gcc

NLopt 2.7.1

https://github.com/stevengj/nlopt/releases

```sh
cd external/nlopt
mkdir build
cd build
cmake -G "MinGW Makefiles" -DCMAKE_TOOLCHAIN_FILE="../cmake/toolchain-x86_64-w64-mingw32.cmake" -DBUILD_SHARED_LIBS=ON -DNLOPT_PYTHON=OFF -DNLOPT_OCTAVE=OFF -DNLOPT_MATLAB=OFF -DNLOPT_GUILE=OFF -DNLOPT_SWIG=OFF ..
mingw32-make
```

toolchain-x86_64-w64-mingw32.cmake

```
set(CMAKE_SYSTEM_NAME Windows)

# specify the compiler paths
set(MINGW_PATH "C:/Program Files/CodeBlocks/MinGW")
set(CMAKE_C_COMPILER "${MINGW_PATH}/bin/gcc.exe")
set(CMAKE_CXX_COMPILER "${MINGW_PATH}/bin/g++.exe")

# where is the target environment
set(CMAKE_FIND_ROOT_PATH "${MINGW_PATH}")

# search for programs in the build host directories
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
# for libraries and headers in the target directories
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE ONLY)

# Set resource compiler
set(CMAKE_RC_COMPILER "${MINGW_PATH}/bin/windres.exe")

# Set make program
set(CMAKE_MAKE_PROGRAM "${MINGW_PATH}/bin/mingw32-make.exe")

# These are needed for compiling
set(CMAKE_AR "${MINGW_PATH}/bin/ar.exe")
set(CMAKE_RANLIB "${MINGW_PATH}/bin/ranlib.exe")
```

# Method 2: Using Pacman
To install the `mingw-w64-x86_64-nlopt` package using `pacman` on MSYS2, follow these steps:

### Step-by-Step Guide to Install `mingw-w64-x86_64-nlopt` with MSYS2

1. **Download and Install MSYS2**:
   - If you don't have MSYS2 installed, download and install it from the [MSYS2 website](https://www.msys2.org/).

2. **Update MSYS2**:
   - Open the MSYS2 MSYS terminal from your Start menu.
   - Update the package database and core system packages:
     ```sh
     pacman -Syu
     ```
   - Close the MSYS2 MSYS terminal and open it again to complete the installation of updated core packages.
   - Run the update command again to ensure that all packages are up to date:
     ```sh
     pacman -Su
     ```

3. **Install `mingw-w64-x86_64-toolchain`**:
   - Open the MSYS2 MinGW 64-bit terminal from your Start menu.
   - Install the base development tools:
     ```sh
     pacman -S --needed base-devel mingw-w64-x86_64-toolchain
     ```

4. **Install `mingw-w64-x86_64-nlopt`**:
   - In the same MSYS2 MinGW 64-bit terminal, install the NLopt package:
     ```sh
     pacman -S mingw-w64-x86_64-nlopt
     ```

5. **Verify Installation**:
   - Verify that the NLopt library and headers are installed in the appropriate directories, typically under `/mingw64/include` for headers and `/mingw64/lib` for libraries.

### Update VSCode Configuration

After installing the `mingw-w64-x86_64-nlopt` package, you need to update your VSCode configuration to use the installed NLopt library.

#### Updated `tasks.json`

```json name=.vscode/tasks.json
{
    "version": "2.0.0",
    "tasks": [
        {
            "label": "Clean NLopt Test",
            "type": "shell",
            "command": "powershell",
            "args": [
                "-Command",
                "if (Test-Path bin/Debug/NLoptTest.exe) { Remove-Item -Path bin/Debug/NLoptTest.exe -Force }",
                "if (Test-Path obj) { Remove-Item -Recurse -Force obj }"
            ],
            "problemMatcher": [],
            "group": {
                "kind": "none",
                "isDefault": true
            }
        },
        {
            "label": "Build NLopt Test",
            "type": "shell",
            "command": "gcc",
            "args": [
                "-g",
                "-w", // Ignore all warnings
                "-o",
                "bin/Debug/NLoptTest.exe",
                "src/nlopt_test.c",
                "-I/mingw64/include", // Include NLopt headers
                "-L/mingw64/lib", // Library path for NLopt
                "-lnlopt", // Link NLopt library
                "-lm"
            ],
            "group": {
                "kind": "build",
                "isDefault": true
            },
            "dependsOn": ["Clean NLopt Test"]
        },
        {
            "label": "Run NLopt Test",
            "type": "shell",
            "command": "bin/Debug/NLoptTest.exe",
            "problemMatcher": [],
            "group": {
                "kind": "test",
                "isDefault": true
            },
            "dependsOn": ["Build NLopt Test"]
        }
    ]
}
```

#### Updated `c_cpp_properties.json`

```json name=.vscode/c_cpp_properties.json
{
  "configurations": [
    {
      "name": "windows-gcc-x64",
      "includePath": [
        "${workspaceFolder}/**",
        "/mingw64/include"
      ],
      "compilerPath": "C:/msys64/mingw64/bin/gcc.exe",
      "cStandard": "${default}",
      "cppStandard": "${default}",
      "intelliSenseMode": "windows-gcc-x64",
      "compilerArgs": [
        ""
      ]
    }
  ],
  "version": 4
}
```

### Steps to Run the Test

1. **Install MSYS2 and Update**:
   - Follow the steps to install MSYS2 and update the package database.

2. **Install NLopt Library**:
   - Use `pacman` to install the `mingw-w64-x86_64-nlopt` package.

3. **Update VSCode Configuration**:
   - Update the `tasks.json` and `c_cpp_properties.json` files to use the installed NLopt library.

4. **Build and Run the Test**:
   - Use VSCode tasks to clean, build, and run the NLopt test.
     - Press `Ctrl+Shift+B` to build the project and select `Build NLopt Test`.
     - Run the `Run NLopt Test` task to execute the test.

By following these steps, you can ensure that the NLopt library is correctly installed and integrated with your VSCode setup using MSYS2 and `pacman`.