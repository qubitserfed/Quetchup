#!/usr/bin/env python3

import subprocess
import os
import sys

def compile_testing():
    """Compile testing.cpp with all dependencies and debugging symbols"""
    
    # Directories
    src_dir = "src"
    build_dir = "build"
    
    # Create build directory if it doesn't exist
    os.makedirs(build_dir, exist_ok=True)
    
    # Get all source files except interface.cpp, testing.cpp, and python_module.cpp
    source_files = []
    for file in os.listdir(src_dir):
        if file.endswith('.cpp') and file not in ['interface.cpp', 'testing.cpp', 'python_module.cpp']:
            source_files.append(os.path.join(src_dir, file))
    
    # Add testing.cpp
    source_files.append(os.path.join(src_dir, 'testing.cpp'))
    
    # Compiler and flags
    compiler = "g++"
    flags = [
        "-std=c++20",
        "-g",           # Debug symbols
        "-O0",          # No optimization for debugging
        "-pthread",
        "-Wall",        # Enable warnings
        "-Wextra"       # Extra warnings
    ]
    
    # Output executable
    output = os.path.join(build_dir, "testing")
    
    # Build command
    cmd = [compiler] + flags + source_files + ["-o", output]
    
    print("Compiling testing executable with debug symbols...")
    print("Command:", " ".join(cmd))
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("Compilation successful!")
        print(f"Executable created: {output}")
        if result.stdout:
            print("Stdout:", result.stdout)
    except subprocess.CalledProcessError as e:
        print("Compilation failed!")
        print("Return code:", e.returncode)
        if e.stdout:
            print("Stdout:", e.stdout)
        if e.stderr:
            print("Stderr:", e.stderr)
        sys.exit(1)

if __name__ == "__main__":
    compile_testing()
