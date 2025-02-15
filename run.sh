#!/bin/bash

# Compile the C++ files
g++ -std=c++17 -o main main.cpp Process.cpp

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful. Running the program..."
    ./main
else
    echo "Compilation failed."
fi