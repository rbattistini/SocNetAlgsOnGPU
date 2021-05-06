#!/bin/sh

alias clang-format="clang-format-10"

# Format all C++ and CUDA files in src.
clang-format -style=file ../src/*.cpp -i
clang-format -style=file ../src/*.cu -i

# Format all header files in include.
clang-format -style=file ../include/*.h -i
clang-format -style=file ../include/*.cuh -i
