# How to run this project

Only for GNU/Linux and MacOS.

1. Download and unzip the SNAP source code from [here](https://snap.stanford.edu/snap/download.html).
2. `cd snap-core`
3. `make libinstall` to obtain the static library `libsnap.a`.
4. build the project with CMake issuing the following commands:

```shell
mkdir build
# use the following command to build if you are in project root
cmake -Bbuild -DCMAKE_BUILD_TYPE=Default

# or use the following command if you are in build directory
# cmake ../ -DCMAKE_BUILD_TYPE=Release
```

5. compile with `make`.
