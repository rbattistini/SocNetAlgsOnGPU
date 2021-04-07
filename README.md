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
# cmake ../ -DCMAKE_BUILD_TYPE=Default
```

5. compile with `make`.

## CUDA best practices and optimizations

- BSR instead of CSR with cuSPARSE
- CUDA Streams
- CUDA Graphs
- `cudaMemcpyAsync()`

Two samples and blog posts on NVIDIA developers.

- cuSPARSE for handling sparse matrices (CSR or COO).
    - https://docs.nvidia.com/cuda/cusparse/index.html#csc-format
    - https://github.com/NVIDIA/CUDALibrarySamples/tree/master/cuSPARSE

- cub
    - https://github.com/NVIDIA/cub

- Thrust
    - https://docs.nvidia.com/cuda/thrust/index.html

- cuGRAPH Rapids
    - https://docs.rapids.ai/api

## Work in progress

device_management -> 
error_handling -> cudaSafeCall(), cudaCheckError()
stream_management -> ...

# Why work with warps?

Using a “per warp” approach has two practical benefits. First, thread in a warp
 can communicate efficiently using warp vote and __shfl (shuffle) operations.
 Second, there is no need for synchronization because threads in a warp work
 synchronously.

Applying the concepts shown later at the thread block level is significantly
slower because of the slower communication through shared memory and the
necessary synchronization.

## Async mem 101

https://stackoverflow.com/questions/13743039/about-cudamemcpyasync-function/13746353#13746353

## Other useful articles and Q&A

- https://developer.nvidia.com/blog/cuda-pro-tip-increase-performance-with-vectorized-memory-access/
- https://developer.nvidia.com/blog/gpu-pro-tip-cuda-7-streams-simplify-concurrency/
- https://developer.nvidia.com/blog/cuda-pro-tip-optimized-filtering-warp-aggregated-atomics/
- https://developer.nvidia.com/blog/voting-and-shuffling-optimize-atomic-operations/
- https://developer.nvidia.com/blog/cuda-pro-tip-write-flexible-kernels-grid-stride-loops/
- https://forums.developer.nvidia.com/t/warp-reduction-in-kernel-with-if-guard/164699/5
- https://developer.nvidia.com/blog/using-cuda-warp-level-primitives/

- https://stackoverflow.com/questions/21986542/is-cudamallocmanaged-slower-than-cudamalloc/21990899#21990899
- https://stackoverflow.com/questions/50639194/shfl-down-and-shfl-down-sync-give-different-results/50647054#50647054
- https://stackoverflow.com/questions/12936986/why-does-cudamalloc-use-pointer-to-pointer/12937162#12937162
- https://stackoverflow.com/questions/17599189/what-is-the-purpose-of-using-multiple-arch-flags-in-nvidias-nvcc-compiler/17599585#17599585
- https://stackoverflow.com/questions/35656294/cuda-how-to-use-arch-and-code-and-sm-vs-compute/35657430#35657430
- https://stackoverflow.com/questions/18020647/cuda-constant-memory-best-practices/18021374#18021374
- https://stackoverflow.com/questions/16483685/get-gpu-memory-usage-programmatically/16483933#16483933
- https://stackoverflow.com/questions/40589814/cuda-runtime-version-vs-cuda-driver-version-whats-the-difference/40590196#40590196
