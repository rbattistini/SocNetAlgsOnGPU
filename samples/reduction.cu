//
// from
// - https://stackoverflow.com/questions/50639194/shfl-down-and-shfl-down-sync-give-different-results/50647054#50647054
// - https://developer.nvidia.com/blog/cooperative-groups/
//
#include <cstdio>
#include <cooperative_groups.h>

namespace cg = cooperative_groups;

__global__ void warp_reduction()
{
    double  temp = 2.0;
    __syncthreads();

    for (int offset = 32/2; offset > 0; offset /= 2){
        temp+=__shfl_down_sync(0xFFFFFFFF, temp, offset,32);
    }
    printf("%d %f %d \n",threadIdx.x ,temp,blockDim.x * gridDim.x);
}

__device__ int reduce_sum(cg::thread_group g, int *temp, int val)
{
    int lane = (int) g.thread_rank();

    // Each iteration halves the number of active threads
    // Each thread adds its partial sum[i] to sum[lane+i]
    for (int i = (int) g.size() / 2; i > 0; i /= 2)
    {
        temp[lane] = val;
        g.sync(); // wait for all threads to store
        if(lane<i) val += temp[lane + i];
        g.sync(); // wait for all threads to load
    }
    return val; // note: only thread 0 will return full sum
}

template <int tile_sz>
__device__ int reduce_sum_tile_shfl(cg::thread_block_tile<tile_sz> g, int val)
{
    // Each iteration halves the number of active threads
    // Each thread adds its partial sum[i] to sum[lane+i]
    for (int i = g.size() / 2; i > 0; i /= 2) {
        val += g.shfl_down(val, i);
    }

    return val; // note: only thread 0 will return full sum
}

template <int tile_sz>
__global__ void sum_kernel_tile_shfl(int *sum, const int *input, int n)
{
    int my_sum = thread_sum(input, n);

    auto tile = cg::tiled_partition<tile_sz>(cg::this_thread_block());
    int tile_sum = reduce_sum_tile_shfl(tile, my_sum);

    if (tile.thread_rank() == 0) atomicAdd(sum, tile_sum);
}

int main(){
    warp_reduction<<<1, 64>>>();
    cudaDeviceSynchronize();
}
