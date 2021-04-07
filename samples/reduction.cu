//
// from https://stackoverflow.com/questions/50639194/shfl-down-and-shfl-down-sync-give-different-results/50647054#50647054
//
#include <stdio.h>

__global__ void shufledown(double* a, double *b,double *c, int N)
{
    double  temp = 2.0;
    __syncthreads();

    for (int offset = 32/2; offset > 0; offset /= 2){
        temp+=__shfl_down_sync(0xFFFFFFFF, temp, offset,32);
    }
    printf("%d %f %d \n",threadIdx.x ,temp,blockDim.x * gridDim.x);
}

int main(){
    double *a = NULL, *b = NULL, *c = NULL;
    shufledown<<<1,64>>>(a, b, c, 0);
    cudaDeviceSynchronize();
}
