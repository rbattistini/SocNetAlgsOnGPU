/****************************************************************************
 *
 * utils.h - basic serial prefix Sum and reduce operators and functions for
 * printing some types of arrays related to internal memory storage of graphs
 *
 * Copyright 2021 (c) 2021 by Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ---------------------------------------------------------------------------
 *
 * TODO implement check_bc()
 *
 ****************************************************************************/

#ifndef SOCNETALGSONGPU_UTILS_H
#define SOCNETALGSONGPU_UTILS_H

/*
 * from http://www.graphics.stanford.edu/~seander/bithacks.html#DetermineIfPowerOf2
 */
int is_power_of2(int n) {
    return (n & (n - 1)) == 0;
}

/**
 * Initialize an array.
 *
 * @param arr pointer to the array to be initialized
 * @param n length of the array to be initialized
 * @param v value used to initialize the array
 */
void fill( int *arr, int n, int v) {
    for (int i = 0; i < n; i++)
        arr[i] = v;
}

/**
 * Compute the prefix sum of an integers array arr in a given array prefixSum.
 *
 * @param arr array of which to compute the prefix sum
 * @param n length of arr
 * @param prefixSum the array where the prefix sum is put. It must be allocated
 * by the caller
 */
void fill_prefix_sum(const int *arr, int n, int *prefixSum) {
    prefixSum[0] = arr[0];
    for (int i = 1; i < n; i++)
        prefixSum[i] = prefixSum[i - 1] + arr[i];
}

/**
 * Compute the partial sum of an integers array arr in a given interval
 * of indexes.
 *
 * @param arr the array of which the sum reduction will be computed
 * @param n number of elements in arr
 * @param start first index of the interval given
 * @param end last index of the interval given
 * @return the result of the sum reduction or -1 if the given interval is
 * invalid.
 */
int reduce_sum(const int *arr, int n, int start, int end) {
    if(end >= n || start < 0)
        return -1;

    int tmp = 0;
    for (int i = start; i <= end; i++)
        tmp += arr[i];
    return tmp;
}

void print_array(const int *arr, int n) {
    printf("[ ");

    for(int i = 0; i < n + 1; i++)
        printf("%d ", arr[i]);

    printf("]\n");
}

void print_edge_list(const int *row_offsets, const int *cols, int nrows) {

    printf("Edge lists for each vertex: \n");

    for(int i = 0; i < nrows; i++) {

        int begin = row_offsets[i];
        int end = row_offsets[i+1];

        for(int j = begin; j < end; j++) {

            if(j == begin)
                printf("%d | %d", i, cols[j]);
            else
                printf(", %d", cols[j]);
        }

        if(begin == end) {
            printf("%d | ", i);
        }

        printf("\n");
    }
}

void check_bc(cudaGraph_t g, const float *bc_cpu, const float *bc_gpu) {
    for(int i = 0; i < g; i++) {
        assert();
    }
}

#endif //SOCNETALGSONGPU_UTILS_H
