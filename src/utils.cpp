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
 ****************************************************************************/

#include "utils.h"

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

/*
 * Properly close a stream with error checking.
 *
 * Taken from:
 * https://stackoverflow.com/questions/4972994/how-to-close-stdout-and-stderr-in-c
 */
int close_stream(FILE *stream) {

    const bool some_pending = (__fpending (stream) != 0);
    const bool prev_fail = (ferror (stream) != 0);
    const bool fclose_fail = (fclose (stream) != 0);

    /* Return an error indication if there was a previous failure or if
       fclose failed, with one exception: ignore an fclose failure if
       there was no previous error, no data remains to be flushed, and
       fclose failed with EBADF.  That can happen when a program like cp
       is invoked like this 'cp a b >&-' (i.e., with standard output
       closed) and doesn't generate any output (hence no previous error
       and nothing to be flushed).  */

    if (prev_fail || (fclose_fail && (some_pending || errno != EBADF))) {

        if (! fclose_fail)
            errno = 0;
        return EOF;
    }

    return 0;
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
        int end = row_offsets[i + 1];

        for(int j = begin; j < end; j++) {

            if(j == begin)
                printf("%d | %d", i, cols[j]);
            else
                printf(", %d", cols[j]);
        }

        /*
         * For isolated vertices.
         */
        if(begin == end) {
            printf("%d | ", i);
        }

        printf("\n");
    }
}
