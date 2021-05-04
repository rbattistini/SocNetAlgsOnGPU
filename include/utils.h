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

#ifndef SOCNETALGSONGPU_UTILS_H
#define SOCNETALGSONGPU_UTILS_H

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cerrno>
#include <stdio_ext.h>

#define array_length(x) (sizeof(x) / sizeof((x)[0]))

void fill( int *arr, int n, int v);

int close_stream(FILE *stream);

void print_array(const int *arr, int n);

void print_edge_list(const int *row_offsets, const int *cols, int nrows);

void bucket_sort(const int *arr, int n);

#endif //SOCNETALGSONGPU_UTILS_H
