/****************************************************************************
 * @file utils.h
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * Basic functions for printing some types of arrays and matrices related to
 * the internal memory storage of graphs and other commonly used utilities.
 *
 * Copyright 2021 (c) 2021 by Riccardo Battistini
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

#pragma once
#ifndef SOCNETALGSONGPU_UTILS_H
#define SOCNETALGSONGPU_UTILS_H

#include <cassert>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <stdio_ext.h>
#include <vector>

inline int max(int a, int b) {
    return a > b ? a : b;
}

inline int min(int a, int b) {
    return a < b ? a : b;
}

int *stlvector_to_array_int(const std::vector<int> &v, int n);

/**
 * Initialize an array.
 *
 * @param arr pointer to the array to be initialized
 * @param n length of the array to be initialized
 * @param v value used to initialize the array
 */
void fill(int *arr, int n, int v);

void fill(bool *arr, int n, int v);

int get_max_idx(const int *arr, int n);

/**
 * Properly close a stream with error checking.
 *
 * Taken from: https://stackoverflow.com/questions/4972994/how-to-close-stdout-and-stderr-in-c
 */
int close_stream(FILE *stream);

void print_array(const int *arr, int n);

void print_array(const float *arr, int n);

void print_edge_list(const int *row_offsets, const int *cols, int nrows);

#endif//SOCNETALGSONGPU_UTILS_H
