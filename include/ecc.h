/****************************************************************************
 * @file ecc.h
 * @author Riccardo Battistini <riccardo.battistini2(at)studio.unibo.it>
 *
 * @brief Functions to compute the diameter of an undirected graph, the
 * eccentricity of its vertices and to get the graph density.
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
 ****************************************************************************/

#pragma once
#ifndef ECC_H
#define ECC_H

#include "graphs.h"
#include "matds.h"
#include <climits>

/**
 * @brief Compute the diameter of the given undirected graph.
 *
 * @note The diameter is computed without using sampling techniques. It is
 * obtained by performing breadth-first searches starting from each
 * vertex and by iterating over the eccentricities' array to find the maximum.
 *
 * @param g input graph in CSR format stored as sparse pattern matrix
 * @return 0 if successful, 1 otherwise
 */
int get_diameter(matrix_pcsr_t *g);

/**
 * @brief Get the eccentricity of each vertex in the given graph.
 *
 * @note Eccentricity of a vertex is the maximum distance from one vertex to
 * all other vertices in the graph.
 *
 * @param g input graph in CSR format stored as sparse pattern matrix
 * @param eccentricity array that stores eccentricity of each vertex
 * @return 0 if successful, 1 otherwise
 */
int get_vertices_eccentricity(matrix_pcsr_t *g, int *eccentricity);

/**
 * @brief Get graph density.
 *
 * @note Graph density indicates how many of all possible connections
 * actually exists.
 *
 * @param nvertices
 * @param nedges
 * @return the density of the graph.
 */
double get_density(double nvertices, double nedges);

#endif//ECC_H
