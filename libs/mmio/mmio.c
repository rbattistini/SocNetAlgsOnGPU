/* 
 * Matrix Market I/O library for ANSI C
 * See http://math.nist.gov/MatrixMarket for details.
 *
 * Last modified by Riccardo Battistini 09/05/21
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "mmio.h"

int mm_read_banner(FILE *f, MM_typecode *matcode) {

    char line[MM_MAX_LINE_LENGTH];
    char banner[MM_MAX_TOKEN_LENGTH];
    char mtx[MM_MAX_TOKEN_LENGTH];
    char crd[MM_MAX_TOKEN_LENGTH];
    char data_type[MM_MAX_TOKEN_LENGTH];
    char storage_scheme[MM_MAX_TOKEN_LENGTH];
    char *p;

    mm_clear_typecode(matcode);

    if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
        return MM_PREMATURE_EOF;

    if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type,
               storage_scheme) != 5)
        return MM_PREMATURE_EOF;

    for (p = mtx;
         *p != '\0'; *p = tolower(*p), p++);  /* convert to lower case */
    for (p = crd; *p != '\0'; *p = tolower(*p), p++);
    for (p = data_type; *p != '\0'; *p = tolower(*p), p++);
    for (p = storage_scheme; *p != '\0'; *p = tolower(*p), p++);

    /* check for banner */
    if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
        return MM_NO_HEADER;

    /* first field should be "mtx" */
    if (strcmp(mtx, MM_MTX_STR) != 0)
        return MM_UNSUPPORTED_TYPE;
    mm_set_matrix(matcode);


    /* second field describes whether this is a sparse matrix (in coordinate
            storage) or a dense array */
    if (strcmp(crd, MM_SPARSE_STR) == 0)
        mm_set_sparse(matcode);
    else if (strcmp(crd, MM_DENSE_STR) == 0)
        mm_set_dense(matcode);
    else
        return MM_UNSUPPORTED_TYPE;


    /* third field */
    if (strcmp(data_type, MM_REAL_STR) == 0)
        mm_set_real(matcode);
    else if (strcmp(data_type, MM_COMPLEX_STR) == 0)
        mm_set_complex(matcode);
    else if (strcmp(data_type, MM_PATTERN_STR) == 0)
        mm_set_pattern(matcode);
    else if (strcmp(data_type, MM_INT_STR) == 0)
        mm_set_integer(matcode);
    else
        return MM_UNSUPPORTED_TYPE;


    /* fourth field */
    if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
        mm_set_general(matcode);
    else if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
        mm_set_symmetric(matcode);
    else if (strcmp(storage_scheme, MM_HERM_STR) == 0)
        mm_set_hermitian(matcode);
    else if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
        mm_set_skew(matcode);
    else
        return MM_UNSUPPORTED_TYPE;


    return 0;
}

char *mm_typecode_to_str(const MM_typecode matcode) {
    char buffer[MM_MAX_LINE_LENGTH];
    const char *types[4];

    /* check for MTX type */
    if (mm_is_matrix(matcode))
        types[0] = MM_MTX_STR;

    /* check for CRD or ARR matrix */
    if (mm_is_sparse(matcode))
        types[1] = MM_SPARSE_STR;
    else if (mm_is_dense(matcode))
        types[1] = MM_DENSE_STR;
    else
        return NULL;

    /* check for element data type */
    if (mm_is_real(matcode))
        types[2] = MM_REAL_STR;
    else if (mm_is_complex(matcode))
        types[2] = MM_COMPLEX_STR;
    else if (mm_is_pattern(matcode))
        types[2] = MM_PATTERN_STR;
    else if (mm_is_integer(matcode))
        types[2] = MM_INT_STR;
    else
        return NULL;


    /* check for symmetry type */
    if (mm_is_general(matcode))
        types[3] = MM_GENERAL_STR;
    else if (mm_is_symmetric(matcode))
        types[3] = MM_SYMM_STR;
    else if (mm_is_hermitian(matcode))
        types[3] = MM_HERM_STR;
    else if (mm_is_skew(matcode))
        types[3] = MM_SKEW_STR;
    else
        return NULL;

    sprintf(buffer, "%s %s %s %s", types[0], types[1], types[2], types[3]);
    return strdup(buffer);

}

int mm_write_banner(FILE *f, MM_typecode matcode) {

    char *str = mm_typecode_to_str(matcode);
    int ret_code;

    ret_code = fprintf(f, "%s %s\n", MatrixMarketBanner, str);
    free(str);
    if (ret_code < 0 )
        return MM_COULD_NOT_WRITE_FILE;
    else
        return 0;
}

int mm_write_mtx_crd_size(FILE *f, int M, int N, int nz) {

    if (fprintf(f, "%d %d %d\n", M, N, nz) < 0)
        return MM_COULD_NOT_WRITE_FILE;
    else
        return 0;
}

int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz ) {

    char line[MM_MAX_LINE_LENGTH];
    int num_items_read;

    /* set return null parameter values, in case we exit with errors */
    *M = *N = *nz = 0;

    /* now continue scanning until you reach the end-of-comments */
    do
    {
        if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL)
            return MM_PREMATURE_EOF;
    }while (line[0] == '%');

    /* line[] is either blank or has M,N, nz */
    if (sscanf(line, "%d %d %d", M, N, nz) == 3)
        return 0;

    else
        do
        {
            num_items_read = fscanf(f, "%d %d %d", M, N, nz);
            if (num_items_read == EOF) return MM_PREMATURE_EOF;
        }
        while (num_items_read != 3);

    return 0;
}
