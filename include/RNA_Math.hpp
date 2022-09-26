#ifndef RNA_MATH_HPP
#define RNA_MATH_HPP

#include "RNACMB.hpp"
#include "Kabsch.hpp"
#define print_gsl_matrix(M) gsl_matrix_print(M, #M);
#define print_gsl_matrix_row(M, i) gsl_matrix_print_row(M, i);

double rmsd_generic(gsl_matrix *A, gsl_matrix *B);
void get_matrix_COM(gsl_matrix *M, double *COM);
void center_matrix(gsl_matrix *M, double *COM);
void translate_matrix(double *__restrict__ tV, gsl_matrix *M, double scalar);
void apply_rotation_matrix(gsl_matrix *R, gsl_matrix *M);
void gsl_matrix_print(gsl_matrix *M, const char* name);
void gsl_matrix_print_row(gsl_matrix *M, size_t i);
void gsl_vector_print(gsl_vector *V); 
int array_min_idx_for_energy(int *Arr, int N);
/* Inline Functions */

inline double distance_vec2vec(gsl_vector *A, gsl_vector *B)
{
    double x, y, z;
    x = gsl_vector_get(B, 0) - gsl_vector_get(A, 0);
    y = gsl_vector_get(B, 1) - gsl_vector_get(A, 1);
    z = gsl_vector_get(B, 2) - gsl_vector_get(A, 2);
    return sqrt(square(x) + square(y) + square(z));
}

inline double distance_mat2vec(gsl_matrix *A, size_t Row, gsl_vector *B)
{
    double x, y, z;
    x = gsl_vector_get(B, 0) - gsl_matrix_get(A, Row, 0);
    y = gsl_vector_get(B, 1) - gsl_matrix_get(A, Row, 1);
    z = gsl_vector_get(B, 2) - gsl_matrix_get(A, Row, 2);
    return sqrt(square(x) + square(y) + square(z));
}

inline double distance_mat2mat(gsl_matrix *A, size_t Row_A, gsl_matrix *B, size_t Row_B)
{
    double x, y, z;
    x = gsl_matrix_get(B, Row_B, 0) - gsl_matrix_get(A, Row_A, 0);
    y = gsl_matrix_get(B, Row_B, 1) - gsl_matrix_get(A, Row_A, 1);
    z = gsl_matrix_get(B, Row_B, 2) - gsl_matrix_get(A, Row_A, 2);
    return sqrt(square(x) + square(y) + square(z));
}

inline void gsl_matrix_row_copy(gsl_matrix *dest, size_t row_d, gsl_matrix *src, size_t row_s)
{
    memmove(&dest->data[IDX_FLAT2D(row_d, 0, MATRIX_DIMENSION2)], &src->data[IDX_FLAT2D(row_s, 0, MATRIX_DIMENSION2)], MATRIX_DIMENSION2 * sizeof(double));
}

#endif