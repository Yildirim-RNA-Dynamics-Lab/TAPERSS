#ifndef RNA_MATH_HPP
#define RNA_MATH_HPP

#include "RNACMB.hpp"
#define print_gsl_matrix(M) gsl_matrix_print(M, #M);

double rmsd_generic(gsl_matrix *A, gsl_matrix *B);
void get_matrix_COM(gsl_matrix *M, double *COM);
void center_matrix(gsl_matrix *M, double *COM);
void translate_matrix(double *__restrict__ tV, gsl_matrix *M, double scalar);
void apply_rotation_matrix(gsl_matrix *R, gsl_matrix *M);
void gsl_matrix_print(gsl_matrix *M, const char* name);
void gsl_vector_print(gsl_vector *V); 

/* Inline Functions */

inline double distance(gsl_vector *A, gsl_vector *B)
{
    double x, y, z;
    x = gsl_vector_get(B, 0) - gsl_vector_get(A, 0);
    y = gsl_vector_get(B, 1) - gsl_vector_get(A, 1);
    z = gsl_vector_get(B, 2) - gsl_vector_get(A, 2);
    return sqrt(square(x) + square(y) + square(z));
}

#endif