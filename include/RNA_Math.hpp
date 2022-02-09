#ifndef RNA_MATH_HPP
#define RNA_MATH_HPP

#include "RNACMB.hpp"

double rmsd_generic(gsl_matrix *A, gsl_matrix *B);
void get_matrix_COM(gsl_matrix *M, double *COM);
void center_matrix(gsl_matrix *M, double *COM);
void translate_matrix(double *__restrict__ tV, gsl_matrix *M, double scalar);
void apply_rotation_matrix(gsl_matrix *R, gsl_matrix *M);
double distance(gsl_vector *A, gsl_vector *B);
void gsl_matrix_print(gsl_matrix *M); 

#endif