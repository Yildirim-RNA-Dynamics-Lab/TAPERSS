#ifndef KABSCH_HPP
#define KABSCH_HPP

#include "RNACMB.hpp"
#include "RNA_Math.hpp"

double get_determinant(gsl_matrix *A, bool inPlace);
gsl_matrix *kabsch_get_rotation_matrix_generic(gsl_matrix *P, gsl_matrix *Q, double * __restrict__ COMP, double * __restrict__ COMQ);

void kabsch_create(size_t memsize);
void kabsch_destroy();
gsl_matrix* kabsch_allocate_work_matrix(gsl_matrix *P);
double get_determinant_3x3fast(gsl_matrix *A);
void kabsch_calculate_rotation_matrix_Nx3fast(gsl_matrix *P, gsl_matrix *Q, gsl_matrix *P_WORK, double *__restrict__ COMP, double *__restrict__ COMQ);
/**
 * Gets pointer to gsl_matrix pointer to rotation matrix. Rotation matrix memory is entirely managed by Kabsch.
 **/
gsl_matrix* kabsch_get_rotation_matrix();

#endif 