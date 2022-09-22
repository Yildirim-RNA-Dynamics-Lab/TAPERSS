#ifndef KABSCH_HPP
#define KABSCH_HPP

#include "RNACMB.hpp"
#include "RNA_Math.hpp"

enum kabsch_matrix {KABSCH_MATRIX_P, KABSCH_MATRIX_Q};

double get_determinant(gsl_matrix *A, bool inPlace);
gsl_matrix *kabsch_get_rotation_matrix_generic(gsl_matrix *P, gsl_matrix *Q, double * __restrict__ COMP, double * __restrict__ COMQ);
void kabsch_create(size_t M, size_t N);
void kabsch_destroy();
template <kabsch_matrix M_REQUEST> gsl_matrix* kabsch_prepare_matrix(size_t M, size_t N, uint16_t *Rows, gsl_matrix *source);
gsl_matrix* kabsch_get_work_matrix(size_t M, size_t N);
gsl_matrix* kabsch_allocate_work_matrix(gsl_matrix *P);
double get_determinant_3x3fast(gsl_matrix *A);
void kabsch_calculate_rotation_matrix_Nx3fast(gsl_matrix *P, gsl_matrix *Q, gsl_matrix *P_WORK, double *__restrict__ COMP, double *__restrict__ COMQ);
/**
 * Gets pointer to gsl_matrix pointer to rotation matrix. Rotation matrix memory is entirely managed by Kabsch.
 * Updated after running: kabsch_calculate_rotation_matrix_Nx3fast(...)
 **/
gsl_matrix* kabsch_get_rotation_matrix();

#endif 