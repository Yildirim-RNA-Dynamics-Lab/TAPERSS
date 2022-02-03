#ifndef KABSCH_HPP
#define KABSCH_HPP

#include "RNACMB.hpp"
#include "RNA_Math.hpp"

double get_determinant(gsl_matrix *A, bool inPlace);
gsl_matrix *kabsch_get_rotation_matrix_generic(gsl_matrix *P, gsl_matrix *Q, double * __restrict__ COMP, double * __restrict__ COMQ);

#endif 