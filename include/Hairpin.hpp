#ifndef HAIRPIN_HPP
#define HAIRPIN_HPP

#include "RNAData.hpp"
#include "DimerLib.hpp"
#include "Kabsch.hpp"
#include "RNACMB.hpp"

gsl_matrix *make_WC_submatrix_gsl(gsl_matrix *A, gsl_matrix *B);
gsl_matrix *make_WC_submatrix(RNAData *A, RNAData *B);
bool is_WC_pair(RNADataArray& sequence, DimerLibArray& WC_Lib, int i, int j, int WC_pair_idx);

#endif