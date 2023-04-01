#ifndef WCPAIR_HPP
#define WCPAIR_HPP

#include "RNAData.hpp"
#include "DimerLib.hpp"
#include "Kabsch.hpp"
#include "RNACMB.hpp"

void overwrite_WC_submatrix_gsl(gsl_matrix *A, gsl_matrix *B, gsl_matrix *WC_pair);
//gsl_matrix *make_WC_submatrix(RNAData *A, RNAData *B); // Deprecated
double WC_check_pair(int WC_pair_idx);
void WC_prepare_structure_matrix(int WC_pair_idx, RNAData* A_data, uint16_t idxA, RNAData* B_data, uint16_t idxB);
double RMSD_WC_pair_gsl(gsl_matrix *WC_matrix, gsl_matrix *Sequence_matrix);
void WC_create(DimerLibArray& WC_Library);
void WC_destroy();

#endif