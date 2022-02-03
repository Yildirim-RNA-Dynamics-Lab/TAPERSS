#ifndef HAIRPIN_HPP
#define HAIRPIN_HPP

#include "RNA_data.hpp"
#include "DimerLib.hpp"
#include "Kabsch.hpp"
#include "RNACMB.hpp"

gsl_matrix *make_WC_submatrix(RNA_data *A, RNA_data *B);
bool is_hairpin(RNA_data_array& sequence, DimerLibArray& WC_Lib);

#endif