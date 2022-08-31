#ifndef STERIC_CLASH_CHECK_HPP
#define STERIC_CLASH_CHECK_HPP

#include "Attaching_Functions.hpp"
#include "RNACMB.hpp"
#include "RNAData.hpp"
#include "RNA_Math.hpp"

bool steric_clash_check(RNADataArray &sequence, RNAData *attach);
void SCC_record_COM_distance(RNADataArray &sequence, RNAData *attach);
attach_status steric_clash_check_COM_tester(RNADataArray &sequence, RNAData *attach);
attach_status steric_clash_check_COM(RNADataArray &sequence, RNAData *attach);
void gsl_vector_print(gsl_vector* V);

#endif