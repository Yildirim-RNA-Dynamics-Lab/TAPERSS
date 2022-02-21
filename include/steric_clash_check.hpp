#ifndef STERIC_CLASH_CHECK_HPP
#define STERIC_CLASH_CHECK_HPP

#include "RNACMB.hpp"
#include "RNA_data.hpp"

bool steric_clash_check(RNA_data_array &sequence, RNA_data *attach);
void SCC_record_COM_distance(RNA_data_array &sequence, RNA_data *attach);

#endif