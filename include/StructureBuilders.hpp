#ifndef STRUCTUREBUILDERS_HPP
#define STRUCTUREBUILDERS_HPP

#include "DimerLib.hpp"
#include "RNA_data.hpp"
#include "output_string.hpp"
#include "Attaching_Functions.hpp"
#include "steric_clash_check.hpp"
#include "Hairpin.hpp"

void create_custom_structure(DimerLibArray& Lib, DimerLibArray& WC_Lib,RNA_data_array& assembled, output_string& o_string, int *indices);
void create_custom_structure_list(DimerLibArray& Lib, DimerLibArray& WC_Lib,RNA_data_array& assembled, output_string& o_string, int num_strs);
void create_custom_structure_list_testing(DimerLibArray &Lib, RNA_data_array &assembled, output_string &o_string, int num_strs);

#endif