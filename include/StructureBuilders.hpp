#ifndef STRUCTUREBUILDERS_HPP
#define STRUCTUREBUILDERS_HPP

#include "DimerLib.hpp"
#include "RNAData.hpp"
#include "output_string.hpp"
#include "Attaching_Functions.hpp"
#include "steric_clash_check.hpp"
#include "Hairpin.hpp"
#include "RNADataArrayInternalLoop.hpp"

void create_custom_structure(DimerLibArray& Lib, DimerLibArray& WC_Lib,RNADataArray& assembled, output_string& o_string, int *indices);
void create_custom_structure_IL(DimerLibArray& Lib, DimerLibArray& WC_Lib,RNADataArrayInternalLoop& assembled, output_string& o_string, int *indices);
void create_custom_structure_list(DimerLibArray& Lib, DimerLibArray& WC_Lib,RNADataArray& assembled, output_string& o_string, int num_strs);
void create_custom_structure_list_testing(DimerLibArray &Lib, RNADataArray &assembled, output_string &o_string, int num_strs);

#endif