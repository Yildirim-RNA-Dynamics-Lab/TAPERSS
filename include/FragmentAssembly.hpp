#ifndef FRAG_ASSEM_HPP
#define FRAG_ASSEM_HPP

#include "RNACMB.hpp"
#include "RNAData.hpp"
#include "RNADataArrayInternalLoop.hpp"
#include "RNA_Math.hpp"
#include "Kabsch.hpp"
#include "steric_clash_check.hpp"
#include "CMB_Manager.hpp"

attach_status rotate(RNAData *reference, RNAData *rotated);
attach_status fragment_assembly(DimerLibArray &Lib, RNADataArray &assembled, CMB_Manager &manager);
attach_status prepare_right(RNAData* to_be_assembled, DimerLibArray &WC_Lib, RNADataArrayInternalLoop &assembled);
attach_status fragment_assembly_IL(DimerLibArray &Lib, DimerLibArray &WC_Lib, RNADataArrayInternalLoop &assembled, CMB_Manager &manager);


#endif
