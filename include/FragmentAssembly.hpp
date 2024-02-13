#ifndef FRAG_ASSEM_HPP
#define FRAG_ASSEM_HPP

#include "TAPERSS.hpp"
#include "RNAData.hpp"
#include "RNADataArrayInternalLoop.hpp"
#include "RNA_Math.hpp"
#include "Kabsch.hpp"
#include "steric_clash_check.hpp"
#include "CMB_Manager.hpp"

AttachStatus attach_5p_to_3p(RNAData *reference, RNAData *rotated, RunInfo& run_info);
template <uint32_t OPTS> AttachStatus fragment_assembly(DimerLibArray &lib, DimerLibArray &wc_lib, RNADataArray &assembled, CMB_Manager &manager, RunInfo& run_info);
template <uint32_t OPTS> AttachStatus prepare_right(RNAData* to_be_assembled, DimerLibArray &WC_Lib, RNADataArray &assembled, RunInfo& run_info);
#endif
