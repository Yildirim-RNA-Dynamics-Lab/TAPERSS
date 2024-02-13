#ifndef STERIC_CLASH_CHECK_HPP
#define STERIC_CLASH_CHECK_HPP

#include "TAPERSS.hpp"
#include "RNA_Math.hpp"
#include "RNAData.hpp"
#include "RNADataArrayInternalLoop.hpp"

AttachStatus cg_scc(RNADataArray& Sequence, RNAData* Attach);
template <bool CountGreater>AttachStatus cg_scc_ds(RNADataArray& Sequence, RNAData* Attach);
AttachStatus cg_scc_ds_check_new_2nd(RNADataArray& Sequence, RNAData* Attach);

#endif
