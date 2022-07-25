#ifndef ATTACH_FUNC_HPP
#define ATTACH_FUNC_HPP

#include "RNACMB.hpp"
#include "RNAData.hpp"
#include "RNA_Math.hpp"
#include "Kabsch.hpp"
#include "steric_clash_check.hpp"

attach_status rotate(RNAData *reference, RNAData *rotated);
attach_status check_attachment(RNADataArray& sequence, RNAData *attach);

#endif
