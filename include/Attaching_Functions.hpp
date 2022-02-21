#ifndef ATTACH_FUNC_HPP
#define ATTACH_FUNC_HPP

#include "RNACMB.hpp"
#include "RNA_data.hpp"
#include "RNA_Math.hpp"
#include "Kabsch.hpp"
#include "steric_clash_check.hpp"

enum attach_status{FAILED, ATTACHED, NOT_CHECKED};

attach_status rotate(RNA_data *reference, RNA_data *rotated);
attach_status check_attachment(RNA_data_array& sequence, RNA_data *attach);

#endif