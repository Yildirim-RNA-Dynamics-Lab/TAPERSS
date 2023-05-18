#include "StructureBuilders.hpp"

void create_custom_structure(DimerLibArray &Lib, RNADataArray &assembled, output_string &o_string, uint32_t *indices)
{
    attach_status status;
    double RMSD;
    assembled.overwrite(0, indices[0], Lib);
    for (int i = 1; i <= assembled.iterator_max; i++)
    {
        assembled.overwrite(i, indices[i], Lib);
        rotate(assembled.current(), assembled[i]);
        status = steric_clash_check_COMFast(assembled, assembled[i]);
        assembled.keep();
        if (status == ATTACHED)
            printf("Successful attachment on: %d\n", i);
        else
            printf("Unsuccessful attachment on: %d\n", i);
    }
    if (GLOBAL_PERFORM_STRUCTCHECK == HAIRPIN)
    {
        WC_prepare_structure_matrix(0, assembled[0], 0, assembled[assembled.iterator_max], 1);
        RMSD = WC_check_pair(0);
        assembled.update_WC_rmsd(RMSD);
    }
    assembled.update_energy();
    o_string.add_string(assembled.to_string());
}

void create_custom_structure_IL(DimerLibArray &Lib, DimerLibArray &WC_Lib, RNADataArrayInternalLoop &assembled, output_string &o_string, uint32_t *indices)
{
    attach_status status;
    double RMSD;
    assembled.overwrite_initialize(0, indices[0], Lib);
    for (int i = 1; i <= assembled.iterator_max; i++)
    {
        //printf("i = %d, indicies: %d\n", i, indices[i]);
        assembled.overwrite(i, indices[i], Lib);
        if (assembled.should_prepare_right(i) == true)
        {
            //printf("i = %d\n", i);
            prepare_right(assembled[i], WC_Lib, assembled);
            status = steric_clash_check_COMFast_IL_1st_right(assembled, assembled[i]);
            assembled.keep();
        }
        else
        {
            rotate(assembled.current(), assembled[i]);
            if(assembled.count > assembled.WC_size_left)
            {
                status = steric_clash_check_COMFast_IL<true>(assembled, assembled[i]);
            }
            else 
            {
                status = steric_clash_check_COMFast_IL<false>(assembled, assembled[i]);
            }
            assembled.keep();
        }
        if (status == ATTACHED)
            printf("Successful attachment on: %d\n", i);
        else
            printf("Unsuccessful attachment on: %d\n", i);
    }
    WC_prepare_structure_matrix(0, assembled[0], 1, assembled[assembled.iterator_max], 1);
    RMSD = WC_check_pair(0);
    assembled.update_WC_rmsd(RMSD);
    assembled.update_energy();
    o_string.add_string(assembled.to_string());
}

template <bool PerformChecks, STRUCTFILTER_TYPE StructCheck>
void create_custom_structure_list(DimerLibArray &Lib, RNADataArray &assembled, output_string &o_string, uint32_t num_strs)
{
    attach_status status = ATTACHED;
    // double RMSD = 0;
    if (PerformChecks == true)
    {
        for (uint32_t i = 0; i < num_strs; i++)
        {
            assembled.overwrite(0, GLOBAL_INPUT_INDICES_LIST[i][0], Lib);
            for (int j = 1; j <= assembled.iterator_max; j++)
            {
                assembled.overwrite(j, GLOBAL_INPUT_INDICES_LIST[i][j], Lib);
                rotate(assembled.current(), assembled[j]);
                status = steric_clash_check_COMFast(assembled, assembled[j]);
                if (status != ATTACHED)
                {
                    // assembled.rollback_by(j + 1);
                    printf("Failed Attach\n");
                    break;
                }
                else
                {
                    // printf("Attached\n");
                    assembled.keep();
                }
            }
            if (status == ATTACHED)
            {
                if constexpr (StructCheck == HAIRPIN)
                {
                    WC_prepare_structure_matrix(0, assembled[0], 0, assembled[assembled.iterator_max], 1);
                    double RMSD = WC_check_pair(0);
                    assembled.update_WC_rmsd(RMSD);
                }
                assembled.update_energy();
                o_string.add_string(assembled.to_string());
            }
            assembled.rollback_by(assembled.iterator);
        }
        assembled.count = assembled.iterator_max + 1;
        assembled.iterator = assembled.iterator_max;
    }
    else
    {
        for (uint32_t i = 0; i < num_strs; i++)
        {
            assembled.overwrite(0, GLOBAL_INPUT_INDICES_LIST[i][0], Lib);
            for (int j = 1; j <= assembled.iterator_max; j++)
            {
                assembled.overwrite(j, GLOBAL_INPUT_INDICES_LIST[i][j], Lib);
                rotate(assembled.current(), assembled[j]);
                assembled.keep();
            }
            if constexpr (StructCheck == HAIRPIN)
            {
                WC_prepare_structure_matrix(0, assembled[0], 0, assembled[assembled.iterator_max], 1);
                double RMSD = WC_check_pair(0);
                assembled.update_WC_rmsd(RMSD);
            }
            assembled.update_energy();
            o_string.add_string(assembled.to_string());
            assembled.rollback_by(assembled.iterator);
        }
        assembled.count = assembled.iterator_max + 1;
        assembled.iterator = assembled.iterator_max;
    }
}
template void create_custom_structure_list<PERFORM_CHECKS_ON_CUSTOM_BUILD, HAIRPIN>(DimerLibArray &Lib, RNADataArray &assembled, output_string &o_string, uint32_t num_strs);
template void create_custom_structure_list<PERFORM_CHECKS_ON_CUSTOM_BUILD, NONE>(DimerLibArray &Lib,  RNADataArray &assembled, output_string &o_string, uint32_t num_strs);
