#include "Combinatorial_Addition.hpp"

bool combinatorial_addition(DimerLibArray &Lib, RNADataArray &assembled, CMB_Manager &manager, output_string &o_string)
{
    int working_position = assembled.iterator + 1; // Position in sequence where new DNT will be attached.
    DimerLib *Library = Lib[working_position];
    //RNAData *base;                 // Already attached base which will have new DNT attached to
    //RNAData *attach;               // DNT which will be attached
    attach_status status = FAILED; // Output from checking functions
    double RMSD;                   // Output of RMSD function.
    //printf("WORKING ON POSITION: %d\n", working_position);
    for (int i = 0; i < Library->count; i++)
    {
        if (Lib.Flags[working_position][i] != NO_FLAG)
        {
            continue;
        }
        //printf("Using Structure: %d\n", i);
        if (assembled.is_empty())
        {
            assembled.overwrite(working_position, i, Lib);
            assembled.keep();
            manager.attach_attempt(working_position, i);
            status = NOT_CHECKED;
            break;
        }
        assembled.overwrite(working_position, i, Lib);
        manager.attach_attempt(working_position, i);
        //assembled.print_index(1);
        if ((status = rotate(assembled.current(), assembled[working_position])) != ATTACHED)
        {
            continue;
        }
        if ((status = steric_clash_check_COMFast(assembled, assembled[working_position])) == ATTACHED)
        {
            break;
        }
        else
        {
            *assembled[working_position]->_flag = NOT_USABLE;
        }
    }

    if (status == FAILED)
    {
        //printf("All failed!\n");
        //assembled.print_index();
        if (manager.is_at_end())
        {
            manager.check_lib_completion();
            if ((manager.get_reset_count()) == manager.last_attempted[0] + 1)
            {
                return true;
            }
            assembled.rollback_by(manager.get_reset_count() - 1);
            Lib.reset_flags(manager.libs_completed);
            manager.clear_attempts();
        }
        else
        {
            assembled.rollback();
        }
        return false;
    }
    else if (status == NOT_CHECKED)
    {
        return false;
    }
    else if (status == ATTACHED)
    {
        assembled.keep();
    }
    if (assembled.is_complete())
    {
        //assembled.print_index(0);
        //DIE;
        if (GLOBAL_PERFORM_STRUCTCHECK == HAIRPIN)
        {
            WC_prepare_structure_matrix(0, assembled[0]->data_matrix, assembled[0]->WC_submatrix_rows[0], assembled[0]->count_per_WC_sub[0],
                                        assembled[assembled.iterator_max]->data_matrix, assembled[assembled.iterator_max]->WC_submatrix_rows[1],
                                        assembled[assembled.iterator_max]->count_per_WC_sub[1]);
            RMSD = WC_check_pair(0);
            //assembled.print_index();
            
            //assembled.print_index();
            if (RMSD <= GLOBAL_WC_RMSD_LIMIT)
            {
                //printf("RMSD WC = %f\n",RMSD);
                assembled.update_WC_rmsd(RMSD);
                assembled.update_energy();
                o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
                manager.hairpins_built++;
                if constexpr (STRUCTURE_BUILD_LIMIT == true)
                {
                    if (manager.hairpins_built == GLOBAL_STRUCTURE_LIMIT_COUNT)
                    {
                        return true;
                    }
                }
            }
        }
        else
        {
            assembled.update_energy();
            o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
        }
        manager.strs_built++;
        if (manager.is_at_end())
        {
            manager.check_lib_completion();
            if ((manager.get_reset_count()) == Lib.count)
            {
                return true;
            }
            assembled.rollback_by(manager.get_reset_count());
            Lib.reset_flags(manager.libs_completed);
            manager.clear_attempts();
        }
        else
        {
            DEBUG(printf("Rolling back b/c: COMPLETED NOT AT LIB END\n"));
            assembled.rollback();
        }
    }
    return false;
}

bool combinatorial_addition_IL(DimerLibArray &Lib, RNADataArrayInternalLoop &assembled, CMB_Manager &manager, output_string &o_string, DimerLibArray &WC_Lib)
{
    int working_position = assembled.iterator + 1; // Position in sequence where new DNT will be attached.
    DimerLib *Library = Lib[working_position];
    //RNAData *base;                 // Already attached base which will have new DNT attached to
    //RNAData *attach;               // DNT which will be attached
    attach_status status = FAILED; // Output from checking functions
    double RMSD;

    DEBUG(printf("Early: checking assembled 0: %ld, working position: %d\n", assembled[0]->id, working_position));
    for (int i = 0; i < Library->count; i++)
    {
        if (Lib.Flags[working_position][i] != NO_FLAG)
        {
            continue;
        }
        if (assembled.is_empty())
        {
            assembled.overwrite(working_position, i, Lib);
            assembled.keep();
            manager.attach_attempt(working_position, i);
            DEBUG(printf("attach: %ld moved @ empty\n", assembled[working_position]->id));
            status = NOT_CHECKED;
            break;
        }
        if (assembled.inLeft_or_inRight(working_position) == true)
        {
            assembled.overwrite(working_position, i, Lib);
            manager.attach_attempt(working_position, i);
            if (assembled.prepare_right(assembled[working_position], WC_Lib) == true)
            {
                status = ATTACHED;
                break;
            }
            else
            {
                continue;
            }
        }
        assembled.overwrite(working_position, i, Lib);
        manager.attach_attempt(working_position, i);
        if ((status = rotate(assembled.current(), assembled[working_position])) != ATTACHED)
        {
            continue;
        }
        // printf("---------Before check attachment\n");
        if ((status = steric_clash_check_COMFast_IL(assembled, assembled[working_position])) == ATTACHED)
        {
            break;
        }
        // DIE;
    }
    if (status == FAILED)
    {
        if (manager.is_at_end())
        {
            manager.check_lib_completion();
            if ((manager.get_reset_count()) == manager.last_attempted[0] + 1)
            {
                return true;
            }
            assembled.rollback_by(manager.get_reset_count() - 1);
            Lib.reset_flags(manager.libs_completed);
            manager.clear_attempts();
        }
        else
        {
            assembled.rollback();
        }
        return false;
    }
    else if (status == NOT_CHECKED)
    {
        return false;
    }
    else if (status == ATTACHED)
    {
        assembled.keep();
        //DEBUG(printf("attach: %ld moved @ Attached @ %d\n", attach->id, working_position));
    }
    if (assembled.is_complete())
    {
        assembled.update_energy();
        if (GLOBAL_PERFORM_STRUCTCHECK == INTERNAL_LOOP)
        {
            WC_prepare_structure_matrix(0, assembled[0]->data_matrix, assembled[0]->WC_submatrix_rows[0], assembled[0]->count_per_WC_sub[0],
                                        assembled[assembled.iterator_max]->data_matrix, assembled[assembled.iterator_max]->WC_submatrix_rows[1],
                                        assembled[assembled.iterator_max]->count_per_WC_sub[1]);
            RMSD = WC_check_pair(0);
            if (RMSD <= GLOBAL_WC_RMSD_LIMIT)
            {
                assembled.update_WC_rmsd(RMSD);
                o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
                manager.internal_loops_built++;
                if (assembled.TMP_END == true)
                {
                    return true;
                }
                if constexpr (STRUCTURE_BUILD_LIMIT == true)
                {
                    if (manager.internal_loops_built == GLOBAL_STRUCTURE_LIMIT_COUNT)
                    {
                        return true;
                    }
                }
            }
            else
            {
                // printf("IS NOT WC PAIR\n");
                // o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
            }
        }
        else
        {
            o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
        }
        manager.strs_built++;
        if (manager.is_at_end())
        {
            manager.check_lib_completion();
            if ((manager.get_reset_count()) == Lib.count)
            {
                return true;
            }
            assembled.rollback_by(manager.get_reset_count());
            Lib.reset_flags(manager.libs_completed);
            manager.clear_attempts();
        }
        else
        {
            DEBUG(printf("Rolling back b/c: COMPLETED NOT AT LIB END\n"));
            assembled.rollback();
            // return true;
        }
    }
    return false;
}
