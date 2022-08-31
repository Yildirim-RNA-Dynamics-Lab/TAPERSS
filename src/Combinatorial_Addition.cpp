#include "Combinatorial_Addition.hpp"


bool combinatorial_addition(DimerLibArray &Lib, RNADataArray &assembled, CMB_Manager &manager, output_string &o_string, DimerLibArray &WC_Lib)
{
    int working_position = assembled.iterator + 1; // Position in sequence where new DNT will be attached.
    DimerLib *Library = Lib[working_position];
    RNAData *base;                // Already attached base which will have new DNT attached to
    RNAData *attach;              // DNT which will be attached
    attach_status status = FAILED; // Output from checking functions

    DEBUG(printf("Early: checking assembled 0: %ld, working position: %d\n", assembled[0]->id, working_position));

    attach = new RNAData(Lib, working_position, 0); // For initialization only
    for (int i = 0; i < Library->count; i++)
    {
        if (Library->flags[i] != NO_FLAG)
        {
            continue;
        }
        if (assembled.is_empty())
        {
            attach->overwrite(Lib, working_position, i);
            assembled.add_move(attach);
            manager.attach_attempt(working_position, i);
            DEBUG(printf("attach: %ld moved @ empty\n", attach->id));
            status = NOT_CHECKED;
            break;
        }
        base = assembled.current();
        attach->overwrite(Lib, working_position, i);
        manager.attach_attempt(working_position, i);
        if ((status = rotate(base, attach)) != ATTACHED)
        {
            continue;
        }
        attach->update_submatrices();
        if ((status = check_attachment(assembled, attach)) == ATTACHED)
        {
            break;
        }
        /*
        else if (status == FAILED_SC)
        {
            assembled.add_move(attach);
            o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
            o_string.add_string((char *)"!!!!!!Unexpected steric clash!!!!!!\n", sizeof("!!!!!!Unexpected steric clash!!!!!!\n"));
            return true;
        }
        */
    }
    if (status == FAILED)
    {
        delete attach;
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
        assembled.add_move(attach);
        DEBUG(printf("attach: %ld moved @ Attached @ %d\n", attach->id, working_position));
    }
    if (assembled.is_complete())
    {
        assembled.update_energy();
        if (GLOBAL_PERFORM_STRUCTCHECK == HAIRPIN)
        {
            if (is_WC_pair(assembled, WC_Lib, 0, assembled.iterator_max, 0))
            {
                o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
                manager.hairpins_built++;
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
        }
    }
    return false;
}

bool combinatorial_addition_IL(DimerLibArray &Lib, RNADataArrayInternalLoop &assembled, CMB_Manager &manager, output_string &o_string, DimerLibArray &WC_Lib)
{
    int working_position = assembled.iterator + 1; // Position in sequence where new DNT will be attached.
    DimerLib *Library = Lib[working_position];
    RNAData *base;                // Already attached base which will have new DNT attached to
    RNAData *attach;              // DNT which will be attached
    attach_status status = FAILED; // Output from checking functions

    DEBUG(printf("Early: checking assembled 0: %ld, working position: %d\n", assembled[0]->id, working_position));

    if(!assembled.is_empty())
    {
        base = assembled.current();
    }
    
    //DIE;
    //printf("Base name: %s\n", base->name);
    attach = new RNAData(Lib, working_position, 0); // For initialization only
    for (int i = 0; i < Library->count; i++)
    {
        if (Library->flags[i] != NO_FLAG)
        {
            continue;
        }
        if (assembled.is_empty())
        {
            attach->overwrite(Lib, working_position, i);
            assembled.add_move(attach);
            manager.attach_attempt(working_position, i);
            DEBUG(printf("attach: %ld moved @ empty\n", attach->id));
            status = NOT_CHECKED;
            break;
        }
        if(assembled.inLeft_or_inRight(working_position) == true) 
        {
            attach->overwrite(Lib, working_position, i);
            manager.attach_attempt(working_position, i);
            if(assembled.prepare_right(attach, WC_Lib) == true) 
            {
                status = ATTACHED;
                break;
            }
            else 
            {
                continue;
            }
            /*
            o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
            o_string.add_string((char *)"!!!!!!Unexpected steric clash!!!!!!\n", sizeof("!!!!!!Unexpected steric clash!!!!!!\n"));
            */
        }
        base = assembled.current();
        //printf("Base name: %s\n", base->name);
        //print_gsl_matrix(base->data_matrix);
        attach->overwrite(Lib, working_position, i);
        //print_gsl_matrix(attach->data_matrix);
        manager.attach_attempt(working_position, i);
        //printf("Attaching attempt: working pos = %d, i = %d\n", working_position, i);
        if ((status = rotate(base, attach)) != ATTACHED)
        {
            continue;
        }
        attach->update_submatrices();
        //printf("---------Before check attachment\n");
        if ((status = check_attachment(assembled, attach)) == ATTACHED)
        {
            break;
        }
        else if (status == FAILED_SC)
        {
            assembled.add_move(attach);
            o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
            o_string.add_string((char *)"!!!!!!Unexpected steric clash!!!!!!\n", sizeof("!!!!!!Unexpected steric clash!!!!!!\n"));
            return true;
        }
        //DIE;
    }
    if (status == FAILED)
    {
        delete attach;
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
        assembled.add_move(attach);
        DEBUG(printf("attach: %ld moved @ Attached @ %d\n", attach->id, working_position));
    }
    if (assembled.is_complete())
    {
        assembled.update_energy();
        if (GLOBAL_PERFORM_STRUCTCHECK == INTERNAL_LOOP)
        {
            if (is_WC_pair(assembled, WC_Lib, 0, assembled.iterator_max, 0))
            {
                //printf("IS WC PAIR\n");
                o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
                manager.internal_loops_built++;
                if(manager.internal_loops_built > 10) 
                {
                    return true;
                }
            }
            else 
            {
                //printf("IS NOT WC PAIR\n");
                //o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
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
            //return true;
        }
    }
    return false;
}
