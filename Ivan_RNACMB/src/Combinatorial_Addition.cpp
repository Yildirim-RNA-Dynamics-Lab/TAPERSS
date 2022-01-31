#include "Combinatorial_Addition.hpp"

bool combinatorial_addition(DimerLibArray &Lib, RNA_data_array &assembled, CMB_Manager &manager, output_string &o_string, DimerLibArray &WC_Lib)
{
    int working_position = assembled.iterator + 1; // Position in sequence where new DNT will be attached.
    // printf("Working postion: %d, iterator:%d, iterator_max:%d\n", working_position, assembled.iterator, assembled.iterator_max);
    DimerLib *Library = Lib[working_position];
    RNA_data *base;                // Already attached base which will have new DNT attached to
    RNA_data *attach;              // DNT which will be attached
    attach_status status = FAILED; // Output from checking functions

    DEBUG(printf("Early: checking assembled 0: %ld, working position: %d\n", assembled[0]->id, working_position));

    attach = new RNA_data(Lib, working_position, 0); // For initialization only
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
        if (GLOBAL_PERFORM_HAIRPIN_CHECK)
        {
            if (is_hairpin(assembled, WC_Lib))
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