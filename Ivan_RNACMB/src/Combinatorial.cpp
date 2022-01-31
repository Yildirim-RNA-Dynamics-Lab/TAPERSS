#include "Combinatorial.hpp"

CMB_Manager::CMB_Manager(DimerLibArray &LA)
{
    count = LA.count;

    count_per_lib = (int *)malloc(sizeof(int) * LA.count);
    attach_attempted = (bool **)malloc(sizeof(bool *) * LA.count);

    for (int i = 0; i < LA.count; i++)
    {
        count_per_lib[i] = LA[i]->count;
        attach_attempted[i] = (bool *)calloc(LA[i]->count, sizeof(bool));
    }

    last_attempted[0] = 0;
    last_attempted[1] = 0;

    strs_built = 0;
    hairpins_built = 0;
    libs_completed = (bool *)calloc(count, sizeof(bool));
}

CMB_Manager::~CMB_Manager()
{
    for (int i = 0; i < count; i++)
    {
        free(attach_attempted[i]);
    }
    free(attach_attempted);
    free(count_per_lib);
    free(libs_completed);
}

void CMB_Manager::attach_attempt(int i, int j)
{
    attach_attempted[i][j] = true;
    last_attempted[0] = i;
    last_attempted[1] = j;
}

bool CMB_Manager::is_at_end()
{
    if (last_attempted[1] == count_per_lib[last_attempted[0]] - 1)
    {
        return true;
    }
    return false;
}

void CMB_Manager::check_lib_completion()
{
    int counter = 0;

    int first_completed = last_attempted[0];

    for (int i = last_attempted[0]; i >= 0; i--)
    {
        libs_completed[i] = false;
        if (attach_attempted[i][count_per_lib[i] - 1] == true)
        {
            libs_completed[i] = true;
            counter++;
            if (i < first_completed)
                first_completed = i;
            // printf("Lib %d complete! Max = %d\n", i, count_per_lib[i] - 1);
        }
        /*else
        {
            printf("Lib %d not complete! Max = %d\n", i, count_per_lib[i] - 1);
        }*/
    }
    if (counter != (count - first_completed))
    {
        bool swap = true;
        for (int i = last_attempted[0]; i >= first_completed; i--)
        {
            if (libs_completed[i] == false)
            {
                swap = false;
            }
            libs_completed[i] = swap;
        }
    }

    if (counter == last_attempted[0] + 1)
        ;
    else
        libs_completed[0] = false;
    return;
}

void CMB_Manager::clear_attempts()
{
    for (int i = 0; i < count; i++)
    {
        if (libs_completed[i] == true)
        {
            // printf("reseting attempts for %d\n", i);
            for (int j = 0; j < count_per_lib[i]; j++)
            {
                attach_attempted[i][j] = false;
            }
        }
    }
}

int CMB_Manager::get_reset_count()
{
    int n_reset = 0;
    for (int i = 0; i < last_attempted[0] + 1; i++)
    {
        if (libs_completed[i] == true)
        {
            n_reset++;
        }
    }
    // printf("n_reset = %d\n", n_reset);
    return n_reset;
}

void CMB_Manager::successful_construction()
{
    strs_built++;
}

void update_energy(RNA_data_array &sequence)
{
    float energy = 0.0;
    gsl_vector_view A, B;
    for (int i = 0; i < sequence.count; i++)
    {
        energy += sequence[i]->energy;
    }
    // printf("Stacking Energy:: %f\n", energy);
    for (int i = 0; i < sequence.count - 1; i++)
    {
        for (int j = sequence.count - 1; j > i + 1; j--)
        {
            if (i != j)
            {
                // printf("testing %d with %d\n", i, j);
                for (int k = 0; k < sequence[i]->count; k++)
                {
                    if (i != 0)
                    {
                        if (sequence[i]->atom_data->dnt_pos[k] != 2)
                        {
                            continue;
                        }
                    }
                    if (sequence[i]->atom_data->charges[k] == NEUTRAL)
                    {
                        continue;
                    }
                    if (sequence[i]->has_interaction[k] == true)
                    {
                        continue;
                    }
                    for (int l = 0; l < sequence[j]->count; l++)
                    {
                        if (sequence[j]->has_interaction[l] == true)
                        {
                            continue;
                        }
                        if (sequence[j]->atom_data->dnt_pos[l] != 2)
                        {
                            continue;
                        }
                        // printf("charge = %d\n", sequence[j]->atom_data->charges[l]);
                        if (sequence[j]->atom_data->charges[l] == NEUTRAL)
                        {
                            continue;
                        }
                        A = gsl_matrix_row(sequence[i]->data_matrix, k);
                        B = gsl_matrix_row(sequence[j]->data_matrix, l);
                        if (distance(&A.vector, &B.vector) < INTERACTION_DISTANCE)
                        {
                            if (sequence[i]->atom_data->charges[k] == sequence[j]->atom_data->charges[l])
                            {
                                // energy += 1;
                            }
                            else
                            {
                                // printf("A: Index: %d, Resid: %d, name: %s, charge = %d\n", i, sequence[i]->atom_data->dnt_pos[k], sequence[i]->atom_data->name[k], sequence[i]->atom_data->charges[k]);
                                // printf("B: Index: %d, Resid: %d, name: %s, charge = %d\n", j, sequence[j]->atom_data->dnt_pos[l], sequence[j]->atom_data->name[l], sequence[j]->atom_data->charges[l]);
                                energy -= 1;
                                sequence[j]->has_interaction[l] = true;
                                sequence[i]->has_interaction[k] = true;
                            }
                        }
                    }
                }
            }
        }
        sequence.reset_interactions();
    }
    // printf("energy = %f\n", energy);
    sequence.structure_energy = energy;
}

attach_status rotate(RNA_data *reference, RNA_data *rotated)
{
    double rmsd_;
    double COMP[] = {0, 0, 0};
    double COMQ[] = {0, 0, 0};
    attach_status status = FAILED;

    gsl_matrix *MODEL = rotated->data_matrix;
    gsl_matrix *P = rotated->get_target_matrix_copy(0);
    gsl_matrix *Q = reference->get_target_matrix_copy(1);
    gsl_matrix *R;

    R = kabsch_get_rotation_matrix_generic(P, Q, COMP, COMQ);
    rmsd_ = rmsd_generic(P, Q);
    if (rmsd_ <= GLOBAL_RMSD_LIMIT)
    {
        status = ATTACHED;
        translate_matrix(COMP, MODEL, -1.0F);
        apply_rotation_matrix(R, MODEL);
        translate_matrix(COMQ, MODEL, 1.0F);
    }
    else
    {
        // printf("RMSD FAIL\n");
    }

    gsl_matrix_free(P);
    gsl_matrix_free(Q);
    gsl_matrix_free(R);
    return status;
}

bool overlap_check(RNA_data_array &sequence, RNA_data *attach)
{
    gsl_vector_view A, B;
    double radius_1, radius_2;
    double dist;
    // printf("Attaching resid: %d\n", sequence.count);
    for (int i = 0; i < sequence.count; i++)
    {
        for (int j = 0; j < sequence[i]->count; j++)
        {
            if (i == sequence.count - 1 && sequence[i]->atom_data->dnt_pos[j] == 2)
            {
                // printf("residue: %d, pos: %s\n", sequence[i]->atom_data->dnt_pos[j], sequence[i]->atom_data->name[j]);
                switch (sequence[i]->atom_data->atom_ids[j])
                {
                case OP1:
                case OP2:
                case O1P:
                case O2P:
                case P:
                case O5p:
                case C5p:
                    break;
                default:
                    continue;
                }
            }
            else if (i != 0 && sequence[i]->atom_data->dnt_pos[j] == 1)
            {
                continue;
            }
            for (int k = 0; k < attach->count; k++)
            {
                if (attach->atom_data->dnt_pos[k] == 1)
                {
                    continue;
                }
                switch (sequence[i]->atom_data->name[j][0])
                {
                case 'C':
                    radius_1 = RADIUS_C;
                    break;
                case 'N':
                    radius_1 = RADIUS_N;
                    break;
                case 'O':
                    radius_1 = RADIUS_O;
                    break;
                case 'P':
                    radius_1 = RADIUS_P;
                    break;
                }
                switch (attach->atom_data->name[k][0])
                {
                case 'C':
                    radius_2 = RADIUS_C;
                    break;
                case 'N':
                    radius_2 = RADIUS_N;
                    break;
                case 'O':
                    radius_2 = RADIUS_O;
                    break;
                case 'P':
                    radius_2 = RADIUS_P;
                    break;
                }
                A = gsl_matrix_row(attach->data_matrix, k);
                B = gsl_matrix_row(sequence[i]->data_matrix, j);
                dist = distance(&A.vector, &B.vector);
                // printf("ref %d: s: %s:%d a: %s:%d\tdist = %f\tvdw = %f\n", i,sequence[i]->atom_data->name[j], sequence[i]->atom_data->dnt_pos[j], attach->atom_data->name[k], attach->atom_data->dnt_pos[k], dist, (radius_1 + radius_2));
                if (dist < (radius_1 + radius_2))
                {
                    // printf("ATT: %f %f %f\n", gsl_vector_get(&A.vector, 0), gsl_vector_get(&A.vector, 1), gsl_vector_get(&A.vector, 2));
                    // printf("REF: %f %f %f\n", gsl_vector_get(&B.vector, 0), gsl_vector_get(&B.vector, 1), gsl_vector_get(&B.vector, 2));

                    // printf("VDW FAIL:: ref %d: s: %s:%d a: %s:%d\tdist = %f\tvdw = %f\n", i, sequence[i]->atom_data->name[j], sequence[i]->atom_data->dnt_pos[j], attach->atom_data->name[k], attach->atom_data->dnt_pos[k], dist, (radius_1 + radius_2));
                    return false;
                }
            }
        }
    }
    // exit(0);
    return true;
}

attach_status check_attachment(RNA_data_array &sequence, RNA_data *attach)
{
    if (!overlap_check(sequence, attach))
    {
        // printf("Overlap detected!\n");
        *attach->_flag = NOT_USABLE;
        return FAILED;
    }
    return ATTACHED;
}

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
        update_energy(assembled);
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
            // printf("lib completed = %d\n", manager.get_reset_count());

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
