#include "RNA_data.hpp"

RNA_data::RNA_data(DimerLibArray &L, int i, int j, bool WC)
{
    atom_data = (L[i]->atom_data);
    data_matrix = gsl_matrix_alloc(L[i]->data_matrices[j]->size1 + 1, L[i]->data_matrices[j]->size2);
    gsl_matrix_memcpy(data_matrix, L[i]->data_matrices[j]);
    energy = L[i]->energy[j];
    name = (char *)malloc(sizeof(char) * 3);
    memcpy(name, L[i]->name, sizeof(char) * 3);
    count = L[i]->atom_data->count;
    _flag = &(L[i]->flags[j]);
    position_in_lib[0] = i;
    position_in_lib[1] = j;
    position_max = L[i]->count - 1;
    id = rna_dat_tracker++;
    is_for_WC = WC;
    WC ? make_WC_submatrices(true) : make_submatrices();
    WC_secondary = false;
    has_interaction = (bool *)calloc(count, sizeof(bool));
}

void RNA_data::overwrite(DimerLibArray &L, int i, int j)
{
    gsl_matrix_memcpy(data_matrix, L[i]->data_matrices[j]);
    energy = L[i]->energy[j];
    memcpy(name, L[i]->name, sizeof(char) * 3);
    count = L[i]->atom_data->count;
    _flag = &(L[i]->flags[j]);
    position_in_lib[0] = i;
    position_in_lib[1] = j;
    is_for_WC ? update_WC_submatrices() : update_submatrices();
}

RNA_data::~RNA_data()
{
    // printf("Deallocationg %ld\n", id);
    gsl_matrix_free(data_matrix);
    free(name);
    free(has_interaction);
    if (is_for_WC)
    {
        gsl_matrix_free(WC_submatrices[0]);
        gsl_matrix_free(WC_submatrices[1]);
        free(WC_submatrix_rows[0]);
        free(WC_submatrix_rows[1]);
        free(WC_submatrices);
    }
    else if (WC_secondary)
    {
        gsl_matrix_free(WC_submatrices[0]);
        gsl_matrix_free(WC_submatrices[1]);
        free(WC_submatrix_rows[0]);
        free(WC_submatrix_rows[1]);
        free(WC_submatrices);
        gsl_matrix_free(submatrices[0]);
        gsl_matrix_free(submatrices[1]);
        free(submatrix_rows[0]);
        free(submatrix_rows[1]);
        free(submatrices);
    }
    else
    {
        gsl_matrix_free(submatrices[0]);
        gsl_matrix_free(submatrices[1]);
        free(submatrix_rows[0]);
        free(submatrix_rows[1]);
        free(submatrices);
    }
}

void RNA_data::make_submatrices()
{
    submatrices = (gsl_matrix **)malloc(sizeof(gsl_matrix *) * 2);

    for (int i = 0; i < 2; i++)
    {
        count_per_sub[i] = get_target(i);
        submatrices[i] = gsl_matrix_alloc(count_per_sub[i], 3);
        submatrix_rows[i] = (int *)malloc(sizeof(int) * count_per_sub[i]);

        int rel_idx = 0;

        for (int j = 0; j < count; j++)
        {
            for (int k = 0; k < count_per_sub[i]; k++)
            {
                if ((target[k] == atom_data->atom_ids[j]) && (atom_data->dnt_pos[j] - 1) == i)
                {
                    gsl_matrix_set(submatrices[i], rel_idx, 0, gsl_matrix_get(data_matrix, j, 0));
                    gsl_matrix_set(submatrices[i], rel_idx, 1, gsl_matrix_get(data_matrix, j, 1));
                    gsl_matrix_set(submatrices[i], rel_idx, 2, gsl_matrix_get(data_matrix, j, 2));
                    submatrix_rows[i][rel_idx] = j;
                    rel_idx++;
                }
            }
        }
        free(target);
    }
}

void RNA_data::update_submatrices()
{
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < count_per_sub[i]; j++)
        {
            gsl_matrix_set(submatrices[i], j, 0, gsl_matrix_get(data_matrix, submatrix_rows[i][j], 0));
            gsl_matrix_set(submatrices[i], j, 1, gsl_matrix_get(data_matrix, submatrix_rows[i][j], 1));
            gsl_matrix_set(submatrices[i], j, 2, gsl_matrix_get(data_matrix, submatrix_rows[i][j], 2));
        }
    }
}

size_t RNA_data::get_target(int res)
{
    constexpr atom_id targetA[] = {N9, C8, N7, C5, C6, N1, C2, N3, C4, C1p, C2p, C3p, C4p, O4p}; // 14
    constexpr atom_id targetC[] = {N1, C2, N3, C4, C5, C6, C1p, C2p, C3p, C4p, O4p};             // 11
    constexpr atom_id targetG[] = {N9, C8, N7, C5, C6, N1, C2, N3, C4, C1p, C2p, C3p, C4p, O4p}; // 14
    constexpr atom_id targetU[] = {N1, C2, N3, C4, C5, C6, C1p, C2p, C3p, C4p, O4p};             // 11

    if (name[res] == 'A')
    {
        target = (atom_id *)malloc(sizeof(targetA));
        memcpy(target, targetA, sizeof(targetA));
        return (sizeof(targetA) / sizeof(atom_id));
    }
    else if (name[res] == 'C')
    {
        target = (atom_id *)malloc(sizeof(targetC));
        memcpy(target, targetC, sizeof(targetC));
        return (sizeof(targetC) / sizeof(atom_id));
    }
    else if (name[res] == 'G')
    {
        target = (atom_id *)malloc(sizeof(targetG));
        memcpy(target, targetG, sizeof(targetG));
        return (sizeof(targetG) / sizeof(atom_id));
    }
    else if (name[res] == 'U')
    {
        target = (atom_id *)malloc(sizeof(targetU));
        memcpy(target, targetU, sizeof(targetU));
        return (sizeof(targetU) / sizeof(atom_id));
    }
    else
    {
        printf("RNA_data: input file has problems, please check at %d\n", __LINE__);
        exit(1);
    }
}

gsl_matrix *RNA_data::get_target_matrix(int res)
{
    return submatrices[res];
}

gsl_matrix *RNA_data::get_target_matrix_copy(int res)
{
    gsl_matrix *rtn = gsl_matrix_alloc(count_per_sub[res], 3);
    gsl_matrix_memcpy(rtn, submatrices[res]);
    return rtn;
}

void RNA_data::make_WC_submatrices(bool first_run)
{
    // printf("%ld: name is %s\n", id, name);
    if (!first_run)
    {
        if (is_for_WC)
        {
            // printf("%ld: Was allocated, already as evidenced by %d\n", id, count_per_WC_sub[0]);
            return;
        }
        else if (WC_secondary)
        {
            update_WC_submatrices();
            return;
        }
        else if (!WC_secondary)
        {
            WC_secondary = true;
        }
    }
    /*else
    {
        // printf("first run for %ld\n", id);
    }*/

    WC_submatrices = (gsl_matrix **)malloc(sizeof(gsl_matrix *) * 2);

    for (int i = 0; i < 2; i++)
    {
        count_per_WC_sub[i] = get_WC_target(i);
        WC_submatrices[i] = gsl_matrix_alloc(count_per_WC_sub[i], 3);
        WC_submatrix_rows[i] = (int *)malloc(sizeof(int) * count_per_WC_sub[i]);

        int rel_idx = 0;

        for (int j = 0; j < count; j++)
        {
            for (unsigned int k = 0; k < count_per_WC_sub[i]; k++)
            {
                if ((target[k] == atom_data->atom_ids[j]) && (atom_data->dnt_pos[j] - 1) == i)
                {
                    gsl_matrix_set(WC_submatrices[i], rel_idx, 0, gsl_matrix_get(data_matrix, j, 0));
                    gsl_matrix_set(WC_submatrices[i], rel_idx, 1, gsl_matrix_get(data_matrix, j, 1));
                    gsl_matrix_set(WC_submatrices[i], rel_idx, 2, gsl_matrix_get(data_matrix, j, 2));
                    WC_submatrix_rows[i][rel_idx] = j;
                    rel_idx++;
                    // printf("%s:%ld\tGetting Atom %s from row %d\n", name, id, atom_data->name[j], j);
                }
            }
        }
        free(target);
    }
}

void RNA_data::update_WC_submatrices()
{
    for (int i = 0; i < 2; i++)
    {
        for (unsigned int j = 0; j < count_per_WC_sub[i]; j++)
        {
            gsl_matrix_set(WC_submatrices[i], j, 0, gsl_matrix_get(data_matrix, WC_submatrix_rows[i][j], 0));
            gsl_matrix_set(WC_submatrices[i], j, 1, gsl_matrix_get(data_matrix, WC_submatrix_rows[i][j], 1));
            gsl_matrix_set(WC_submatrices[i], j, 2, gsl_matrix_get(data_matrix, WC_submatrix_rows[i][j], 2));
        }
    }
}

size_t RNA_data::get_WC_target(int res)
{
    constexpr atom_id targetA[] = {N9, C8, N7, C5, C6, N1, C2, N3, C4, N6 /*, C1p, C2p, C3p, C4p, O4p*/};
    constexpr atom_id targetC[] = {N1, C2, N3, C4, C5, C6, N4, O2 /*, C1p, C2p, C3p, C4p, O4p*/};
    constexpr atom_id targetG[] = {N9, C8, N7, C5, C6, N1, C2, N3, C4, N2, O6 /*, C1p, C2p, C3p, C4p, O4p*/};
    constexpr atom_id targetU[] = {N1, C2, N3, C4, C5, C6, O4, O2 /*, C1p, C2p, C3p, C4p, O4p*/};

    if (name[res] == 'A')
    {
        target = (atom_id *)malloc(sizeof(targetA));
        memcpy(target, targetA, sizeof(targetA));
        return (sizeof(targetA) / sizeof(atom_id));
    }
    else if (name[res] == 'C')
    {
        target = (atom_id *)malloc(sizeof(targetC));
        memcpy(target, targetC, sizeof(targetC));
        return (sizeof(targetC) / sizeof(atom_id));
    }
    else if (name[res] == 'G')
    {
        target = (atom_id *)malloc(sizeof(targetG));
        memcpy(target, targetG, sizeof(targetG));
        return (sizeof(targetG) / sizeof(atom_id));
    }
    else if (name[res] == 'U')
    {
        target = (atom_id *)malloc(sizeof(targetU));
        memcpy(target, targetU, sizeof(targetU));
        return (sizeof(targetU) / sizeof(atom_id));
    }
    else
    {
        printf("RNA_data: input file has problems, please check at %d\n", __LINE__);
        exit(1);
    }
}

gsl_matrix *RNA_data::get_WC_target_matrix(int res)
{
    // printf("size of WC_submatrix = %ld\n", WC_submatrices[res]->size1);
    return WC_submatrices[res];
}

gsl_matrix *RNA_data::get_WC_target_matrix_copy(int res)
{
    gsl_matrix *rtn = gsl_matrix_alloc(count_per_WC_sub[res], 3);
    gsl_matrix_memcpy(rtn, WC_submatrices[res]);
    return rtn;
}

void RNA_data::print()
{
    for (unsigned int i = 0; i < data_matrix->size1; i++)
    {
        atom_data->print_at(i);
        for (unsigned int j = 0; j < data_matrix->size2; j++)
            printf("%8.3f", gsl_matrix_get(data_matrix, i, j));
        putchar('\n');
    }
}

void RNA_data::print_target(int res)
{
    for (int i = 0; i < count_per_sub[res]; i++)
    {
        atom_data->print_at(submatrix_rows[res][i]);
        for (unsigned int j = 0; j < data_matrix->size2; j++)
            printf("%8.3f", gsl_matrix_get(data_matrix, submatrix_rows[res][i], j));
        putchar('\n');
    }
}

void RNA_data::print_WC_target(int res)
{
    for (unsigned int i = 0; i < count_per_WC_sub[res]; i++)
    {
        atom_data->print_at(WC_submatrix_rows[res][i]);
        for (unsigned int j = 0; j < data_matrix->size2; j++)
            printf("%8.3f", gsl_matrix_get(data_matrix, WC_submatrix_rows[res][i], j));
        putchar('\n');
    }
}

void RNA_data::reset_interactions()
{
    for (int i = 0; i < count; i++)
    {
        has_interaction[i] = false;
    }
}

void RNA_data::print(int res)
{
    for (unsigned int i = 0; i < data_matrix->size1; i++)
    {
        if (atom_data->dnt_pos[i] != (res + 1))
            continue;
        atom_data->print_at(i);
        for (unsigned int j = 0; j < data_matrix->size2; j++)
            printf("%8.3f", gsl_matrix_get(data_matrix, i, j));
        putchar('\n');
    }
}

int RNA_data::to_string(char *s, int buffer_size, int string_index)
{
    for (unsigned int i = 0; i < data_matrix->size1; i++)
    {
        string_index += snprintf(&s[string_index], buffer_size - string_index, "%-6s%5d %4s %3c  %4d    ", "ATOM", atom_data->index[i], atom_data->name[i], atom_data->residue[i], atom_data->dnt_pos[i]);
        for (unsigned int j = 0; j < data_matrix->size2; j++)
        {
            string_index += snprintf(&s[string_index], buffer_size - string_index, "%8.3f", gsl_matrix_get(data_matrix, i, j));
        }
        string_index += snprintf(&s[string_index], buffer_size - string_index, "\n");
    }
    return string_index;
}

void RNA_data::print_offset(int res, int position)
{
    for (unsigned int i = 0; i < data_matrix->size1; i++)
    {
        if (atom_data->dnt_pos[i] != (res + 1))
            continue;
        printf("%-4s %-4d %-4c %-4d ", atom_data->name[i], atom_data->index[i], atom_data->residue[i], atom_data->dnt_pos[i] + position);
        for (unsigned int j = 0; j < data_matrix->size2; j++)
            printf("%-3.2f ", gsl_matrix_get(data_matrix, i, j));
        putchar('\n');
    }
}

int RNA_data::to_string_offset(int res, int position, char *s, int buffer_size, int string_index, int *idx_offset)
{
    for (unsigned int i = 0; i < data_matrix->size1; i++)
    {
        if (atom_data->dnt_pos[i] != (res + 1))
            continue;
        *idx_offset += 1;
        string_index += snprintf(&s[string_index], buffer_size - string_index, "%-6s%5d %4s %3c  %4d    ", "ATOM", *idx_offset, atom_data->name[i], atom_data->residue[i], atom_data->dnt_pos[i] + position);
        for (unsigned int j = 0; j < data_matrix->size2; j++)
            string_index += snprintf(&s[string_index], buffer_size - string_index, "%8.3f", gsl_matrix_get(data_matrix, i, j));
        string_index += snprintf(&s[string_index], buffer_size - string_index, "\n");
    }
    return string_index;
}

RNA_data_array::RNA_data_array(int size)
{
    iterator_max = size - 1;
    sequence = (RNA_data **)malloc(sizeof(RNA_data *) * size);
    iterator = -1;
    count = 0;
    structure_energy = 0;
}
RNA_data_array::~RNA_data_array()
{
    // printf("iterator is @ %d\n", iterator);
    for (int i = 0; i < count; i++)
    {
        delete sequence[i];
    }
    free(sequence);
    if (string_initialized)
        free(string_out);
}
RNA_data *RNA_data_array::operator[](int i)
{
    return sequence[i];
}
void RNA_data_array::add_copy(RNA_data *A)
{
    *sequence[++iterator] = *A;
    count++;
}
void RNA_data_array::add_move(RNA_data *A)
{
    sequence[++iterator] = A;
    count++;
    *A->_flag = USED;
}
RNA_data *RNA_data_array::current()
{
    return sequence[iterator];
}

bool RNA_data_array::is_complete()
{
    return iterator == iterator_max ? true : false;
}

void RNA_data_array::rollback()
{
    DEBUG(printf("attach: %ld deleted @ rollback @ %d\n", sequence[iterator]->id, iterator));
    delete sequence[iterator];
    iterator--;
    count--;
}

void RNA_data_array::safe_rollback() // unused
{
    if (sequence[iterator]->position_in_lib[1] == sequence[iterator]->position_max)
        ;
    else
        delete sequence[iterator];
    iterator--;
    count--;
}

void RNA_data_array::rollback_by(int amount)
{
    for (int i = iterator; i > (iterator - (amount + 1)); i--)
    {
        DEBUG(printf("attach: %ld deleted @ rollback_by @ %d\n", sequence[i]->id, i));
        delete sequence[i];
    }
    iterator -= (amount + 1);
    count -= (amount + 1);
    // printf("moving to pos: %d\n", iterator);
}

bool RNA_data_array::is_empty()
{
    if (count == 0)
        return true;
    return false;
}

void RNA_data_array::update_WC_rmsd(float rmsd_val)
{
    WC_rmsd = rmsd_val;
}

void RNA_data_array::reset_interactions()
{
    for (int i = 0; i < count; i++)
    {
        sequence[i]->reset_interactions();
    }
}

void RNA_data_array::update_energy()
{
    float energy_ = 0.0;
    gsl_vector_view A, B;
    for (int i = 0; i < count; i++)
    {
        energy_ += sequence[i]->energy;
    }
    for (int i = 0; i < count - 1; i++)
    {
        for (int j = count - 1; j > i + 1; j--)
        {
            if (i != j)
            {
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
                                energy_ -= 1;
                                sequence[j]->has_interaction[l] = true;
                                sequence[i]->has_interaction[k] = true;
                            }
                        }
                    }
                }
            }
        }
        reset_interactions();
    }
    structure_energy = energy_;
}

void RNA_data_array::printall()
{
    sequence[0]->print();
    for (int i = 1; i < count; i++)
    {
        sequence[i]->print_offset(1, i);
    }
}

void RNA_data_array::initialize_string()
{
    get_atom_sum();
    // printf("atomsum = %d\n", atom_sum);

    string_buffer = 54 * atom_sum;
    string_out = (char *)malloc(sizeof(char) * string_buffer);
    string_index = 0;
    model_count = 0;

    string_initialized = true;
}

int RNA_data_array::get_atom_sum()
{
    if (string_initialized)
        return atom_sum;
    atom_sum = 0;
    for (int i = 0; i < count; i++)
    {
        atom_sum += sequence[i]->count;
    }
    return atom_sum;
}

int RNA_data_array::out_string_header()
{
    if (GLOBAL_WRITE_COORDINATES)
    {
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "MODEL %d\n", ++model_count);
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "REMARK ");
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%s ", GLOBAL_INPUT_SEQUENCE);
        for (int i = 0; i < count; i++)
        {
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%d", sequence[i]->position_in_lib[1]);
            if (i != count - 1)
            {
                string_index += snprintf(&string_out[string_index], string_buffer - string_index, "-");
            }
        }
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, " %f ", structure_energy);
        if (GLOBAL_PERFORM_HAIRPIN_CHECK)
        {
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%f", WC_rmsd);
        }
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\n");
    }
    else
    {
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%s ", GLOBAL_INPUT_SEQUENCE);
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#INDEX ");
        for (int i = 0; i < count; i++)
        {
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%d", sequence[i]->position_in_lib[1]);
            if (i != count - 1)
            {
                string_index += snprintf(&string_out[string_index], string_buffer - string_index, "-");
            }
        }
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\t");
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#ENERGY ");
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%f\t", structure_energy);
        if (GLOBAL_PERFORM_HAIRPIN_CHECK)
        {
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#WC-RMSD ");
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%f", WC_rmsd);
        }
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\n");
    }
    return string_index;
}

char *RNA_data_array::to_string()
{
    int idx_offset;

    if (string_initialized)
    {
        string_index = 0;
    }
    else
    {
        initialize_string();
    }

    string_index = out_string_header();

    if (string_print_coordinates == false)
        return string_out;

    string_index = sequence[0]->to_string(string_out, string_buffer, string_index);
    idx_offset = sequence[0]->count;

    for (int i = 1; i < count; i++)
    {
        string_index = sequence[i]->to_string_offset(1, i, string_out, string_buffer, string_index, &idx_offset);
    }

    string_index += snprintf(&string_out[string_index], string_buffer - string_index, "ENDMDL\n");

    return string_out;
}

int *RNA_data_array::get_index()
{
    int *ar = (int *)malloc(sizeof(int) * count);
    for (int i = 0; i < count; i++)
    {
        ar[i] = sequence[i]->position_in_lib[1];
    }
    return ar;
}
