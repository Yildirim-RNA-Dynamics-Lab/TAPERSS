#include "RNAData.hpp"
#include "RNADataArrayInternalLoop.hpp"

/**
 * @brief Construct a new rna data::rna data object
 * 
 * @param L Dimer Library Array
 * @param i Index of DimerLib (which DNT will be used)
 * @param j Index of model in library
 * @param WC Is this RNA_data for watson crick data
 */
RNAData::RNAData()
{

}

void RNAData::initialize(DimerLibArray &L, int idx, int idx_L, gsl_block *MemBlock, size_t *offset_matrix, uint16_t *ArrayMemBlock, size_t *offset_array)
{
    int rel_idx1 = 0;
    int rel_idx2 = 0;
    int rel_idx3 = 0;
    int rel_idx4 = 0;

    atom_data = (L[idx]->atom_data);
    memcpy(name, L[idx]->name, sizeof(char) * 3);
    energy = L[idx]->energy[idx_L];
    count = L[idx]->atom_data->count;
    _flag = &(L.Flags[idx][idx_L]);
    position_in_lib[0] = idx;
    position_in_lib[1] = idx_L;
    position_max = L[idx]->count - 1;
    ResBoundaries[0] = 0;
    ResBoundaries[1] = atom_data->count_per_res[0];
    ResBoundaries[2] = atom_data->count_per_res[0];
    ResBoundaries[3] = atom_data->count;

    id = rna_dat_tracker++;

    //printf("RNAData at %d,%d with ID: %lu Allocated\n", idx, idx_L, id);

    COM_Radii[0] = L[idx]->radii[0][idx_L];
    COM_Radii[1] = L[idx]->radii[1][idx_L];
    
    data_matrix = gsl_matrix_alloc_from_block(MemBlock, *offset_matrix, L[idx]->data_matrices[idx_L]->size1, L[idx]->data_matrices[idx_L]->size2, MATRIX_DIMENSION2);
    gsl_matrix_memcpy(data_matrix, L[idx]->data_matrices[idx_L]);

    *offset_matrix += L[idx]->data_matrices[idx_L]->size1 * L[idx]->data_matrices[idx_L]->size2;

    count_per_sub[0] = get_target(name[0], &target1);

    submatrix_rows[0] = &ArrayMemBlock[*offset_array];
    *offset_array += count_per_sub[0];

    count_per_sub[1] = get_target(name[1], &target2);

    submatrix_rows[1] = &ArrayMemBlock[*offset_array];
    *offset_array += count_per_sub[1];


    //printf("IN RNAData: Array Offset: %lu\n", *offset_array);

    for (unsigned int j = 0; j < count; j++)
    {
        for (int k = 0; k < count_per_sub[0]; k++)
        {
            if ((target1[k] == atom_data->atom_ids[j]) && (atom_data->dnt_pos[j] - 1) == 0)
            {
                //gsl_matrix_set(submatrices[0], rel_idx1, 0, gsl_matrix_get(data_matrix, j, 0));
                //gsl_matrix_set(submatrices[0], rel_idx1, 1, gsl_matrix_get(data_matrix, j, 1));
                //gsl_matrix_set(submatrices[0], rel_idx1, 2, gsl_matrix_get(data_matrix, j, 2));
                //printf("lib: %d, j = %d, k = %d\n", idx, j, target1[k]);
                submatrix_rows[0][k] = j;
                rel_idx1++;
            }
        }
        for(int k = 0; k < count_per_sub[1]; k++)
        {
            if ((target2[k] == atom_data->atom_ids[j]) && (atom_data->dnt_pos[j] - 1) == 1)
            {
                //gsl_matrix_set(submatrices[1], rel_idx2, 0, gsl_matrix_get(data_matrix, j, 0));
                //gsl_matrix_set(submatrices[1], rel_idx2, 1, gsl_matrix_get(data_matrix, j, 1));
                //gsl_matrix_set(submatrices[1], rel_idx2, 2, gsl_matrix_get(data_matrix, j, 2));
                submatrix_rows[1][k] = j;
                rel_idx2++;
            }
        }
    }

    count_per_WC_sub[0] = get_WC_target(name[0], &WC_target1);
    WC_submatrix_rows[0] = &ArrayMemBlock[*offset_array];
    *offset_array += count_per_WC_sub[0];

    count_per_WC_sub[1] = get_WC_target(name[1], &WC_target2);
    WC_submatrix_rows[1] = &ArrayMemBlock[*offset_array];
    *offset_array += count_per_WC_sub[1];

    //printf("IN RNAData: Array Offset: %lu\n", *offset_array);

    rel_idx1 = 0;
    rel_idx2 = 0;

    for (unsigned int j = 0; j < count; j++)
    {
        for (unsigned int k = 0; k < count_per_WC_sub[0]; k++)
        {
            if ((WC_target1[k] == atom_data->atom_ids[j]) && (atom_data->dnt_pos[j] - 1) == 0)
            {
                WC_submatrix_rows[0][k] = j;
                rel_idx1++;
            }
        }
        for (unsigned int k = 0; k < count_per_WC_sub[1]; k++)
        {
            if ((WC_target2[k] == atom_data->atom_ids[j]) && (atom_data->dnt_pos[j] - 1) == 1)
            {
                WC_submatrix_rows[1][k] = j;
                rel_idx2++;
            }
        }
    }


    //printf("IN RNAData: BFORE STERIC Array Offset: %lu\n", *offset_array);

    StericIndices[0] = &ArrayMemBlock[*offset_array];

    rel_idx1 = 0;
    rel_idx2 = 0;

    for(int i = 0; i < L[idx]->atom_data->count; i++)
    {
        if(L[idx]->atom_data->dnt_pos[i] == 2)
        {
            continue;
        }
        if( L[idx]->atom_data->atom_ids[i] == P   || L[idx]->atom_data->atom_ids[i] == OP1 || 
            L[idx]->atom_data->atom_ids[i] == OP2 || L[idx]->atom_data->atom_ids[i] == C5p || 
            L[idx]->atom_data->atom_ids[i] == O5p)
        {
            StericIndices[0][rel_idx1++] = i;
        }
    }

    count_per_Steric[0] = rel_idx1;
    *offset_array += rel_idx1;
    StericIndices[1] = &ArrayMemBlock[*offset_array];


    for(int i = 0; i < L[idx]->atom_data->count; i++)
    {
        if(L[idx]->atom_data->dnt_pos[i] == 1)
        {
            continue;
        }
        if( L[idx]->atom_data->atom_ids[i] == P   || L[idx]->atom_data->atom_ids[i] == OP1 || 
            L[idx]->atom_data->atom_ids[i] == OP2 || L[idx]->atom_data->atom_ids[i] == C5p || 
            L[idx]->atom_data->atom_ids[i] == O5p)
        {
            StericIndices[1][rel_idx2++] = i;
        }
    }

    count_per_Steric[1] = rel_idx2;
    *offset_array += rel_idx2;

    EnergyIndices[0] = &ArrayMemBlock[*offset_array];
    //printf("IN RNAData: Array Offset: %lu\n", *offset_array);

    rel_idx1 = 0;
    rel_idx2 = 0;

    for(int i = 0; i < L[idx]->atom_data->count; i++)
    {
        if( L[idx]->atom_data->charges[i] != atom_charge::POSITIVE)
        {
            continue;
        }
        else if(L[idx]->atom_data->dnt_pos[i] == 1)
        {
            EnergyIndices[0][rel_idx1++] = i;
        }
    }

    count_per_Energy[0] = rel_idx1;
    *offset_array += rel_idx1;
    EnergyIndices[1] = &ArrayMemBlock[*offset_array];

    //printf("IN RNAData: Array Offset: %lu\n", *offset_array);

    for(int i = 0; i < L[idx]->atom_data->count; i++)
    {
        if( L[idx]->atom_data->charges[i] != atom_charge::NEGATIVE)
        {
            continue;
        }
        else if(L[idx]->atom_data->dnt_pos[i] == 1)
        {
            EnergyIndices[1][rel_idx2++] = i;
        }
    }
    count_per_Energy[1] = rel_idx2;
    *offset_array += rel_idx2;
    EnergyIndices[2] = &ArrayMemBlock[*offset_array];

    for(int i = 0; i < L[idx]->atom_data->count; i++)
    {
        if( L[idx]->atom_data->charges[i] != atom_charge::POSITIVE)
        {
            continue;
        }
        else if(L[idx]->atom_data->dnt_pos[i] == 2)
        {
            EnergyIndices[2][rel_idx3++] = i;
        }
    }

    count_per_Energy[2] = rel_idx3;
    *offset_array += rel_idx3;
    EnergyIndices[3] = &ArrayMemBlock[*offset_array];

    //printf("IN RNAData: Array Offset: %lu\n", *offset_array);

    for(int i = 0; i < L[idx]->atom_data->count; i++)
    {
        if( L[idx]->atom_data->charges[i] != atom_charge::NEGATIVE)
        {
            continue;
        }
        else if(L[idx]->atom_data->dnt_pos[i] == 2)
        {
            EnergyIndices[3][rel_idx4++] = i;
        }
    }
    count_per_Energy[3] = rel_idx4;
    *offset_array += rel_idx4;
    //printf("IN RNAData: Array Offset: %lu\n", *offset_array);
}

void RNAData::overwrite(DimerLibArray &L, int i, int j)
{
    gsl_matrix_memcpy(data_matrix, L[i]->data_matrices[j]);
    energy = L[i]->energy[j];
    //memcpy(name, L[i]->name, sizeof(char) * 3);
    //count = L[i]->atom_data->count;
    _flag = &(L.Flags[i][j]);
    position_in_lib[0] = i;
    position_in_lib[1] = j;    
    //COM_Radii[0] = L[i]->radii[0][j];
    //COM_Radii[1] = L[i]->radii[1][j];
    //is_for_WC ? update_WC_submatrices() : update_submatrices();
}

RNAData::~RNAData()
{
    // printf("Deallocationg %ld\n", id);
    gsl_matrix_free(data_matrix);
    /*if (is_for_WC)
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
    }*/
}

void RNAData::destroy()
{
    gsl_matrix_free(data_matrix);
}

/** Deprecated **/
//void RNAData::make_submatrices()
//{
    /*submatrices = (gsl_matrix **)malloc(sizeof(gsl_matrix *) * 2);

    for (int i = 0; i < 2; i++)
    {
        count_per_sub[i] = get_target(name[i], &target);
        submatrices[i] = gsl_matrix_alloc(count_per_sub[i], 3);
        submatrix_rows[i] = (uint16_t *)malloc(sizeof(size_t) * count_per_sub[i]);

        int rel_idx = 0;

        for (unsigned int j = 0; j < count; j++)
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
    }*/
//}

/*void RNAData::update_submatrices()
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
}*/

size_t get_target(char res, const atom_id** dest)
{
    if (res == 'A')
    {
        *dest = targetA;        
        return (sizeof(targetA) / sizeof(atom_id));
    }
    else if (res == 'C')
    {
        *dest = targetC;
        return (sizeof(targetC) / sizeof(atom_id));
    }
    else if (res == 'G')
    {
        *dest = targetG;        
        return (sizeof(targetG) / sizeof(atom_id));
    }
    else if (res == 'U')
    {
        *dest = targetU;
        return (sizeof(targetU) / sizeof(atom_id));
    }
    else
    {
        printf("RNA_data: input file has problems, please check at %d\n", __LINE__);
        exit(1);
    }
}

/*gsl_matrix *RNAData::get_target_matrix(int res)
{
    //return submatrices[res];
}*/

/*gsl_matrix *RNAData::get_target_matrix_copy(int res)
{
    gsl_matrix *rtn = gsl_matrix_alloc(count_per_sub[res], 3);
    gsl_matrix_memcpy(rtn, submatrices[res]);
    return rtn;
}*/

/*void RNAData::make_WC_submatrices(bool first_run)
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
    else
    {
        // printf("first run for %ld\n", id);
    }
    
    WC_submatrices = (gsl_matrix **)malloc(sizeof(gsl_matrix *) * 2);

    for (int i = 0; i < 2; i++)
    {
        count_per_WC_sub[i] = get_WC_target(name[i], &WC_target);
        WC_submatrices[i] = gsl_matrix_alloc(count_per_WC_sub[i], 3);
        WC_submatrix_rows[i] = (uint16_t *)malloc(sizeof(size_t) * count_per_WC_sub[i]);

        int rel_idx = 0;

        for (unsigned int j = 0; j < count; j++)
        {
            for (unsigned int k = 0; k < count_per_WC_sub[i]; k++)
            {
                if ((WC_target[k] == atom_data->atom_ids[j]) && (atom_data->dnt_pos[j] - 1) == i)
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
    }
}*/

/*void RNAData::update_WC_submatrices()
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
}*/

size_t get_WC_target(char res, const atom_id** dest)
{
    if (res == 'A')
    {
        *dest = WCtargetA;        
        return (sizeof(WCtargetA) / sizeof(atom_id));
    }
    else if (res == 'C')
    {
        *dest = WCtargetC;
        return (sizeof(WCtargetC) / sizeof(atom_id));
    }
    else if (res == 'G')
    {
        *dest = WCtargetG;
        return (sizeof(WCtargetG) / sizeof(atom_id));
    }
    else if (res == 'U')
    {
        *dest = WCtargetU;
        return (sizeof(WCtargetU) / sizeof(atom_id));
    }
    else
    {
        printf("RNA_data: input file has problems, please check at %d\n", __LINE__);
        exit(1);
    }
}

/*gsl_matrix *RNAData::get_WC_target_matrix(int res)
{
    //printf("size of WC_submatrix = %ld\n", WC_submatrices[res]->size1);
    return WC_submatrices[res];
}*/

/*gsl_matrix *RNAData::get_WC_target_matrix_copy(int res)
{
    gsl_matrix *rtn = gsl_matrix_alloc(count_per_WC_sub[res], 3);
    gsl_matrix_memcpy(rtn, WC_submatrices[res]);
    return rtn;
}*/

int RNAData::get_residue_COM_index(int res)
{
    return (count + res);
}

void RNAData::print()
{
    for (unsigned int i = 0; i < data_matrix->size1; i++)
    {
        atom_data->print_at(i);
        for (unsigned int j = 0; j < data_matrix->size2; j++)
            printf("%8.3f", gsl_matrix_get(data_matrix, i, j));
        putchar('\n');
    }
}

void RNAData::print_target(int res)
{
    for (int i = 0; i < count_per_sub[res]; i++)
    {
        atom_data->print_at(submatrix_rows[res][i]);
        for (unsigned int j = 0; j < data_matrix->size2; j++)
            printf("%8.3f", gsl_matrix_get(data_matrix, submatrix_rows[res][i], j));
        putchar('\n');
    }
}

void RNAData::print_WC_target(int res)
{
    for (unsigned int i = 0; i < count_per_WC_sub[res]; i++)
    {
        atom_data->print_at(WC_submatrix_rows[res][i]);
        for (unsigned int j = 0; j < data_matrix->size2; j++)
            printf("%8.3f", gsl_matrix_get(data_matrix, WC_submatrix_rows[res][i], j));
        putchar('\n');
    }
}

void RNAData::print(int res)
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

int RNAData::to_string(char *s, int buffer_size, int string_index)
{
    for (unsigned int i = 0; i < count; i++)
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

void RNAData::print_offset(int res, int position)
{
    for (unsigned int i = 0; i < count; i++)
    {
        if (atom_data->dnt_pos[i] != (res + 1))
            continue;
        printf("%-4s %-4d %-4c %-4d ", atom_data->name[i], atom_data->index[i], atom_data->residue[i], atom_data->dnt_pos[i] + position);
        for (unsigned int j = 0; j < data_matrix->size2; j++)
            printf("%-3.2f ", gsl_matrix_get(data_matrix, i, j));
        putchar('\n');
    }
}

int RNAData::to_string_offset(int res, int position, char *s, int buffer_size, int string_index, int *idx_offset)
{
    for (unsigned int i = 0; i < count; i++)
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

void RNADataArray::initialize(size_t size, DimerLibArray& Library)
{
    int PLAC = Library.PositiveAtomCount;
    int NLAC = Library.NegativeAtomCount;
    //printf("PLAC = %d, NLAC = %d\n", PLAC, NLAC);
    iterator_max = size - 1;
    sequence = (RNAData *)malloc(sizeof(RNAData) * size);
    uint64_t matrix_memsize = 0;
    uint64_t array_memsize = 0;
    for(size_t i = 0; i < size; i++)
    {
        matrix_memsize   += calculate_matrix_memory_needed(Library, i);
        array_memsize    += calculate_index_array_memory_needed(Library, i);
    }
    //printf("Array Memsize : %lu\nMatrix Memsize: %lu\n", array_memsize / sizeof(uint16_t), matrix_memsize);
    matrix_memsize += (PLAC * NLAC); // Interaction Table Matrix
    //printf("Array Memsize : %lu\nMatrix Memsize: %lu\n", array_memsize / sizeof(uint16_t), matrix_memsize);
    MatrixMemBlock = gsl_block_alloc(matrix_memsize); //Allocate Block of memory to pass to individual RNAData's. This is to maintain contiguous memory.
    ArrayMemBlock = (uint16_t *)malloc(array_memsize); // Same as ^.

    size_t MatrixOffset = 0;
    size_t ArrayOffset = 0;

    for(size_t i = 0; i < size; i++)
    {
        //printf("BEFORE: Array Offset : %lu\nMatrix Offset: %lu\n", ArrayOffset, MatrixOffset);
        sequence[i].initialize(Library, i, 0, MatrixMemBlock, &MatrixOffset, ArrayMemBlock, &ArrayOffset);
        
    }

    //printf("AFTER: Array Offset : %lu\nMatrix Offset: %lu\n", ArrayOffset, MatrixOffset);
    *sequence[0]._flag = USED;
    iterator = 0;  //Start with first already DNT included
    count = 1;
    structure_energy = 0;
    COMS = gsl_matrix_alloc_from_block(MatrixMemBlock, MatrixOffset, 2 * size, MATRIX_DIMENSION2, MATRIX_DIMENSION2);
    MatrixOffset += (2 * size * MATRIX_DIMENSION2);
    //printf("AFTER2: Array Offset : %lu\nMatrix Offset: %lu\n", ArrayOffset, MatrixOffset);
    InteractionTable = gsl_matrix_alloc_from_block(MatrixMemBlock, MatrixOffset, PLAC, NLAC, NLAC);
    gsl_matrix_set_zero(InteractionTable); // Set all to 0
    PositiveAtomMap = Library.PositiveAtomMap;
    NegativeAtomMap = Library.NegativeAtomMap;
    TableRowCount = PLAC;
    TableColCount = NLAC;
    LAC = Library.LargestAtomCount;
    Radii = (double *)malloc(2 * size * sizeof(double));
    PassedCOMCheck = (bool*)malloc((size + 1) * sizeof(bool)); // in normal building of structure only first DNT will have valid 1st residue.
                                                               // in all other DNTs, only second residue is needed.
    overwrite_initialize(0,0, Library);
}

void RNADataArray::overwrite(size_t LibIdx, size_t IdxInLib, DimerLibArray &Library)
{
    sequence[LibIdx].overwrite(Library, LibIdx, IdxInLib);
    //gsl_matrix_row_copy(COMS, LibIdx * 2, Library[LibIdx]->data_matrices[IdxInLib], Library[LibIdx]->atom_data->count + 0);
    //gsl_matrix_row_copy(COMS, (LibIdx * 2) + 1, Library[LibIdx]->data_matrices[IdxInLib], Library[LibIdx]->atom_data->count + 1);
    //gsl_matrix_row_copy(COMS, LibIdx * 2, Library[LibIdx]->data_matrices[IdxInLib], Library[LibIdx]->atom_data->count + 0);
    //gsl_matrix_row_copy(COMS, (LibIdx * 2) + 1, Library[LibIdx]->data_matrices[IdxInLib], Library[LibIdx]->atom_data->count + 1);
    Radii[LibIdx * 2] = Library[LibIdx]->radii[0][IdxInLib];
    Radii[LibIdx * 2 + 1] = Library[LibIdx]->radii[1][IdxInLib];
}

void RNADataArray::overwrite_initialize(size_t LibIdx, size_t IdxInLib, DimerLibArray &Library)
{
    sequence[LibIdx].overwrite(Library, LibIdx, IdxInLib);
    gsl_matrix_row_copy(COMS, LibIdx * 2, Library[LibIdx]->data_matrices[IdxInLib], Library[LibIdx]->atom_data->count + 0);
    gsl_matrix_row_copy(COMS, (LibIdx * 2) + 1, Library[LibIdx]->data_matrices[IdxInLib], Library[LibIdx]->atom_data->count + 1);
    //gsl_matrix_row_copy(COMS, LibIdx * 2, Library[LibIdx]->data_matrices[IdxInLib], Library[LibIdx]->atom_data->count + 0);
    //gsl_matrix_row_copy(COMS, (LibIdx * 2) + 1, Library[LibIdx]->data_matrices[IdxInLib], Library[LibIdx]->atom_data->count + 1);
    Radii[LibIdx * 2] = Library[LibIdx]->radii[0][IdxInLib];
    Radii[LibIdx * 2 + 1] = Library[LibIdx]->radii[1][IdxInLib];
}

RNADataArray::RNADataArray(){}

RNADataArray::~RNADataArray()
{
    // printf("iterator is @ %d\n", iterator);
    for(int i = 0; i < count; i++)
    {
        sequence[i].destroy();
    }
    free(sequence);
    gsl_matrix_free(COMS);
    free(Radii);
    free(PassedCOMCheck);
    if (string_initialized)
        free(string_out);
    gsl_matrix_free(InteractionTable);
    free(PositiveAtomMap);
    free(NegativeAtomMap);
    gsl_block_free(MatrixMemBlock);
    free(ArrayMemBlock);
}

uint_fast64_t RNADataArray::calculate_matrix_memory_needed(DimerLibArray& L, int idx)
{
    uint_fast64_t memsize = 0;
    size_t data_mat_nrows = L[idx]->data_matrices[0]->size1;

    memsize += data_mat_nrows * MATRIX_DIMENSION2;
    memsize += 2 * MATRIX_DIMENSION2; //For COM Matrix

    /* Size for Submatrices */
    //memsize += get_target(L[idx]->name[0], &dummy) * MATRIX_DIMENSION2; //Submatrix for Resid 1
    //memsize += get_target(L[idx]->name[1], &dummy) * MATRIX_DIMENSION2; //Submatrix for Resid 2

    /* Size for WCSubmatrices 
    memsize += get_target(L[idx]->name[0], &dummy) * MATRIX_DIMENSION2;
    memsize += get_target(L[idx]->name[1], &dummy) * MATRIX_DIMENSION2;
    */

    return memsize;
}

uint_fast64_t RNADataArray::calculate_index_array_memory_needed(DimerLibArray& L, int idx)
{
    uint_fast64_t memsize = 0;
    const atom_id *dummy = NULL;

    /* Size for Submatrices Rows*/
    memsize += get_target(L[idx]->name[0], &dummy);
    memsize += get_target(L[idx]->name[1], &dummy);

    //printf("IN CALC: Array Offset: %lu\n", memsize);

    /* Size for WCSubmatrices Rows */
    memsize += get_WC_target(L[idx]->name[0], &dummy);
    memsize += get_WC_target(L[idx]->name[1], &dummy);

    //printf("IN CALC: Array Offset: %lu\n", memsize);

    /* Size for Steric Rows */
    //int counter = 0;
    for(int i = 0; i < L[idx]->atom_data->count; i++)
    {
        if( L[idx]->atom_data->atom_ids[i] == P   || L[idx]->atom_data->atom_ids[i] == OP1 || 
            L[idx]->atom_data->atom_ids[i] == OP2 || L[idx]->atom_data->atom_ids[i] == C5p ||
            L[idx]->atom_data->atom_ids[i] == O5p)
        {
            memsize++;
            //printf("CALC STERIC: %d\n", counter++);
        }
    }

    //printf("IN CALC: Array Offset: %lu\n", memsize);

    /* Size for Charged Rows */
    for(int i = 0; i < L[idx]->atom_data->count; i++)
    {
        if( L[idx]->atom_data->charges[i] != atom_charge::NEUTRAL)
        {
            memsize++;
        }
    }

    //printf("IN CALC: Array Offset: %lu\n", memsize);
    memsize *= sizeof(uint16_t);

    return memsize;
}

RNAData *RNADataArray::operator[](int i)
{
    return &sequence[i];
}

/* Unused and should not be used for any reason (will remove at some point)*/
void RNADataArray::add_copy(RNAData *A)
{
    sequence[++iterator] = *A;
    count++;
}

void RNADataArray::add_move(RNAData *A)
{
    iterator++;
    gsl_matrix_row_copy(COMS, iterator * 2, A->data_matrix, A->get_residue_COM_index(0));
    gsl_matrix_row_copy(COMS, iterator * 2 + 1, A->data_matrix, A->get_residue_COM_index(1));
    Radii[iterator * 2] = A->COM_Radii[0];
    Radii[iterator * 2 + 1] = A->COM_Radii[1];
    sequence[iterator] = *A;
    count++;
    *A->_flag = USED;
}


RNAData *RNADataArray::current()
{
    return &sequence[iterator];
}

bool RNADataArray::is_complete()
{
    return (iterator == iterator_max);
}

void RNADataArray::keep()
{
    iterator++;
    gsl_matrix_row_copy(COMS, iterator * 2, sequence[iterator].data_matrix, sequence[iterator].get_residue_COM_index(0));
    gsl_matrix_row_copy(COMS, (iterator * 2) + 1, sequence[iterator].data_matrix, sequence[iterator].get_residue_COM_index(1));
    count++;
    *sequence[iterator]._flag = USED;
}

void RNADataArray::rollback()
{
    //DEBUG(printf("attach: %ld deleted @ rollback @ %d\n", sequence[iterator]->id, iterator));
    //delete sequence[iterator];
    iterator--;
    count--;
}

void RNADataArray::safe_rollback() // unused
{
    /*if (sequence[iterator]->position_in_lib[1] == sequence[iterator]->position_max)
        ;
    else
        delete sequence[iterator];*/
    iterator--;
    count--;
}

void RNADataArray::rollback_by(int amount)
{
    //print_index(0);
    /*for (int i = iterator; i > (iterator - (amount + 1)); i--)
    {
        DEBUG(printf("attach: %ld deleted @ rollback_by @ %d\n", sequence[i]->id, i));
        delete sequence[i];
    }*/
    iterator -= (amount + 1);
    count -= (amount + 1);
    //printf("Count: %d, iterator: %d\n", count, iterator);
    //printf("moving to pos: %d\n", iterator);
    //print_index(0);
}

bool RNADataArray::is_empty()
{
    if (count == 0)
        return true;
    return false;
}

void RNADataArray::update_WC_rmsd(float rmsd_val)
{
    WC_rmsd1_6 = rmsd_val;
}

void RNADataArray::update_energy()
{
    float energy_ = 0.0;
    int Interactions = 0;


    for (int i = 0; i < count; i++)
    {
        energy_ += sequence[i].energy;
    }

    for(int i = count - 1; i > 1; i--)
    {
        Interactions += FindInteraction(0, 0, &sequence[0], 1, i, &sequence[i],
                        InteractionTable, PositiveAtomMap, NegativeAtomMap, LAC);
        Interactions += FindInteraction(1, 0, &sequence[0], 1, i, &sequence[i],
                        InteractionTable, PositiveAtomMap, NegativeAtomMap, LAC);
    }

    for (int i = 1; i <= count - 3; i++)
    {
        for (int j = count - 1; j > i + 1; j--)
        {
            Interactions += FindInteraction(1, i, &sequence[i], 1, j, &sequence[j],
                        InteractionTable, PositiveAtomMap, NegativeAtomMap, LAC);
        }
    }
    
    if(Interactions > 1)
    {
        Interactions = HK_GetMaxMatching(InteractionTable);
    }
    energy_ -= Interactions;
    structure_energy = energy_;

    gsl_matrix_set_zero(InteractionTable);
}

int FindInteraction(int A_ResId, int A_Idx, RNAData *AData, int B_ResId, int B_Idx, RNAData *BData, gsl_matrix *AdjMatrix, int *PAdjMap, int*NAdjMap, size_t DIM2)
{
    gsl_matrix *A = AData->data_matrix;
    gsl_matrix *B = BData->data_matrix;
    uint16_t *A_P_Rows = AData->EnergyIndices[A_ResId * 2];
    uint16_t *A_N_Rows = AData->EnergyIndices[A_ResId * 2 + 1];
    uint16_t *B_P_Rows = BData->EnergyIndices[B_ResId * 2];
    uint16_t *B_N_Rows = BData->EnergyIndices[B_ResId * 2 + 1];
    size_t A_P_NRow = AData->count_per_Energy[A_ResId * 2];
    size_t A_N_NRow = AData->count_per_Energy[A_ResId * 2 + 1];
    size_t B_P_NRow = BData->count_per_Energy[B_ResId * 2];
    size_t B_N_NRow = BData->count_per_Energy[B_ResId * 2 + 1];
    size_t IDX1, IDX2;
    int Interactions = 0;
    for(size_t i = 0; i < A_P_NRow; i++)
    {
        for(size_t j = 0; j < B_N_NRow; j++)
        {
            if(distance_mat2mat(A, A_P_Rows[i], B, B_N_Rows[j]) < INTERACTION_DISTANCE)
            {
                IDX1 = IDX_FLAT2D(A_Idx, A_P_Rows[i], DIM2);
                IDX2 = IDX_FLAT2D(B_Idx, B_N_Rows[j], DIM2);
                //printf("(%d,%d):(%d,%d) == %d:%d\n", A_Idx, A_P_Rows[i], B_Idx, B_N_Rows[j],PAdjMap[IDX1], NAdjMap[IDX2]);
                gsl_matrix_set(AdjMatrix, PAdjMap[IDX1], NAdjMap[IDX2], 1);
                Interactions++;
            }
        }
    }
    for(size_t i = 0; i < A_N_NRow; i++)
    {
        for(size_t j = 0; j < B_P_NRow; j++)
        {
            if(distance_mat2mat(A, A_N_Rows[i], B, B_P_Rows[j]) < INTERACTION_DISTANCE)
            {
                IDX1 = IDX_FLAT2D(A_Idx, A_N_Rows[i], DIM2);
                IDX2 = IDX_FLAT2D(B_Idx, B_P_Rows[j], DIM2);
                gsl_matrix_set(AdjMatrix, PAdjMap[IDX2], NAdjMap[IDX1], 1);
                //printf("(%d,%d):(%d,%d) == %d:%d\n", A_Idx, A_N_Rows[i], B_Idx, B_P_Rows[j],PAdjMap[IDX2], NAdjMap[IDX1]);
                Interactions++;
            }
        }
    }
    return Interactions;
}

void RNADataArray::printall()
{
    sequence[0].print();
    for (int i = 1; i < count; i++)
    {
        sequence[i].print_offset(1, i);
    }
}

void RNADataArray::initialize_string()
{
    get_atom_sum();
    // printf("atomsum = %d\n", atom_sum);

    string_buffer = 54 * atom_sum;
    string_out = (char *)malloc(sizeof(char) * string_buffer);
    string_index = 0;
    model_count = 0;

    string_initialized = true;
}

int RNADataArray::get_atom_sum()
{
    if (string_initialized)
        return atom_sum;
    atom_sum = 0;
    for (int i = 0; i < count; i++)
    {
        atom_sum += sequence[i].count;
    }
    return atom_sum;
}
int RNADataArray::out_string_header_coord()
{
    string_index += snprintf(&string_out[string_index], string_buffer - string_index, "MODEL %d\n", ++model_count);
    string_index += snprintf(&string_out[string_index], string_buffer - string_index, "REMARK ");
    string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%s ", GLOBAL_INPUT_SEQUENCE);
    for (int i = 0; i < count; i++)
    {
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%d", sequence[i].position_in_lib[1]);
        if (i != count - 1)
        {
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "-");
        }
    }
    string_index += snprintf(&string_out[string_index], string_buffer - string_index, " %f ", structure_energy);
    if (GLOBAL_PERFORM_STRUCTCHECK == HAIRPIN)
    {
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%f", WC_rmsd1_6);
    }
    string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\n");
    return string_index;
}
int RNADataArray::out_string_header()
{
    if (GLOBAL_WRITE_COORDINATES)
    {
        string_index = out_string_header_coord();
    }
    else
    {
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%s ", GLOBAL_INPUT_SEQUENCE);
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#INDEX ");
        for (int i = 0; i < count; i++)
        {
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%d", sequence[i].position_in_lib[1]);
            if (i != count - 1)
            {
                string_index += snprintf(&string_out[string_index], string_buffer - string_index, "-");
            }
        }
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\t");
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#ENERGY ");
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%f\t", structure_energy);
        if (GLOBAL_PERFORM_STRUCTCHECK == HAIRPIN)
        {
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#WC-RMSD ");
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%f", WC_rmsd1_6);
        }
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\n");
    }
    return string_index;
}

char *RNADataArray::to_string()
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

    string_index = sequence[0].to_string(string_out, string_buffer, string_index);
    idx_offset = sequence[0].count;

    for (int i = 1; i < count; i++)
    {
        string_index = sequence[i].to_string_offset(1, i, string_out, string_buffer, string_index, &idx_offset);
    }

    string_index += snprintf(&string_out[string_index], string_buffer - string_index, "ENDMDL\n");

    return string_out;
}

int* RNADataArray::get_index()
{
    int *ar = (int *)malloc(sizeof(int) * count);
    for (int i = 0; i < count; i++)
    {
        ar[i] = sequence[i].position_in_lib[1];
    }
    return ar;
}

void RNADataArray::print_index(int offset)
{
    for (int i = 0; i < count + offset; i++)
    {
        if(i == (count - 1 + offset))
        {
            printf("%d\n", sequence[i].position_in_lib[1]);
            break;   
        }
        printf("%d-", sequence[i].position_in_lib[1]);
    }
}
