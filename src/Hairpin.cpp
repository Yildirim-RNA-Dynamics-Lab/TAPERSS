#include "Hairpin.hpp"

gsl_block  *WC_memblock;
gsl_matrix **WC_reference_matrices;
gsl_matrix **WC_structure_matrices;
gsl_matrix **WC_work;
size_t     N_matrices;

void overwrite_WC_submatrix_gsl(gsl_matrix *A, gsl_matrix *B, gsl_matrix *WC_matrix)
{
    int iterator = 0;
    for (unsigned int i = 0; i < A->size1; i++)
    {
        gsl_matrix_set(WC_matrix, iterator, 0, gsl_matrix_get(A, i, 0));
        gsl_matrix_set(WC_matrix, iterator, 1, gsl_matrix_get(A, i, 1));
        gsl_matrix_set(WC_matrix, iterator, 2, gsl_matrix_get(A, i, 2));
        printf("A: Size:%lu %f, %f, %f\n", A->size1, gsl_matrix_get(WC_matrix, iterator, 0), gsl_matrix_get(WC_matrix, iterator, 1), gsl_matrix_get(WC_matrix, iterator, 2));
        iterator++;
    }

    for (unsigned int i = 0; i < B->size1; i++)
    {
        printf("Iterator = %d\n", iterator);
        printf("Size of WC: %lu\n", WC_matrix->size1);
        gsl_matrix_set(WC_matrix, iterator, 0, gsl_matrix_get(B, i, 0));
        gsl_matrix_set(WC_matrix, iterator, 1, gsl_matrix_get(B, i, 1));
        gsl_matrix_set(WC_matrix, iterator, 2, gsl_matrix_get(B, i, 2));
        printf("B: Size:%lu %f, %f, %f\n", B->size1,gsl_matrix_get(WC_matrix, iterator, 0), gsl_matrix_get(WC_matrix, iterator, 1), gsl_matrix_get(WC_matrix, iterator, 2));
        iterator++;
    }
}

/*gsl_matrix *make_WC_submatrix(RNAData *A, RNAData *B)
{
    gsl_matrix *A_target;
    gsl_matrix *B_target;
    gsl_matrix *WC_matrix;
    
    if (A->id != B->id)
    {
        A->make_WC_submatrices();
        B->make_WC_submatrices();
    }
    else
    {
        A->make_WC_submatrices();
    }

    A_target = A->get_WC_target_matrix(0);
    B_target = B->get_WC_target_matrix(1);
    WC_matrix = gsl_matrix_alloc((A_target->size1 + B_target->size1), 3);
    overwrite_WC_submatrix_gsl(A_target, B_target, WC_matrix);
    //printf("A: ID = %ld, B: ID = %ld\n", A->id, B->id);
    //printf("WC MAT SIZE = %ld\n", (A->WC_submatrices[0]->size1));
    return WC_matrix;
}*/

double RMSD_WC_pair_gsl(gsl_matrix *WC_matrix, gsl_matrix *Sequence_matrix)
{
    double COMP[] = {0, 0, 0};
    double COMQ[] = {0, 0, 0};
    gsl_matrix* Work_matrix;
    Work_matrix = kabsch_allocate_work_matrix(WC_matrix);
    kabsch_calculate_rotation_matrix_Nx3fast(WC_matrix, Sequence_matrix, Work_matrix, COMP, COMQ);
    gsl_matrix_free(Work_matrix);
    return rmsd_generic(Sequence_matrix, WC_matrix);
}

double WC_check_pair(int WC_pair_idx)
{
    gsl_matrix *sequence_matrix = WC_structure_matrices[WC_pair_idx];
    gsl_matrix *WC_model_matrix = WC_reference_matrices[WC_pair_idx];
    gsl_matrix *Work_matrix = WC_work[WC_pair_idx];
    double COMP[] = {0, 0, 0}, COMQ[] = {0, 0, 0};
    kabsch_calculate_rotation_matrix_Nx3fast(WC_model_matrix, sequence_matrix, Work_matrix, COMP, COMQ);
    return (rmsd_generic(sequence_matrix, WC_model_matrix));
}

void WC_prepare_structure_matrix(int WC_pair_idx, gsl_matrix *A, uint16_t *A_rows, size_t A_row_count, gsl_matrix *B, uint16_t *B_rows, size_t B_row_count)
{
    int rel_idx = 0;
    for(size_t i = 0; i < A_row_count; i++)
    {
        gsl_matrix_row_copy(WC_structure_matrices[WC_pair_idx], rel_idx, A, A_rows[i]);
        rel_idx++;
    }
    for(size_t i = 0; i < B_row_count; i++)
    {
        gsl_matrix_row_copy(WC_structure_matrices[WC_pair_idx], rel_idx, B, B_rows[i]);
        rel_idx++;
    }
}

void WC_create(DimerLibArray &WC_Library) 
{
    size_t memsize = 0;
    const atom_id *target1 = nullptr;
    const atom_id *target2 = nullptr;
    size_t offset = 0;
    size_t matrixsize = 0;
    size_t target1size = 0;
    size_t target2size = 0;
    N_matrices = WC_Library.count;
    for(uint64_t i = 0; i < WC_Library.count; i++)
    {
        memsize += get_WC_target(WC_Library[i]->name[0], &target1) * MATRIX_DIMENSION2;
        memsize += get_WC_target(WC_Library[i]->name[1], &target2) * MATRIX_DIMENSION2;
    }
    WC_memblock = gsl_block_alloc(memsize * 3); // Times 3 because space is needed for Work matrices AND matrices from built structure
    WC_reference_matrices = (gsl_matrix**)malloc(sizeof(gsl_matrix*) * WC_Library.count);
    WC_structure_matrices = (gsl_matrix**)malloc(sizeof(gsl_matrix*) * WC_Library.count);
    WC_work = (gsl_matrix**)malloc(sizeof(gsl_matrix*) * WC_Library.count);
    for(uint64_t i = 0; i < WC_Library.count; i++)
    {
        matrixsize = 0;
        matrixsize = get_WC_target(WC_Library[i]->name[0], &target1);
        matrixsize += get_WC_target(WC_Library[i]->name[1], &target2);
        WC_reference_matrices[i] = gsl_matrix_alloc_from_block(WC_memblock, offset, matrixsize, MATRIX_DIMENSION2, MATRIX_DIMENSION2);
        offset += matrixsize * MATRIX_DIMENSION2;
        WC_structure_matrices[i] = gsl_matrix_alloc_from_block(WC_memblock, offset, matrixsize, MATRIX_DIMENSION2, MATRIX_DIMENSION2);
        offset += matrixsize * MATRIX_DIMENSION2;
        WC_work[i] = gsl_matrix_alloc_from_block(WC_memblock, offset, MATRIX_DIMENSION2, matrixsize, matrixsize);
        offset += matrixsize * MATRIX_DIMENSION2;

        target1size = get_WC_target(WC_Library[i]->name[0], &target1);
        target2size = get_WC_target(WC_Library[i]->name[1], &target2);
        size_t rel_idx = 0;

        for (uint64_t j = 0; j < WC_Library[i]->atom_data->count; j++)
        {
            for (unsigned int k = 0; k < target1size; k++)
            {
                if ((target1[k] == WC_Library[i]->atom_data->atom_ids[j]) && (uint8_t)(WC_Library[i]->atom_data->dnt_pos[j] - 1) == 0)
                {
                    gsl_matrix_set(WC_reference_matrices[i], rel_idx, 0, gsl_matrix_get(WC_Library[i]->data_matrices[0], j, 0));
                    gsl_matrix_set(WC_reference_matrices[i], rel_idx, 1, gsl_matrix_get(WC_Library[i]->data_matrices[0], j, 1));
                    gsl_matrix_set(WC_reference_matrices[i], rel_idx, 2, gsl_matrix_get(WC_Library[i]->data_matrices[0], j, 2));
                    rel_idx++;
                    // printf("%s:%ld\tGetting Atom %s from row %d\n", name, id, atom_data->name[j], j);
                }
            }
            for (unsigned int k = 0; k < target2size; k++)
            {
                if ((target2[k] == WC_Library[i]->atom_data->atom_ids[j]) && (uint8_t)(WC_Library[i]->atom_data->dnt_pos[j] - 1) == 1)
                {
                    gsl_matrix_set(WC_reference_matrices[i], rel_idx, 0, gsl_matrix_get(WC_Library[i]->data_matrices[0], j, 0));
                    gsl_matrix_set(WC_reference_matrices[i], rel_idx, 1, gsl_matrix_get(WC_Library[i]->data_matrices[0], j, 1));
                    gsl_matrix_set(WC_reference_matrices[i], rel_idx, 2, gsl_matrix_get(WC_Library[i]->data_matrices[0], j, 2));
                    rel_idx++;
                    // printf("%s:%ld\tGetting Atom %s from row %d\n", name, id, atom_data->name[j], j);
                }
            }
        }
    }
}

void WC_destroy()
{
    for(size_t i = 0; i < N_matrices; i++)
    {
        gsl_matrix_free(WC_reference_matrices[i]);
        gsl_matrix_free(WC_structure_matrices[i]);
        gsl_matrix_free(WC_work[i]);
    }
    free(WC_reference_matrices);
    free(WC_structure_matrices);
    free(WC_work);
    gsl_block_free(WC_memblock);
}