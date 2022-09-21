#include "Hairpin.hpp"


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

gsl_matrix *make_WC_submatrix(RNAData *A, RNAData *B)
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
}

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

bool is_WC_pair(RNADataArray &sequence, DimerLibArray &WC_Lib, int i, int j, int WC_pair_idx)
{
    //printf("int i = %d, int j = %d\n", i, j);
    //printf("Library index: %d\n", sequence[i]->position_in_lib[1]);
    //printf("Library index: %d\n", sequence[j]->position_in_lib[1]);
    //RNAData *WC_pair = new RNAData(WC_Lib, WC_pair_idx, 0, true);
    //WC_pair->print_WC_target(0);
    //WC_pair->print_WC_target(1);
    //sequence[i]->print_WC_target(0);
    //sequence[j]->print_WC_target(1);
    gsl_matrix *sequence_matrix;
    gsl_matrix *WC_model_matrix = nullptr;
    gsl_matrix *Work_matrix;
    double _rmsd = 0;
    double COMP[] = {0, 0, 0}, COMQ[] = {0, 0, 0};
    bool is_WC_pair = false;
    sequence_matrix = make_WC_submatrix(sequence[i], sequence[j]);
    //WC_model_matrix = make_WC_submatrix(WC_pair, WC_pair);
    Work_matrix = kabsch_allocate_work_matrix(WC_model_matrix);
    kabsch_calculate_rotation_matrix_Nx3fast(WC_model_matrix, sequence_matrix, Work_matrix, COMP, COMQ);
    _rmsd = rmsd_generic(sequence_matrix, WC_model_matrix);
    //printf("WC RMSD: %f VS %f\n", _rmsd, GLOBAL_WC_RMSD_LIMIT);
    if (_rmsd <= GLOBAL_WC_RMSD_LIMIT)
    {
        is_WC_pair = true;
        //printf("RMSD: %f\n", _rmsd);
    }
    sequence.update_WC_rmsd(_rmsd);
    gsl_matrix_free(Work_matrix);
    gsl_matrix_free(sequence_matrix);
    gsl_matrix_free(WC_model_matrix);
    //delete WC_pair;
    return is_WC_pair;
}
