#include "Hairpin.hpp"

gsl_matrix *make_WC_submatrix(RNA_data *A, RNA_data *B)
{
    gsl_matrix *A_target;
    gsl_matrix *B_target;
    gsl_matrix *WC_matrix;
    int iterator = 0;
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

    // printf("A: ID = %ld, B: ID = %ld\n", A->id, B->id);
    // printf("WC MAT SIZE = %ld\n", (A->WC_submatrices[0]->size1));

    WC_matrix = gsl_matrix_alloc((A_target->size1 + B_target->size1), 3);

    for (unsigned int i = 0; i < A_target->size1; i++)
    {
        gsl_matrix_set(WC_matrix, iterator, 0, gsl_matrix_get(A_target, i, 0));
        gsl_matrix_set(WC_matrix, iterator, 1, gsl_matrix_get(A_target, i, 1));
        gsl_matrix_set(WC_matrix, iterator, 2, gsl_matrix_get(A_target, i, 2));
        // printf("A: atom id: %c %d %s\t %f, %f, %f\n", A->atom_data->residue[A->WC_submatrix_rows[0][i]], A->atom_data->dnt_pos[A->WC_submatrix_rows[0][i]], A->atom_data->name[A->WC_submatrix_rows[0][i]], gsl_matrix_get(WC_matrix, iterator, 0), gsl_matrix_get(WC_matrix, iterator, 1), gsl_matrix_get(WC_matrix, iterator, 2));
        iterator++;
    }

    for (unsigned int i = 0; i < B_target->size1; i++)
    {
        gsl_matrix_set(WC_matrix, iterator, 0, gsl_matrix_get(B_target, i, 0));
        gsl_matrix_set(WC_matrix, iterator, 1, gsl_matrix_get(B_target, i, 1));
        gsl_matrix_set(WC_matrix, iterator, 2, gsl_matrix_get(B_target, i, 2));
        // printf("B: atom id: %c %d %s\t %f, %f, %f\n", B->atom_data->residue[B->WC_submatrix_rows[1][i]], B->atom_data->dnt_pos[B->WC_submatrix_rows[1][i]], B->atom_data->name[B->WC_submatrix_rows[1][i]], gsl_matrix_get(WC_matrix, iterator, 0), gsl_matrix_get(WC_matrix, iterator, 1), gsl_matrix_get(WC_matrix, iterator, 2));
        iterator++;
    }
    return WC_matrix;
}

bool is_hairpin(RNA_data_array &sequence, DimerLibArray &WC_Lib)
{
    RNA_data *WC_pair = new RNA_data(WC_Lib, 0, 0, true);
    gsl_matrix *sequence_matrix;
    gsl_matrix *WC_model_matrix;
    gsl_matrix *R;
    double _rmsd = 0;
    double COMP[] = {0, 0, 0}, COMQ[] = {0, 0, 0};
    bool is_h_pin = false;
    sequence_matrix = make_WC_submatrix(sequence[0], sequence[sequence.iterator_max]);
    WC_model_matrix = make_WC_submatrix(WC_pair, WC_pair);
    R = kabsch_get_rotation_matrix_generic(WC_model_matrix, sequence_matrix, COMP, COMQ);
    _rmsd = rmsd_generic(sequence_matrix, WC_model_matrix);
    if (_rmsd <= GLOBAL_WC_RMSD_LIMIT)
    {
        is_h_pin = true;
        sequence.update_WC_rmsd(_rmsd);
    }
    gsl_matrix_free(R);
    gsl_matrix_free(sequence_matrix);
    gsl_matrix_free(WC_model_matrix);
    delete WC_pair;
    return is_h_pin;
}
