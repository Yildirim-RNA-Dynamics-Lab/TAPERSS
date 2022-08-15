#include "Attaching_Functions.hpp"

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
        *rotated->_flag = NOT_USABLE;
        // printf("RMSD FAIL\n");
    }

    gsl_matrix_free(P);
    gsl_matrix_free(Q);
    gsl_matrix_free(R);
    return status;
}

attach_status check_attachment(RNA_data_array &sequence, RNA_data *attach)
{
    attach_status status;
    if ((status = steric_clash_check_COM(sequence, attach)) != ATTACHED)
    {
        *attach->_flag = NOT_USABLE;
        return status;
    }
    return ATTACHED;
}