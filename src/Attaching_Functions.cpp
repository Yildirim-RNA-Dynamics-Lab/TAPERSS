#include "Attaching_Functions.hpp"

attach_status rotate(RNAData *reference, RNAData *rotated)
{
    double rmsd_;
    double COMP[] = {0, 0, 0};
    double COMQ[] = {0, 0, 0};
    attach_status status = FAILED;

    gsl_matrix *MODEL = rotated->data_matrix;
    gsl_matrix *P = kabsch_prepare_matrix<KABSCH_MATRIX_P>(rotated->count_per_sub[0], MATRIX_DIMENSION2, rotated->submatrix_rows[0], rotated->data_matrix);
    gsl_matrix *Q = kabsch_prepare_matrix<KABSCH_MATRIX_Q>(reference->count_per_sub[1], MATRIX_DIMENSION2, reference->submatrix_rows[1], reference->data_matrix);
    gsl_matrix *P_WORK = kabsch_get_work_matrix(MATRIX_DIMENSION2, rotated->count_per_sub[0]);
    gsl_matrix *R;
    
    
    //printf("################Residue name: %s\n", reference->name);
    //print_gsl_matrix(reference->data_matrix);
    print_gsl_matrix(P);
    print_gsl_matrix(Q);
    //printf("################Residue name: %s\n", rotated->name);

    kabsch_calculate_rotation_matrix_Nx3fast(P, Q, P_WORK, COMP, COMQ);
    //print_gsl_matrix(P);
    //print_gsl_matrix(Q);
    R = kabsch_get_rotation_matrix();
    //print_gsl_matrix(R);
    rmsd_ = rmsd_generic(P, Q);
    //print_gsl_matrix(P);
    //print_gsl_matrix(Q);
    printf("RMSD: %f\n", rmsd_);
    if (rmsd_ <= GLOBAL_RMSD_LIMIT)
    {
        status = ATTACHED;
        translate_matrix(COMP, MODEL, -1.0F);
        apply_rotation_matrix(R, MODEL);
        //print_gsl_matrix(R);
        translate_matrix(COMQ, MODEL, 1.0F);
        *rotated->_flag = USED;
    }
    else
    {
        *rotated->_flag = NOT_USABLE;
    }

    return status;
}

attach_status check_attachment(RNADataArray &sequence, RNAData *attach)
{
    attach_status status;
    if ((status = steric_clash_check_COM(sequence, attach)) != ATTACHED)
    {
        //SCC_record_COM_distance(sequence, attach);
        *attach->_flag = NOT_USABLE;
        return status;
    }
    return ATTACHED;
}