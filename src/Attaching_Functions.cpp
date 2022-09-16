#include "Attaching_Functions.hpp"

attach_status rotate(RNAData *reference, RNAData *rotated)
{
    double rmsd_;
    double COMP[] = {0, 0, 0};
    double COMQ[] = {0, 0, 0};
    attach_status status = FAILED;

    gsl_matrix *MODEL = rotated->data_matrix;
    gsl_matrix *P = rotated->get_target_matrix_copy(0);
    gsl_matrix *Q = reference->get_target_matrix_copy(1);
    gsl_matrix *P_WORK = kabsch_allocate_work_matrix(P);
    gsl_matrix *R;
    
    /*
    printf("################Residue name: %s\n", reference->name);
    print_gsl_matrix(reference->data_matrix);
    print_gsl_matrix(Q);
    printf("################Residue name: %s\n", rotated->name);
    print_gsl_matrix(P);
    */
    kabsch_calculate_rotation_matrix_Nx3fast(P, Q, P_WORK, COMP, COMQ);
    R = kabsch_get_rotation_matrix();
    print_gsl_matrix(R);
    rmsd_ = rmsd_generic(P, Q);
    if (rmsd_ <= GLOBAL_RMSD_LIMIT)
    {
        status = ATTACHED;
        translate_matrix(COMP, MODEL, -1.0F);
        apply_rotation_matrix(R, MODEL);
        print_gsl_matrix(R);
        translate_matrix(COMQ, MODEL, 1.0F);
        *rotated->_flag = USED;
    }
    else
    {
        *rotated->_flag = NOT_USABLE;
        // printf("RMSD FAIL\n");
        *rotated->_flag = NOT_USABLE;
    }

    gsl_matrix_free(P);
    gsl_matrix_free(Q);
    gsl_matrix_free(P_WORK);
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