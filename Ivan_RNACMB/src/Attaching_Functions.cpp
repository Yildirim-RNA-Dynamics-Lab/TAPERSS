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
    for (int i = 0; i < sequence.count; i++)
    {
        for (int j = 0; j < sequence[i]->count; j++)
        {
            if (i == sequence.count - 1 && sequence[i]->atom_data->dnt_pos[j] == 2)
            {
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
                if (dist < (radius_1 + radius_2))
                {
                    return false;
                }
            }
        }
    }
    return true;
}

attach_status check_attachment(RNA_data_array &sequence, RNA_data *attach)
{
    if (!overlap_check(sequence, attach))
    {
        *attach->_flag = NOT_USABLE;
        return FAILED;
    }
    return ATTACHED;
}