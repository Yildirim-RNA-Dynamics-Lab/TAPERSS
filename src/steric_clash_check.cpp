#include "steric_clash_check.hpp"

bool steric_clash_check(RNA_data_array &sequence, RNA_data *attach)
{
    gsl_vector_view A, B;
    double radius_1, radius_2;
    double dist;
    for (int i = 0; i < sequence.count; i++)
    {
        for (unsigned int j = 0; j < sequence[i]->count; j++)
        {
            if (i == sequence.count - 1 && sequence[i]->atom_data->dnt_pos[j] == 2)
            {
                switch (sequence[i]->atom_data->atom_ids[j])
                {
                case OP1:
                case OP2:
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
            for (unsigned int k = 0; k < attach->count; k++)
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

attach_status steric_clash_check_COM(RNA_data_array &sequence, RNA_data *attach)
{
    gsl_vector_view A, B;
    double radius_1, radius_2;
    double dist;    
    unsigned int iterator_start;


    gsl_vector_view attach_COM, sequence_COM;
    double SCC_dist[2];
    double attach_radius = attach->COM_Radii[1];
    double sequence_radius[2];
    int this_index = sequence.count;

    attach_COM = gsl_matrix_row(attach->data_matrix, attach->get_residue_COM_index(1));    

    for (int i = 0; i < sequence.count; i++)
    {
        if(i == 0)        
        {
            sequence_COM = gsl_matrix_row(sequence[i]->data_matrix, sequence[i]->get_residue_COM_index(0));
            SCC_dist[0] = distance(&attach_COM.vector, &sequence_COM.vector);
            iterator_start = 0;
            sequence_radius[0] = sequence[i]->COM_Radii[0];
        } 
        else
        {
            iterator_start = (unsigned int)sequence[i]->atom_data->count_per_res[0];
        }

        sequence_COM = gsl_matrix_row(sequence[i]->data_matrix, sequence[i]->get_residue_COM_index(1));
        SCC_dist[1] = distance(&attach_COM.vector, &sequence_COM.vector);                    
        sequence_radius[1] = sequence[i]->COM_Radii[1];

        for (unsigned int j = iterator_start; j < sequence[i]->count; j++)
        {
            if (i == sequence.count - 1 && sequence[i]->atom_data->dnt_pos[j] == 2)
            {
                switch (sequence[i]->atom_data->atom_ids[j])
                {
                case OP1:
                case OP2:
                case P:
                case O5p:
                case C5p:
                    break;
                default:
                    continue;
                }
            }
            for (unsigned int k = (unsigned int)attach->atom_data->count_per_res[0]; k < attach->count; k++)
            {
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
                    if(sequence[i]->atom_data->dnt_pos[j] == 1)
                    {   
                        if(SCC_dist[0] >= (sequence_radius[0] + attach_radius))
                        {
                            printf("FAILED: %d:%c->%d:%c: %f\n", i, sequence[i]->name[0], this_index, attach->name[1], SCC_dist[0]);
                            printf("%s :: %s = %f\n", sequence[i]->atom_data->name[j], attach->atom_data->name[j], dist);
                            print_vector(&sequence_COM.vector);
                            print_vector(&attach_COM.vector);
                            return FAILED_SC;
                        }
                    }
                    else
                    {
                        if(SCC_dist[1] >= (sequence_radius[1] + attach_radius))
                        {
                            printf("FAILED: %d:%c->%d:%c: %f\n", i, sequence[i]->name[1], this_index, attach->name[1], SCC_dist[1]);
                            printf("%s :: %s = %f\n", sequence[i]->atom_data->name[j], attach->atom_data->name[j], dist);
                            print_vector(&sequence_COM.vector);
                            print_vector(&attach_COM.vector);
                            return FAILED_SC;                            
                        }
                    }
                    return FAILED;
                }
            }
        }
    }
    return ATTACHED;
}

void print_vector(gsl_vector* V)
{
    printf("%f %f %f\n", gsl_vector_get(V, 0), gsl_vector_get(V, 1), gsl_vector_get(V, 2));
}

void SCC_record_COM_distance(RNA_data_array &sequence, RNA_data *attach)
{
    gsl_vector_view A, B, D;
    double dist;
    int this_index = sequence.count;

    //C = gsl_matrix_row(attach->data_matrix, attach->get_residue_COM_index(0));
    D = gsl_matrix_row(attach->data_matrix, attach->get_residue_COM_index(1));

    for (int i = 0; i < sequence.count; i++)
    {
        if(i == 0)        
        {
            A = gsl_matrix_row(sequence[i]->data_matrix, sequence[i]->get_residue_COM_index(0));
            //dist = distance(&A.vector, &C.vector);
            //printf("%d:%c->%d:%c: %f\n", i, sequence[i]->name[0], this_index, attach->name[0], dist);
            dist = distance(&A.vector, &D.vector);
            printf("%d:%c->%d:%c: %f\n", i, sequence[i]->name[0], this_index, attach->name[1], dist);
        }        
        B = gsl_matrix_row(sequence[i]->data_matrix, sequence[i]->get_residue_COM_index(1));
        //dist = distance(&B.vector, &C.vector);
        //printf("%d:%c->%d:%c: %f\n", i, sequence[i]->name[1], this_index, attach->name[0], dist);
        dist = distance(&B.vector, &D.vector);
        printf("%d:%c->%d:%c: %f\n", i, sequence[i]->name[1], this_index, attach->name[1], dist);
    }
}