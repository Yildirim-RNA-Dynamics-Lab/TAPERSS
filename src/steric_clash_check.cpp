#include "steric_clash_check.hpp"

bool steric_clash_check(RNADataArray &sequence, RNAData *attach)
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

attach_status steric_clash_check_COM_tester(RNADataArray& __restrict__ sequence, RNAData* __restrict__ attach)
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
            SCC_dist[0] = 0;
            iterator_start = sequence[i]->atom_data->count_per_res[0];
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
            for (unsigned int k = attach->atom_data->count_per_res[0]; k < attach->count; k++)
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
                        if(SCC_dist[0] > (sequence_radius[0] + attach_radius))
                        {
                            printf("FAILED: %d:%c->%d:%c: %f\n", i, sequence[i]->name[0], this_index, attach->name[1], SCC_dist[0]);
                            printf("New Residue Radius == %f, Model Residue Radius == %f\n", attach_radius, sequence_radius[0]);
                            printf("%s :: %s = %f\n", sequence[i]->atom_data->name[j], attach->atom_data->name[k], dist);
                            gsl_vector_print(&sequence_COM.vector);
                            gsl_vector_print(&attach_COM.vector);
                            return FAILED_SC;
                        }
                       // printf("New Residue Radius == %f, Model Residue Radius == %f\n", attach_radius, sequence_radius[0]);
                    }
                    else
                    {
                        if(SCC_dist[1] > (sequence_radius[1] + attach_radius))
                        {
                            printf("FAILED: %d:%c->%d:%c: %f\n", i, sequence[i]->name[1], this_index, attach->name[1], SCC_dist[1]);
                            printf("New Residue Radius == %f, Model Residue Radius == %f\n", attach_radius, sequence_radius[1]);
                            printf("%s :: %s = %f\n", sequence[i]->atom_data->name[j], attach->atom_data->name[k], dist);
                            gsl_vector_print(&sequence_COM.vector);
                            gsl_vector_print(&attach_COM.vector);
                            return FAILED_SC;                            
                        }
                        //printf("New Residue Radius == %f, Model Residue Radius == %f\n", attach_radius, sequence_radius[1]);
                    }
                    return FAILED;
                }
            }
        }
    }
    return ATTACHED;
}



attach_status steric_clash_check_COM(RNA_data_array& __restrict__ sequence, RNA_data* __restrict__ attach)
{
    gsl_vector_view A, B;
    double radius_1, radius_2;
    double dist;    
    unsigned int iterator_start;


    gsl_vector_view attach_COM, sequence_COM;
    double SCC_dist[2];
    double attach_radius = attach->COM_Radii[1];
    double sequence_radius[2];
    //int this_index = sequence.count;

    attach_COM = gsl_matrix_row(attach->data_matrix, attach->get_residue_COM_index(1));    

    for (int i = 0; i < sequence.count; i++)
    {
        //printf("steric clash index i: %d\n", i);
        if(i == 0)        
        {
            steric_clash_checks_attempted++;
            sequence_COM = gsl_matrix_row(sequence[i]->data_matrix, sequence[i]->get_residue_COM_index(0));
            SCC_dist[0] = distance(&attach_COM.vector, &sequence_COM.vector);
            iterator_start = 0;
            sequence_radius[0] = sequence[i]->COM_Radii[0];
            //printf("%d::Attach: \n", this_index);
            //gsl_vector_print(&attach_COM.vector);
            //printf("Radius: %f\n", attach_radius);

            //printf("0::Sequence: \n");
            //gsl_vector_print(&sequence_COM.vector);
            //printf("Radius: %f\n", sequence_radius[0]);
    

            if(SCC_dist[0] > (sequence_radius[0] + attach_radius))
            {
                steric_clash_checks_skipped++;
                continue;
            }
        } 
        else
        {
            SCC_dist[0] = 0;
            iterator_start = sequence[i]->atom_data->count_per_res[0];
        }

        steric_clash_checks_attempted++;

        sequence_COM = gsl_matrix_row(sequence[i]->data_matrix, sequence[i]->get_residue_COM_index(1));
        SCC_dist[1] = distance(&attach_COM.vector, &sequence_COM.vector);                    
        sequence_radius[1] = sequence[i]->COM_Radii[1];

        if(SCC_dist[1] > (sequence_radius[1] + attach_radius))
        {
            steric_clash_checks_skipped++;
            continue;
        }

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
            for (unsigned int k = attach->atom_data->count_per_res[0]; k < attach->count; k++)
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
                    return FAILED;
                }
            }
        }
    }
    return ATTACHED;
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