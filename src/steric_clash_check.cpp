#include "steric_clash_check.hpp"

bool steric_clash_check_fast_all(gsl_matrix *A, size_t A_start, size_t A_end, gsl_matrix *B, size_t B_start, size_t B_end)
{
    for(size_t i = A_start; i < A_end; i++)
    {
        for(size_t j = B_start; j < B_end; j++)
        {
            if(sqr_distance_mat2mat(A, i, B, j) < square(VDW_RADIUS))
            {
                return false;
            }
        }
    }
    return true;
}

bool steric_clash_check_fast_last(gsl_matrix *A, uint16_t *A_Rows, size_t A_Count, gsl_matrix *B, size_t B_start, size_t B_end)
{
    for(size_t i = 0; i < A_Count; i++)
    {
        for(size_t j = B_start; j < B_end; j++)
        {
            if(sqr_distance_mat2mat(A, A_Rows[i], B, j) < square(VDW_RADIUS))
            {
                return false;
            }
        }
    }
    return true;
}

void steric_clash_checkCOM(gsl_matrix* M_COMS, double* M_Radii, int Count, gsl_matrix* A, size_t A_index, double A_Radius, bool* PassArray)
{    
    PassArray[0] = (sqr_distance_mat2mat(M_COMS,  0, A, A_index) > square((M_Radii[0] + A_Radius)));
    PassArray[1] = (sqr_distance_mat2mat(M_COMS, 1, A, A_index) > square((M_Radii[1] + A_Radius)));

    for(int i = 1; i < Count; i++)
    {
        PassArray[i + 1] = (sqr_distance_mat2mat(M_COMS, i * 2 + 1 , A, A_index) > square((M_Radii[i * 2 + 1] + A_Radius)));
    }
}

attach_status steric_clash_check_COMFast(RNADataArray& Sequence, RNAData* Attach)
{  
    steric_clash_checkCOM(Sequence.COMS, Sequence.Radii, Sequence.count, Attach->data_matrix, Attach->get_residue_COM_index(1), Attach->COM_Radii[1], 
                          Sequence.PassedCOMCheck);
    if(Sequence.PassedCOMCheck[0] == false)
    {
        if(steric_clash_check_fast_all(Sequence[0]->data_matrix, Sequence[0]->ResBoundaries[0], Sequence[0]->ResBoundaries[1],
        Attach->data_matrix, Attach->ResBoundaries[2], Attach->ResBoundaries[3]) == false)
        {
            return attach_status::FAILED;
        }
    }
    if(Sequence.PassedCOMCheck[Sequence.count] == false)
    {
        if(steric_clash_check_fast_last(Sequence[Sequence.count - 1]->data_matrix, Sequence[Sequence.count - 1]->StericIndices[1], 
                                        Sequence[Sequence.count - 1]->count_per_Steric[1], Attach->data_matrix, 
                                        Attach->ResBoundaries[2], Attach->ResBoundaries[3]) == false)
        {
            return attach_status::FAILED;
        }
                    
    }
    for(int i = 0; i < Sequence.count - 1; i++)
    {
        if(Sequence.PassedCOMCheck[i + 1] == false)
        {
            if(steric_clash_check_fast_all(Sequence[i]->data_matrix, Sequence[i]->ResBoundaries[2], Sequence[i]->ResBoundaries[3],
              Attach->data_matrix, Attach->ResBoundaries[2], Attach->ResBoundaries[3]) == false)
            {
                return attach_status::FAILED;
            }
        }
    }
    return attach_status::ATTACHED;
}

template <bool countGreater>void steric_clash_checkCOM_IL(gsl_matrix* M_COMS, double* M_Radii, int Count, gsl_matrix* A, size_t A_index, double A_Radius, bool* PassArray, int size1)
{    
    PassArray[0] = (sqr_distance_mat2mat(M_COMS,  0, A, A_index) > square((M_Radii[0] + A_Radius)));
    PassArray[1] = (sqr_distance_mat2mat(M_COMS, 1, A, A_index) > square((M_Radii[1] + A_Radius)));
    if constexpr (countGreater == true) 
    {
        for(int i = 1; i < size1; i++)
        {
            PassArray[i + 1] = (sqr_distance_mat2mat(M_COMS, i * 2 + 1 , A, A_index) > square((M_Radii[i * 2 + 1] + A_Radius)));
        }
    
        PassArray[size1] = (sqr_distance_mat2mat(M_COMS, size1 * 2, A, A_index) > square((M_Radii[size1] + A_Radius)));
        PassArray[size1 + 1] = (sqr_distance_mat2mat(M_COMS, size1 * 2 + 1, A, A_index) > square((M_Radii[size1 + 1] + A_Radius)));

        for(int i = size1 + 2; i < Count; i++)
        {
            PassArray[i + 1] = (sqr_distance_mat2mat(M_COMS, i * 2 + 1, A, A_index) > square((M_Radii[i * 2 + 1] + A_Radius)));
        }
        return;
    }   
    if constexpr (countGreater == false)
    {
        for(int i = 1; i < Count; i++)
        {
            PassArray[i + 1] = (sqr_distance_mat2mat(M_COMS, i * 2 + 1 , A, A_index) > square((M_Radii[i * 2 + 1] + A_Radius)));
        }
        return;
    }
}
template void steric_clash_checkCOM_IL<true>(gsl_matrix* M_COMS, double* M_Radii, int Count, gsl_matrix* A, size_t A_index, double A_Radius, bool* PassArray, int size1);
template void steric_clash_checkCOM_IL<false>(gsl_matrix* M_COMS, double* M_Radii, int Count, gsl_matrix* A, size_t A_index, double A_Radius, bool* PassArray, int size1);

/**
 * @brief Performs a steric clash check on the newly started right side strand after proper alignment. Checks both 5' and 3' residues of the new nucleotide
 * using COM check first, then a further check on any failed COM checks.
 * @param Sequence - RNADataArrayInternalLoop reference.
 * @param Attach  - RNAData (new nucleotides) to be attached
 * @return attach_status - FAILED if steric clash is detected, ATTACHED otherwise
 */
attach_status steric_clash_check_COMFast_IL_1st_right(RNADataArrayInternalLoop& Sequence, RNAData* Attach)
{
    steric_clash_checkCOM_IL<false>(Sequence.COMS, Sequence.Radii, Sequence.count, Attach->data_matrix, Attach->get_residue_COM_index(0), Attach->COM_Radii[0], 
    Sequence.PassedCOMCheck, Sequence.WC_size_left);
    if(Sequence.PassedCOMCheck[0] == false)
    {
        if(steric_clash_check_fast_all(Sequence[0]->data_matrix, Sequence[0]->ResBoundaries[0], Sequence[0]->ResBoundaries[1],
        Attach->data_matrix, Attach->ResBoundaries[0], Attach->ResBoundaries[1]) == false)
        {
            return attach_status::FAILED;
        }
    }
    for(int i = 0; i < Sequence.count; i++)
    {
        if(Sequence.PassedCOMCheck[i + 1] == false)
        {
            if(steric_clash_check_fast_all(Sequence[i]->data_matrix, Sequence[i]->ResBoundaries[2], Sequence[i]->ResBoundaries[3],
              Attach->data_matrix, Attach->ResBoundaries[0], Attach->ResBoundaries[1]) == false)
            {
                return attach_status::FAILED;
            }
        }
    }
    steric_clash_checkCOM_IL<false>(Sequence.COMS, Sequence.Radii, Sequence.count, Attach->data_matrix, Attach->get_residue_COM_index(1), Attach->COM_Radii[1], 
    Sequence.PassedCOMCheck, Sequence.WC_size_left);
    if(Sequence.PassedCOMCheck[0] == false)
    {
        if(steric_clash_check_fast_all(Sequence[0]->data_matrix, Sequence[0]->ResBoundaries[0], Sequence[0]->ResBoundaries[1],
        Attach->data_matrix, Attach->ResBoundaries[2], Attach->ResBoundaries[3]) == false)
        {
            return attach_status::FAILED;
        }
    }
    for(int i = 0; i < Sequence.count; i++)
    {
        if(Sequence.PassedCOMCheck[i + 1] == false)
        {
            if(steric_clash_check_fast_all(Sequence[i]->data_matrix, Sequence[i]->ResBoundaries[2], Sequence[i]->ResBoundaries[3],
              Attach->data_matrix, Attach->ResBoundaries[2], Attach->ResBoundaries[3]) == false)
            {
                return attach_status::FAILED;
            }
        }
    }
    return attach_status::ATTACHED;
}

/**
 * @brief Performs steric clash check for addition of a new (non 5') nucleotide to internal loop. First the COM distances are calculated, 
 * then any which fail COM check have a full atom-to-atom check. When performing a steric clash check on the previous nucleotide (i.e. 4th being added
 * residue is checked with 3rd) only the phosphate groups are checked with steric_clash_check_fast_last(...).
 * @tparam CountGreater - template bool indicating whether or not the nucleotide being added is on the left side or right side. If true, nt is on the right
 * @param Sequence - RNADataArrayInternalLoop reference.
 * @param Attach  - RNAData (new nucleotide) to be attached
 * @return attach_status - FAILED if steric clash is detected, ATTACHED otherwise
 */
template <bool CountGreater>attach_status steric_clash_check_COMFast_IL(RNADataArrayInternalLoop& Sequence, RNAData* Attach)
{
    steric_clash_checkCOM_IL<CountGreater>(Sequence.COMS, Sequence.Radii, Sequence.count, Attach->data_matrix, Attach->get_residue_COM_index(1), Attach->COM_Radii[1], 
    Sequence.PassedCOMCheck, Sequence.WC_size_left);
    if(Sequence.PassedCOMCheck[0] == false)
    {
        if(steric_clash_check_fast_all(Sequence[0]->data_matrix, Sequence[0]->ResBoundaries[0], Sequence[0]->ResBoundaries[1],
        Attach->data_matrix, Attach->ResBoundaries[2], Attach->ResBoundaries[3]) == false)
        {
            return attach_status::FAILED;
        }
    }
    if(Sequence.PassedCOMCheck[Sequence.count] == false)
    {
        if(steric_clash_check_fast_last(Sequence[Sequence.count - 1]->data_matrix, Sequence[Sequence.count - 1]->StericIndices[1], 
                                        Sequence[Sequence.count - 1]->count_per_Steric[1], Attach->data_matrix, 
                                        Attach->ResBoundaries[2], Attach->ResBoundaries[3]) == false)
        {
            return attach_status::FAILED;
        }
                    
    }
    if constexpr(CountGreater == true) 
    {
        if(Sequence.PassedCOMCheck[Sequence.WC_size_left] == false)
        {
            if(steric_clash_check_fast_all(Sequence[Sequence.WC_size_left]->data_matrix, Sequence[Sequence.WC_size_left]->ResBoundaries[0], Sequence[Sequence.WC_size_left]->ResBoundaries[1],
            Attach->data_matrix, Attach->ResBoundaries[2], Attach->ResBoundaries[3]) == false)
            {
                return attach_status::FAILED;
            }
        }
        for(int i = 0; i < Sequence.count - 1; i++)
        {
            if(Sequence.PassedCOMCheck[i + 1] == false)
            {
                if(steric_clash_check_fast_all(Sequence[i]->data_matrix, Sequence[i]->ResBoundaries[2], Sequence[i]->ResBoundaries[3],
                Attach->data_matrix, Attach->ResBoundaries[2], Attach->ResBoundaries[3]) == false)
                {
                    return attach_status::FAILED;
                }
            }
        }
    }
    else if constexpr(CountGreater == false)
    {
        for(int i = 0; i < Sequence.count - 1; i++)
        {
            if(Sequence.PassedCOMCheck[i + 1] == false)
            {
                if(steric_clash_check_fast_all(Sequence[i]->data_matrix, Sequence[i]->ResBoundaries[2], Sequence[i]->ResBoundaries[3],
                Attach->data_matrix, Attach->ResBoundaries[2], Attach->ResBoundaries[3]) == false)
                {
                    return attach_status::FAILED;
                }
            }
        }
    }
    return attach_status::ATTACHED;
}

template attach_status steric_clash_check_COMFast_IL<true>(RNADataArrayInternalLoop& Sequence, RNAData* Attach);
template attach_status steric_clash_check_COMFast_IL<false>(RNADataArrayInternalLoop& Sequence, RNAData* Attach);

/* LEGACY, UNUSED FUNCTIONS BELOW...*/

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
                dist = distance_vec2vec(&A.vector, &B.vector);
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
            SCC_dist[0] = distance_vec2vec(&attach_COM.vector, &sequence_COM.vector);
            iterator_start = 0;
            sequence_radius[0] = sequence[i]->COM_Radii[0];
        } 
        else
        {
            SCC_dist[0] = 0;
            iterator_start = sequence[i]->atom_data->count_per_res[0];
        }

        sequence_COM = gsl_matrix_row(sequence[i]->data_matrix, sequence[i]->get_residue_COM_index(1));
        SCC_dist[1] = distance_vec2vec(&attach_COM.vector, &sequence_COM.vector);                    
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
                dist = distance_vec2vec(&A.vector, &B.vector);
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

attach_status steric_clash_check_COM(RNADataArray& sequence, RNAData* __restrict__ attach)
{
    gsl_vector_view A, B;
    double radius_1, radius_2;
    double dist;    
    unsigned int iterator_start;


    gsl_vector_view attach_COM, sequence_COM;
    double SCC_dist[2];
    double attach_radius = attach->COM_Radii[1];
    double sequence_radius[2];
    bool passed_check[2] = {true, true};
    //int this_index = sequence.count;

    attach_COM = gsl_matrix_row(attach->data_matrix, attach->get_residue_COM_index(1));    

    for (int i = 0; i < sequence.count; i++)
    {
        passed_check[0] = true;
        passed_check[1] = false;
        //printf("_______________________\n");
        //printf("steric clash index i: %d\n", i);
        if(i == 0)        
        {
            printf("Resid = %d\n", i);
            passed_check[0] = false;
            steric_clash_checks_attempted++;
            sequence_COM = gsl_matrix_row(sequence[i]->data_matrix, sequence[i]->get_residue_COM_index(0));
            SCC_dist[0] = distance_vec2vec(&attach_COM.vector, &sequence_COM.vector);
            iterator_start = 0;
            sequence_radius[0] = sequence[i]->COM_Radii[0];
            //printf("%d::Attach: \n", this_index);
            //gsl_vector_print(&attach_COM.vector);
            //printf("Radius: %f\n", attach_radius);

            //printf("0::Sequence: \n");
            //gsl_vector_print(&sequence_COM.vector);
            //printf("Radius: %f\n", sequence_radius[0]);

            printf("Distance 0: %f\n", SCC_dist[0]);
            printf("Limit 0: %f\n", sequence_radius[0] + attach_radius);
    

            if(SCC_dist[0] > (sequence_radius[0] + attach_radius))
            {                
                printf("PASSED CHECK\n");
                steric_clash_checks_skipped++;
                passed_check[0] = true;
            }
        } 
        else
        {
            SCC_dist[0] = 0;
            iterator_start = sequence[i]->atom_data->count_per_res[0];
        }

        printf("Resid = %d\n", i + 1);
        steric_clash_checks_attempted++;

        sequence_COM = gsl_matrix_row(sequence[i]->data_matrix, sequence[i]->get_residue_COM_index(1));
        SCC_dist[1] = distance_vec2vec(&attach_COM.vector, &sequence_COM.vector);                    
        sequence_radius[1] = sequence[i]->COM_Radii[1];

        printf("Distance 1: %f\n", SCC_dist[1]);
        printf("Limit 1: %f\n", sequence_radius[1] + attach_radius);        

        if(SCC_dist[1] > (sequence_radius[1] + attach_radius))
        {
            steric_clash_checks_skipped++;
            printf("PASSED CHECK\n");
            passed_check[1] = true;
        }

        //printf("_______________________\n");

        if(passed_check[0] == true && passed_check[1] == true)
        {
            //printf("Continuing...\n");
            continue;
        }
        if(passed_check[0] == true && passed_check[1] == false)
        {
            iterator_start = sequence[i]->atom_data->count_per_res[0];
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
                dist = distance_vec2vec(&A.vector, &B.vector);
                if (dist < (radius_1 + radius_2))
                {                    
                    return FAILED;
                }
            }
        }
    }
    return ATTACHED;
}

void SCC_record_COM_distance(RNADataArray &sequence, RNAData *attach)
{
    gsl_vector_view A, B, C, D;
    double dist;
    double limit;
    int this_index = sequence.count;

    C = gsl_matrix_row(attach->data_matrix, attach->get_residue_COM_index(0));
    D = gsl_matrix_row(attach->data_matrix, attach->get_residue_COM_index(1));

    for (int i = 0; i < sequence.count; i++)
    {
        if(i == 0)        
        {
            A = gsl_matrix_row(sequence[i]->data_matrix, sequence[i]->get_residue_COM_index(0));
            dist = distance_vec2vec(&A.vector, &C.vector);
            limit = sequence[i]->COM_Radii[0] + attach->COM_Radii[0];
            printf("%d:%c(%f)->%d:%c(%f): %f vs %f\n", i, sequence[i]->name[0], sequence[i]->COM_Radii[0], this_index, attach->name[0], attach->COM_Radii[0],dist, limit);
            dist = distance_vec2vec(&A.vector, &D.vector);
            limit = sequence[i]->COM_Radii[0] + attach->COM_Radii[1];
            printf("%d:%c(%f)->%d:%c(%f): %f vs %f\n", i, sequence[i]->name[0], sequence[i]->COM_Radii[0], this_index, attach->name[0], attach->COM_Radii[1],dist, limit);
        }        
        B = gsl_matrix_row(sequence[i]->data_matrix, sequence[i]->get_residue_COM_index(1));
        dist = distance_vec2vec(&B.vector, &C.vector);
        limit = sequence[i]->COM_Radii[0] + attach->COM_Radii[0];
        printf("%d:%c(%f)->%d:%c(%f): %f vs %f\n", i, sequence[i]->name[0], sequence[i]->COM_Radii[0], this_index, attach->name[0], attach->COM_Radii[0],dist, limit);          
        dist = distance_vec2vec(&B.vector, &D.vector);
        limit = sequence[i]->COM_Radii[0] + attach->COM_Radii[1];
        printf("%d:%c(%f)->%d:%c(%f): %f vs %f\n", i, sequence[i]->name[0], sequence[i]->COM_Radii[0], this_index, attach->name[0], attach->COM_Radii[1],dist, limit);
    }
}
