#ifndef RNADATAINTERNALLOOP_HPP
#define RNADATAINTERNALLOOP_HPP

#include "RNAData.hpp"
#include "RNACMB.hpp"
#include "Kabsch.hpp"
#include "Hairpin.hpp"

enum {left = 0, right = 1};
static double COMZERO[] = {0, 0, 0};

struct RNADataArrayInternalLoop : public RNADataArray
{
    float WC_rmsd3_4 = -1.0F;
    int iterator_max_1;
    //RNAData** sequence;
    RNADataArrayInternalLoop() : RNADataArray() {};

    void initialize(int size1, int size2, int* AtomMap, int LAC)
    {
        //printf("~~~~~~~~~~~Size 1: %d\n", size1);
        iterator_max_1 = size1 - 1;
        iterator_max = size1 + size2 - 1;
        sequence = (RNAData **)malloc(sizeof(RNAData *) * (size1 + size2));
        iterator = -1;
        count = 0;
        structure_energy = 0;
        InteractionTable = gsl_matrix_alloc(LAC, LAC);
        gsl_matrix_set_zero(InteractionTable); // Set diagonals to 1.
        InteractionTableSum = (int *)calloc(LAC, sizeof(int));
        InteractionTableMap = AtomMap;
        TableRowCount = LAC;
    }

    bool inLeft_or_inRight(int working_pos) //return true if in left strand
    {
        if(working_pos == iterator_max_1 + 1)
        {
            return true;
        }
        else
        {
            return false;
        }

    }

    void add_move(RNAData *A)
    {
        sequence[++iterator] = A;
        count++;
        *A->_flag = USED;
    }

    bool prepare_right(RNAData* to_be_assembled, DimerLibArray &WC_Lib) 
    {
        bool rmsd_pass = true;
        RNAData *WC_pair = new RNAData(WC_Lib, 1, 0, false);
        RNAData *assembled_ref = sequence[iterator_max_1];
        gsl_matrix *R1, *R2;
        double COMP[] = {0, 0, 0}, COMQ[] = {0, 0, 0};
        
        gsl_matrix *MODEL = WC_pair->data_matrix;
        gsl_matrix *P1 = WC_pair->get_target_matrix_copy(0);
        //assembled_ref->make_submatrices();
        gsl_matrix *Q1 = assembled_ref->get_target_matrix_copy(1);
        /*
        printf("################Residue name: %s\n", WC_pair->name);
        printf("################Residue name: %s\n", assembled_ref->name);
        print_gsl_matrix(P1);
        print_gsl_matrix(Q1);
        */

        R1 = kabsch_get_rotation_matrix_generic(P1, Q1, COMP, COMQ);
        

        translate_matrix(COMP, MODEL, -1.0F);
        apply_rotation_matrix(R1, MODEL);
        translate_matrix(COMQ, MODEL, 1.0F);

        //WC_pair->make_submatrices();
        WC_pair->update_submatrices();
        //gsl_matrix_free(P1);
        //gsl_matrix_free(Q1);
        gsl_matrix_free(R1);

        memcpy(COMP, COMZERO, sizeof(COMZERO));
        memcpy(COMQ, COMZERO, sizeof(COMZERO));

        MODEL = to_be_assembled->data_matrix;
        //to_be_assembled->make_submatrices();
        gsl_matrix *P2 = to_be_assembled->get_target_matrix_copy(0);
        gsl_matrix *Q2 = WC_pair->get_target_matrix_copy(1);

        /*
        printf("################Residue name: %s\n", to_be_assembled->name);
        printf("################Residue name: %s\n", WC_pair->name);
        print_gsl_matrix(P2);
        print_gsl_matrix(Q2);
        */ 
        
        R2 = kabsch_get_rotation_matrix_generic(P2, Q2, COMP, COMQ);

        translate_matrix(COMP, MODEL, -1.0F);
        apply_rotation_matrix(R2, MODEL);
        translate_matrix(COMQ, MODEL, 1.0F);

        memcpy(COMP, COMZERO, sizeof(COMZERO));
        memcpy(COMQ, COMZERO, sizeof(COMZERO));

        gsl_matrix *sequence_matrix = gsl_matrix_alloc(Q1->size1 + P2->size2, MATRIX_DIMENSION2);
        gsl_matrix *WC_model_matrix = gsl_matrix_alloc(P1->size1 + Q2->size2, MATRIX_DIMENSION2);
        gsl_matrix *work_matrix;
        gsl_matrix *R;
        double _rmsd;

        overwrite_WC_submatrix_gsl(Q1, P2, sequence_matrix);
        overwrite_WC_submatrix_gsl(P1, Q2, WC_model_matrix);
        work_matrix = kabsch_allocate_work_matrix(sequence_matrix);
        kabsch_calculate_rotation_matrix_Nx3fast(WC_model_matrix, sequence_matrix, work_matrix, COMP, COMQ);
        R = kabsch_get_rotation_matrix();
        apply_rotation_matrix(R, WC_model_matrix);
        _rmsd = rmsd_generic(sequence_matrix, WC_model_matrix);

        if(_rmsd > GLOBAL_WC_RMSD_LIMIT)
        {
            rmsd_pass = false;
        }
        //printf("WC RMSD Resid 3 & 4: %f\n", _rmsd);
        WC_rmsd3_4 = (float) _rmsd;

        gsl_matrix_free(P2);
        gsl_matrix_free(Q2);
        gsl_matrix_free(R2);
        
        delete WC_pair;
        to_be_assembled->update_submatrices();
        //add_move(to_be_assembled);
        gsl_matrix_free(P1);
        gsl_matrix_free(Q1);
        gsl_matrix_free(sequence_matrix);
        gsl_matrix_free(WC_model_matrix);
        gsl_matrix_free(work_matrix);
        gsl_matrix_free(R);

        return rmsd_pass;
    }

    char* to_string()
    {
        int idx_offset;

        if (string_initialized)
        {
            string_index = 0;
        }
        else
        {
            initialize_string();
        }

        string_index = out_string_header();

        if (string_print_coordinates == false)
            return string_out;

        string_index = sequence[0]->to_string(string_out, string_buffer, string_index);
        idx_offset = sequence[0]->count;

        for (int i = 1, j = 1; j < count; i++, j++)
        {
            if(j == iterator_max_1 + 1)
            {
                i++;
                string_index = sequence[j]->to_string_offset(0, i, string_out, string_buffer, string_index, &idx_offset);
                string_index = sequence[j]->to_string_offset(1, i, string_out, string_buffer, string_index, &idx_offset);
            }
            else
            {
                string_index = sequence[j]->to_string_offset(1, i, string_out, string_buffer, string_index, &idx_offset);
            }
        }

        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "ENDMDL\n");

        return string_out;
    }
    int out_string_header_coord()
    {
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "MODEL %d\n", ++model_count);
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "REMARK ");
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%s ", GLOBAL_INPUT_SEQUENCE);
        for (int i = 0; i < count; i++)
        {
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%d", sequence[i]->position_in_lib[1]);
            if (i != count - 1)
            {
                string_index += snprintf(&string_out[string_index], string_buffer - string_index, "-");
            }
        }
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, " %f ", structure_energy);
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%f ", WC_rmsd1_6);
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%f", WC_rmsd3_4);
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\n");
        return string_index;
    }
    int out_string_header()
    {
        if (GLOBAL_WRITE_COORDINATES)
        {
            string_index = out_string_header_coord();
        }
        else
        {
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%s ", GLOBAL_INPUT_SEQUENCE);
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#INDEX ");
            for (int i = 0; i < count; i++)
            {
                string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%d", sequence[i]->position_in_lib[1]);
                if (i != count - 1)
                {
                    string_index += snprintf(&string_out[string_index], string_buffer - string_index, "-");
                }
            }
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\t");
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#ENERGY ");
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%f\t", structure_energy);
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#WC-RMSD ");
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%f", WC_rmsd1_6);
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\n");
        }
        return string_index;
    }
};

#endif