#ifndef RNADATAINTERNALLOOP_HPP
#define RNADATAINTERNALLOOP_HPP

#include "RNAData.hpp"
#include "RNACMB.hpp"
#include "Kabsch.hpp"

enum {left = 0, right = 1};

struct RNADataArrayInternalLoop : public RNADataArray
{
    int iterator_max_1;
    //RNAData** sequence;
    RNADataArrayInternalLoop() : RNADataArray() {};

    void initialize(int size1, int size2)
    {
        printf("~~~~~~~~~~~Size 1: %d\n", size1);
        iterator_max_1 = size1 - 1;
        iterator_max = size1 + size2 - 1;
        sequence = (RNAData **)malloc(sizeof(RNAData *) * (size1 + size2));
        iterator = -1;
        count = 0;
        structure_energy = 0;
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

    void prepare_right(RNAData* to_be_assembled, DimerLibArray &WC_Lib) 
    {
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
        gsl_matrix_free(P1);
        gsl_matrix_free(Q1);
        gsl_matrix_free(R1);

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

        gsl_matrix_free(P2);
        gsl_matrix_free(Q2);
        gsl_matrix_free(R2);
        
        delete WC_pair;
        to_be_assembled->update_submatrices();
        //add_move(to_be_assembled);
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
};

#endif