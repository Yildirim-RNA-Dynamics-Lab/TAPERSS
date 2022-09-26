#ifndef RNADATAINTERNALLOOP_HPP
#define RNADATAINTERNALLOOP_HPP

#include "RNAData.hpp"
#include "RNACMB.hpp"
#include "Kabsch.hpp"
#include "Hairpin.hpp"

enum {left = 0, right = 1};

struct RNADataArrayInternalLoop : public RNADataArray
{
    float WC_rmsd3_4 = -1.0F;
    int iterator_max_1;
    //RNAData** sequence;
    RNADataArrayInternalLoop() : RNADataArray() {};
    int WC_size_left;

    void initialize(int size1, int size2, DimerLibArray& Library)
    {
        int LAC = Library.LargestAtomCount;
        WC_size_left = size1;
        //printf("~~~~~~~~~~~Size 1: %d\n", size1);
        iterator_max_1 = size1 - 1;
        iterator_max = size1 + size2 - 1;
        
        sequence = (RNAData **)malloc(sizeof(RNAData *) * (size1 + size2));
        uint64_t matrix_memsize = 0;
        uint64_t array_memsize = 0;
        for(int i = 0; i < size1 + size2; i++)
        {
            matrix_memsize   += calculate_matrix_memory_needed(Library, i);
            array_memsize    += calculate_index_array_memory_needed(Library, i);
        }
        printf("Array Memsize : %lu\nMatrix Memsize: %lu\n", array_memsize / sizeof(uint16_t), matrix_memsize);
        MatrixMemBlock = gsl_block_alloc(matrix_memsize); //Allocate Block of memory to pass to individual RNAData's. This is to maintain contiguous memory.
        ArrayMemBlock = (uint16_t *)malloc(array_memsize); // Same as ^.

        size_t MatrixOffset = 0;
        size_t ArrayOffset = 0;

        for(int i = 0; i < size1 + size2; i++)
        {
            //printf("BEFORE: Array Offset : %lu\nMatrix Offset: %lu\n", ArrayOffset, MatrixOffset);
            sequence[i] = new RNAData();
            sequence[i]->initialize(Library, i, 0, MatrixMemBlock, &MatrixOffset, ArrayMemBlock, &ArrayOffset);
            printf("AFTER: Array Offset : %lu\nMatrix Offset: %lu\n", ArrayOffset, MatrixOffset);
        }

        iterator = -1;
        count = 0;
        structure_energy = 0;
        InteractionTable = gsl_matrix_alloc(LAC, LAC);
        gsl_matrix_set_zero(InteractionTable); // Set diagonals to 1.
        //InteractionTableSum = (int *)calloc(LAC, sizeof(int));
        //PositiveAtomMap = Library.AtomMap;
        TableRowCount = LAC;
        COMS = gsl_matrix_alloc(2 * (iterator_max + 1), MATRIX_DIMENSION2);
        Radii = (double *)malloc(2 * (iterator_max + 1) * sizeof(double));
        PassedCOMCheck = (bool*)malloc((iterator_max + 1 + 2) * sizeof(gsl_vector*)); // +2 B/C 1st residue of left and 1st resiude of right need to be considered.
        memset(PassedCOMCheck, false, iterator_max + 1);
        for(int i = 0; i < iterator_max + 1; i++) 
        {
            printf("PassCOMCheck Boolean: %s\n", PassedCOMCheck[i] == true ? "true" : "false");
        }
    }

    bool inLeft_or_inRight(int working_pos) //return true if in left strand
    {
        return (working_pos == iterator_max_1 + 1);
    }

    void add_move(RNAData *A)
    {
        iterator++;
        gsl_matrix_row_copy(COMS, iterator * 2, A->data_matrix, A->get_residue_COM_index(0));
        gsl_matrix_row_copy(COMS, iterator * 2 + 1, A->data_matrix, A->get_residue_COM_index(1));
        Radii[iterator * 2] = A->COM_Radii[0];
        Radii[iterator * 2 + 1] = A->COM_Radii[1];
        sequence[iterator] = A;
        count++;
        *A->_flag = USED;
    }

    bool prepare_right(RNAData* to_be_assembled, DimerLibArray &WC_Lib) 
    {
        bool rmsd_pass = true;
        RNAData *assembled_ref = sequence[iterator_max_1];
        gsl_matrix *R1, *R2;
        double COMP[] = {0, 0, 0}, COMQ[] = {0, 0, 0};
        
        gsl_matrix *MODEL = WC_Lib[1]->data_matrices[0];
        gsl_matrix *P1 = kabsch_prepare_matrix<KABSCH_MATRIX_P>(assembled_ref->count_per_sub[1], MATRIX_DIMENSION2, assembled_ref->submatrix_rows[1], MODEL);
        gsl_matrix *Q1 = kabsch_prepare_matrix<KABSCH_MATRIX_Q>(assembled_ref->count_per_sub[1], MATRIX_DIMENSION2, assembled_ref->submatrix_rows[1], assembled_ref->data_matrix);
        gsl_matrix *P1WORK = kabsch_get_work_matrix(MATRIX_DIMENSION2, assembled_ref->count_per_sub[1]);
        kabsch_calculate_rotation_matrix_Nx3fast(P1, Q1, P1WORK, COMP, COMQ);
        R1 = kabsch_get_rotation_matrix();
        
        translate_matrix(COMP, MODEL, -1.0F);
        apply_rotation_matrix(R1, MODEL);
        translate_matrix(COMQ, MODEL, 1.0F);

        memset(COMP, 0, sizeof(COMP));
        memset(COMQ, 0, sizeof(COMQ));

        MODEL = to_be_assembled->data_matrix;
        gsl_matrix *P2 = kabsch_prepare_matrix<KABSCH_MATRIX_P>(to_be_assembled->count_per_sub[0], MATRIX_DIMENSION2, to_be_assembled->submatrix_rows[0], MODEL);
        gsl_matrix *Q2 = kabsch_prepare_matrix<KABSCH_MATRIX_Q>(assembled_ref->count_per_sub[1], MATRIX_DIMENSION2, assembled_ref->submatrix_rows[1], assembled_ref->data_matrix);
        gsl_matrix *P2WORK = kabsch_get_work_matrix(MATRIX_DIMENSION2, to_be_assembled->count_per_sub[0]);
        kabsch_calculate_rotation_matrix_Nx3fast(P2, Q2, P2WORK, COMP, COMQ);

        R2 = kabsch_get_rotation_matrix();

        translate_matrix(COMP, MODEL, -1.0F);
        apply_rotation_matrix(R2, MODEL);
        translate_matrix(COMQ, MODEL, 1.0F);

        memset(COMP, 0, sizeof(COMP));
        memset(COMQ, 0, sizeof(COMQ));

        WC_prepare_structure_matrix(1, assembled_ref->data_matrix, assembled_ref->WC_submatrix_rows[0], assembled_ref->count_per_WC_sub[0], 
                                       to_be_assembled->data_matrix, to_be_assembled->WC_submatrix_rows[0], to_be_assembled->count_per_WC_sub[0]);
        double RMSD = WC_check_pair(1);
        if(RMSD > GLOBAL_WC_RMSD_LIMIT)
        {
            rmsd_pass = false;
        }
        WC_rmsd3_4 = (float) RMSD;
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