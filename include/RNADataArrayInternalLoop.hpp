#ifndef RNADATAINTERNALLOOP_HPP
#define RNADATAINTERNALLOOP_HPP

#include "RNACMB.hpp"
#include "RNAData.hpp"
#include "Kabsch.hpp"
#include "WatsonCrickPair.hpp"

enum {left = 0, right = 1};

struct RNADataArrayInternalLoop : public RNADataArray
{
    double WC_rmsd3_4 = -1.0F;
    int iterator_max_1;
    int WC_size_left;

    uint16_t *WC_rows[2];
    size_t WC_rowsize[2];
    void initialize(int size1, int size2, DimerLibArray& Library, DimerLibArray& WC_Library)
    {
        const atom_id* WC_target1;
        const atom_id* WC_target2;

        size_t target1Size = get_target(WC_Library[1]->name[0], &WC_target1);
        size_t target2Size = get_target(WC_Library[1]->name[1], &WC_target2);

        WC_rowsize[0] = target1Size;
        WC_rowsize[1] = target2Size;

        int PLAC = Library.PositiveAtomCount;
        int NLAC = Library.NegativeAtomCount;
        //printf("PLAC = %d, NLAC = %d\n", PLAC, NLAC);
        WC_size_left = size1;
        iterator_max_1 = size1 - 1;
        iterator_max = size1 + size2 - 1;
        size_t size = size1 + size2;
        printf("size1 = %d, size2 = %d\n", size1, size2);
        sequence = (RNAData *)malloc(sizeof(RNAData) * (size));
        uint64_t matrix_memsize = 0;
        uint64_t array_memsize = 0;
        for(size_t i = 0; i < size; i++)
        {
            matrix_memsize   += calculate_matrix_memory_needed(Library, i);
            array_memsize    += calculate_index_array_memory_needed(Library, i);
        }
        //printf("Array Memsize : %lu\nMatrix Memsize: %lu\n", array_memsize / sizeof(uint16_t), matrix_memsize);
        array_memsize += sizeof(uint16_t) * (target1Size + target2Size);
        matrix_memsize += (PLAC * NLAC); // Interaction Table Matrix
        //printf("Array Memsize : %lu\nMatrix Memsize: %lu\n", array_memsize / sizeof(uint16_t), matrix_memsize);
        MatrixMemBlock = gsl_block_alloc(matrix_memsize); //Allocate Block of memory to pass to individual RNAData's. This is to maintain contiguous memory.
        ArrayMemBlock = (uint16_t *)malloc(array_memsize); // Same as ^.

        size_t MatrixOffset = 0;
        size_t ArrayOffset = 0;

        for(size_t i = 0; i < size; i++)
        {
            //printf("BEFORE: Array Offset : %lu\nMatrix Offset: %lu\n", ArrayOffset, MatrixOffset);
            sequence[i].initialize(Library, i, 0, MatrixMemBlock, &MatrixOffset, ArrayMemBlock, &ArrayOffset);
            
        }

        //printf("AFTER: Array Offset : %lu\nMatrix Offset: %lu\n", ArrayOffset, MatrixOffset);
        *sequence[0]._flag = USED;
        iterator = 0;  //Start with first already DNT included
        count = 1;
        structure_energy = 0;
        COMS = gsl_matrix_alloc_from_block(MatrixMemBlock, MatrixOffset, 2 * size, MATRIX_DIMENSION2, MATRIX_DIMENSION2);
        MatrixOffset += (2 * size * MATRIX_DIMENSION2);
        WC_rows[0] = &ArrayMemBlock[ArrayOffset];
        ArrayOffset += target1Size;
        //printf("AFTER2: Array Offset : %lu\nMatrix Offset: %lu\n", ArrayOffset, MatrixOffset);
        //printf("Target 2 size: %lu\n", target2Size);
        //printf("Array memsize: %lu\n", array_memsize);
        WC_rows[1] = &ArrayMemBlock[ArrayOffset];
        for(int i = 0; i < WC_Library[1]->atom_data->count; i++)
        {
            for(size_t j = 0; j < target1Size; j++) 
            {
                if(WC_Library[1]->atom_data->atom_ids[i] == WC_target1[j] && WC_Library[1]->atom_data->dnt_pos[i] == 1)
                {
                    WC_rows[0][j] = (uint16_t)i;
                }
            }
            for(size_t j = 0; j < target2Size; j++) 
            {
                if(WC_Library[1]->atom_data->atom_ids[i] == WC_target2[j] && WC_Library[1]->atom_data->dnt_pos[i] == 2)
                {
                    WC_rows[1][j] = (uint16_t)i;
                }
            }
        }
        //printf("AFTER2: Array Offset : %lu\nMatrix Offset: %lu\n", ArrayOffset, MatrixOffset);
        InteractionTable = gsl_matrix_alloc_from_block(MatrixMemBlock, MatrixOffset, PLAC, NLAC, NLAC);
        gsl_matrix_set_zero(InteractionTable); // Set all to 0
        PositiveAtomMap = Library.PositiveAtomMap;
        NegativeAtomMap = Library.NegativeAtomMap;
        TableRowCount = PLAC;
        TableColCount = NLAC;
        LAC = Library.LargestAtomCount;
        Radii = (double *)malloc(2 * size * sizeof(double));
        PassedCOMCheck = (bool*)malloc((size + 2) * sizeof(bool)); // in normal building of structure only first DNT will have valid 1st residue.
                                                                        // in all other DNTs, only second residue is needed.
        overwrite_initialize(0, 0, Library);
    }

    bool should_prepare_right(int working_pos) //return true if in left strand
    {
        return (working_pos == iterator_max_1 + 1);
    }

    /*
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
    */

    void print_target(gsl_matrix* A, size_t size, uint16_t* rows, atom_info* a_data)
    {
        for (unsigned int i = 0; i < size; i++)
        {
            a_data->print_at(rows[i]);
            for (unsigned int j = 0; j < A->size2; j++)
                printf("%8.3f", gsl_matrix_get(A, rows[i], j));
            putchar('\n');
        }
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

        string_index = sequence[0].to_string(string_out, string_buffer, string_index);
        idx_offset = sequence[0].count;

        for (int i = 1, j = 1; j < count; i++, j++)
        {
            if(j == iterator_max_1 + 1)
            {
                i++;
                string_index += snprintf(&string_out[string_index], string_buffer - string_index, "TER\n");
                string_index = sequence[j].to_string_offset(0, i, string_out, string_buffer, string_index, &idx_offset);
                string_index = sequence[j].to_string_offset(1, i, string_out, string_buffer, string_index, &idx_offset);
            }
            else
            {
                string_index = sequence[j].to_string_offset(1, i, string_out, string_buffer, string_index, &idx_offset);
            }
        }
        
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "TER\n");
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
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%d", sequence[i].position_in_lib[1]);
            if (i != count - 1)
            {
                string_index += snprintf(&string_out[string_index], string_buffer - string_index, "-");
            }
        }
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, " %f ", structure_energy);
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%f ", WC_rmsd0_N);
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
                string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%d", sequence[i].position_in_lib[1]);
                if (i != count - 1)
                {
                    string_index += snprintf(&string_out[string_index], string_buffer - string_index, "-");
                }
            }
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\t");
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#ENERGY ");
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%f\t", structure_energy);
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#WC-RMSD ");
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%f", WC_rmsd0_N);
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\n");
        }
        return string_index;
    }
};

#endif
