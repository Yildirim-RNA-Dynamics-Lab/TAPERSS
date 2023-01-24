#ifndef RNA_DATA_HPP
#define RNA_DATA_HPP

#include "RNACMB.hpp"
#include "Atom_info.hpp"
#include "DimerLib.hpp"
#include "RNA_Math.hpp"
#include "HopcroftKarp.hpp"

const atom_id targetA[] = {N9, C8, N7, C5, C6, N1, C2, N3, C4, C1p, C2p, C3p, C4p, O4p}; // 14
const atom_id targetC[] = {N1, C2, N3, C4, C5, C6, C1p, C2p, C3p, C4p, O4p};             // 11
const atom_id targetG[] = {N9, C8, N7, C5, C6, N1, C2, N3, C4, C1p, C2p, C3p, C4p, O4p}; // 14
const atom_id targetU[] = {N1, C2, N3, C4, C5, C6, C1p, C2p, C3p, C4p, O4p};             // 11

const atom_id WCtargetA[] = {N9, C8, N7, C5, C6, N1, C2, N3, /*C4, N6, C1p, C2p, C3p, C4p, O4p*/};
const atom_id WCtargetC[] = {N1, C2, N3, C4, C5, C6, N4, O2, /*C1p, C2p, C3p, C4p, O4p*/};
const atom_id WCtargetG[] = {N9, C8, N7, C5, C6, N1, C2, N3, /*C4, N2, O6, C1p, C2p, C3p, C4p, O4p*/};
const atom_id WCtargetU[] = {N1, C2, N3, C4, C5, C6, O4, O2, /*C1p, C2p, C3p, C4p, O4p*/};

struct RNAData
{
    gsl_matrix   *data_matrix;
    //gsl_matrix   **submatrices; // Submatrix for each residue corresponding to the target atoms. Ordered 5' to 3'.
    //gsl_matrix   **WC_submatrices;
    uint16_t     *StericIndices[2];
    uint16_t     *EnergyIndices[4];  // 0,1: Positive and Negative for Res 1. 2,3: Positive and Negative for Res 2
    uint16_t     *submatrix_rows[2];
    uint16_t     *WC_submatrix_rows[2];
    atom_info    *atom_data;
    float energy;
    char name[3];
    flag *_flag;
    unsigned int count; // Number of atoms in structure

    int position_in_lib[2];
    int count_per_sub[2];
    int count_per_Steric[2];
    size_t ResBoundaries[4]; // 0,1: Start End of Res 1 Row Idx. 2,3: Start End of Res 2 Row Idx
    size_t count_per_Energy[4];
    size_t count_per_WC_sub[2];
    

    double COM_Radii[2];

    int position_max;
    const atom_id *target1;
    const atom_id *target2;
    const atom_id *WC_target1;
    const atom_id *WC_target2;

    long int id;

    RNAData();
    void initialize(DimerLibArray &L, int idx, int idx_L, gsl_block *MemBlock, size_t *offset_matrix, uint16_t *ArrayMemBlock, size_t *offset_array);
    void overwrite(DimerLibArray &L, int i, int j);
    void destroy();
    ~RNAData();

    /* Submatrices are now managed by the functions that used them for contiguous memory and better memory management */
    //void make_submatrices();
    //void update_submatrices();
    //gsl_matrix *get_target_matrix(int res);
    //gsl_matrix *get_target_matrix_copy(int res);
    //void make_WC_submatrices(bool first_run = false);
    //void update_WC_submatrices();
    //gsl_matrix *get_WC_target_matrix(int res);
    //gsl_matrix *get_WC_target_matrix_copy(int res);

    int get_residue_COM_index(int res);
    void print();
    void print_target(int res);
    void print_WC_target(int res);
    void reset_interactions();
    void print(int res);
    int to_string(char *s, int buffer_size, int string_index);
    void print_offset(int res, int position);
    int to_string_offset(int res, int position, char *s, int buffer_size, int string_index, int *idx_offset);
};

struct RNADataArray
{
    RNAData    *sequence;
    gsl_matrix *InteractionTable;
    gsl_matrix *COMS;
    gsl_block  *MatrixMemBlock;
    uint16_t   *ArrayMemBlock;

    double *Radii;
    bool *PassedCOMCheck;
    int *PositiveAtomMap;
    int *NegativeAtomMap;
    uint32_t LAC;
    int TableRowCount;
    int TableColCount;
    int count;
    int iterator;
    int iterator_max;
    float structure_energy;
    float WC_rmsd1_6 = -1.0F;

    char *string_out;
    int string_buffer;
    int string_index;
    int model_count;
    bool string_initialized = false;
    bool string_print_coordinates = GLOBAL_WRITE_COORDINATES;
    bool TMP_END = false;

    int atom_sum;

    RNADataArray();
    void initialize(size_t size, DimerLibArray &Library);
    void overwrite(size_t LibIdx, size_t IdxInLib, DimerLibArray &Library);
    void overwrite_initialize(size_t LibIdx, size_t IdxInLib, DimerLibArray &Library);
    ~RNADataArray();
    uint_fast64_t calculate_matrix_memory_needed(DimerLibArray &L, int i);
    uint_fast64_t calculate_index_array_memory_needed(DimerLibArray &L, int i);
    RNAData *operator[](int i);
    void add_copy(RNAData *A);
    void add_move(RNAData *A);
    RNAData *current();
    bool is_complete();
    void keep();
    void rollback();
    void safe_rollback(); // unused
    void rollback_by(int amount);
    bool is_empty();
    void update_WC_rmsd(float rmsd_val);
    void update_energy();
    void printall();
    void initialize_string();
    int get_atom_sum();
    int out_string_header_coord();
    int out_string_header();
    char *to_string();
    int *get_index();
    void print_index(int offset);
};

size_t get_target(char res, const atom_id **dest);
size_t get_WC_target(char res, const atom_id **dest);
int FindInteraction(int A_ResId, int A_Idx, RNAData *AData, int B_ResId, int B_Idx, RNAData *BData, 
                    gsl_matrix *AdjMatrix, int *PAdjMap, int*NAdjMap, size_t DIM2);

#endif