#ifndef DIMERLIB_HPP
#define DIMERLIB_HPP

#include "Atom_info.hpp"
#include "RNA_Math.hpp" //Temporary for testing...
#include "TAPERSS.hpp"

enum flag{NO_FLAG, USED, NOT_USABLE};

struct DimerLib
{
    gsl_block*        LibraryMemblock;
    gsl_matrix**      data_matrices;
    atom_info*        atom_data;
    float*            energy;
    char*             name;
    int               count; //Number of structures in Library
    double*           radii[2];

    void initialize(int n, int a_n, int e_n);
    void destroy();
    void save_lib(gsl_matrix** d_m, float* e, char* n, double** r);
    void clear_flags();
};

struct DimerLibArray
{    
    DimerLib** library;
    flag** Flags;
    int* PositiveAtomMap;           //Deallocation is handled by RNAData.
    int* NegativeAtomMap;           //Deallocation is handled by RNAData.
    bool* is_duplicate;
    uint32_t count = 0;          //Number of "Libraries" in array
    uint32_t full_structure_element_sum = 0; //Number of elements (atoms + COM)  which will be in a full structure.
    int iterator;
    bool was_initialized = false;
    uint32_t PositiveAtomCount;
    uint32_t NegativeAtomCount;
    uint32_t LargestAtomCount = 0;

    void initialize(int s);
		void destroy();
    DimerLib* operator[](int i);
    void alloc_lib(size_t n, size_t a_n, size_t e_n);
    void add_to_atom_info(char *N, int i, char r, int p);
    void map_duplicate(size_t org, size_t dupli);
    void add_lib(gsl_matrix** d_m, float* e, char* n, double** r);
    void get_charged_atom_map();
    void initialize_flags();
    void reset_flags(bool *reset);
    void print_dimer_info(int i);
    void print_matrix(int i, int j);
    gsl_matrix *get_matrix(int s, int index);
};


void get_model_count(FILE* fp, int* i);
void calculate_dnt_COM(gsl_matrix *A, atom_info *A_info);
void load_libs(RunInfo& run_info, DimerLibArray& rtn_lib_array, bool for_WC);

#endif
