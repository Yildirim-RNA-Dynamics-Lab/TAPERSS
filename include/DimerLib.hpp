#ifndef DIMERLIB_HPP
#define DIMERLIB_HPP

#include "Atom_info.hpp"
#include "RNA_Math.hpp" //Temporary for testing...
#include "RNACMB.hpp"

enum flag{NO_FLAG, USED, NOT_USABLE};

struct DimerLib
{
    gsl_matrix**      data_matrices;
    atom_info*        atom_data;
    float*            energy;
    char*             name;
    int               count; //Number of structures in Library
    flag*             flags; 

    DimerLib(int n, int a_n, int e_n);
    ~DimerLib();
    void save_lib(gsl_matrix** d_m, float* e, char* n);
    void clear_flags();
};

struct DimerLibArray
{    
    DimerLib** library;
    int count = 0;          //Number of "Libraries" in array
    int iterator;
    bool was_initialized = false;

    DimerLibArray();
    DimerLibArray(int s);
    ~DimerLibArray();
    void initialize(int s);
    DimerLib* operator[](int i);
    void alloc_lib(int n, int a_n, int e_n);
    void add_to_atom_info(char *N, int i, char r, int p);
    void add_lib(gsl_matrix** d_m, float* e, char* n);
    void reset_flags(bool *reset);
    void print_dimer_info(int i);
    void print_matrix(int i, int j);
    gsl_matrix *get_matrix(int s, int index);
};


void get_model_count(FILE* fp, int* i);
void load_libs(char **LibNames, int N_diNts, DimerLibArray& RTN, bool for_WC = false);

#endif