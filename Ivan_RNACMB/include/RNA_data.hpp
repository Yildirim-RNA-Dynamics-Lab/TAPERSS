#ifndef RNA_DATA_HPP
#define RNA_DATA_HPP

#include "RNACMB.hpp"
#include "Atom_info.hpp"
#include "DimerLib.hpp"

struct RNA_data
{
    gsl_matrix *data_matrix;
    gsl_matrix **submatrices; //Submatrix for each residue corresponding to the target atoms. Ordered 5' to 3'.
    gsl_matrix **WC_submatrices;
    float energy;
    atom_info *atom_data;
    char *name;
    flag *_flag;
    int count; //Number of atoms in structure
    bool *has_interaction;
    
    int *submatrix_rows[2];
    int *WC_submatrix_rows[2];

    int position_in_lib[2];
    int count_per_sub[2];    
    int sub_starts_at[2];
    size_t count_per_WC_sub[2];

    int position_max;
    atom_id *target;
    
    bool is_for_WC;
    bool WC_secondary;
    
    long int id;

    RNA_data(DimerLibArray& L, int i, int j, bool WC = false);
    void overwrite(DimerLibArray& L, int i, int j);
    ~RNA_data();
    void make_submatrices();
    void update_submatrices();
    size_t get_target(int res);
    gsl_matrix* get_target_matrix(int res);
    gsl_matrix* get_target_matrix_copy(int res);
    void make_WC_submatrices(bool first_run = false);
    void update_WC_submatrices();
    size_t get_WC_target(int res);
    gsl_matrix* get_WC_target_matrix(int res);
    gsl_matrix* get_WC_target_matrix_copy(int res);
    void print();
    void print_target(int res);
    void print_WC_target(int res);
    void reset_interactions();
    void print(int res);
    int to_string(char *s, int buffer_size, int string_index);
    void print_offset(int res, int position);
    int to_string_offset(int res, int position, char *s, int buffer_size, int string_index, int *idx_offset);
};

struct RNA_data_array
{
    RNA_data **sequence;
    gsl_matrix *WC_submatrix;
    int count;
    int iterator;
    int iterator_max;
    float structure_energy;
    float WC_rmsd = -1.0F;

    char *string_out;
    int string_buffer;
    int string_index;
    int model_count;
    bool string_initialized = false;
    bool string_print_coordinates = GLOBAL_WRITE_COORDINATES;

    int atom_sum;

    RNA_data_array(int size);
    ~RNA_data_array();
    RNA_data* operator[](int i);
    void add_copy(RNA_data* A);
    void add_move(RNA_data* A);
    RNA_data* current();
    bool is_complete();
    void rollback();
    void safe_rollback();    //unused
    void rollback_by(int amount);
    bool is_empty();
    void update_WC_rmsd(float rmsd_val);
    void reset_interactions();
    void printall();
    void initialize_string();
    int get_atom_sum();
    int out_string_header();
    char* to_string();
    int *get_index();
};

#endif