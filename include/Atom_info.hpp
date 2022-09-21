#ifndef ATOM_INFO_HPP
#define ATOM_INFO_HPP

#include "RNACMB.hpp"

enum atom_id{N1, N2, N3, N4, N6, N7, N9, C2, C4, C5, C6, C8, C1p, C2p, C3p, C4p, C5p, O2, O4, O6, O2p, O3p, O4p, O5p, O1P, O2P, OP1 = O1P, OP2 = O2P, P, NOTATOM};
enum atom_charge{POSITIVE = 1, NEGATIVE = -1, NEUTRAL = 0}; //simplistic charge model for rudimentary energy calculations.

struct atom_info
{
    char**        name;
    int*          index;
    char*         residue;
    int*          dnt_pos;
    atom_id*      atom_ids;
    atom_charge*  charges;
    uint64_t      count;            // Number of atoms per structure
    uint_fast32_t count_per_res[2]; // Number of atoms per residue 1 and 2
    uint_fast32_t iterator;
    uint_fast32_t num_charged_atoms = 0;

    
    atom_info();
    atom_info(int n);
    atom_info(char** n, int* in, char* r, int* d, int c);
    void initialize_atom_info(int n);
    ~atom_info();
    void add_atom(char *N, int i, char r, int p);
    void clear();
    atom_info(const atom_info &A);
    void move(atom_info& A);
    int get_idx_of_atom(atom_id atom_name, int residx);
    void print_at(int line);
};


#endif