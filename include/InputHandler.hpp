#ifndef INPUTHANDLER_HPP
#define INPUTHANDLER_HPP

#include "RNACMB.hpp"

#define LINESIZE 1000

int get_diNt_names(char* sequence, char** rtn, int *duplicates, int N_diNts);
void get_WC_partner(char* sequence, char** rtn, int Nt1, int Nt2);
void get_index_int(char *index, int *indices);
int read_input_index_file(char *file_name, int n_DiNts);
void read_input_file(char *file_name);
void input_handler(int argc, char *ARGV[]);
void free_libs(char **A, int s);

#endif