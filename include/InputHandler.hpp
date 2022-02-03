#ifndef INPUTHANDLER_HPP
#define INPUTHANDLER_HPP

#include "RNACMB.hpp"

char** get_diNt_names(char* sequence, int *N_diNts);
char** get_WC_partner(char* sequence, int *N_WC);
void get_index_int(char *index, int *indices);
int read_input_index_file(char *file_name, int n_DiNts);
void read_input_file(char *file_name);
void input_handler(int argc, char *ARGV[]);
void free_libs(char **A, int s);

#endif