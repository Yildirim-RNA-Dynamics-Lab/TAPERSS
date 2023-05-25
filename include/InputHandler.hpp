#ifndef INPUTHANDLER_HPP
#define INPUTHANDLER_HPP

#include "RNACMB.hpp"


void get_DNMP_lib_names(char* sequence, char** rtn, int N_diNts);
int find_duplicates(char** rtn, int *duplicates, int N_diNts);
void get_WC_partner(char* sequence, char** rtn, int Nt1, int Nt2);
uint32_t get_index_int(char *index, uint32_t *indices);
int read_input_index_file(char *file_name, int n_DiNts);
void read_input_file(char *file_name);
void input_handler(int argc, char *ARGV[], RunInfo& run_info);
void free_libs(char **A, int s);

#endif
