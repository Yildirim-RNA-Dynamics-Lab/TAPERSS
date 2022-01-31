#ifndef OUTPUT_STRING_HPP
#define OUTPUT_STRING_HPP

#include "RNACMB.hpp"

struct output_string
{
    FILE *output_file;
    char **string_storage;
    int iterator;
    int num_allocated;
    int max_string;

    bool all_init = false;

    output_string(const char *F, int max_s, char *action);
    ~output_string();
    void add_string(char *s, int max_atoms);
};

#endif