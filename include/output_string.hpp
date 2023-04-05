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
	int max_size;

	output_string(const char *F, int max_s, char *action, char* string_prototype);
	~output_string();
	void add_string(char *s);
	void overwrite_string(char *src, int dest_idx);
};

#endif
