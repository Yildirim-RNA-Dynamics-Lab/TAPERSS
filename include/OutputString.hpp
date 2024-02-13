#ifndef OUTPUT_STRING_HPP
#define OUTPUT_STRING_HPP

#include "TAPERSS.hpp"

struct OutputString
{
	FILE *output_file;
	char **string_storage;
	int iterator;
	int num_allocated;
	int max_string;
	int max_size;
	bool type_insert;

	void initialize(RunInfo& run_info, char* string_prototype);
	void destroy();
	void add_string(char *s);
	void insert_string(char *src, int dest_idx);
	void fix_model_counts();
};

#endif
