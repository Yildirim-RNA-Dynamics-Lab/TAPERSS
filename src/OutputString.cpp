#include "OutputString.hpp"

void OutputString::initialize(RunInfo& run_info, char* string_prototype)
{
	output_file = fopen(run_info.output_file, "w");
	if(output_file == NULL) {
		fprintf(stderr, "Cannot open output file: %s!\n", run_info.output_file);
		exit(2);
	}
	max_size = strlen(string_prototype); 
	max_string = (uint32_t)(run_info.memory_limit / (max_size + sizeof(char *)));
	string_storage = (char **)malloc(sizeof(char *) * max_string);
	string_storage[0] = (char *)malloc(sizeof(char) * max_size * max_string + 1); //+1 b/c end of packed string needs to have '\0' included.
	string_storage[0][max_size * max_string] = '\0';
	for(int i = 0, j = 0; i < max_string; i++, j +=max_size)
	{
		string_storage[i] = &string_storage[0][j];
	}
	iterator = 0;
	num_allocated = 0;
}

void OutputString::destroy()
{
	fprintf(output_file, "%s", string_storage[0]);
	free(string_storage[0]);
	free(string_storage);
	fclose(output_file);
}

void OutputString::add_string(char *s)
{
	if (iterator == max_string - 1)
	{
		memcpy(string_storage[iterator], s, max_size);
		iterator = 0;
		fprintf(output_file, "%s", string_storage[0]);
		return;
	}
	else
	{
		memcpy(string_storage[iterator], s, max_size);
		iterator++;
	}
	return;
}

void OutputString::insert_string(char *src, int dest_idx)
{
	int shift_size = 0;
	if(dest_idx < iterator)
	{
		shift_size = iterator - dest_idx;
		memmove(string_storage[dest_idx + 1], string_storage[dest_idx], max_size * (shift_size));
	}
	memmove(string_storage[dest_idx], src, max_size);
	if(dest_idx > iterator)
	{
		iterator = dest_idx;
	}
}
