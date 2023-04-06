#include "output_string.hpp"

output_string::output_string(const char *F, int max_s, char *action, char* string_prototype)
{
	output_file = fopen(F, action);
	string_storage = (char **)malloc(sizeof(char *) * max_s);
	max_size = strlen(string_prototype); 
	max_string = max_s;
	string_storage[0] = (char *)malloc(sizeof(char) * max_size * max_string + 1); //+1 b/c end of packed string needs to have '\0' included.
	string_storage[0][max_size * max_string] = '\0';
	for(int i = 0, j = 0; i < max_string; i++, j +=max_size)
	{
		string_storage[i] = &string_storage[0][j];
	}
	iterator = 0;
	num_allocated = 0;
}

output_string::~output_string()
{
	fprintf(output_file, "%s", string_storage[0]);
	free(string_storage[0]);
	free(string_storage);
	fclose(output_file);
}

void output_string::add_string(char *s)
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

void output_string::overwrite_string(char *src, int dest_idx)
{
	int shift_size = 0;
	if(dest_idx < iterator)
	{
		shift_size = iterator - dest_idx;
		//printf("Dest: %d, Max: %d, iter: %d\n", dest_idx + 1, max_string - 1, iterator);
		memmove(string_storage[dest_idx + 1], string_storage[dest_idx], max_size * (shift_size));
	}
	memmove(string_storage[dest_idx], src, max_size);
	if(dest_idx > iterator)
	{
		iterator = dest_idx;
	}
}
