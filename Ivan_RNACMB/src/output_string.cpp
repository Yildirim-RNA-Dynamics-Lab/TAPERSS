#include "output_string.hpp"

output_string::output_string(const char *F, int max_s, char *action)
{
    output_file = fopen(F, action);
    string_storage = (char **)malloc(sizeof(char *) * max_s);
    iterator = 0;
    max_string = max_s;
    num_allocated = 0;
}
output_string::~output_string()
{
    for (int i = 0; i < max_string; i++)
    {
        if (all_init)
        {
            if (i < iterator && iterator < max_string)
            {
                fprintf(output_file, "%s", string_storage[i]);
            }
            free(string_storage[i]);
        }
        else
        {
            if (i < iterator && iterator < max_string)
            {
                fprintf(output_file, "%s", string_storage[i]);
                free(string_storage[i]);
            }
        }
    }
    free(string_storage);
    fclose(output_file);
}

void output_string::add_string(char *s, int max_atoms)
{
    if (iterator == max_string - 1)
    {
        if (all_init)
        {
            strcpy(string_storage[iterator], s);
            iterator = 0;
        }
        else
        {
            all_init = true;
            int size = strlen(s);
            string_storage[iterator] = (char *)malloc(sizeof(char) * (size + (max_atoms * 3)));
            num_allocated++;
            strcpy(string_storage[iterator], s);
            iterator = 0;
        }
        for (int i = 0; i < max_string; i++)
        {
            fprintf(output_file, "%s", string_storage[i]);
        }
        return;
    }
    else
    {
        if (all_init)
        {
            strcpy(string_storage[iterator], s);
            iterator++;
        }
        else
        {
            int size = strlen(s);
            string_storage[iterator] = (char *)malloc(sizeof(char) * (size + (3 * max_atoms)));
            num_allocated++;
            strcpy(string_storage[iterator], s);
            iterator++;
        }
    }
    return;
}
