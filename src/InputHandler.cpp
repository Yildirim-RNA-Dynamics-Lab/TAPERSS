#include "InputHandler.hpp"

int get_diNt_names(char *sequence, char** rtn, int *duplicates, int N_diNts)
{
    int name_position = 0;
    int num_duplicates = 0;
    //*N_diNts = strlen(sequence) - 1;
    //char **rtn = (char **)malloc(sizeof(char *) * *N_diNts);
    while (LIBRARY_FILENAME_PROTOTYPE[name_position] != 'X')
    {
        name_position++;
    };
    //name_position -= 1;

    for (int i = 0; i < N_diNts; i++)
    {
        rtn[i] = (char *)malloc(sizeof(LIBRARY_FILENAME_PROTOTYPE));
        memcpy(rtn[i], LIBRARY_FILENAME_PROTOTYPE, sizeof(LIBRARY_FILENAME_PROTOTYPE));
        memcpy(&rtn[i][name_position], sequence++, sizeof(char) * 2);
    }


    for (int i = N_diNts - 1; i > -1; i--)
    {
        for(int j = 0; j < i; j++)
        {
            if(!strcmp(rtn[i], rtn[j]))
            {
                duplicates[i] = j;
                num_duplicates++;
            }
        }
    }
    return num_duplicates;
}

void get_WC_partner(char *sequence, char** rtn, int Nt1, int Nt2)
{
    char WC_pair[3];

    int name_position = 0;
    while (WATSON_CRICK_LIBRARY_PROTOTYPE[name_position] != 'X')
    {
        name_position++;
    }

    WC_pair[0] = sequence[Nt1];
    WC_pair[1] = sequence[Nt2];
    WC_pair[2] = '\0';

    if (WC_pair[0] == 'A')
    {
        if (WC_pair[1] != 'C' && WC_pair[1] != 'U')
        {
            printf("Not potential hairpin, exiting... ( this is a temporary fix while program is dumb )\n");
            printf("Pair Attempted: %s\n", WC_pair);
            exit(2);
        }
    }
    else if (WC_pair[0] == 'U')
    {
        if (WC_pair[1] != 'A' && WC_pair[1] != 'G')
        {
            printf("Not potential hairpin, exiting... ( this is a temporary fix while program is dumb )\n");
            printf("Pair Attempted: %s\n", WC_pair);
            exit(2);
        }
    }
    else if (WC_pair[0] == 'C')
    {
        if (WC_pair[1] != 'A' && WC_pair[1] != 'G')
        {
            printf("Not potential hairpin, exiting... ( this is a temporary fix while program is dumb )\n");
            printf("Pair Attempted: %s\n", WC_pair);
            exit(2);
        }
    }
    else if (WC_pair[0] == 'G')
    {
        if (WC_pair[1] != 'C' && WC_pair[1] != 'U')
        {
            printf("Not potential hairpin, exiting... ( this is a temporary fix while program is dumb. Instead it should look through whole sequence to find h-bond pairs and then check if hairpin is possible)\n");
            printf("Pair Attempted: %s\n", WC_pair);
            exit(2);
        }
    }

    *rtn = (char *)malloc(sizeof(WATSON_CRICK_LIBRARY_PROTOTYPE));
    memcpy(*rtn, WATSON_CRICK_LIBRARY_PROTOTYPE, sizeof(WATSON_CRICK_LIBRARY_PROTOTYPE));
    memcpy(&rtn[0][name_position], WC_pair, sizeof(char) * 2);
}

void get_index_int(char *index, uint32_t *indices)
{
    char buffer[5];
    int count = 0, buf_c = 0;
    for (char c = *index; c != '\0'; c = *(++index))
    {
        if (c == '-')
        {
            buffer[buf_c] = '\0';
            indices[count] = atoi(buffer);
            buf_c = 0;
            count++;
            continue;
        }
        buffer[buf_c++] = *index;
    }
    buffer[buf_c] = '\0';
    indices[count] = atoi(buffer);
}

int read_input_index_file(char *file_name, int n_DiNts)
{
    FILE *input = fopen(file_name, "r");
    char line[GLOBAL_STANDARD_STRING_LENGTH];
    char indices[GLOBAL_STANDARD_STRING_LENGTH];
    int str_count = 0;

    while (fgets(line, sizeof(line), input)) // Get line count
    {
        char *header;
        header = strtok(line, " ");
        if (!strcmp(header, GLOBAL_INPUT_SEQUENCE))
            str_count++;
    }
    rewind(input);

    GLOBAL_INPUT_INDICES_LIST = (uint32_t **)malloc(str_count * sizeof(int *));

    for (int i = 0; i < str_count; i++)
    {
        GLOBAL_INPUT_INDICES_LIST[i] = (uint32_t *)malloc(n_DiNts * sizeof(int));
        if(fgets(line, sizeof(line), input) != NULL)
        {
            char *str1;
            char *header = strtok(line, " ");
            if (!strcmp(header, GLOBAL_INPUT_SEQUENCE))
            {
                str1 = strtok(NULL, " ");
                str1 = strtok(NULL, " ");
                str1[strcspn(str1, "\t")] = '\0'; // equivalent to chomp() from perl
                strcpy(indices, str1);
            }
            get_index_int(indices, GLOBAL_INPUT_INDICES_LIST[i]);
        }
    }
    fclose(input);
    return str_count;
}

void read_input_file(char *file_name)
{
    FILE *input = fopen(file_name, "r");
    if(input == nullptr)
    {
        printf("Could not open file: %s\n", file_name);
        exit(3);
    }
    GLOBAL_OUTPUT_FILE[0] = '\0';
<<<<<<< HEAD
    char line[GLOBAL_STANDARD_STRING_LENGTH];
=======
    char line[LINESIZE];
>>>>>>> cmb_optimization
    while (fgets(line, sizeof(line), input))
    {
        char *str1;
        char *header = strtok(line, " ");
        if (!strcmp(header, "SEQUENCE"))
        {
            str1 = strtok(NULL, " ");
            str1 = strtok(NULL, " ");
            str1[strcspn(str1, "\n")] = '\0'; // equivalent to chomp() from perl
            strcpy(GLOBAL_INPUT_SEQUENCE, str1);
            // printf("SEQ: %s\n", GLOBAL_INPUT_SEQUENCE);
            continue;
        }
        if (!strcmp(header, "OUTPUT"))
        {
            str1 = strtok(NULL, " ");
            str1 = strtok(NULL, " ");
            str1[strcspn(str1, "\n")] = '\0';
            strcpy(GLOBAL_OUTPUT_FILE, str1);
            continue;
        }
        if (!strcmp(header, "OUTPUTACTION"))
        {
            str1 = strtok(NULL, " ");
            str1 = strtok(NULL, " ");
            GLOBAL_IO_ACTION[0] = *str1;
            continue;
        }
        if (!strcmp(header, "TYPE"))
        {
            str1 = strtok(NULL, " ");
            str1 = strtok(NULL, " ");
            str1[strcspn(str1, "\n")] = '\0';
            if (!strcasecmp(str1, "cmb"))
            {
                GLOBAL_RUN_COMBINATORIAL = true;
            }
            else if (!strcasecmp(str1, "str"))
            {
                GLOBAL_RUN_BUILD_STRUCTURE = true;
            }
            if (!strcasecmp(str1, "strlist"))
            {
                GLOBAL_RUN_BUILD_STRUCTURE_LIST = true;
            }
            if (!strcasecmp(str1, "strlisttest"))
            {
                GLOBAL_RUN_BUILD_STRUCTURE_LIST_TESTING = true;
            }

            continue;
        }
        if (!strcmp(header, "STRUCTCHECK"))
        {
            str1 = strtok(NULL, " ");
            str1 = strtok(NULL, " ");
            str1[strcspn(str1, "\n")] = '\0';
            if (!strcasecmp(str1, "HAIRPIN"))
            {
                GLOBAL_PERFORM_STRUCTCHECK = HAIRPIN;
            }
            if(!strcasecmp(str1, "INTERNAL_LOOP"))
            {
                GLOBAL_PERFORM_STRUCTCHECK = INTERNAL_LOOP;
            }
            if (!strcasecmp(str1, "NONE"))
            {
                GLOBAL_PERFORM_STRUCTCHECK = NONE;
            }
        }
        if (!strcmp(header, "INPUTFILE"))
        {
            str1 = strtok(NULL, " ");
            str1 = strtok(NULL, " ");
            str1[strcspn(str1, "\n")] = '\0';
            strcpy(GLOBAL_INPUT_FILE, str1);
            continue;
        }
        if (!strcmp(header, "INDICES"))
        {
            str1 = strtok(NULL, " ");
            str1 = strtok(NULL, " ");
            str1[strcspn(str1, "\n")] = '\0';
            strcpy(GLOBAL_INPUT_INDICES, str1);
            // printf("%s\n", GLOBAL_INPUT_INDICES);
            continue;
        }
        if (!strcmp(header, "RMSD"))
        {
            str1 = strtok(NULL, " ");
            str1 = strtok(NULL, " ");
            GLOBAL_RMSD_LIMIT = atof(str1);
            continue;
        }
        if (!strcmp(header, "WCRMSD"))
        {
            str1 = strtok(NULL, " ");
            str1 = strtok(NULL, " ");
            GLOBAL_WC_RMSD_LIMIT = atof(str1);
            continue;
        }
        if (!strcmp(header, "LIBRARY"))
        {
            str1 = strtok(NULL, " ");
            str1 = strtok(NULL, " ");
            str1[strcspn(str1, "\n")] = '\0';
            strcpy(LIBRARY_FILENAME_PROTOTYPE, str1);
            // printf("%s\n", LIBRARY_FILENAME_PROTOTYPE);
            continue;
        }
        if (!strcmp(header, "WCLIBRARY"))
        {
            str1 = strtok(NULL, " ");
            str1 = strtok(NULL, " ");
            str1[strcspn(str1, "\n")] = '\0';
            strcpy(WATSON_CRICK_LIBRARY_PROTOTYPE, str1);
            // printf("%s\n", WATSON_CRICK_LIBRARY_PROTOTYPE);
            continue;
        }
        if (!strcmp(header, "WRITE_COORDINATES"))
        {
            str1 = strtok(NULL, " ");
            str1 = strtok(NULL, " ");
            str1[strcspn(str1, "\n")] = '\0';
            if (!strcasecmp(str1, "TRUE"))
            {
                GLOBAL_WRITE_COORDINATES = true;
            }
            else if (!strcasecmp(str1, "FALSE"))
            {
                GLOBAL_WRITE_COORDINATES = false;
            }
            continue;
        }
    }
    fclose(input);
}

void input_handler(int argc, char *ARGV[])
{
    if (argc == 1)
    {
        printf("No Input\n");
        exit(1);
    }
    for (int i = 0; i < argc; i++)
    {
        if (ARGV[i][0] == '-')
        {
            switch (ARGV[i][1])
            {
            case 'i':
                read_input_file(ARGV[i + 1]);
                break;
            case 'n':
                printf("%s\n", ARGV[i + 1]); 
                strcpy(GLOBAL_INPUT_FILE, ARGV[i + 1]);
                printf("%s\n", GLOBAL_INPUT_FILE);
                break;
            case 'o':
                // printf("%s\n", ARGV[i + 1]);
                strcpy(GLOBAL_OUTPUT_FILE, ARGV[i + 1]);
                break;
            case 's':
                // printf("%s\n", ARGV[i + 1]);
                strcpy(GLOBAL_INPUT_SEQUENCE, ARGV[i + 1]);
                break;
            case 't':
                if (!strcasecmp(ARGV[i + 1], "cmb"))
                {
                    GLOBAL_RUN_COMBINATORIAL = true;
                    // printf("Performing Combinatorial Run\n");
                }
                else if (!strcasecmp(ARGV[i + 1], "str"))
                {
                    GLOBAL_RUN_BUILD_STRUCTURE = true;
                    // printf("Performing Single Structure Build\n");
                }
                break;
            case 'x':
                strcpy(GLOBAL_INPUT_INDICES, ARGV[i + 1]);
                break;
            case 'c':
                GLOBAL_WRITE_COORDINATES = true;
                break;
            case 'r':
                GLOBAL_RMSD_LIMIT = atof(ARGV[i + 1]);
                printf("RMSD CUTOFF: %f\n", GLOBAL_RMSD_LIMIT);
                break;
            case 'a':
                GLOBAL_IO_ACTION[0] = 'a';
                printf("Appending to output...\n");
                break;
            }
        }
    }
}

void free_libs(char **A, int s)
{
    for (int i = 0; i < s; i++)
        free(A[i]);
    free(A);
}