#include "InputHandler.hpp"

char **get_diNt_names(char *sequence, int *N_diNts)
{
    int name_position = -1;
    *N_diNts = strlen(sequence) - 1;
    char **rtn = (char **)malloc(sizeof(char *) * *N_diNts);
    while (LIBRARY_FILENAME_PROTOTYPE[name_position++] != 'X')
        ;
    name_position -= 1;
    for (int i = 0; i < *N_diNts; i++)
    {
        rtn[i] = (char *)malloc(sizeof(char) * sizeof(LIBRARY_FILENAME_PROTOTYPE));
        memcpy(rtn[i], LIBRARY_FILENAME_PROTOTYPE, sizeof(LIBRARY_FILENAME_PROTOTYPE));
        memcpy(&rtn[i][name_position], sequence++, sizeof(char) * 2);
    }
    return rtn;
}

char **get_WC_partner(char *sequence, int *N_WC)
{
    char WC_pair[3];

    int name_position = 0;
    while (WATSON_CRICK_LIBRARY_PROTOTYPE[name_position++] != 'X')
        ;
    name_position -= 1;

    int seq_count = 0;
    while (sequence[seq_count++] != '\0')
        ;
    seq_count -= 2;

    WC_pair[0] = sequence[0];
    WC_pair[1] = sequence[seq_count];
    WC_pair[2] = '\0';
    *N_WC = 1;

    if (WC_pair[0] == 'A')
    {
        if (WC_pair[1] != 'C' && WC_pair[1] != 'U')
        {
            printf("Not potential hairpin, exiting... ( this is a temporary fix while program is dumb )\n");
            exit(2);
        }
    }
    else if (WC_pair[0] == 'U')
    {
        if (WC_pair[1] != 'A' && WC_pair[1] != 'G')
        {
            printf("Not potential hairpin, exiting... ( this is a temporary fix while program is dumb )\n");
            exit(2);
        }
    }
    else if (WC_pair[0] == 'C')
    {
        if (WC_pair[1] != 'A' && WC_pair[1] != 'G')
        {
            printf("Not potential hairpin, exiting... ( this is a temporary fix while program is dumb )\n");
            exit(2);
        }
    }
    else if (WC_pair[0] == 'G')
    {
        if (WC_pair[1] != 'C' && WC_pair[1] != 'U')
        {
            printf("Not potential hairpin, exiting... ( this is a temporary fix while program is dumb. Instead it should look through whole sequence to find h-bond pairs and then check if hairpin is possible)\n");
            exit(2);
        }
    }

    char **rtn = (char **)malloc(sizeof(char *) * *N_WC);
    rtn[0] = (char *)malloc(sizeof(WATSON_CRICK_LIBRARY_PROTOTYPE));

    memcpy(rtn[0], WATSON_CRICK_LIBRARY_PROTOTYPE, sizeof(WATSON_CRICK_LIBRARY_PROTOTYPE));
    memcpy(&rtn[0][name_position], WC_pair, sizeof(char) * 2);
    return rtn;
}

void get_index_int(char *index, int *indices)
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
    char line[100];
    char indices[100];
    int str_count = 0;

    while (fgets(line, sizeof(line), input)) // Get line count
    {
        char *header;
        header = strtok(line, " ");
        if (!strcmp(header, GLOBAL_INPUT_SEQUENCE))
            str_count++;
    }
    rewind(input);

    GLOBAL_INPUT_INDICES_LIST = (int **)malloc(str_count * sizeof(int *));

    for (int i = 0; i < str_count; i++)
    {
        GLOBAL_INPUT_INDICES_LIST[i] = (int *)malloc(n_DiNts * sizeof(int));
        fgets(line, sizeof(line), input);
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

    return str_count;
}

void read_input_file(char *file_name)
{
    FILE *input = fopen(file_name, "r");
    char line[100];
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
                GLOBAL_WRITE_COORDINATES = true;
            }
            if (!strcasecmp(str1, "strlist"))
            {
                GLOBAL_RUN_BUILD_STRUCTURE_LIST = true;
                GLOBAL_WRITE_COORDINATES = true;
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
                GLOBAL_PERFORM_HAIRPIN_CHECK = true;
            }
            if (!strcasecmp(str1, "NONE"))
            {
                GLOBAL_PERFORM_HAIRPIN_CHECK = false;
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