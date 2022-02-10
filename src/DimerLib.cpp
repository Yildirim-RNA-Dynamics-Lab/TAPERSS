#include "DimerLib.hpp"
/**
 * @brief Construct a DimerLib
 * 
 * @param n Number of structures in a library
 * @param a_n Number of elements in each structure (atoms + extra)
 * @param e_n Number of extra elements in a_n
 */
DimerLib::DimerLib(int n, int a_n, int e_n)
{
    data_matrices = (gsl_matrix **)malloc(sizeof(gsl_matrix *) * n);
    name = (char *)malloc(sizeof(char) * 3);
    flags = (flag *)calloc(n, sizeof(flag));
    energy = (float *)malloc(sizeof(float) * n);
    count = n;
    atom_data = new atom_info(a_n - e_n);
}
DimerLib::~DimerLib()
{
    // printf("DimerLib %s Destructor Called\n", name);
    for (int i = 0; i < count; i++)
    {
        gsl_matrix_free(data_matrices[i]);
    }
    free(energy);
    free(data_matrices);
    free(flags);
    free(name);
    delete atom_data;
}

void DimerLib::save_lib(gsl_matrix **d_m, float *e, char *n)
{
    for (int i = 0; i < count; i++)
    {
        data_matrices[i] = d_m[i];
    }
    memcpy(energy, e, sizeof(float) * count);
    memcpy(name, n, sizeof(char) * 3);
}

void DimerLib::clear_flags()
{
    for (int i = 0; i < count; i++)
    {
        flags[i] = NO_FLAG;
    }
}

DimerLibArray::DimerLibArray()
{
    // Do nothing
}

DimerLibArray::DimerLibArray(int s)
{
    library = (DimerLib **)malloc(sizeof(DimerLib *) * s);
    count = s;
    iterator = 0;
    was_initialized = true;
}

DimerLibArray::~DimerLibArray()
{
    if (was_initialized)
    {
        for (int i = 0; i < count; i++)
        {
            delete library[i];
        }
        free(library);
    }
}

void DimerLibArray::initialize(int s)
{
    library = (DimerLib **)malloc(sizeof(DimerLib *) * s);
    count = s;
    iterator = 0;
    was_initialized = true;
}

DimerLib *DimerLibArray::operator[](int i)
{
    return library[i];
}

void DimerLibArray::alloc_lib(int n, int a_n, int e_n)
{
    library[iterator] = new DimerLib(n, a_n, e_n);
}

void DimerLibArray::DimerLibArray::add_to_atom_info(char *N, int i, char r, int p)
{
    library[iterator]->atom_data->add_atom(N, i, r, p);
}

void DimerLibArray::add_lib(gsl_matrix **d_m, float *e, char *n)
{
    library[iterator]->save_lib(d_m, e, n);
    iterator++;
}

void DimerLibArray::reset_flags(bool *reset)
{
    for (int i = 0; i < count; i++)
    {
        if (reset[i] == true)
        {
            library[i]->clear_flags();
            // printf("flags for %d reset\n", i);
        }
    }
}

void DimerLibArray::print_dimer_info(int i)
{
    printf("DNT: %s, # Models: %d, # Atoms Per Model: %d\n", library[i]->name, library[i]->count, library[i]->atom_data->count);
}

void DimerLibArray::print_matrix(int i, int j)
{
    for (unsigned int k = 0; k < library[i]->data_matrices[j]->size1; k++)
    {
        for (int l = 0; l < 3; l++)
        {
            l == 0 ? printf("%s %3.2f ", library[i]->atom_data->name[k], gsl_matrix_get(library[i]->data_matrices[j], k, l)) : printf("%3.2f ", gsl_matrix_get(library[i]->data_matrices[j], k, l));
        }
        printf("\n");
    }
}

gsl_matrix *DimerLibArray::get_matrix(int s, int index)
{
    return library[s]->data_matrices[index];
}

void get_model_count(FILE *fp, int *i)
{
    char line[100];
    while (fgets(line, sizeof(line), fp))
    {
        char *header;
        header = strtok(line, " ");
        if (!strcmp(header, "MODEL"))
            i[0]++;
        if (!strcmp(header, "ATOM"))
            i[1]++;
    }
    i[1] /= i[0];
    rewind(fp);
}

void calculate_dnt_COM(gsl_matrix *A, atom_info *A_info)
{
    int n_r1 = A_info->count_per_res[0];
    int n_r2 = A_info->count_per_res[1];

    gsl_matrix *A_1 = gsl_matrix_alloc(n_r1, 3);
    gsl_matrix *A_2 = gsl_matrix_alloc(n_r2, 3);

    double COMA_1[3] = {0, 0, 0};
    double COMA_2[3] = {0, 0, 0};

    for(int i = 0; i < n_r1; i++)
    {
        gsl_matrix_set(A_1, i, 0, gsl_matrix_get(A, i, 0));
        gsl_matrix_set(A_1, i, 1, gsl_matrix_get(A, i, 1));
        gsl_matrix_set(A_1, i, 2, gsl_matrix_get(A, i, 2));
    }

    for(int i = n_r1; i < n_r1 + n_r2; i++)
    {
        gsl_matrix_set(A_2, i - n_r1, 0, gsl_matrix_get(A, i, 0));
        gsl_matrix_set(A_2, i - n_r1, 1, gsl_matrix_get(A, i, 1));
        gsl_matrix_set(A_2, i - n_r1, 2, gsl_matrix_get(A, i, 2));
    }

    get_matrix_COM(A_1, COMA_1);
    get_matrix_COM(A_2, COMA_2);

    gsl_matrix_set(A, n_r1 + n_r2, 0, COMA_1[0]);
    gsl_matrix_set(A, n_r1 + n_r2, 1, COMA_1[1]);
    gsl_matrix_set(A, n_r1 + n_r2, 2, COMA_1[2]);

    gsl_matrix_set(A, n_r1 + n_r2 + 1, 0, COMA_2[0]);
    gsl_matrix_set(A, n_r1 + n_r2 + 1, 1, COMA_2[1]);
    gsl_matrix_set(A, n_r1 + n_r2 + 1, 2, COMA_2[2]);

    gsl_matrix_free(A_1);
    gsl_matrix_free(A_2);
}

void load_libs(char **LibNames, int N_diNts, DimerLibArray &RTN, bool for_WC)
{
    enum
    {
        model_count = 0,
        atom_count = 1
    };
    float *energies;
    gsl_matrix **data_mats;
    int model_info[2];

    char tmp[3];
    char line[100];

    RTN.initialize(N_diNts);

    for (int i = 0; i < N_diNts; i++)
    {
        bool first_itr = true;
        model_info[model_count] = model_info[atom_count] = 0;

        FILE *LibFile = fopen(LibNames[i], "r");
        if (LibFile == NULL)
        {
            printf("Cannot open library file: %s\n", LibNames[i]);
            exit(3);
        }
        printf("%d: %s ", i, LibNames[i]);

        get_model_count(LibFile, model_info);
        printf("Models: %d, Atoms per model: %d\n", model_info[model_count], model_info[atom_count]);

        model_info[atom_count] += 2; //Plus 2 B/C COM for each DNT will be included in data matrix
        
        RTN.alloc_lib(model_info[model_count], model_info[atom_count], 2); 
        data_mats = (gsl_matrix **)malloc(sizeof(gsl_matrix *) * model_info[model_count]);
        
        for (int i = 0; i < model_info[model_count]; i++)
        {
            data_mats[i] = gsl_matrix_alloc(model_info[atom_count], 3);
        }

        energies = (float *)malloc(sizeof(float) * model_info[model_count]);
        
        int iterator = 0;
        int row = 0;
        while (fgets(line, sizeof(line), LibFile))
        {
            char line_origin[100];
            char *header;

            strncpy(line_origin, line, 99);
            header = strtok(line, " ");

            if (!strcmp(header, "ENERGY"))
            {
                char *str1;
                str1 = strtok(NULL, " ");
                str1 = strtok(NULL, " ");
                energies[iterator] = atof(str1);
            }
            if ((!strcmp(header, "ATOM")))
            {
                // printf("Iterator = %d\n", iterator);
                char *index, *name, *residue, *position, *X, *Y, *Z;
                index = strtok(NULL, " ");
                name = strtok(NULL, " ");
                residue = strtok(NULL, " ");
                position = strtok(NULL, " ");
                X = strtok(NULL, " ");
                Y = strtok(NULL, " ");
                Z = strtok(NULL, " ");

                // printf("name = %s\n", name);

                if (first_itr)
                    RTN.add_to_atom_info(name, atoi(index), *residue, atoi(position));

                gsl_matrix_set(data_mats[iterator], row, 0, atof(X));
                gsl_matrix_set(data_mats[iterator], row, 1, atof(Y));
                gsl_matrix_set(data_mats[iterator], row, 2, atof(Z));
                row++;
            }
            else if (!strcmp(header, "ENDMDL\n"))
            {
                if (first_itr)
                    first_itr = false;

                calculate_dnt_COM(data_mats[iterator], RTN[RTN.iterator]->atom_data);
                iterator++;
                row = 0;
            }
        }
        fclose(LibFile);

        const char *default_name;
        for_WC ? default_name = WATSON_CRICK_LIBRARY_PROTOTYPE : default_name = LIBRARY_FILENAME_PROTOTYPE;
        int name_position = -1;
        while (default_name[name_position++] != 'X')
            ;
        name_position -= 1;
        tmp[0] = for_WC ? LibNames[i][name_position] : LibNames[i][name_position];
        tmp[1] = for_WC ? LibNames[i][name_position + 1] : LibNames[i][name_position + 1];
        tmp[2] = '\0';

        RTN.add_lib(data_mats, energies, tmp);
        free(data_mats);
        free(energies);
    }
}
