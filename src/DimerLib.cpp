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
    LibraryMemblock = gsl_block_alloc(n * a_n * MATRIX_DIMENSION2);
    data_matrices = (gsl_matrix **)malloc(sizeof(gsl_matrix *) * n);
    for(int i = 0, offset = 0; i < n; i++, offset += a_n * MATRIX_DIMENSION2)
    {
        data_matrices[i] = gsl_matrix_alloc_from_block(LibraryMemblock, offset, a_n, MATRIX_DIMENSION2, MATRIX_DIMENSION2);
    }

    name = (char *)malloc(sizeof(char) * 3);
    energy = (float *)malloc(sizeof(float) * n);
    radii[0] = (double *)malloc(sizeof(double) *n);
    radii[1] = (double *)malloc(sizeof(double) *n);
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
    gsl_block_free(LibraryMemblock);
    free(energy);
    free(data_matrices);
    free(name);
    free(radii[0]);
    free(radii[1]);
    delete atom_data;
}

void DimerLib::save_lib(gsl_matrix **d_m, float *e, char *n, double** r)
{
    for (int i = 0; i < count; i++)
    {
        data_matrices[i] = d_m[i];
    }
    memcpy(energy, e, sizeof(float) * count);
    memcpy(name, n, sizeof(char) * 3);
    memcpy(radii[0], r[0], sizeof(double) * count);
    memcpy(radii[1], r[1], sizeof(double) * count);
    //for(int i = 0; i < count; i++)
    //{
     //   printf("%d  %s: 1:%f 2:%f\n", i, name, radii[0][i], radii[1][i]);
    //}
    //putchar('\n');
}

void DimerLib::clear_flags()
{
}

DimerLibArray::DimerLibArray()
{
    // Do nothing
}

DimerLibArray::DimerLibArray(int s)
{
    library = (DimerLib **)malloc(sizeof(DimerLib *) * s);
    Flags = (flag**)malloc(sizeof(flag) * s);
    count = s;
    iterator = 0;
    was_initialized = true;
}

DimerLibArray::~DimerLibArray()
{
    if (was_initialized)
    {
        for (uint64_t i = 0; i < count; i++)
        {
            if(is_duplicate[i] != true)
            {
                delete library[i];
            }            
        }
        free(library);
        free(is_duplicate);
    }
}

void DimerLibArray::initialize(int s)
{
    //printf("DIMERLIBARRAY:INIT s = %d\n", s);
    library = (DimerLib **)malloc(sizeof(DimerLib *) * s);
    is_duplicate = (bool*)malloc(sizeof(bool) * s);
    memset(is_duplicate, false, sizeof(bool) * s);
    Flags = (flag **)malloc(sizeof(flag *) * s);
    count = s;
    iterator = 0;
    was_initialized = true;
}

DimerLib *DimerLibArray::operator[](int i)
{
    //printf("DIMERLIBARRAY:ACCESS i = %d\n", i);
    return library[i];
}

void DimerLibArray::map_duplicate(size_t dupli, size_t orig)
{
    library[dupli] = library[orig];
    is_duplicate[dupli] = true;
    iterator++;
}

void DimerLibArray::alloc_lib(size_t n, size_t a_n, size_t e_n)
{
    library[iterator] = new DimerLib(n, a_n, e_n);
    full_structure_element_sum += a_n;
    if((a_n - e_n) > LargestAtomCount)
    {
        LargestAtomCount = (a_n - e_n);
    }
}

void DimerLibArray::DimerLibArray::add_to_atom_info(char *N, int i, char r, int p)
{
    library[iterator]->atom_data->add_atom(N, i, r, p);
}

void DimerLibArray::add_lib(gsl_matrix **d_m, float *e, char *n, double** r)
{
    library[iterator]->save_lib(d_m, e, n, r);
    iterator++;
}

void DimerLibArray::get_charged_atom_map()
{
    int LargestAtomCount = 0;
    int MapTracker = 0;
    for(uint64_t i = 0; i < count; i++)
    {
        if(library[i]->atom_data->count > LargestAtomCount)
        {
            LargestAtomCount = library[i]->atom_data->count;
        }        
        
    }
    AtomMap = (int*)calloc(LargestAtomCount * count, sizeof(int));
    for(uint64_t i = 0; i < count; i++)
    {
        for(int j = 0; j < LargestAtomCount; j++)
        {
            if(j < library[i]->atom_data->count)
            {
                if(library[i]->atom_data->charges[j] != atom_charge::NEUTRAL)
                {
                    AtomMap[IDX_FLAT2D(i,j,LargestAtomCount)] = MapTracker;
                    MapTracker++;
                }
                else
                {
                    AtomMap[IDX_FLAT2D(i,j,LargestAtomCount)] = -1;
                }
                //printf("Lib: %d, Row: %d, Index: %d, Atom: %s -> %d\n", i, j, IDX_FLAT2D(i, j, LargestAtomCount), library[i]->atom_data->name[j], AtomMap[IDX_FLAT2D(i,j,LargestAtomCount)]);
            }
            else
            {
                //printf("Lib: %d, Row: %d, Index: %d, Atom: %s\n", i, j, IDX_FLAT2D(i, j, LargestAtomCount), "DUMMY");
            }           
        }
    }
    ChargedAtomCount = LargestAtomCount;
}

void DimerLibArray::initialize_flags()
{
    int TotalSum = 0;
    int offset = 0;
    for (uint64_t i = 0; i < count; i++)
    {        
        TotalSum += library[i]->count;
    }
    Flags[0] = (flag*)malloc(sizeof(flag) * TotalSum);
    memset(&Flags[0][0], flag::NO_FLAG, sizeof(flag) * TotalSum);
    for (uint64_t i = 1; i < count; i++)
    {
        offset += library[i - 1]->count;
        Flags[i] = &Flags[0][offset];
    }
}

void DimerLibArray::reset_flags(bool *reset)
{
    for (uint64_t i = 0; i < count; i++)
    {
        if (reset[i] == true)
        {
            for(int j = 0; j < library[i]->count; j++)
            {
                Flags[i][j] = flag::NO_FLAG;
            }
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

    //printf("1: %f, %f, %f\n", COMA_1[0], COMA_1[1], COMA_1[2]);
    //printf("2: %f, %f, %f\n", COMA_2[0], COMA_2[1], COMA_2[2]);
    //exit(0);

    gsl_matrix_free(A_1);
    gsl_matrix_free(A_2);
}

/**
 * @brief Calculates Steric Clash COM (SCC) radius for SC Optimization. First the distance between the phosphate oxygens (OP1 and OP2) and the COM is calculated.
 * Then the same is done but for the functional groups of the residue (N6, N2, O6, O2, N4, O4). The greatest distance between COM and either functional group or
 * residue used to set the radius of the COM sphere. (Radius of the COM sphere = distance + 1)
 * 
 * @param A Matrix of DNT
 * @param A_info Atom Info of DNT
 * @param radius Radius for residue
 * @param res_id Index of residue 
 */
void calculate_SCC_radii(gsl_matrix *A, atom_info *A_info, double* radius, int res_id)
{
    double dists_OP[] = {0, 0};
    double dists_FG[] = {0, 0};

    int FG_idxs[] = {-1, -1};

    char res_name = A_info->residue[res_id == 0 ? 0 : A_info->count_per_res[0]];

    gsl_vector_view OP_vec;
    gsl_vector_view FG_vec;
    gsl_vector_view COM_vec;

    *radius = 0;

    COM_vec = gsl_matrix_row(A, A->size1 - (2 - res_id));

    OP_vec = gsl_matrix_row(A, A_info->get_idx_of_atom(OP1, 1));
    dists_OP[0] = distance_vec2vec(&COM_vec.vector, &OP_vec.vector);
    
    OP_vec = gsl_matrix_row(A, A_info->get_idx_of_atom(OP2, 1));
    dists_OP[1] = distance_vec2vec(&COM_vec.vector, &OP_vec.vector);

    if(dists_OP[1] > dists_OP[0])
    {
        dists_OP[0] = dists_OP[1]; //I only care about the largest value
    }

    switch (res_name)
    {
    case 'A':
        FG_idxs[0] = A_info->get_idx_of_atom(N6, res_id);
        break;
    case 'C':
        FG_idxs[0] = A_info->get_idx_of_atom(N4, res_id);
        FG_idxs[1] = A_info->get_idx_of_atom(O2, res_id);
        break;
    case 'G':
        FG_idxs[0] = A_info->get_idx_of_atom(N2, res_id);
        FG_idxs[1] = A_info->get_idx_of_atom(O4, res_id);
        break;
    case 'U':
        FG_idxs[0] = A_info->get_idx_of_atom(O4, res_id);
        FG_idxs[1] = A_info->get_idx_of_atom(O2, res_id);
        break;
    }

    FG_vec = gsl_matrix_row(A, FG_idxs[0]);
    dists_FG[0] = distance_vec2vec(&COM_vec.vector, &FG_vec.vector);
    if(FG_idxs[1] != -1)
    {
        FG_vec = gsl_matrix_row(A, FG_idxs[1]);
        dists_FG[1] = distance_vec2vec(&COM_vec.vector, &FG_vec.vector);
    }

    if(dists_FG[1] > dists_FG[0])
    {
        dists_FG[0] = dists_FG[1]; //I only care about the largest value
    }

    if(dists_FG[0] > dists_OP[0])
    {
        *radius = dists_FG[0] + 1;
    }
    else
    {
        *radius = dists_OP[0] + 1;
    }
}

void check_if_all_in_sphere(gsl_matrix *A, atom_info *A_info, double* radius, int res_id)
{
    gsl_vector_view V, COM;
    double dist;

    COM = gsl_matrix_row(A, A->size1 - (2 - res_id));

    for(unsigned int i = (res_id == 0 ? 0 : A_info->count_per_res[0]); i < (res_id == 0 ? A_info->count_per_res[0] : A_info->count_per_res[1]);  i++)
    {
        V = gsl_matrix_row(A, i);
        dist = distance_vec2vec(&COM.vector, &V.vector);
        if(dist > *radius)
        {
            printf("%s is outside COM sphere!\n", A_info->name[i]);
            *radius = dist + 1;
        }
    }
}

void load_libs(char **LibNames, int N_diNts, DimerLibArray &RTN, int* duplicate_record, bool for_WC)
{
    enum { model_count = 0, atom_count = 1 };
    float *energies;
    gsl_matrix **data_mats;
    int model_info[2];
    double* _radii[2];
    char line[100];

    RTN.initialize(N_diNts);

    for (int i = 0; i < N_diNts; i++)
    {
        if(!for_WC && duplicate_record[i] != -1)
        {
            RTN.map_duplicate(i, duplicate_record[i]);
            printf("%d: %s (Points to library %d)\n", i, LibNames[i], duplicate_record[i]);
            continue;
        }
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
        
        data_mats = RTN[i]->data_matrices;
        energies  = RTN[i]->energy;
        _radii[0] = RTN[i]->radii[0];
        _radii[1] = RTN[i]->radii[1];
        
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
                {
                    RTN.add_to_atom_info(name, atoi(index), *residue, atoi(position));
                }

                gsl_matrix_set(data_mats[iterator], row, 0, atof(X));
                gsl_matrix_set(data_mats[iterator], row, 1, atof(Y));
                gsl_matrix_set(data_mats[iterator], row, 2, atof(Z));
                row++;
            }
            else if (!strcmp(header, "ENDMDL\n"))
            {
                if (first_itr)
                    first_itr = false;
                
                if(!for_WC)
                {
                    calculate_dnt_COM(data_mats[iterator], RTN[RTN.iterator]->atom_data);
                    calculate_SCC_radii(data_mats[iterator], RTN[RTN.iterator]->atom_data, &_radii[0][iterator], 0);
                    calculate_SCC_radii(data_mats[iterator], RTN[RTN.iterator]->atom_data, &_radii[1][iterator], 1);
                    check_if_all_in_sphere(data_mats[iterator], RTN[RTN.iterator]->atom_data, &_radii[0][iterator], 0);
                    check_if_all_in_sphere(data_mats[iterator], RTN[RTN.iterator]->atom_data, &_radii[1][iterator], 1);
                }                
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
        RTN[i]->name[0] = for_WC ? LibNames[i][name_position] : LibNames[i][name_position];
        RTN[i]->name[1] = for_WC ? LibNames[i][name_position + 1] : LibNames[i][name_position + 1];
        RTN[i]->name[2] = '\0';

        for(int i = 0; i < model_info[model_count]; i++)
        {
            //printf("%d  %s: 1:%f 2:%f\n", i, tmp, _radii[0][i], _radii[1][i]);
        }
        //putchar('\n'); 
        RTN.iterator++;    
    }
    RTN.initialize_flags();
}
