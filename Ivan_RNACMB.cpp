//First attempt for Combinatorial building RNA  **03/5/21**
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>
#include <ctime>
#include <vector>
using namespace std;

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <sys/stat.h>
#include <sys/types.h>

#define LIBRARY_FILENAME_STANDARD "XX_library_combined.txt"
#define idx(i,j,s) (i * s + j)
#define square(a) (a * a)

#define DEBUG_SWITCH false
#define DEBUG(a) if(DEBUG_SWITCH == true) a;

#define RADIUS_C 1.00
#define RADIUS_N 1.00
#define RADIUS_P 1.00
#define RADIUS_O 1.00

#define RMSD_LIMIT 0.5

#define MAX_STRINGS 100

long int rna_dat_tracker = 0;

enum flag{NO_FLAG, USED, NOT_USABLE};

struct atom_info
{
    char**  name;
    int*    index;
    char*   residue;
    int*    dnt_pos;

    int     count;      //Number of atoms per structure
    int     iterator;

    atom_info()
    {
        iterator = 0;
    }

    atom_info(int n)
    {
        name    = (char**)malloc(sizeof(char*) * n);
        index   = (int*)  malloc(sizeof(int)   * n);
        residue = (char*) malloc(sizeof(char)  * n);
        dnt_pos = (int*)  malloc(sizeof(int)   * n);
        count   = n;
        for(int i = 0; i < n; i++)
        {
            name[i] = (char*)malloc(sizeof(char) * 5);
        }
        iterator = 0;
    }

    atom_info(char** n, int* in, char* r, int* d, int c)
    {
        name = (char**)malloc(sizeof(char*) * c);
        for(int i = 0; i < c; i++)
        {
            name[i] = (char*)malloc(sizeof(char) * 5);
            strcpy(name[i], n[i]);
        }

        index       =   (int*) malloc(sizeof(int)  * c);
        residue     =   (char*)malloc(sizeof(char) * c);
        dnt_pos     =   (int*) malloc(sizeof(int)  * c);
        
        memcpy(index,   in, sizeof(int)  * c);
        memcpy(residue, r,  sizeof(char) * c);
        memcpy(dnt_pos, d,  sizeof(int)  * c);
        
        count = c;
    }

    void initialize_atom_info(int n)
    {
        name    = (char**)malloc(sizeof(char*) * n);
        index   = (int*)  malloc(sizeof(int)   * n);
        residue = (char*) malloc(sizeof(char)  * n);
        dnt_pos = (int*)  malloc(sizeof(int)   * n);
        count   = n;
        for(int i = 0; i < n; i++)
        {
            name[i] = (char*)malloc(sizeof(char) * 5);
        }
        iterator = 0;
    }

    ~atom_info()
    {
        //printf("Destructor Called!\n");
        for(int i = 0; i < count; i++)
        {
            free(name[i]);
        }
        free(index);
        free(residue);
        free(dnt_pos);
        free(name);
    }

    void add_atom(char *N, int i, char r, int p)
    {  
        strcpy(name[iterator], N);
        index  [iterator] = i;
        residue[iterator] = r;
        dnt_pos[iterator] = p;
        iterator++;
    }

    void clear()
    {
        iterator = 0;
    }

    atom_info(const atom_info &A)
    {
        name = (char**)malloc(sizeof(char*) * A.count);
        for(int i = 0; i < A.count; i++)
        {
            name[i] = (char*)malloc(sizeof(char) * 5);
            strcpy(name[i], A.name[i]);
        }

        index       =   (int*) malloc(sizeof(int)  * A.count);
        residue     =   (char*)malloc(sizeof(char) * A.count);
        dnt_pos     =   (int* )malloc(sizeof(int)  * A.count);
        
        memcpy(index,   A.index,    sizeof(int)  * A.count);
        memcpy(residue, A.residue,  sizeof(char) * A.count);
        memcpy(dnt_pos, A.dnt_pos,  sizeof(int)  * A.count);
        
        count = A.count;
    }

    atom_info& operator=(const atom_info& A)
    {
        name = (char**)malloc(sizeof(char*) * A.count);
        for(int i = 0; i < A.count; i++)
        {
            name[i] = (char*)malloc(sizeof(char) * 5);
            strcpy(name[i], A.name[i]);
        }

        index       =   (int*) malloc(sizeof(int)  * A.count);
        residue     =   (char*)malloc(sizeof(char) * A.count);
        dnt_pos     =   (int* )malloc(sizeof(int)  * A.count);
        
        memcpy(index,   A.index,    sizeof(int)  * A.count);
        memcpy(residue, A.residue,  sizeof(char) * A.count);
        memcpy(dnt_pos, A.dnt_pos,  sizeof(int)  * A.count);
        
        count = A.count;
        //printf("Operator= Called!\n");
    }

    atom_info(atom_info &&A) noexcept
    {
        name    = A.name;
        index   = A.index;
        residue = A.residue;
        dnt_pos = A.dnt_pos;
        count   = A.count;
        //printf("Move Called\n");
    }

    void move(atom_info& A)
    {
        name    = A.name;
        index   = A.index;
        residue = A.residue;
        dnt_pos = A.dnt_pos;
        count   = A.count;
    }

    void print_at(int line)
    {
        printf("%-4s %-4d %-4c %-4d ", name[line], index[line], residue[line], dnt_pos[line]);
    }
};

struct DimerLib
{
    gsl_matrix**      data_matrices;
    atom_info*        atom_data;
    double*           energy;
    char*             name;
    int               count; //Number of structures in Library
    
    flag*              flags;

    DimerLib(int n, int a_n)
    {
        data_matrices = (gsl_matrix**)malloc(sizeof(gsl_matrix*) * n);
        name          = (char*)       malloc(sizeof(char)        * 3);
        flags         = (flag*)calloc(n, sizeof(flag));
        count         = n;
        atom_data     = new atom_info(a_n);
    }
    ~DimerLib()
    {
        //printf("DimerLib %s Destructor Called\n", name);
        for(int i = 0; i < count; i++)
        {
            gsl_matrix_free(data_matrices[i]);
        }
        free(energy);
        free(data_matrices);
        free(flags);
        free(name);
        delete atom_data;
    }
    
    void save_lib(gsl_matrix** d_m, double* e, char* n)
    {
        for(int i = 0; i < count; i++)
        {
            data_matrices[i] = d_m[i];
        }
        energy = e;
        memcpy(name, n, sizeof(char) * 3);
    }

    void clear_flags()
    {
        for(int i = 0; i < count; i++)
        {
            flags[i] = NO_FLAG;
        }
    }
};

struct DimerLibArray
{    
    DimerLib** library;
    int count;          //Number of "Libraries" in array

    int iterator;

    DimerLibArray(int s)
    {   
        library = (DimerLib**)malloc(sizeof(DimerLib*) * s);
        count = s;
        iterator = 0;
    }

    ~DimerLibArray()
    {
        for(int i = 0; i < count; i++)
            delete library[i];
        free(library);
    }

    DimerLib* operator[](int i)
    {
        return library[i];
    }
    
    void alloc_lib(int n, int a_n)
    {
        library[iterator] = new DimerLib(n, a_n);
    }

    void add_to_atom_info(char *N, int i, char r, int p)
    {
        library[iterator]->atom_data->add_atom(N, i, r, p);
    }

    void add_lib(gsl_matrix** d_m, double* e, char* n)
    {
        library[iterator]->save_lib(d_m, e, n);
        iterator++;
    }

    void reset_flags(bool *reset)
    {
        for(int i = 0; i < count; i++)
        {
            if(reset[i] == true)
            {
                library[i]->clear_flags();
                printf("flags for %d reset\n", i);
            }
        }
    }

    void print_dimer_info(int i)
    {
        printf("DNT: %s, # Models: %d, # Atoms Per Model: %d\n", library[i]->name, library[i]->count, library[i]->atom_data->count);
    }

    void print_matrix(int i, int j)
    {
        for(int k = 0; k < library[i]->data_matrices[j]->size1; k++)
        {
            for(int l = 0; l < 3; l++)
            {
                l == 0 ? printf("%s %3.2f ", library[i]->atom_data->name[k], gsl_matrix_get(library[i]->data_matrices[j], k, l)) : printf("%3.2f ", gsl_matrix_get(library[i]->data_matrices[j], k, l));
            }
            printf("\n");
        }
    }
    
    gsl_matrix *get_matrix(int s, int index)
    {
        return library[s]->data_matrices[index];
    }
};

struct CMB_Manager
{
    int *count_per_lib;         //Number of structures in each library
    bool **attach_attempted;    //Tracks if all structures in each respective library has been tested
    int last_attempted[2];

    bool *libs_completed;
    
    int strs_built;

    int count;                  //For deallocation

    CMB_Manager(DimerLibArray& LA)
    {
        count = LA.count;

        count_per_lib = (int *)malloc(sizeof(int) * LA.count);
        attach_attempted = (bool **)malloc(sizeof(bool *) * LA.count);

        for(int i = 0; i < LA.count; i++)
        {
            count_per_lib[i] = LA[i]->count;
            attach_attempted[i] = (bool *)calloc(LA[i]->count, sizeof(bool));
        }

        last_attempted[0] = 0;
        last_attempted[1] = 0;

        strs_built = 0;
        libs_completed = (bool *)calloc(count, sizeof(bool));
    }

    ~CMB_Manager()
    {
        for(int i = 0; i < count; i++)
        {
            free(attach_attempted[i]);
        }
        free(attach_attempted);
        free(count_per_lib);
        free(libs_completed);
    }

    void attach_attempt(int i, int j)
    {
        attach_attempted[i][j] = true;
        last_attempted[0] = i;
        last_attempted[1] = j;
    }

    bool is_at_end()
    {
        if(last_attempted[1] == count_per_lib[last_attempted[0]] - 1)
        {
            return true;
        }
        return false;
    }

    void check_lib_completion()
    {
        int counter = 0;

        int first_completed = last_attempted[0];

        for(int i = last_attempted[0]; i >= 0; i--)
        {
            libs_completed[i] = false;
            if(attach_attempted[i][count_per_lib[i] - 1] == true)
            {
                libs_completed[i] = true;
                counter++;
                if(i < first_completed)
                    first_completed = i;
                printf("Lib %d complete! Max = %d\n", i, count_per_lib[i] - 1);
            }
            else
            {
                printf("Lib %d not complete! Max = %d\n", i, count_per_lib[i] - 1);
            }
        }
        if(counter != (count - first_completed))
        {
            bool swap = true;
            for(int i = last_attempted[0]; i >= first_completed; i--)
            {
                if(libs_completed[i] == false)
                {
                    swap = false;
                }
                libs_completed[i] = swap;
            }
        }
        
        if(counter == last_attempted[0] + 1)
            ;
        else
            libs_completed[0] = false;
        return;
    }

    void clear_attempts()
    {
        for(int i = 0; i < count; i++)
        {
            if(libs_completed[i] == true)
            {
                printf("reseting attempts for %d\n", i);
                for(int j = 0; j < count_per_lib[i]; j++)
                {
                    attach_attempted[i][j] = false;
                }
            }
        }
    }

    int get_reset_count()
    {
        int n_reset = 0;
        for(int i = 0; i < last_attempted[0] + 1; i++)
        {
            if(libs_completed[i] == true)
            {
                n_reset++;
            }
        }
        //printf("n_reset = %d\n", n_reset);
        return n_reset;
    }

    void successful_construction()
    {
        strs_built++;
    }

};

struct output_string
{
    FILE *output_file;
    char **string_storage;
    int iterator;
    int num_allocated;
    int max_string;

    bool all_init = false;

    output_string(const char *F, int max_s)
    {
        
        output_file = fopen(F, "w");
        string_storage = (char **)malloc(sizeof(char *) * max_s);
        iterator = 0;
        max_string = max_s;
        num_allocated = 0;
    }
    ~output_string()
    {
        for(int i = 0; i < max_string; i++)
        {
            if(iterator < max_string)
            {
                fprintf(output_file, "%s", string_storage[i]);
            }
            free(string_storage[i]);
        }
        free(string_storage);
        fclose(output_file);
    }

    void add_string(char *s, int max_atoms)
    {
        if(iterator == max_string - 1)
        {
            if(all_init)
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
                //printf("itr = %d\n", iterator);
                iterator = 0;
            }
            for(int i = 0; i < max_string; i++)
            {
                fprintf(output_file, "%s", string_storage[i]);
            }
            return;
        }
        else
        {
            if(all_init)
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
};

struct RNA_data
{
    gsl_matrix *data_matrix;
    gsl_matrix **submatrices; //Submatrix for each residue corresponding to the target atoms. Ordered 5' to 3'.
    float energy;
    atom_info *atom_data;
    char *name;
    flag *_flag;
    int count; //Number of atoms in structure
    
    int position_in_lib[2];
    int count_per_sub[2];
    int sub_starts_at[2];

    int position_max;

    char** target;
    
    long int id;

    RNA_data(DimerLibArray& L, int i, int j)
    {
        atom_data = (L[i]->atom_data);
        data_matrix = gsl_matrix_alloc(L[i]->data_matrices[j]->size1, L[i]->data_matrices[j]->size2);
        gsl_matrix_memcpy(data_matrix, L[i]->data_matrices[j]);
        energy = L[i]->energy[j];
        name = (char *)malloc(sizeof(char) * 3);
        memcpy(name, L[i]->name, sizeof(char) * 3);
        count = L[i]->atom_data->count;
        _flag = &(L[i]->flags[j]);
        position_in_lib[0] = i;
        position_in_lib[1] = j;
        position_max = L[i]->count - 1;
        make_submatrices();
        id = rna_dat_tracker++;
        DEBUG(printf("ID: %ld allocated\n", id));
    }

    void overwrite(DimerLibArray& L, int i, int j)
    {
        gsl_matrix_memcpy(data_matrix, L[i]->data_matrices[j]);
        energy = L[i]->energy[j];
        memcpy(name, L[i]->name, sizeof(char) * 3);
        count = L[i]->atom_data->count;
        _flag = &(L[i]->flags[j]);
        position_in_lib[0] = i;
        position_in_lib[1] = j;
        update_submatrices();
        //printf("ID: %d %d allocated\n", position_in_lib[0], position_in_lib[1]);
    }

    ~RNA_data()
    {
        DEBUG(printf("ID: %ld deallocated\n", id));
        
        gsl_matrix_free(data_matrix);
        free(name);
        gsl_matrix_free(submatrices[0]);
        gsl_matrix_free(submatrices[1]);
        free(submatrices);
    }
    
    void make_submatrices()
    {
        submatrices = (gsl_matrix**)malloc(sizeof(gsl_matrix*) * 2);

        for(int i = 0; i < 2; i++)
        {
            get_target(i);
            if(name[i] == 'A' || name[i] == 'G')
                count_per_sub[i] = 14; //Number of atoms in target for A & G
            else
                count_per_sub[i] = 11; //Number of atoms in target for C & U
            
            submatrices[i] = gsl_matrix_alloc(count_per_sub[i], 3);

            int rel_idx = 0;

            for(int j = 0; j < count; j++)
            {
                for(int k = 0; k < count_per_sub[i]; k++)
                {
                    if(!strcmp(target[k], atom_data->name[j]) && (atom_data->dnt_pos[j] - 1) == i)
                    {
                        //printf("target: %s  name: %s \n", target[k], atom_data->name[j]);
                        //printf("Rel_idx: %d\n", rel_idx);
                        //printf("dnt_pos: %d\n", atom_data->dnt_pos[j]);
                        gsl_matrix_set(submatrices[i], rel_idx, 0, gsl_matrix_get(data_matrix, j, 0));
                        gsl_matrix_set(submatrices[i], rel_idx, 1, gsl_matrix_get(data_matrix, j, 1));
                        gsl_matrix_set(submatrices[i], rel_idx, 2, gsl_matrix_get(data_matrix, j, 2));
                        rel_idx++;
                    }
                }
            }
            free_target(count_per_sub[i]);
        }
    }

    void update_submatrices() /* Need to be rewritten for efficiency */
    {
        gsl_matrix_free(submatrices[0]);
        gsl_matrix_free(submatrices[1]);
        free(submatrices);
        make_submatrices();
    }

    void free_target(int size)
    {
        for(int z = 0; z < size; z++)
            free(target[z]);
        free(target);
    }

    gsl_matrix* get_target_matrix(int res)
    {
        return submatrices[res];
    }

    gsl_matrix* get_target_matrix_copy(int res)
    {
        gsl_matrix *rtn = gsl_matrix_alloc(count_per_sub[res], 3);
        gsl_matrix_memcpy(rtn, submatrices[res]);
        return rtn;
    }

    void get_target(int res)
    {
        char targetA[][5] = {"N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4", "C1'", "C2'", "C3'", "C4'", "O4'" }; // 14
        char targetC[][5] = {"N1", "C2", "N3", "C4", "C5", "C6", "C1'", "C2'", "C3'", "C4'", "O4'"}; // 11
        char targetG[][5] = {"N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4", "C1'", "C2'", "C3'", "C4'", "O4'"}; // 14
        char targetU[][5] = {"N1", "C2", "N3", "C4", "C5", "C6", "C1'", "C2'", "C3'", "C4'", "O4'"}; // 11

        if( name[res] == 'A' ){
            target = (char **)malloc(sizeof(char *) * 14);
            for(int i = 0; i < 14; i++)
            {
                target[i] = (char *)malloc(sizeof(char) * 5);
                strcpy(target[i], targetA[i]);
            }
        }
        else if( name[res] == 'C' ){
            target = (char **)malloc(sizeof(char *) * 11);
            for(int i = 0; i < 11; i++)
            {
                target[i] = (char *)malloc(sizeof(char) * 5);
                strcpy(target[i], targetC[i]);
            }
        }
        else if( name[res] == 'G' ){
            target = (char **)malloc(sizeof(char *) * 14);
            for(int i = 0; i < 14; i++)
            {
                target[i] = (char *)malloc(sizeof(char) * 5);
                strcpy(target[i], targetG[i]);
            }
        }
        else if( name[res] == 'U' ){
            target = (char **)malloc(sizeof(char *) * 11);
            for(int i = 0; i < 11; i++)
            {
                target[i] = (char *)malloc(sizeof(char) * 5);
                strcpy(target[i], targetU[i]);
            }
        }
        else{
            printf("input file has problems, please check\n"); exit(0);
        }
    }

    RNA_data& operator=(RNA_data& A)
    {
        name = (char *)malloc(sizeof(char) * 3);
        strcpy(name, A.name);
        data_matrix = gsl_matrix_alloc(A.data_matrix->size1, A.data_matrix->size2);
        gsl_matrix_memcpy(data_matrix, A.data_matrix);
        atom_data = A.atom_data;
        
        energy = A.energy;
        _flag = A._flag;

        position_in_lib[0] = A.position_in_lib[0];
        position_in_lib[1] = A.position_in_lib[1];

        count = A.count;
    }

    void print()
    {
        for(int i = 0; i < data_matrix->size1; i++)
        {
            atom_data->print_at(i);
            for(int j = 0; j < data_matrix->size2; j++)
                printf("%-3.2f ", gsl_matrix_get(data_matrix,i,j));
            putchar('\n');
        }
    }
    void print(int res)
    {
        for(int i = 0; i < data_matrix->size1; i++)
        {
            if(atom_data->dnt_pos[i] != (res + 1))
                continue;
            atom_data->print_at(i);
            for(int j = 0; j < data_matrix->size2; j++)
                printf("%-3.2f ", gsl_matrix_get(data_matrix,i,j));
            putchar('\n');
        }
    }

    int to_string(char *s, int buffer_size, int string_index)
    {
        for(int i = 0; i < data_matrix->size1; i++)
        {
            string_index += snprintf(&s[string_index], buffer_size - string_index, "%-6s%5d %4s %3c  %4d    ", "ATOM",  atom_data->index[i], atom_data->name[i], atom_data->residue[i], atom_data->dnt_pos[i]);
            for(int j = 0; j < data_matrix->size2; j++)
            {
                string_index += snprintf(&s[string_index], buffer_size - string_index, "%8.3f", gsl_matrix_get(data_matrix,i,j));
            }
            string_index += snprintf(&s[string_index], buffer_size - string_index, "\n");
        }
        return string_index;
    }

    void print_offset(int res, int position)
    {
        for(int i = 0; i < data_matrix->size1; i++)
        {
            if(atom_data->dnt_pos[i] != (res + 1))
                continue;
            printf("%-4s %-4d %-4c %-4d ", atom_data->name[i], atom_data->index[i], atom_data->residue[i], atom_data->dnt_pos[i] + position);
            for(int j = 0; j < data_matrix->size2; j++)
                printf("%-3.2f ", gsl_matrix_get(data_matrix,i,j));
            putchar('\n');
        }
    }

    int to_string_offset(int res, int position, char *s, int buffer_size, int string_index, int *idx_offset)
    {
        for(int i = 0; i < data_matrix->size1; i++)
        {
            if(atom_data->dnt_pos[i] != (res + 1))
                continue;
            *idx_offset += 1;
            string_index += snprintf(&s[string_index], buffer_size - string_index, "%-6s%5d %4s %3c  %4d    ", "ATOM",  *idx_offset, atom_data->name[i], atom_data->residue[i], atom_data->dnt_pos[i] + position);
            for(int j = 0; j < data_matrix->size2; j++)
                string_index += snprintf(&s[string_index], buffer_size - string_index, "%8.3f", gsl_matrix_get(data_matrix,i,j));
            string_index += snprintf(&s[string_index], buffer_size - string_index, "\n");
        }
        return string_index;
    }
};

struct RNA_data_array
{
    RNA_data **sequence;
    int count;
    int iterator;
    int iterator_max;
    
    char *string_out;
    int string_buffer;
    int string_index;
    bool string_initialized = false;

    int atom_sum;

    RNA_data_array(int size)
    {
        iterator_max = size - 1;
        sequence = (RNA_data**)malloc(sizeof(RNA_data*) * size);
        iterator = -1;
        count = 0;
    }
    ~RNA_data_array()
    {
        //printf("iterator is @ %d\n", iterator);
        for(int i = 0; i < count; i++)
        {
            delete sequence[i];
        }
        free(sequence);
        free(string_out);
    }
    RNA_data* operator[](int i)
    {
        return sequence[i];
    }
    void add_copy(RNA_data* A)
    {
        *sequence[++iterator] = *A;
        count++;
    }
    void add_move(RNA_data* A)
    {
        sequence[++iterator] = A;
        count++;
        *A->_flag = USED;
    }
    RNA_data* current()
    {
        return sequence[iterator];
    }

    bool is_complete()
    {
        return iterator == iterator_max ? true : false;
    }

    void rollback()
    {
        DEBUG(printf("attach: %ld deleted @ rollback @ %d\n", sequence[iterator]->id, iterator));
        delete sequence[iterator];
        iterator--;
        count--;
    }

    void safe_rollback()
    {
        if(sequence[iterator]->position_in_lib[1] == sequence[iterator]->position_max)
            ;
        else
            delete sequence[iterator];
        iterator--;
        count--;
    }

    void rollback_by(int amount)
    {
        for(int i = iterator; i > (iterator - (amount + 1)); i--)
        {
            DEBUG(printf("attach: %ld deleted @ rollback_by @ %d\n", sequence[i]->id, i));
            delete sequence[i];
        }
        iterator -= (amount + 1);
        count -= (amount + 1);
        //printf("moving to pos: %d\n", iterator);
    }

    bool is_empty()
    {
        if(count == 0)
            return true;
        return false;
    }

    void printall()
    {
        sequence[0]->print();
        for(int i = 1; i < count; i++)
        {
            sequence[i]->print_offset(1, i);
        }
    }

    void initialize_string()
    {
        get_atom_sum();
        printf("atomsum = %d\n", atom_sum);

        string_buffer = 54 * atom_sum;
        string_out = (char *)malloc(sizeof(char) * string_buffer);
        string_index = 0;

        string_initialized = true;
    }

    int get_atom_sum()
    {
        if(string_initialized) return atom_sum;
        atom_sum = 0;
        for(int i = 0; i < count; i++)
        {
            atom_sum += sequence[i]->count;
        }
        return atom_sum;
    }

    int out_string_header()
    {
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#INDEX ");
        for(int i = 0; i < count; i++)
        {
            string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%d ", sequence[i]->position_in_lib[1]);
        }
        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\n");
        return string_index;
    }

    char* to_string()
    {

        int idx_offset;

        if(string_initialized)
        {
            string_index = 0;
        }
        else
        {
            initialize_string();
        }

        string_index = out_string_header();
        string_index = sequence[0]->to_string(string_out, string_buffer, string_index);
        idx_offset = sequence[0]->count;

        for(int i = 1; i < count; i++)
        {
            string_index = sequence[i]->to_string_offset(1, i, string_out, string_buffer, string_index, &idx_offset);
            //printf("str_idx = %d\n", string_index);
        }

        string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\n");

        return string_out;
    }

    int *get_index()
    {
        int *ar = (int *)malloc(sizeof(int) * count);
        for(int i = 0; i < count; i++)
        {
            ar[i] = sequence[i]->position_in_lib[1];
        }
        return ar;
    }
};

char** get_diNt_names(char* sequence, int *N_diNts)
{
    const char default_name[] = LIBRARY_FILENAME_STANDARD;
    *N_diNts = strlen(sequence) - 1;
    char **rtn = (char **)malloc(sizeof(char *) * *N_diNts);
    for(int i = 0; i < *N_diNts; i++)
    {
        rtn[i] = (char *)malloc(sizeof(char) * sizeof(default_name));
        memcpy(rtn[i], default_name, sizeof(default_name));
        memcpy(rtn[i], sequence++, sizeof(char) * 2);
    }
    return rtn;
}

void gsl_matrix_print(gsl_matrix *M)
{
    for(int i = 0; i < M->size1; i++){
        for(int j = 0; j < M->size2; j++)
            printf("%3.2f ", gsl_matrix_get(M,i,j));
        putchar('\n');
    }

}

gsl_matrix *make_gsl_matrix(vector<double> data)
{
    gsl_matrix *rtn = gsl_matrix_alloc(data.size()/3, 3);
    for(int i = 0; i < data.size()/3; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            gsl_matrix_set(rtn, i, j, data[idx(i,j,3)]);
            DEBUG(printf("%3.2f ", gsl_matrix_get(rtn, i, j)));
        }
        DEBUG(printf("\n"));
    }
    return rtn;
}

void get_model_count(FILE* fp, int* i)
{
    char line[100];
    while(fgets(line, sizeof(line), fp))
    {
        char *header;
        header = strtok(line, " ");
        if(!strcmp(header, "MODEL"))
            i[0]++;
        if(!strcmp(header,"ATOM"))
            i[1]++;
    }
    i[1] /= i[0];
    rewind(fp);
}

DimerLibArray load_libs(char **LibNames, int N_diNts)
{
    enum {model_count = 0, atom_count = 1};
    DimerLibArray RTN(N_diNts);
    gsl_matrix*   matrix;
    double*       energies;
    gsl_matrix**  data_mats;
    int*          model_id;
    int           model_info[2];

    char          tmp[3];
    char          line[100];

    for(int i = 0; i < N_diNts; i++)
    {
        bool first_itr = true;
        model_info[model_count] = model_info[atom_count] = 0;
        
        FILE *LibFile = fopen(LibNames[i], "r");
        if(LibFile == NULL)
        {
            printf("Cannot open library file: %s\n", LibNames[i]);
            exit(3);
        }
        printf("%d: %s ", i, LibNames[i]);


        get_model_count(LibFile, model_info);
        printf("Models: %d, Atoms per model: %d\n", model_info[model_count], model_info[atom_count]);
        
        RTN.alloc_lib(model_info[model_count], model_info[atom_count]);
        data_mats = (gsl_matrix**)malloc(sizeof(gsl_matrix*) * model_info[model_count]);
        for(int i = 0; i < model_info[model_count]; i++)
            data_mats[i] = gsl_matrix_alloc(model_info[atom_count], 3);
        energies = (double*)malloc(sizeof(double) * model_info[model_count]);

        int iterator = 0;
        int row = 0;
        while(fgets(line, sizeof(line), LibFile))
        {
            char line_origin[100];
            char *header;

            strncpy(line_origin, line, 99);
            header = strtok(line, " ");

            if(!strcmp(header, "MODEL") ){
                char *str1;
                str1 = strtok(NULL, " ");
            }

            if(!strcmp(header, "ENERGY") )
            {
                char *str1;
                str1 = strtok(NULL, " ");
                energies[iterator] = atof(str1);
                //printf("Iterator = %d\n", iterator);
            }         
            if( (!strcmp(header, "ATOM")))
            {
                char *index, *name, *residue, *position, *X, *Y, *Z;
                index = strtok(NULL, " ");
                name = strtok(NULL, " ");
                residue = strtok(NULL, " ");
                position = strtok(NULL, " ");
                X = strtok(NULL, " ");
                Y = strtok(NULL, " ");
                Z = strtok(NULL, " ");

                if(first_itr) RTN.add_to_atom_info(name, atoi(index), *residue, atoi(position));

                gsl_matrix_set(data_mats[iterator], row, 0, atof(X));
                gsl_matrix_set(data_mats[iterator], row, 1, atof(Y));
                gsl_matrix_set(data_mats[iterator], row, 2, atof(Z));
                row++;
            }
            else if(!strcmp(header, "ENDMDL\n") )
            {
                if(first_itr) first_itr = false;
                iterator++;
                row = 0;
            }
        }
        fclose(LibFile);
        tmp[0] = LibNames[i][0];
        tmp[1] = LibNames[i][1];
        tmp[2] = '\0';
        RTN.add_lib(data_mats, energies, tmp);
        free(data_mats);
    }
    return RTN;
}

inline void free_libs(char **A, int s)
{
    for(int i = 0; i < s; i++)
        free(A[i]);
    free(A);
}

double cross_prod(double u1, double u2, double u3, double v1, double v2, double v3, double *uvi, double *uvj, double *uvk)
{
    *uvi = u2*v3 - v2*u3;    
    *uvj = v1*u3 - u1*v3;
    *uvk = u1*v2 - v1*u2;
}

double det_get(gsl_matrix *A, int inPlace) {

    /*
    inPlace = 1 => A is replaced with the LU decomposed copy.
    inPlace = 0 => A is retained, and a copy is used for LU.
    */

    double det;
    int signum;
    gsl_permutation *p = gsl_permutation_alloc(A->size1);
    gsl_matrix *tmpA;

    if (inPlace)
       tmpA = A;
    else {
      gsl_matrix *tmpA = gsl_matrix_alloc(A->size1, A->size2);
      gsl_matrix_memcpy(tmpA , A);
    }


    gsl_linalg_LU_decomp(tmpA , p , &signum);
    det = gsl_linalg_LU_det(tmpA , signum);
    gsl_permutation_free(p);
    if (! inPlace)
       gsl_matrix_free(tmpA);


    return det;
    
}

double* center_matrix(gsl_matrix *M)
{
    double* COM = (double *)calloc(3, sizeof(double));
    for(int i = 0; i < M->size1; i++)
    {
        COM[0] += gsl_matrix_get(M, i, 0);
        COM[1] += gsl_matrix_get(M, i, 1);
        COM[2] += gsl_matrix_get(M, i, 2);
    }

    COM[0] /= M->size1;
    COM[1] /= M->size1;
    COM[2] /= M->size1;
    
    for(int i = 0; i < M->size1; i++)
    {
        gsl_matrix_set(M, i, 0, (gsl_matrix_get(M, i, 0) - COM[0]));
        gsl_matrix_set(M, i, 1, (gsl_matrix_get(M, i, 1) - COM[1]));
        gsl_matrix_set(M, i, 2, (gsl_matrix_get(M, i, 2) - COM[2]));
    }
    return COM;
}

double rotate (RNA_data *pdb_1, RNA_data *pdb_2)
{
    gsl_matrix *MODEL = pdb_2->data_matrix;
    gsl_matrix *MODEL_TEMP = gsl_matrix_alloc(MODEL->size2, MODEL->size1);

    const int dimension2 = MODEL->size2;

    gsl_matrix *P = pdb_2->get_target_matrix_copy(0);
    gsl_matrix *Q = pdb_1->get_target_matrix_copy(1);

    gsl_matrix *H = gsl_matrix_alloc(dimension2, dimension2);
    gsl_matrix *V = gsl_matrix_alloc(dimension2, dimension2);

    gsl_vector *S = gsl_vector_alloc(dimension2);     
    gsl_vector *work = gsl_vector_alloc(dimension2);     

    gsl_matrix *VUt = gsl_matrix_alloc(dimension2, dimension2);
    gsl_matrix *DIA = gsl_matrix_alloc(dimension2, dimension2);

    gsl_matrix *TEMP = gsl_matrix_alloc(dimension2, dimension2);
    gsl_matrix *R = gsl_matrix_alloc(dimension2, dimension2);

    //gsl_matrix_print(P);
//
    double *COMP, *COMQ;
    COMP = center_matrix(P);
    
    for(int i = 0; i < MODEL->size1; i++)
    {
        gsl_matrix_set(MODEL, i, 0, (gsl_matrix_get(MODEL, i, 0) - COMP[0]));
        gsl_matrix_set(MODEL, i, 1, (gsl_matrix_get(MODEL, i, 1) - COMP[1]));
        gsl_matrix_set(MODEL, i, 2, (gsl_matrix_get(MODEL, i, 2) - COMP[2]));
    }    

    COMQ = center_matrix(Q);

    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, P, Q, 0.0, H);    

    gsl_linalg_SV_decomp (H, V, S, work);

    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, V, H, 0.0, VUt);

    double det = det_get(VUt, 1);
    
    gsl_matrix_set_identity( DIA );
    gsl_matrix_set (DIA, dimension2 - 1, dimension2 - 1, det);

    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, DIA, H, 0.0, TEMP);

    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, TEMP, 0.0, R);    
    //gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, TEMP, H, 0.0, R);    
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, R, MODEL, 0.0, MODEL_TEMP);            
    
    gsl_matrix_transpose_memcpy(MODEL, MODEL_TEMP);

    for(int i = 0; i < MODEL->size1; i++)
    {
        gsl_matrix_set(MODEL, i, 0, (gsl_matrix_get(MODEL, i, 0) + COMQ[0]));
        gsl_matrix_set(MODEL, i, 1, (gsl_matrix_get(MODEL, i, 1) + COMQ[1]));
        gsl_matrix_set(MODEL, i, 2, (gsl_matrix_get(MODEL, i, 2) + COMQ[2]));
    }

    gsl_matrix_free(P);
    gsl_matrix_free(Q);

    gsl_matrix_free(H);
    gsl_matrix_free(V);

    gsl_vector_free(S);     
    gsl_vector_free(work);     

    gsl_matrix_free(VUt);
    gsl_matrix_free(DIA);

    gsl_matrix_free(TEMP);
    gsl_matrix_free(R);
    gsl_matrix_free(MODEL_TEMP);

    free(COMP);
    free(COMQ);
}

double rmsd (RNA_data *pdb_1, RNA_data *pdb_2)
{   
    double tempx, tempy, tempz;
    double temprmsd = 0.0;
    int atom_count = 0;
    
    gsl_matrix *A, *B;

    A = pdb_1->get_target_matrix(1);
    B = pdb_2->get_target_matrix(0);

    for(int i = 0; i < B->size1; i++)
    {
        tempx = pow((gsl_matrix_get(A, i, 0) - gsl_matrix_get(B, i, 0)), 2);
        tempy = pow((gsl_matrix_get(A, i, 1) - gsl_matrix_get(B, i, 1)), 2);
        tempz = pow((gsl_matrix_get(A, i, 2) - gsl_matrix_get(B, i, 2)), 2);
        temprmsd += tempx + tempy + tempz;
        //printf("sum = %f\n", temprmsd)       ;
        atom_count++;    
    }
    return sqrt( temprmsd/(double)atom_count ); 
}

double distance(gsl_vector *A, gsl_vector *B)
{
    double x,y,z;
    x = gsl_vector_get(B, 0) - gsl_vector_get(A, 0);
    y = gsl_vector_get(B, 1) - gsl_vector_get(A, 1);
    z = gsl_vector_get(B, 2) - gsl_vector_get(A, 2);
    return sqrt(square(x) + square(y) + square(z));
}

bool overlap_check(RNA_data_array& sequence, RNA_data *attach)
{
    RNA_data *base = sequence.current();
    gsl_vector_view A, B;
    double radius_1, radius_2;
    double dist;
    for(int i = 0; i < sequence.count; i++)
    {
        for(int j = 0; j < sequence[i]->count; j++)
        {
            if(i == sequence.count - 1 && sequence[i]->atom_data->dnt_pos[j] == 2)
            {
                //printf("residue: %d, pos: %s\n", sequence[i]->atom_data->dnt_pos[j], sequence[i]->atom_data->name[j]);
                continue;
            }
            for(int k = 0; k < attach->count; k++)
            {
                if(attach->atom_data->dnt_pos[k] == 1)
                {
                    continue;
                }
                switch(sequence[i]->atom_data->name[j][0]) //Need to change target to hash table/linked list to quickly check if current atoms is in target without massive performance loss.
                {
                    case 'C': radius_1 = RADIUS_C; break; 
                    case 'N': radius_1 = RADIUS_N; break; 
                    case 'O': radius_1 = RADIUS_O; break; 
                    case 'P': radius_1 = RADIUS_P; break; 
                }
                switch(attach->atom_data->name[k][0])
                {
                    case 'C': radius_2 = RADIUS_C; break; 
                    case 'N': radius_2 = RADIUS_N; break; 
                    case 'O': radius_2 = RADIUS_O; break; 
                    case 'P': radius_2 = RADIUS_P; break; 
                }
                A = gsl_matrix_row(attach->data_matrix, k);
                B = gsl_matrix_row(sequence[i]->data_matrix, j);
                dist = distance(&A.vector, &B.vector);
                //printf("s: %s a: %s\tdist = %f\tvdw = %f\n", sequence[i]->atom_data->name[j], attach->atom_data->name[k], dist, (radius_1 + radius_2));
                if(dist < (radius_1 + radius_2))
                {
                    return false;
                }
            }
        }
    }
    return true;
}

enum attach_status{FAILED, ATTACHED, NOT_CHECKED};

attach_status check_attachment(RNA_data_array& sequence, RNA_data *base, RNA_data *attach)
{
    double rmsd_v;
    rmsd_v = rmsd(base, attach);
    //printf("RMSD = %f\n", rmsd_v);
    if(rmsd_v > RMSD_LIMIT)
    {
        *attach->_flag = NOT_USABLE;
        return FAILED;
    }
    if(!overlap_check(sequence, attach))
    {
        //printf("Overlap detected!\n");
        *attach->_flag = NOT_USABLE;
        return FAILED;
    }
    return ATTACHED;
}

bool combinatorial_addition(DimerLibArray& Lib, RNA_data_array& assembled, CMB_Manager& manager, output_string& o_string)
{
    int working_position = assembled.iterator + 1; //Position in sequence where new DNT will be attached.
    //printf("Working postion: %d, iterator:%d, iterator_max:%d\n", working_position, assembled.iterator, assembled.iterator_max);
    DimerLib *Library = Lib[working_position];
    RNA_data *base;         //Already attached base which will have new DNT attached to
    RNA_data *attach;       //DNT which will be attached
    attach_status status = FAILED;   //Output from checking functions

    DEBUG(printf("Early: checking assembled 0: %ld, working position: %d\n", assembled[0]->id, working_position));

    attach = new RNA_data(Lib, working_position, 0); // For initialization only
    int *idxs = assembled.get_index();
    
    printf("New Iteration\n");
    for(int i = 0; i < Library->count; i++)
    {   
        if(Library->flags[i] != NO_FLAG)
        {
            continue;
        }
        for(int j = 0; j < assembled.count; j++)
        {
            printf("%d ",idxs[j]);
        }
        printf("%d\n", i);
        //printf("count for %s library: %d\n", Library->name, i);
        if(assembled.is_empty())
        {
            printf("sequence is empty\n");
            attach->overwrite(Lib, working_position, i);
            assembled.add_move(attach);
            manager.attach_attempt(working_position, i);
            DEBUG(printf("attach: %ld moved @ empty\n", attach->id));
            status = NOT_CHECKED;
            break;
        }
        base = assembled.current();
        attach->overwrite(Lib, working_position, i);
        rotate(base, attach);
        attach->update_submatrices();
        manager.attach_attempt(working_position, i);
        if((status = check_attachment(assembled, base, attach)) == ATTACHED)
        {
            break;
        }
    }
    free(idxs);
    if(status == FAILED)
    {
        printf("ALL FAILED\n");
        DEBUG(printf("attach: %ld deleted @ All failed\n", attach->id));
        delete attach;
        //printf("flagging index %d,%d\n", base->position_in_lib[0], base->position_in_lib[1]);
        DEBUG(printf("All failed: checking assembled 0: %ld\n", assembled[0]->id));
        DEBUG(printf("assembled count = %d\n", assembled.count));
        DEBUG(printf("last added = %d, max in lib: %d\n", Library->flags[Library->count - 1], Library->count));
        if(manager.is_at_end())
        {
            manager.check_lib_completion();
            if((manager.get_reset_count()) == manager.last_attempted[0] + 1)
            {
                return true;
            }
            printf("lib completed = %d; n_reset = %d\n", manager.last_attempted[0] + 1, manager.get_reset_count());
            
            assembled.rollback_by(manager.get_reset_count() - 1);
            Lib.reset_flags(manager.libs_completed);
            manager.clear_attempts();
        }
        else
        {
            assembled.rollback();
        }
        /*printf("Rolling back b/c: FAILED\n");
        //printf("last attempt: %d %d\n", manager.last_attempted[0], manager.last_attempted[1]);
        assembled.rollback();
        manager.check_lib_completion();
        manager.clear_attempts();*/
        return false;
        //printf("assembled iterator: %d\n", assembled.iterator);
    }
    else if(status == NOT_CHECKED)
    {
        //printf("Replaced position 0\n");
        return false;
    }
    else if(status == ATTACHED)
    {
        //printf("addition successful\n");
        assembled.add_move(attach);
        DEBUG(printf("attach: %ld moved @ Attached @ %d\n", attach->id, working_position));
    }
    if(assembled.is_complete())
    {
        o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
        manager.strs_built++;
        if(manager.is_at_end())
        {
            manager.check_lib_completion();
            if((manager.get_reset_count()) == Lib.count)
            {
                return true;
            }
            //printf("lib completed = %d\n", manager.get_reset_count());
            
            assembled.rollback_by(manager.get_reset_count());
            Lib.reset_flags(manager.libs_completed);
            manager.clear_attempts();
        }
        else
        {
            DEBUG(printf("Rolling back b/c: COMPLETED NOT AT LIB END\n"));
            assembled.rollback();
        }
    }
    /* Use CMB_manager to decide next addition by moving iterator back by one if not at end of lib, or back by
    however much once end of each respective library is reached. On each successful addition, save structure to
    string. On each 1000 ( or TBD) write out data to file and overwrite previous strings. Flags could be complex*/
    return false;
}

int main()
{
    int N_diNts = 0;
    char sequence[] = "AUGCG";
    char **Libs2Load = get_diNt_names(sequence, &N_diNts);
    DimerLibArray Library = load_libs(Libs2Load, N_diNts);
    RNA_data_array RNA(N_diNts);
    RNA.add_move(new RNA_data(Library, 0, 0));
    
    CMB_Manager manager(Library);

    output_string output_s("testing_IVANCMB.txt", MAX_STRINGS);

    while(!combinatorial_addition(Library, RNA, manager, output_s));

    /* *** MAKE A HASH FOR DINUCLEOTIDE LIBRARY FOR SIMPLICITY & EFFICIENCY *** */
    /* ALSO MAKE HASH FOR GET_TARGET FOR EFFICIENCY */
    /* CHECK FOR COLLISIONS WITH NON NEIGHBORING NTs TO FLAG AS UNUSABLE TO IMPROVE EFFICIENCY */

    /* PDB FORMAT: printf("%-6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f\n", "ATOM", index , name, residue, residue ID, X, Y, Z)*/


    //RNA.to_string();

    free_libs(Libs2Load, N_diNts);

    printf("# of Structures Built: %d\n", manager.strs_built);
}
