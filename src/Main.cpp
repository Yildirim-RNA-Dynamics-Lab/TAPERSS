// First attempt for Combinatorial building RNA  **03/5/21**
#include "RNACMB.hpp"
#include "Atom_info.hpp"
#include "DimerLib.hpp"
#include "RNAData.hpp"
#include "Combinatorial_Addition.hpp"
#include "output_string.hpp"
#include "Hairpin.hpp"
#include "HBondDetector.hpp"
#include "InputHandler.hpp"
#include "StructureBuilders.hpp"
#include "RNADataArrayInternalLoop.hpp"
#include <unistd.h>

using namespace std;

unsigned int get_x_location(char *sequence)
{
    char *ptr = sequence;
    unsigned int size = strlen(sequence);
    unsigned int n = 0;
    for (unsigned int i = 0; i < size; i++)
    {
        //putchar(ptr[i]);
        if (ptr[i] != 'x')
        {
            n++;
        }
        else
        {
            break;
        }
    }
    //putchar('\n');
    //printf("n = %d\n", n);
    return n;
}

char **get_diNt_wrapper(char *sequence, int *N, int *leftStrand, int *rightStrand, int** duplicate_record)
{
    char **rtn = nullptr;
    int* pointer = duplicate_record[0];
    unsigned int x_loc = get_x_location(sequence);
    int left = 0, right = 0;
    //printf("%s\n", sequence);
    switch (GLOBAL_PERFORM_STRUCTCHECK)
    {
    case HAIRPIN:
        *N = strlen(sequence) - 1;
        pointer = (int *)malloc(sizeof(int) * *N);
        memset(pointer, -1, *N * sizeof(int));
        rtn = (char **)malloc(sizeof(char *) * *N);
        get_diNt_names(sequence, rtn, pointer, *N);
        break;
    case INTERNAL_LOOP:
        *N = strlen(sequence) - 3; // e.g. A G A X U G U, will have 7 - 3 = 4 diN
        pointer = (int *)malloc(sizeof(int) * *N);
        memset(pointer, -1, *N * sizeof(int));
        rtn = (char **)malloc(sizeof(char *) * *N);
        for (unsigned int i = 0; i < strlen(sequence); i++)
        {
            if (i < x_loc)
            {
                left++;
            }
            else if (i == x_loc)
            {
                continue;
            }
            else
            {
                right++;
            }
        }
        left--, right--;
        *leftStrand = left;
        *rightStrand = right;
        get_diNt_names(sequence, rtn, pointer, left);
        get_diNt_names(&sequence[x_loc + 1], &rtn[left], &pointer[left], right);
        break;
    case NONE:
        // get_diNt_names;
        *N = strlen(sequence) - 1;
        pointer = (int *)malloc(sizeof(int) * *N);
        memset(pointer, -1, *N * sizeof(int));
        rtn = (char **)malloc(sizeof(char *) * *N);
        get_diNt_names(sequence, rtn, pointer, *N);
        break;
    }

    duplicate_record[0] = pointer;
    return rtn;
}


/* I should add a duplicate checker for this function as well */
char **get_WC_wrapper(char *sequence, int *N)
{
    char **rtn = nullptr;
    unsigned int x_loc = get_x_location(sequence);
    int length;
    switch (GLOBAL_PERFORM_STRUCTCHECK)
    {
    case HAIRPIN:
        length = strlen(sequence) - 1;
        rtn = (char **)malloc(sizeof(char *) * length);
        get_WC_partner(sequence, rtn, 0, length);
        *N = 1;
        break;
    case INTERNAL_LOOP:
        *N = 2;
        rtn = (char **)malloc(sizeof(char *) * *N);
        get_WC_partner(sequence, rtn, 0, strlen(sequence) - 1);
        get_WC_partner(sequence, &rtn[1], x_loc - 1, x_loc + 1);
        break;
    case NONE:
        break;
    }
    return rtn;
}

template <typename T>
void Run()
{
    int N_diNts = 0;
    int N_WC = 0;
    int *DuplicateRecord = nullptr;

    int leftStrand = 0;
    int rightStrand = 0;

    char **Libs2Load = nullptr;
    char **WCLibs2Load = nullptr;
    T RNA;

    DimerLibArray Library;
    DimerLibArray WC_Library;

    clock_t Start, End;
    double TimeUsed;

    Start = clock();

    Libs2Load = get_diNt_wrapper(GLOBAL_INPUT_SEQUENCE, &N_diNts, &leftStrand, &rightStrand, &DuplicateRecord);
    load_libs(Libs2Load, N_diNts, Library, DuplicateRecord);
    Library.get_charged_atom_map();
    if constexpr (is_same<T, RNADataArray>::value)
    {
        WCLibs2Load = get_WC_wrapper(GLOBAL_INPUT_SEQUENCE, &N_WC);
        load_libs(WCLibs2Load, N_WC, WC_Library, nullptr, true);
        RNA.initialize(N_diNts, Library);
    }
    if constexpr (is_same<T, RNADataArrayInternalLoop>::value)
    {
        WCLibs2Load = get_WC_wrapper(GLOBAL_INPUT_SEQUENCE, &N_WC);
        load_libs(WCLibs2Load, N_WC, WC_Library, nullptr, true);
        RNA.initialize(leftStrand, rightStrand, Library, WC_Library);
    }

    //RNA.add_move(new RNAData(Library, 0, 0));

    CMB_Manager manager(Library);

    if (GLOBAL_OUTPUT_FILE[0] == '\0')
    {
        printf("No output file specified...\n");
        exit(1);
    }
    output_string output_s(GLOBAL_OUTPUT_FILE, GLOBAL_MAX_STRINGS, GLOBAL_IO_ACTION);
    kabsch_create(Library.LargestAtomCount, MATRIX_DIMENSION2);
    WC_create(WC_Library);
    HK_create(Library.PositiveAtomCount, Library.NegativeAtomCount);

    End = clock();

    TimeUsed = ((double)(End - Start)) / CLOCKS_PER_SEC;

    printf("Time to load: %fms\n", TimeUsed * 1000);

    if (GLOBAL_RUN_COMBINATORIAL)
    {
        printf("Performing combinatorial run on sequence: %s\n", GLOBAL_INPUT_SEQUENCE);
        Start = clock();
        if constexpr (is_same<T, RNADataArray>::value)
        {
            while (!combinatorial_addition(Library, RNA, manager, output_s))
            {
                ;
            }  
        }
        if constexpr (is_same<T, RNADataArrayInternalLoop>::value)
        {
            while (!combinatorial_addition_IL(Library, RNA, manager, output_s, WC_Library))
                ;
        }
        End = clock();
        TimeUsed = ((double)(End - Start)) / CLOCKS_PER_SEC;
        printf("# of Structures Built: %ld\n", manager.strs_built);
        if (GLOBAL_PERFORM_STRUCTCHECK == HAIRPIN)
        {
            printf("# of Hairpins Built: %ld\n", manager.hairpins_built);
        }
        if (GLOBAL_PERFORM_STRUCTCHECK == INTERNAL_LOOP)
        {
            printf("# of Internal loops Built: %ld\n", manager.internal_loops_built);
        }
        //printf("# of Steric Clash Checks Attempted: %ld\n", steric_clash_checks_attempted);
        //printf("# of Steric Clash Checks Skipped: %ld\n", steric_clash_checks_skipped);
        //printf("%% of Steric Clash Checks Skipped: %f\n", (float)steric_clash_checks_skipped / (float)steric_clash_checks_attempted * 100);
        printf("Calculation Time: ");
        printf("%dh:%dm:%ds\n", (int)(TimeUsed / (60 * 60)), ((int)TimeUsed % (60 * 60))/60, ((int)TimeUsed % (60 * 60)) % 60);
    }

    if (GLOBAL_RUN_BUILD_STRUCTURE)
    {
        printf("Creating Single Structure of Sequence: %s, with Indices: %s\n", GLOBAL_INPUT_SEQUENCE, GLOBAL_INPUT_INDICES);
        int *indices = (int *)alloca(N_diNts * sizeof(int));
        get_index_int(GLOBAL_INPUT_INDICES, indices);
        if constexpr (is_same<T, RNADataArray>::value)
        {
            create_custom_structure(Library, RNA, output_s, indices);
        }
        if constexpr (is_same<T, RNADataArrayInternalLoop>::value)
        {
            create_custom_structure_IL(Library, WC_Library, RNA, output_s, indices);
        }
    }
    if (GLOBAL_RUN_BUILD_STRUCTURE_LIST)
    {
        printf("Creating structures from index list...\n");
        int num_strs = read_input_index_file(GLOBAL_INPUT_FILE, N_diNts);
        printf("Done reading input...\n");
        if (GLOBAL_PERFORM_STRUCTCHECK == HAIRPIN)
        {
            create_custom_structure_list<PERFORM_CHECKS_ON_CUSTOM_BUILD, HAIRPIN>(Library, RNA, output_s, num_strs); // Need to define for IL at somepoint
        }
        for (int i = 0; i < num_strs; i++)
        {
            free(GLOBAL_INPUT_INDICES_LIST[i]);
        }
        free(GLOBAL_INPUT_INDICES_LIST);
    }
    if (GLOBAL_RUN_BUILD_STRUCTURE_LIST_TESTING)
    {
        printf("Reading input... ");
        int num_strs = read_input_index_file(GLOBAL_INPUT_FILE, N_diNts);
        printf("Done!\n");

        printf("Creating structures from index list for testing...\n");
        //create_custom_structure_list_testing(Library, RNA, output_s, num_strs);
        for (int i = 0; i < num_strs; i++)
        {
            free(GLOBAL_INPUT_INDICES_LIST[i]);
        }
        free(GLOBAL_INPUT_INDICES_LIST);
    }

    free_libs(Libs2Load, N_diNts);
    if (GLOBAL_PERFORM_STRUCTCHECK == HAIRPIN || GLOBAL_PERFORM_STRUCTCHECK == INTERNAL_LOOP)
    {
        free_libs(WCLibs2Load, N_WC);
    }
    kabsch_destroy();
    WC_destroy();
    HK_destroy();
    free(DuplicateRecord);
}

int main(int argc, char *ARGV[])
{
    input_handler(argc, ARGV);
    switch (GLOBAL_PERFORM_STRUCTCHECK)
    {
    case INTERNAL_LOOP:
        Run<RNADataArrayInternalLoop>();
        break;
    default:
        Run<RNADataArray>();
        break;
    }
    /* PDB FORMAT: printf("%-6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f\n", "ATOM", index , name, residue, residue ID, X, Y, Z)*/
}
