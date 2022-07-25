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

using namespace std;

auto n = [](char* sequence) -> unsigned int
{
    char *ptr = sequence;
    unsigned int n = 0;
    while( *(ptr++) != 'x')
    {
       n++;
    }
    return n;
};

char **get_diNt_wrapper(char *sequence, int *N, int *leftStrand, int *rightStrand)
{
    char** rtn = nullptr;
    
    unsigned int x_loc = n(sequence);
    int left = 0, right = 0;
    switch(GLOBAL_PERFORM_STRUCTCHECK)
    {
        case HAIRPIN:
            *N = strlen(sequence) - 1;
            rtn = (char **)malloc(sizeof(char *) * *N);
            get_diNt_names(sequence, rtn, *N);
            break;
        case INTERNAL_LOOP:
            *N = strlen(sequence) - 3; //e.g. A G A X U G U, will have 7 - 3 = 4 diN 
            rtn = (char **)malloc(sizeof(char *) * *N);
            for(unsigned int i = 0; i < strlen(sequence); i++) 
            {
                if(i < x_loc)
                {
                    left++;
                }
                else if(i == x_loc)
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
            get_diNt_names(sequence, rtn, left);
            get_diNt_names(&sequence[x_loc + 1], &rtn[left], right);
            break;
        case NONE:
            //get_diNt_names;
            *N = strlen(sequence) - 1;
            rtn = (char **)malloc(sizeof(char *) * *N);
            get_diNt_names(sequence, rtn, *N);
            break;
    }
    //*N = strlen(sequence) - 1;
    return rtn;
}

char** get_WC_wrapper(char *sequence, int *N)
{
    char** rtn = nullptr;
    unsigned int x_loc = n(sequence);
    switch(GLOBAL_PERFORM_STRUCTCHECK)
    {
        case HAIRPIN:
            *N = strlen(sequence) - 1;
            rtn = (char **)malloc(sizeof(char *) * *N);
            get_WC_partner(sequence, rtn, *N);
            *N = 1;
            break;
        case INTERNAL_LOOP:
            *N = 2;
            rtn = (char **)malloc(sizeof(char *) * *N);
            get_WC_partner(sequence, rtn, strlen(sequence));
            get_WC_partner(&sequence[x_loc - 1], &rtn[1], 3);
            break;
        case NONE:
            break;
    }
    return rtn;
}

template <typename T> void Run()
{
    int N_diNts = 0;
    int N_WC = 0;

    int leftStrand = 0;
    int rightStrand = 0;
    
    char **Libs2Load = nullptr;
    char **WCLibs2Load = nullptr;
    
    T RNA;

    DimerLibArray Library;
    DimerLibArray WC_Library;

    Libs2Load = get_diNt_wrapper(GLOBAL_INPUT_SEQUENCE, &N_diNts, &leftStrand, &rightStrand);
    load_libs(Libs2Load, N_diNts, Library);
    if constexpr (is_same<T,RNADataArray>::value)
    {
        WCLibs2Load = get_WC_wrapper(GLOBAL_INPUT_SEQUENCE, &N_WC);
        load_libs(WCLibs2Load, N_WC, WC_Library, true);
        RNA.initialize(N_diNts);
    }
    if constexpr (is_same<T,RNADataArrayInternalLoop>::value)
    {
        WCLibs2Load = get_WC_wrapper(GLOBAL_INPUT_SEQUENCE, &N_WC);
        load_libs(WCLibs2Load, N_WC, WC_Library, true);
        RNA.initialize(leftStrand, rightStrand);
    }
    
    RNA.add_move(new RNAData(Library, 0, 0));


    CMB_Manager manager(Library);

    if (GLOBAL_OUTPUT_FILE == NULL)
    {
        printf("No output file specified...\n");
        exit(1);
    }
    output_string output_s(GLOBAL_OUTPUT_FILE, MAX_STRINGS, GLOBAL_IO_ACTION);

    if (GLOBAL_RUN_COMBINATORIAL)
    {
        printf("Performing combinatorial run on sequence: %s\n", GLOBAL_INPUT_SEQUENCE);

        if constexpr (is_same<T,RNADataArray>::value)
        {
           while (!combinatorial_addition(Library, RNA, manager, output_s, WC_Library))
            ;
        }
        if constexpr (is_same<T,RNADataArrayInternalLoop>::value)
        {
            while (!combinatorial_addition_IL(Library, RNA, manager, output_s, WC_Library))
            ;
        }
      
        printf("# of Structures Built: %ld\n", manager.strs_built);
        if (GLOBAL_PERFORM_STRUCTCHECK == HAIRPIN)
        {
            printf("# of Hairpins Built: %ld\n", manager.hairpins_built);
        }
        printf("# of Steric Clash Checks Attempted: %ld\n", steric_clash_checks_attempted);
        printf("# of Steric Clash Checks Skipped: %ld\n", steric_clash_checks_skipped);
        printf("%% of Steric Clash Checks Skipped: %f\n", (float)steric_clash_checks_skipped/(float)steric_clash_checks_attempted * 100);
    }
    if (GLOBAL_RUN_BUILD_STRUCTURE)
    {
        printf("Creating Single Structure of Sequence: %s, with Indices: %s\n", GLOBAL_INPUT_SEQUENCE, GLOBAL_INPUT_INDICES);
        int *indices = (int *)alloca(N_diNts * sizeof(int));
        get_index_int(GLOBAL_INPUT_INDICES, indices);
        create_custom_structure(Library, WC_Library, RNA, output_s, indices);
    }
    if (GLOBAL_RUN_BUILD_STRUCTURE_LIST)
    {
        printf("Creating structures from index list...\n");
        int num_strs = read_input_index_file(GLOBAL_INPUT_FILE, N_diNts);
        printf("Done reading input...\n");
        create_custom_structure_list(Library, WC_Library, RNA, output_s, num_strs);
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
        create_custom_structure_list_testing(Library, RNA, output_s, num_strs);
        for (int i = 0; i < num_strs; i++)
        {
            free(GLOBAL_INPUT_INDICES_LIST[i]);
        }
        free(GLOBAL_INPUT_INDICES_LIST);
    }

    free_libs(Libs2Load, N_diNts);
    if (GLOBAL_PERFORM_STRUCTCHECK == HAIRPIN)
    {
        free_libs(WCLibs2Load, N_WC);
    }

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
