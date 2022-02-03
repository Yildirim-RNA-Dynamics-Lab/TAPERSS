// First attempt for Combinatorial building RNA  **03/5/21**
#include "RNACMB.hpp"
#include "Atom_info.hpp"
#include "DimerLib.hpp"
#include "RNA_data.hpp"
#include "Combinatorial_Addition.hpp"
#include "output_string.hpp"
#include "Hairpin.hpp"
#include "InputHandler.hpp"
#include "StructureBuilders.hpp"
using namespace std;

/* Hello */

/* Set Global Defaults */
long int rna_dat_tracker = 0;

char GLOBAL_INPUT_SEQUENCE[100];
char GLOBAL_INPUT_INDICES[100];
char GLOBAL_INPUT_FILE[100];
char GLOBAL_OUTPUT_FILE[100];
char GLOBAL_IO_ACTION[5] = "w";

bool GLOBAL_RUN_COMBINATORIAL = false;
bool GLOBAL_RUN_BUILD_STRUCTURE = false;
bool GLOBAL_RUN_BUILD_STRUCTURE_LIST = false;
bool GLOBAL_WRITE_COORDINATES = false;
bool GLOBAL_PERFORM_HAIRPIN_CHECK = false;

double GLOBAL_RMSD_LIMIT = DEFAULT_RMSD_LIMIT;
double GLOBAL_WC_RMSD_LIMIT = DEFAULT_WC_RMSD_LIMIT;
int **GLOBAL_INPUT_INDICES_LIST;
char LIBRARY_FILENAME_PROTOTYPE[100];     // = "../AGAAAU_test/XX_library_combined.txt";
char WATSON_CRICK_LIBRARY_PROTOTYPE[100]; // = "../AGAAAU_test/WC_XX_library.txt";

int main(int argc, char *ARGV[])
{
    int N_diNts = 0;
    int N_WC = 0;

    char **Libs2Load = nullptr;
    char **WCLibs2Load = nullptr;

    DimerLibArray Library;
    DimerLibArray WC_Library;

    input_handler(argc, ARGV);
    Libs2Load = get_diNt_names(GLOBAL_INPUT_SEQUENCE, &N_diNts);
    load_libs(Libs2Load, N_diNts, Library);

    if (GLOBAL_PERFORM_HAIRPIN_CHECK)
    {
        WCLibs2Load = get_WC_partner(GLOBAL_INPUT_SEQUENCE, &N_WC);
        load_libs(WCLibs2Load, N_WC, WC_Library, true);
    }

    RNA_data_array RNA(N_diNts);
    RNA.add_move(new RNA_data(Library, 0, 0));

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
        while (!combinatorial_addition(Library, RNA, manager, output_s, WC_Library))
            ;
        printf("# of Structures Built: %d\n", manager.strs_built);
        if (GLOBAL_PERFORM_HAIRPIN_CHECK)
        {
            printf("# of Hairpins Built: %d\n", manager.hairpins_built);
        }
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

    /* PDB FORMAT: printf("%-6s%5d %4s %3s  %4d    %8.3f%8.3f%8.3f\n", "ATOM", index , name, residue, residue ID, X, Y, Z)*/

    free_libs(Libs2Load, N_diNts);
    if (GLOBAL_PERFORM_HAIRPIN_CHECK)
    {
        free_libs(WCLibs2Load, N_WC);
    }
}
