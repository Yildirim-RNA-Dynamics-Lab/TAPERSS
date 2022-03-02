#include "RNACMB.hpp"

/* Set Global Defaults */
long int rna_dat_tracker = 0;
long unsigned int steric_clash_checks_attempted = 0;
long unsigned int steric_clash_checks_skipped = 0;

char GLOBAL_INPUT_SEQUENCE[100];
char GLOBAL_INPUT_INDICES[100];
char GLOBAL_INPUT_FILE[100];
char GLOBAL_OUTPUT_FILE[100];
char GLOBAL_IO_ACTION[5] = "w";

bool GLOBAL_RUN_COMBINATORIAL = false;
bool GLOBAL_RUN_BUILD_STRUCTURE = false;
bool GLOBAL_RUN_BUILD_STRUCTURE_LIST = false;
bool GLOBAL_RUN_BUILD_STRUCTURE_LIST_TESTING = false;
bool GLOBAL_WRITE_COORDINATES = false;
bool GLOBAL_PERFORM_HAIRPIN_CHECK = false;

double GLOBAL_RMSD_LIMIT = DEFAULT_RMSD_LIMIT;
double GLOBAL_WC_RMSD_LIMIT = DEFAULT_WC_RMSD_LIMIT;
double GLOBAL_SCC_LIMIT = DEFAULT_STERIC_CLASH_COM_DISTANCE_LIMIT;
int **GLOBAL_INPUT_INDICES_LIST;
char LIBRARY_FILENAME_PROTOTYPE[100];     // = "../AGAAAU_test/XX_library_combined.txt";
char WATSON_CRICK_LIBRARY_PROTOTYPE[100]; // = "../AGAAAU_test/WC_XX_library.txt";
