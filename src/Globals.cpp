#include "RNACMB.hpp"

/* Set Global Defaults */
uint64_t rna_dat_tracker = 0;
uint64_t steric_clash_checks_attempted = 0;
uint64_t steric_clash_checks_skipped = 0;

char GLOBAL_INPUT_SEQUENCE[GLOBAL_STANDARD_STRING_LENGTH];
char GLOBAL_INPUT_INDICES[GLOBAL_STANDARD_STRING_LENGTH];
char GLOBAL_INPUT_FILE[GLOBAL_STANDARD_STRING_LENGTH];
char GLOBAL_OUTPUT_FILE[GLOBAL_STANDARD_STRING_LENGTH];
char GLOBAL_IO_ACTION[5] = "w";

bool GLOBAL_RUN_COMBINATORIAL = false;
bool GLOBAL_RUN_BUILD_STRUCTURE = false;
bool GLOBAL_RUN_BUILD_STRUCTURE_LIST = false;
bool GLOBAL_RUN_BUILD_STRUCTURE_LIST_TESTING = false;
bool GLOBAL_WRITE_COORDINATES = false;
uint32_t GLOBAL_STRUCTURE_LIMIT_COUNT = 1000;
//bool GLOBAL_PERFORM_HAIRPIN_CHECK = false;
//bool GLOBAL_PERFORM_HBOND_CHECK = false;
STRUCTFILTER_TYPE GLOBAL_PERFORM_STRUCTCHECK = NONE;

double GLOBAL_RMSD_LIMIT = DEFAULT_RMSD_LIMIT;
double GLOBAL_WC_RMSD_LIMIT = DEFAULT_WC_RMSD_LIMIT;
double GLOBAL_SCC_LIMIT = DEFAULT_STERIC_CLASH_COM_DISTANCE_LIMIT;
uint32_t **GLOBAL_INPUT_INDICES_LIST;
char LIBRARY_FILENAME_PROTOTYPE[GLOBAL_STANDARD_STRING_LENGTH];     // = "../XX_library_combined.txt";
char WATSON_CRICK_LIBRARY_PROTOTYPE[GLOBAL_STANDARD_STRING_LENGTH]; // = "../AGAAAU_test/WC_XX_library.txt";
