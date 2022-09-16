#ifndef RNACMB_H
#define RNACMB_H
//Includes
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <type_traits>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <cstdint>

/* Semi-useful preprocessor functions */
#define square(a) (a * a)
#define IDX_FLAT2D(row, col, n) (row * n + col)
#define DEBUG_SWITCH false
#define DEBUG(a) if(DEBUG_SWITCH == true) {a;}
#define DIE exit(2)

enum STRUCTCHECK_TYPE{HAIRPIN, INTERNAL_LOOP, NONE};
enum attach_status{FAILED, ATTACHED, NOT_CHECKED, FAILED_SC};

/* Compile-time evaluated constants */
constexpr float RADIUS_C = 0.75;
constexpr float RADIUS_N = 0.75;
constexpr float RADIUS_P = 0.75;
constexpr float RADIUS_O = 0.75;
constexpr float INTERACTION_DISTANCE = 3.6;
constexpr float DEFAULT_RMSD_LIMIT = 0.5;
constexpr float DEFAULT_WC_RMSD_LIMIT = 2.5;
constexpr float DEFAULT_STERIC_CLASH_COM_DISTANCE_LIMIT = 11.2;
constexpr bool STRUCTURE_BUILD_LIMIT = true;
constexpr bool PERFORM_CHECKS_ON_CUSTOM_BUILD = false;
constexpr uint_fast8_t MATRIX_DIMENSION2 = 3; //DIM2 = 3, for X,Y,Z

constexpr int MAX_STRINGS = 1000;

/* Global Variables */
extern uint_fast64_t rna_dat_tracker;
extern uint_fast64_t steric_clash_checks_attempted;
extern uint_fast64_t steric_clash_checks_skipped;

extern char GLOBAL_INPUT_SEQUENCE[100];
extern char GLOBAL_INPUT_INDICES[100];
extern char GLOBAL_INPUT_FILE[100];
extern char GLOBAL_OUTPUT_FILE[100];
extern char GLOBAL_IO_ACTION[5];

extern bool GLOBAL_RUN_COMBINATORIAL;
extern bool GLOBAL_RUN_BUILD_STRUCTURE;
extern bool GLOBAL_RUN_BUILD_STRUCTURE_LIST;
extern bool GLOBAL_RUN_BUILD_STRUCTURE_LIST_TESTING; /* This will probably be a temporary setting */
extern bool GLOBAL_WRITE_COORDINATES;
//extern bool GLOBAL_PERFORM_HAIRPIN_CHECK;
//extern bool GLOBAL_PERFORM_HBOND_CHECK;
extern STRUCTCHECK_TYPE GLOBAL_PERFORM_STRUCTCHECK;

extern double GLOBAL_RMSD_LIMIT;
extern double GLOBAL_WC_RMSD_LIMIT;
extern double GLOBAL_SCC_LIMIT;
extern int **GLOBAL_INPUT_INDICES_LIST;
extern uint_fast64_t GLOBAL_STRUCTURE_LIMIT_COUNT;

extern char LIBRARY_FILENAME_PROTOTYPE[100];        // = "../AGAAAU_test/XX_library_combined.txt";
extern char WATSON_CRICK_LIBRARY_PROTOTYPE[100];    // = "../AGAAAU_test/WC_XX_library.txt";

#endif