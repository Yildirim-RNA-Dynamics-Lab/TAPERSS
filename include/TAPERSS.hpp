#ifndef RNACMB_H
#define RNACMB_H
//Includes
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <type_traits>
#include <cstdint>
#include <cctype>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

/* Semi-useful preprocessor functions */
#define square(a) ((a) * (a))
#define IDX_FLAT2D(row, col, n_col) ((row) * (n_col) + (col))
#define DEBUG_SWITCH false
#define DEBUG(a) if(DEBUG_SWITCH == true) {a;}
#define DEBUG_PRINT(...) do { printf("%s:%d: ", __FILE__,__LINE__); printf(__VA_ARGS__);  } while (0)
#define DIE do { printf("DIED at line %d in file %s!\n", __LINE__, __FILE__); exit(2);} while (0)

enum STRUCTFILTER_TYPE{HAIRPIN, INTERNAL_LOOP, NONE};
enum AttachStatus{FAILED, ATTACHED};

enum RunType{combinatorial, build_from_index, build_from_index_list, runtype_undef};
enum StrType{single_strand, double_strand, strtype_undef};
enum RunOpts{
	build_limit_by_energy=0x00000001, 
	blind_build_limit=0x00000002, 
	write_coordinates=0x00000004, 
	use_structure_filter=0x00000008, 
	strtype_ds=0x00000010,
	str_filter_uses_ds_closing_bp=0x00000020,
};
/* Compile-time evaluated constants */
constexpr float VDW_RADIUS = 1.5;
constexpr float RADIUS_C   = 0.75;
constexpr float RADIUS_N   = 0.75;
constexpr float RADIUS_P   = 0.75;
constexpr float RADIUS_O   = 0.75;
constexpr float INTERACTION_DISTANCE = 3.6;
constexpr float DEFAULT_RMSD_LIMIT = 0.5;
constexpr float DEFAULT_WC_RMSD_LIMIT = 0.5;
constexpr char  DEFAULT_OUTPUT_FILENAME[] = "RNACMB_Output.dat";
constexpr bool STRUCTURE_BUILD_LIMIT = false;
constexpr bool PERFORM_CHECKS_ON_CUSTOM_BUILD = true;
constexpr uint8_t MATRIX_DIMENSION2 = 3; //DIM2 = 3, for X,Y,Z
constexpr size_t DEFAULT_MEMORY_LIMIT = 5 * 1048576; //5MB

constexpr int GLOBAL_MAX_STRINGS = 1000;
constexpr int GLOBAL_STANDARD_STRING_LENGTH = 2048; //2kB

struct Timer {
	clock_t start,end;
	double time_used;
	void start_timer() {
		start = clock();
	}
	void stop_timer() {
		end = clock();
	}
	void print() {
		time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
		uint h = time_used / (60 * 60);
		uint m = (time_used - (h * 60 * 60)) / 60;
		uint s = time_used - ((h * 60 * 60) + (m * 60));
		uint ms = (time_used - ((h * 60 * 60) + (m * 60) + s)) * 1000;
		printf("%dh:%dm:%ds:%dms\n", h, m, s, ms);
	}
};

template <typename T> struct IndexPair {
	T idx1;
	T idx2;
	T operator[](uint i) {
		//Pointer cast allows this to be branchless, just don't index greater than 1
		T* ptr_hack = (T*)&(this->idx1);
		return ptr_hack[i];
	}
};

struct RunInfo {
	StrType							structure_type = StrType::strtype_undef;
	RunType							run_type = RunType::runtype_undef;
	uint32_t						run_options = 0;
	double							rmsd_limit = 0;
	double							wc_rmsd_limit = 0;
	size_t							memory_limit = 0;
	size_t							build_limit = 0;
	size_t							n_fragments = 0;
	size_t							n_wc_pairs = 0;
	size_t							ds_strand1_n_frags = 0;
	size_t							ds_strand2_n_frags = 0;
	uint64_t						n_total_structs_built = 0;
	uint64_t						n_filter_structs_built = 0;
	char								sequence[GLOBAL_STANDARD_STRING_LENGTH] = {'\0'};
	char								library_prototype[GLOBAL_STANDARD_STRING_LENGTH] = {'\0'};
	char								wc_library_prototype[GLOBAL_STANDARD_STRING_LENGTH] = {'\0'};
	char								output_file[GLOBAL_STANDARD_STRING_LENGTH] = {'\0'};
	char								input_file[GLOBAL_STANDARD_STRING_LENGTH] = {'\0'};
	char								index_list_file[GLOBAL_STANDARD_STRING_LENGTH] = {'\0'};
	char								dot_bracket[GLOBAL_STANDARD_STRING_LENGTH] = {'\0'};
	char								parallel_lib_len[GLOBAL_STANDARD_STRING_LENGTH] = {'\0'};
	char								parallel_lib_idx[GLOBAL_STANDARD_STRING_LENGTH] = {'\0'};
	uint32_t*						index = nullptr;
	char**							fragment_lib_list = nullptr;
	char**							wc_lib_list = nullptr;
	size_t*							lib_duplicate_record = nullptr;
	size_t*							wc_lib_duplicate_record = nullptr;
	IndexPair<size_t>*	wc_pair_list = nullptr;
	IndexPair<size_t>*	frag_lib_bounds = nullptr;
	Timer								init_timer;
	Timer								run_timer;
	void print_stats() {
		printf("Total number of structures built: %lu\n", n_total_structs_built);
		if(run_options & RunOpts::use_structure_filter || run_options & RunOpts::build_limit_by_energy || run_options & RunOpts::blind_build_limit) {
			printf("Number of structures saved: %lu\n", n_filter_structs_built);
		}
	}
};
#endif
