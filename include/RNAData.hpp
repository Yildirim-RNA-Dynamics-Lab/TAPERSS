#ifndef RNA_DATA_HPP
#define RNA_DATA_HPP

#include "RNACMB.hpp"
#include "Atom_info.hpp"
#include "DimerLib.hpp"
#include "RNA_Math.hpp"
#include "HopcroftKarp.hpp"

constexpr atom_id targetA[] = {N9, C8, N7, C5, C6, N1, C2, N3, C4, C1p, C2p, C3p, C4p, O4p}; // 14
constexpr atom_id targetC[] = {N1, C2, N3, C4, C5, C6, C1p, C2p, C3p, C4p, O4p};             // 11
constexpr atom_id targetG[] = {N9, C8, N7, C5, C6, N1, C2, N3, C4, C1p, C2p, C3p, C4p, O4p}; // 14
constexpr atom_id targetU[] = {N1, C2, N3, C4, C5, C6, C1p, C2p, C3p, C4p, O4p};             // 11

//constexpr atom_id WCtargetA[] = {N9, C8, N7, C5, C6, N1, C2, N3, C4, N6, /* C1p, C2p, C3p, C4p, O4p */};
//constexpr atom_id WCtargetC[] = {N1, C2, N3, C4, C5, C6, N4, O2, /* C1p, C2p, C3p, C4p, O4p */};
//constexpr atom_id WCtargetG[] = {N9, C8, N7, C5, C6, N1, C2, N3, C4, N2, O6,  /* C1p, C2p, C3p, C4p, O4p */};
//constexpr atom_id WCtargetU[] = {N2, C2, N3, C4, C5, C6, O4, O2, /* C1p, C2p, C3p, C4p, O4p */};

constexpr atom_id WCtargetA[] = {N9, C8, N7, C5, C6, N1, C2, N3, /*C4, N6, C1p, C2p, C3p, C4p, O4p*/};
constexpr atom_id WCtargetC[] = {N1, C2, N3, C4, C5, C6, N4, O2, /*C1p, C2p, C3p, C4p, O4p*/};
constexpr atom_id WCtargetG[] = {N9, C8, N7, C5, C6, N1, C2, N3, /*C4, N2, O6, C1p, C2p, C3p, C4p, O4p*/};
constexpr atom_id WCtargetU[] = {N1, C2, N3, C4, C5, C6, O4, O2, /*C1p, C2p, C3p, C4p, O4p*/};

struct RNAData
{
	gsl_matrix* data_matrix;        // gsl_matrix containing XYZ coordinates for each atom
	uint16_t *StericIndices[2];     // Phosphate backbone indices for SCC with neighboring nucleotide.
	uint16_t *EnergyIndices[4];     // [0],[1]: Positive and Negative atom indices for Res 1. [2],[3]: Positive and Negative for Res 2
	uint16_t *submatrix_rows[2];    // Atom indices for residue [0] and [1] which correspond to the respective target for Kabsch overlap (defined above)
	uint16_t *WC_submatrix_rows[2]; // Same as ^ but for WatsonCrick overlap
	atom_info *atom_data;           // Contains (non-repeating) info about each atom excluding coordinate info.
	float energy;                   // Energy of model (from library)
	char name[3];                   // DNMP name (ex. "AA\0") (this might be removed..)
	flag *_flag;                    // Pointer to flag to update library, not necessary for combinatorial but will be useful for other sampling methods
	uint32_t count;                 // Number of atoms in structure

	int position_in_lib[2];         // What library and model data originates from
	int count_per_sub[2];           // How many atoms are in each "submatrix" (submatrix_rows)
	int count_per_Steric[2];        // How many atoms are in each StericIndices
	size_t ResBoundaries[4];        // [0]= Start of residue 1 row indices, [1]= End of Res 1 row idx. [2],[3]: Start, End of Res 2 row idx
	size_t count_per_Energy[4];     // How many atoms are in each EnergyIndices
	size_t count_per_WC_sub[2];     // How many atoms are in each WC_submatrix

	double COM_Radii[2];            // Radius for each COM for course grained SCC

	const atom_id *target1;         // Pointer to target for residue 1
	const atom_id *target2;         // Same as ^ for residue 2
	const atom_id *WC_target1;      // Pointer to WC target for residue 1
	const atom_id *WC_target2;      // Same as ^ for residue 2

	long int id;                    // id to keep track of allocations and deallocations, used for debugging only.

	//RNAData();
	void initialize(DimerLibArray &L, int idx, int idx_L, gsl_block *MemBlock, size_t *offset_matrix, uint16_t *ArrayMemBlock, size_t *offset_array);
	void overwrite(DimerLibArray &L, int i, int j);
	void destroy();
	//~RNAData();

	int get_residue_COM_index(int res);
	void print();
	void print_target(int res);
	void print_WC_target(int res);
	void reset_interactions();
	void print(int res);
	int to_string(char *s, int buffer_size, int string_index);
	void print_offset(int res, int position);
	int to_string_offset(int res, int position, char *s, int buffer_size, int string_index, uint32_t& idx_offset);
	int generate_string_prototype(int res, int position, char *s, int buffer_size, int string_index, uint32_t& idx_offset);
	int overwrite_string_prototype(int res, char *s, int buffer_size, int string_index, uint32_t& idx_offset);
};

struct DoubleStrandData
{
	int iterator_max_Strand1;
	int strand1_size;
	uint16_t *wc_rows[2];
	size_t wc_rowsize[2];
};

struct RNADataArray
{
	RNAData *sequence;
	gsl_matrix *InteractionTable;
	gsl_matrix *COMS;
	gsl_block *MatrixMemBlock;
	uint16_t *ArrayMemBlock;
	StrType structure_type;

	double *Radii;
	bool *PassedCOMCheck;
	int *PositiveAtomMap;
	int *NegativeAtomMap;
	uint32_t LAC;
	int TableRowCount;
	int TableColCount;
	int count;
	int iterator;
	int iterator_max;
	float structure_energy;
	DoubleStrandData ds_data;
	IndexPair<size_t>* wc_pairs;
	IndexPair<size_t>* wc_pair_res;
	float* wc_rmsd;
	size_t n_wc_pairs;

	char *string_out;
	int string_buffer;
	int string_index;
	int header_pos_record[5];
	int header_length;
	int model_count;
	bool string_initialized = false;
	bool string_print_coordinates = GLOBAL_WRITE_COORDINATES;
	bool TMP_END = false;

	int atom_sum;

	void initialize(RunInfo& run_info, DimerLibArray &library, DimerLibArray& wc_library);
	void initialize_ds(RunInfo& run_info, DimerLibArray &wc_library, size_t& arr_offset);
	void overwrite(size_t LibIdx, size_t IdxInLib, DimerLibArray &Library);
	void overwrite_initialize(size_t LibIdx, size_t IdxInLib, DimerLibArray &Library);
	void destroy();
	uint64_t calculate_matrix_memory_needed(DimerLibArray &L, int i);
	uint64_t calculate_index_array_memory_needed(DimerLibArray &L, int i);
	RNAData *operator[](int i);
	void add_copy(RNAData *A);
	void add_move(RNAData *A);
	RNAData *current();
	inline bool is_complete();
	void keep();
	void rollback();
	void safe_rollback(); // unused
	void rollback_by(int amount);
	bool is_empty();
	bool should_prepare_right(int working_pos); 
	void update_WC_rmsd(float rmsd_val, size_t idx);
	void update_energy();
	void printall();
	void initialize_string();
	int get_atom_sum();
	int out_string_header_coord(RunInfo& run_info);
	int out_string_header(RunInfo& run_info);
	void generate_string_prototype(RunInfo& run_info); 
	void gen_string_prototype_ds(uint32_t& idx_offset);
	void gen_string_prototype_ss(uint32_t& idx_offset);
	char* overwrite_string_prototype(RunInfo& run_info); 
	void overwrite_string_prototype_ds(uint32_t& idx_offset);
	void overwrite_string_prototype_ss(uint32_t& idx_offset);
	int *get_index();
	void print_index(int offset);
	void print_index();
};

size_t get_target(char res, const atom_id **dest);
size_t get_WC_target(char res, const atom_id **dest);
int FindInteraction(int A_ResId, int A_Idx, RNAData *AData, int B_ResId, int B_Idx, RNAData *BData,
									gsl_matrix *AdjMatrix, int *PAdjMap, int *NAdjMap, size_t DIM2);

inline bool RNADataArray::is_complete()
{
	return (iterator == iterator_max);
}

#endif
