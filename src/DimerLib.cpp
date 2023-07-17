#include "DimerLib.hpp"
#include <cstdlib>
/**
 * @brief Construct a DimerLib. Uses the number of atoms + extra elements, namely the steric clash centroids.
 * 
 * @param n Number of structures in a library
 * @param a_n Number of elements in each structure (atoms + extra)
 * @param e_n Number of extra elements in a_n
 */
void DimerLib::initialize(int n, int a_n, int e_n)
{
	LibraryMemblock = gsl_block_alloc(n * a_n * MATRIX_DIMENSION2);
	data_matrices = (gsl_matrix **)malloc(sizeof(gsl_matrix *) * n);
	for(int i = 0, offset = 0; i < n; i++, offset += a_n * MATRIX_DIMENSION2) {
		data_matrices[i] = gsl_matrix_alloc_from_block(LibraryMemblock, offset, a_n, MATRIX_DIMENSION2, MATRIX_DIMENSION2);
	}
	name = (char *)malloc(sizeof(char) * 3);
	energy = (float *)malloc(sizeof(float) * n);
	radii[0] = (double *)malloc(sizeof(double) *n);
	radii[1] = (double *)malloc(sizeof(double) *n);
	count = n;
	atom_data = new atom_info(a_n - e_n);
}

void DimerLib::destroy()
{
	for (int i = 0; i < count; i++) {
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
	for (int i = 0; i < count; i++) {
		data_matrices[i] = d_m[i];
	}
	memcpy(energy, e, sizeof(float) * count);
	memcpy(name, n, sizeof(char) * 3);
	memcpy(radii[0], r[0], sizeof(double) * count);
	memcpy(radii[1], r[1], sizeof(double) * count);
}

void DimerLib::clear_flags()
{
}

void DimerLibArray::initialize(int s)
{
	library = (DimerLib **)malloc(sizeof(DimerLib *) * s);
	is_duplicate = (bool*)malloc(sizeof(bool) * s);
	memset(is_duplicate, false, sizeof(bool) * s);
	Flags = (flag **)malloc(sizeof(flag *) * s);
	count = s;
	iterator = 0;
}

void DimerLibArray::destroy()
{
	for (size_t i = 0; i < count; i++) {
		if(is_duplicate[i] != true) {
			library[i]->destroy();
		}            
	}
	free(library);
	free(Flags[0]);
	free(Flags);
	free(is_duplicate);
}

DimerLib *DimerLibArray::operator[](int i)
{
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
	library[iterator] = (DimerLib *)malloc(sizeof(DimerLib));
	library[iterator]->initialize(n, a_n, e_n);
	full_structure_element_sum += a_n;
	if((a_n - e_n) > LargestAtomCount) {
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
	int LargestPositiveCount = 0, TmpPos = 0;
	int LargestNegativeCount = 0, TmpNeg = 0;
	int PosMapTracker = 0;
	int NegMapTracker = 0;
	for(uint32_t i = 0; i < count; i++) {
		TmpPos = 0;
		TmpNeg = 0;
		for(uint32_t j = 0; j < library[i]->atom_data->count; j++) {
			if(library[i]->atom_data->charges[j] == atom_charge::POSITIVE) {
				TmpPos++;
			}
			if(library[i]->atom_data->charges[j] == atom_charge::NEGATIVE) {
				TmpNeg++;
			}
		}
		if(TmpPos > LargestPositiveCount) {
			LargestPositiveCount = TmpPos;
		}       
		if(TmpNeg > LargestNegativeCount) {
			LargestNegativeCount = TmpNeg;
		}   
	}

	PositiveAtomMap = (int*)malloc(LargestAtomCount * count * sizeof(int));
	NegativeAtomMap = (int*)malloc(LargestAtomCount * count * sizeof(int));
	for(uint32_t i = 0; i < count; i++) {
		for(uint32_t j = 0; j < LargestAtomCount; j++) {
			if(j < library[i]->atom_data->count) {
				if(i != 0 && library[i]->atom_data->dnt_pos[j] != 2) {
					continue;
				}
				if(library[i]->atom_data->charges[j] == atom_charge::POSITIVE) {
					PositiveAtomMap[IDX_FLAT2D(i,j,LargestAtomCount)] = PosMapTracker;
					PosMapTracker++;
				} else if(library[i]->atom_data->charges[j] == atom_charge::NEGATIVE) {
					NegativeAtomMap[IDX_FLAT2D(i,j,LargestAtomCount)] = NegMapTracker;
					NegMapTracker++;
				} else {
					PositiveAtomMap[IDX_FLAT2D(i,j,LargestAtomCount)] = -1;
					NegativeAtomMap[IDX_FLAT2D(i,j,LargestAtomCount)] = -1;
				}
			}
		}
	}
	PositiveAtomCount = PosMapTracker;
	NegativeAtomCount = NegMapTracker;
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
			//printf("flags for %lu reset\n", i);
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
	char line[GLOBAL_STANDARD_STRING_LENGTH];
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
/**
 * @brief Calculates the centroid of each residue for the given gsl_matrix (coordinates) with its associated atom_info. This for use with our coarse
 * grain steric clash check method. It assumes the given matrix has enough room for 2 additional XYZ coordinates which correspond to the 
 * 2 residue centroids. Centroids are calculated using the entire nucleotide for both residues, meaning the phophate group belonging to residue 2
 * is also included for residue 1.
 * 
 * @param A gsl_matrix containing XYZ coordinates for fragment
 * @param A_info atom_info for the gsl_matrix coordinates
 */
void calculate_dnt_COM(gsl_matrix *A, atom_info *A_info)
{
	int n_r1 = A_info->count_per_res[0]; 
	int n_r2 = A_info->count_per_res[1];

	gsl_matrix *A_1 = gsl_matrix_alloc(n_r1 + 3, 3); //+3 for inclusion of phosphate group
	gsl_matrix *A_2 = gsl_matrix_alloc(n_r2, 3);

	constexpr atom_id phosphate_atoms[] = {P, OP1, OP2};

	double COMA_1[3] = {0, 0, 0};
	double COMA_2[3] = {0, 0, 0};

	for(int i = 0; i < n_r1; i++)
	{
		gsl_matrix_set(A_1, i, 0, gsl_matrix_get(A, i, 0));
		gsl_matrix_set(A_1, i, 1, gsl_matrix_get(A, i, 1));
		gsl_matrix_set(A_1, i, 2, gsl_matrix_get(A, i, 2));
	}
	for(int i = n_r1; i < n_r1 + 3; i++)
	{
		gsl_matrix_set(A_1, i , 0, gsl_matrix_get(A, A_info->get_idx_of_atom(phosphate_atoms[i - n_r1], 1), 0));
		gsl_matrix_set(A_1, i , 1, gsl_matrix_get(A, A_info->get_idx_of_atom(phosphate_atoms[i - n_r1], 1), 1));
		gsl_matrix_set(A_1, i , 2, gsl_matrix_get(A, A_info->get_idx_of_atom(phosphate_atoms[i - n_r1], 1), 2));
	}
	for(int i = n_r1; i < n_r1 + n_r2; i++)
	{
		gsl_matrix_set(A_2, i - n_r1, 0, gsl_matrix_get(A, i, 0));
		gsl_matrix_set(A_2, i - n_r1, 1, gsl_matrix_get(A, i, 1));
		gsl_matrix_set(A_2, i - n_r1, 2, gsl_matrix_get(A, i, 2));
	}

	get_matrix_COM(A_1, COMA_1);
	get_matrix_COM(A_2, COMA_2);

	gsl_matrix_set(A, A_info->count, 0, COMA_1[0]);
	gsl_matrix_set(A, A_info->count, 1, COMA_1[1]);
	gsl_matrix_set(A, A_info->count, 2, COMA_1[2]);

	gsl_matrix_set(A, A_info->count + 1, 0, COMA_2[0]);
	gsl_matrix_set(A, A_info->count + 1, 1, COMA_2[1]);
	gsl_matrix_set(A, A_info->count + 1, 2, COMA_2[2]);

	gsl_matrix_free(A_1);
	gsl_matrix_free(A_2);
}

/**
 * @brief Calculates Steric Clash COM (SCC) radius for SC Optimization. First the distance between the phosphate oxygens (OP1 and OP2) and the COM is calculated.
 * Then the same is done but for the functional groups of the residue (N6, N2, O6, O2, N4, O4). The greatest distance between COM and either functional group or
 * residue used to set the radius of the COM sphere. (Radius of the COM sphere = distance + 1)
 * 
 * @param A gsl_matrix containing XYZ coordinates for fragment
 * @param A_info Atom Info of fragment
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

/**
 * @brief This function ensures that all of the atoms lie within the sphere created by calculate_SCC_radii, if they are not 
 * it prints a warning and readjusts the size of the sphere to then include everything. calculate_SCC_radii method does
 * seem to work, but because this function is only run during initialization, I don't really care if its inefficient. 
 * 
 * @param A gsl_matrix containing XYZ coordinates for fragment
 * @param A_info atom_info of fragment
 * @param radius radius of sphere created by calculate_SCC_radii (modified as needed)
 * @param res_id residue index for radius
 */
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
			fprintf(stderr, "Warning: atom: %s in residue: %c is outside COM sphere!\n", A_info->name[i], A_info->residue[i]);
			*radius = dist + 1;
		}
	}
}

void prepare_run_info_lib_bounds(RunInfo& run_info, uint32_t lib_len, uint32_t idx) 
{
	uint32_t split_idx = run_info.frag_lib_bounds[idx].idx1;
	uint32_t split_len = run_info.frag_lib_bounds[idx].idx2;
	if(split_len == 0) {
		run_info.frag_lib_bounds[idx].idx1 = 0;
		run_info.frag_lib_bounds[idx].idx2 = lib_len;
		return;
	}
	uint32_t splits = ceil((float)lib_len/(float)split_len);
	run_info.frag_lib_bounds[idx].idx1 = split_len * split_idx;
	if(split_idx == splits - 1) {
		run_info.frag_lib_bounds[idx].idx2 = lib_len;
	} else {
		run_info.frag_lib_bounds[idx].idx2 = (split_len * (split_idx + 1));
	}
	//DEBUG_PRINT("Lib %d, BOUNDS %lu - %lu\n", idx, run_info.frag_lib_bounds[idx].idx1, run_info.frag_lib_bounds[idx].idx2);
}

/**
 * @breif Reads all files given in lib_files, assuming they are library files, and will load them into the given DimerLibArray
 * This function works for both WatsonCrick libraries and normal DNMP libraries, determined by bool for_WC.
 *
 * @param lib_files array of strings containing file location/name to read from.
 * @param n_libs number of DNMP (fragment) libraries to load.
 * @param &rtn_lib_array reference to  uninitialized DimerLibArray
 * @param duplicate_record array which maps any repeat libries to the index of the first occurance, e.g. there are two AA libraries, 
 * one at index 0, the next at index 4. In this array index [4] = 0. This is to stop any duplicate loading of libraries.
 * @param for_WC boolean value. if true: these libraries are for watson crick pairs, if false: these are for regular DNMP libs.
 **/
void load_libs(RunInfo& run_info, DimerLibArray &rtn_lib_array, bool for_WC)
{
	enum { model_count = 0, atom_count = 1 };
	float *energies;
	gsl_matrix **data_mats;
	int model_info[2];
	double* _radii[2];
	char line[GLOBAL_STANDARD_STRING_LENGTH];
	size_t* duplicate_record = for_WC ? run_info.wc_lib_duplicate_record : run_info.lib_duplicate_record;
	uint32_t n_libs = for_WC ? run_info.n_wc_pairs : run_info.n_fragments;	
	char** lib_files = for_WC ? run_info.wc_lib_list : run_info.fragment_lib_list;

	rtn_lib_array.initialize(n_libs);

	printf("Loading %s libraries...\n", for_WC ? "watson-crick pair" : "fragment");

	for (uint32_t i = 0; i < n_libs; i++) {
		model_info[model_count] = model_info[atom_count] = 0;
		FILE *lib_file = fopen(lib_files[i], "r");
		if (lib_file == NULL) {
			printf("Cannot open library file: %s\n", lib_files[i]);
			exit(3);
		}
		get_model_count(lib_file, model_info);
		if(!for_WC) {prepare_run_info_lib_bounds(run_info, model_info[model_count], i);}

		if(duplicate_record[i] != i) {
			rtn_lib_array.map_duplicate(i, duplicate_record[i]);
			printf("%d: %s (Points to library %lu)\n", i, lib_files[i], duplicate_record[i]);
			continue;
		}

		bool first_itr = true;
		printf("%d: %s ", i, lib_files[i]);
		printf("Models: %d, Atoms per model: %d\n", model_info[model_count], model_info[atom_count]);

		if(for_WC) {
			rtn_lib_array.alloc_lib(model_info[model_count], model_info[atom_count], 0); 
		} else {
			model_info[atom_count] += 2; //Plus 2 B/C COM for each DNT will be included in data matrix
			rtn_lib_array.alloc_lib(model_info[model_count], model_info[atom_count], 2); 
		}

		data_mats = rtn_lib_array[i]->data_matrices;
		energies  = rtn_lib_array[i]->energy;
		_radii[0] = rtn_lib_array[i]->radii[0];
		_radii[1] = rtn_lib_array[i]->radii[1];

		int iterator = 0;
		int row = 0;
		while (fgets(line, sizeof(line), lib_file)) {
			char line_origin[GLOBAL_STANDARD_STRING_LENGTH];
			char *header;

			strncpy(line_origin, line, GLOBAL_STANDARD_STRING_LENGTH);
			header = strtok(line, " ");

			if (!strcmp(header, "ENERGY")) {
				char *str1;
				str1 = strtok(NULL, " ");
				str1 = strtok(NULL, " ");
				energies[iterator] = atof(str1);
			}
			if ((!strcmp(header, "ATOM"))) {
				char *index, *name, *residue, *position, *X, *Y, *Z;
				index = strtok(NULL, " ");
				name = strtok(NULL, " ");
				residue = strtok(NULL, " ");
				position = strtok(NULL, " ");
				X = strtok(NULL, " ");
				Y = strtok(NULL, " ");
				Z = strtok(NULL, " ");

				if (first_itr) {
					rtn_lib_array.add_to_atom_info(name, atoi(index), *residue, atoi(position));
				}

				gsl_matrix_set(data_mats[iterator], row, 0, atof(X));
				gsl_matrix_set(data_mats[iterator], row, 1, atof(Y));
				gsl_matrix_set(data_mats[iterator], row, 2, atof(Z));
				row++;
			}
			else if (!strcmp(header, "ENDMDL\n")) {
				if (first_itr)
					first_itr = false;

				if(!for_WC) {
					atom_info* a_info = rtn_lib_array[rtn_lib_array.iterator]->atom_data;
					calculate_dnt_COM(data_mats[iterator], a_info);
					calculate_SCC_radii(data_mats[iterator], a_info, &_radii[0][iterator], 0);
					calculate_SCC_radii(data_mats[iterator], a_info, &_radii[1][iterator], 1);
					check_if_all_in_sphere(data_mats[iterator], a_info, &_radii[0][iterator], 0);
					check_if_all_in_sphere(data_mats[iterator], a_info, &_radii[1][iterator], 1);
				}     
				iterator++;
				row = 0;
			}
		}
		fclose(lib_file);

		const char *default_name = for_WC ? run_info.wc_library_prototype : run_info.library_prototype;
		size_t name_position = strcspn(default_name, "X");
		rtn_lib_array[i]->name[0] = for_WC ? lib_files[i][name_position] : lib_files[i][name_position];
		rtn_lib_array[i]->name[1] = for_WC ? lib_files[i][name_position + 1] : lib_files[i][name_position + 1];
		rtn_lib_array[i]->name[2] = '\0';
		rtn_lib_array.iterator++;    
	}
	rtn_lib_array.initialize_flags();
}
