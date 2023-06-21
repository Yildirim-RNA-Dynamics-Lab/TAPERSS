#include "RNAData.hpp"
#include "RNACMB.hpp"
#include "RNADataArrayInternalLoop.hpp"
#include <sys/types.h>

void RNAData::initialize(DimerLibArray &L, int idx, int idx_L, gsl_block *MemBlock, size_t *offset_matrix, uint16_t *ArrayMemBlock, size_t *offset_array)
{
	int rel_idx1 = 0;
	int rel_idx2 = 0;
	int rel_idx3 = 0;
	int rel_idx4 = 0;

	atom_data = (L[idx]->atom_data);
	memcpy(name, L[idx]->name, sizeof(char) * 3);
	energy = L[idx]->energy[idx_L];
	count = L[idx]->atom_data->count;
	_flag = &(L.Flags[idx][idx_L]);
	position_in_lib[0] = idx;
	position_in_lib[1] = idx_L;
	ResBoundaries[0] = 0;
	ResBoundaries[1] = atom_data->count_per_res[0];
	ResBoundaries[2] = atom_data->count_per_res[0];
	ResBoundaries[3] = atom_data->count;

	id = rna_dat_tracker++;

	COM_Radii[0] = L[idx]->radii[0][idx_L];
	COM_Radii[1] = L[idx]->radii[1][idx_L];
	
	data_matrix = gsl_matrix_alloc_from_block(MemBlock, *offset_matrix, L[idx]->data_matrices[idx_L]->size1, L[idx]->data_matrices[idx_L]->size2, MATRIX_DIMENSION2);
	gsl_matrix_memcpy(data_matrix, L[idx]->data_matrices[idx_L]);

	*offset_matrix += L[idx]->data_matrices[idx_L]->size1 * L[idx]->data_matrices[idx_L]->size2;

	count_per_sub[0] = get_target(name[0], &target1);

	submatrix_rows[0] = &ArrayMemBlock[*offset_array];
	*offset_array += count_per_sub[0];

	count_per_sub[1] = get_target(name[1], &target2);

	submatrix_rows[1] = &ArrayMemBlock[*offset_array];
	*offset_array += count_per_sub[1];

	for (unsigned int j = 0; j < count; j++)
	{
		for (int k = 0; k < count_per_sub[0]; k++)
		{
			if ((target1[k] == atom_data->atom_ids[j]) && (atom_data->dnt_pos[j] - 1) == 0)
			{
				submatrix_rows[0][k] = j;
				rel_idx1++;
			}
		}
		for(int k = 0; k < count_per_sub[1]; k++)
		{
			if ((target2[k] == atom_data->atom_ids[j]) && (atom_data->dnt_pos[j] - 1) == 1)
			{
				submatrix_rows[1][k] = j;
				rel_idx2++;
			}
		}
	}
	count_per_WC_sub[0] = get_WC_target(name[0], &WC_target1);
	WC_submatrix_rows[0] = &ArrayMemBlock[*offset_array];
	*offset_array += count_per_WC_sub[0];

	count_per_WC_sub[1] = get_WC_target(name[1], &WC_target2);
	WC_submatrix_rows[1] = &ArrayMemBlock[*offset_array];
	*offset_array += count_per_WC_sub[1];

	rel_idx1 = 0;
	rel_idx2 = 0;

	for (unsigned int j = 0; j < count; j++)
	{
		for (unsigned int k = 0; k < count_per_WC_sub[0]; k++)
		{
			if ((WC_target1[k] == atom_data->atom_ids[j]) && (atom_data->dnt_pos[j] - 1) == 0)
			{
				WC_submatrix_rows[0][k] = j;
				rel_idx1++;
			}
		}
		for (unsigned int k = 0; k < count_per_WC_sub[1]; k++)
		{
			if ((WC_target2[k] == atom_data->atom_ids[j]) && (atom_data->dnt_pos[j] - 1) == 1)
			{
				WC_submatrix_rows[1][k] = j;
				rel_idx2++;
			}
		}
	}

	StericIndices[0] = &ArrayMemBlock[*offset_array];

	rel_idx1 = 0;
	rel_idx2 = 0;

	for(int i = 0; i < L[idx]->atom_data->count; i++)
	{
		if(L[idx]->atom_data->dnt_pos[i] == 2)
		{
			continue;
		}
		if( L[idx]->atom_data->atom_ids[i] == P   || L[idx]->atom_data->atom_ids[i] == OP1 || 
				L[idx]->atom_data->atom_ids[i] == OP2 || L[idx]->atom_data->atom_ids[i] == C5p || 
				L[idx]->atom_data->atom_ids[i] == O5p)
		{
			StericIndices[0][rel_idx1++] = i;
		}
	}

	count_per_Steric[0] = rel_idx1;
	*offset_array += rel_idx1;
	StericIndices[1] = &ArrayMemBlock[*offset_array];


	for(int i = 0; i < L[idx]->atom_data->count; i++)
	{
		if(L[idx]->atom_data->dnt_pos[i] == 1)
		{
			continue;
		}
		if( L[idx]->atom_data->atom_ids[i] == P   || L[idx]->atom_data->atom_ids[i] == OP1 || 
				L[idx]->atom_data->atom_ids[i] == OP2 || L[idx]->atom_data->atom_ids[i] == C5p || 
				L[idx]->atom_data->atom_ids[i] == O5p)
		{
			StericIndices[1][rel_idx2++] = i;
		}
	}

	count_per_Steric[1] = rel_idx2;
	*offset_array += rel_idx2;

	EnergyIndices[0] = &ArrayMemBlock[*offset_array];

	rel_idx1 = 0;
	rel_idx2 = 0;

	for(int i = 0; i < L[idx]->atom_data->count; i++)
	{
		if( L[idx]->atom_data->charges[i] != atom_charge::POSITIVE)
		{
			continue;
		}
		else if(L[idx]->atom_data->dnt_pos[i] == 1)
		{
			EnergyIndices[0][rel_idx1++] = i;
		}
	}

	count_per_Energy[0] = rel_idx1;
	*offset_array += rel_idx1;
	EnergyIndices[1] = &ArrayMemBlock[*offset_array];

	for(int i = 0; i < L[idx]->atom_data->count; i++)
	{
		if( L[idx]->atom_data->charges[i] != atom_charge::NEGATIVE)
		{
			continue;
		}
		else if(L[idx]->atom_data->dnt_pos[i] == 1)
		{
			EnergyIndices[1][rel_idx2++] = i;
		}
	}

	count_per_Energy[1] = rel_idx2;
	*offset_array += rel_idx2;
	EnergyIndices[2] = &ArrayMemBlock[*offset_array];

	for(int i = 0; i < L[idx]->atom_data->count; i++)
	{
		if( L[idx]->atom_data->charges[i] != atom_charge::POSITIVE)
		{
			continue;
		}
		else if(L[idx]->atom_data->dnt_pos[i] == 2)
		{
			EnergyIndices[2][rel_idx3++] = i;
		}
	}

	count_per_Energy[2] = rel_idx3;
	*offset_array += rel_idx3;
	EnergyIndices[3] = &ArrayMemBlock[*offset_array];

	for(int i = 0; i < L[idx]->atom_data->count; i++)
	{
		if( L[idx]->atom_data->charges[i] != atom_charge::NEGATIVE)
		{
			continue;
		}
		else if(L[idx]->atom_data->dnt_pos[i] == 2)
		{
			EnergyIndices[3][rel_idx4++] = i;
		}
	}
	count_per_Energy[3] = rel_idx4;
	*offset_array += rel_idx4;
}

void RNAData::overwrite(DimerLibArray &L, int i, int j)
{
	gsl_matrix_memcpy(data_matrix, L[i]->data_matrices[j]);
	energy = L[i]->energy[j];
	_flag = &(L.Flags[i][j]);
	position_in_lib[0] = i;
	position_in_lib[1] = j;    
}

void RNAData::destroy()
{
	gsl_matrix_free(data_matrix);
}

size_t get_target(char res, const atom_id** dest)
{
	if (res == 'A')
	{
		*dest = targetA;        
		return (sizeof(targetA) / sizeof(atom_id));
	}
	else if (res == 'C')
	{
		*dest = targetC;
		return (sizeof(targetC) / sizeof(atom_id));
	}
	else if (res == 'G')
	{
		*dest = targetG;        
		return (sizeof(targetG) / sizeof(atom_id));
	}
	else if (res == 'U')
	{
		*dest = targetU;
		return (sizeof(targetU) / sizeof(atom_id));
	}
	else
	{
		printf("RNA_data: input file has problems, please check at %d\n", __LINE__);
		exit(1);
	}
}

size_t get_WC_target(char res, const atom_id** dest)
{
	if (res == 'A')
	{
		*dest = WCtargetA;        
		return (sizeof(WCtargetA) / sizeof(atom_id));
	}
	else if (res == 'C')
	{
		*dest = WCtargetC;
		return (sizeof(WCtargetC) / sizeof(atom_id));
	}
	else if (res == 'G')
	{
		*dest = WCtargetG;
		return (sizeof(WCtargetG) / sizeof(atom_id));
	}
	else if (res == 'U')
	{
		*dest = WCtargetU;
		return (sizeof(WCtargetU) / sizeof(atom_id));
	}
	else
	{
		printf("RNA_data: input file has problems, please check at %d\n", __LINE__);
		exit(1);
	}
}

int RNAData::get_residue_COM_index(int res)
{
	return (count + res);
}

void RNAData::print()
{
	for (unsigned int i = 0; i < data_matrix->size1; i++)
	{
		atom_data->print_at(i);
		for (unsigned int j = 0; j < data_matrix->size2; j++)
		{
			printf("%8.3f", gsl_matrix_get(data_matrix, i, j));
		}
		putchar('\n');
	}
}

void RNAData::print_target(int res)
{
	for (int i = 0; i < count_per_sub[res]; i++)
	{
		atom_data->print_at(submatrix_rows[res][i]);
		for (unsigned int j = 0; j < data_matrix->size2; j++)
		{
			printf("%8.3f", gsl_matrix_get(data_matrix, submatrix_rows[res][i], j));
		}
		putchar('\n');
	}
}

void RNAData::print_WC_target(int res)
{
	for (unsigned int i = 0; i < count_per_WC_sub[res]; i++)
	{
			atom_data->print_at(WC_submatrix_rows[res][i]);
			for (unsigned int j = 0; j < data_matrix->size2; j++)
			{
				printf("%8.3f", gsl_matrix_get(data_matrix, WC_submatrix_rows[res][i], j));
			}
			putchar('\n');
	}
}

void RNAData::print(int res)
{
	for (unsigned int i = 0; i < data_matrix->size1; i++)
	{
		if (atom_data->dnt_pos[i] != (res + 1))
		{
			continue;
		}
		atom_data->print_at(i);
		for (unsigned int j = 0; j < data_matrix->size2; j++)
		{
			printf("%8.3f", gsl_matrix_get(data_matrix, i, j));
		}
		putchar('\n');
	}
}

int RNAData::to_string(char *s, int buffer_size, int string_index)
{
	for (unsigned int i = 0; i < count; i++)
	{
		string_index += snprintf(&s[string_index], buffer_size - string_index, "%-6s%5d %4s %3c  %4d    ", "ATOM", atom_data->index[i], atom_data->name[i], atom_data->residue[i], atom_data->dnt_pos[i]);
		for (unsigned int j = 0; j < data_matrix->size2; j++)
		{
			string_index += snprintf(&s[string_index], buffer_size - string_index, "%8.3f", gsl_matrix_get(data_matrix, i, j));
		}
		string_index += snprintf(&s[string_index], buffer_size - string_index, "\n");
	}
	return string_index;
}

void RNAData::print_offset(int res, int position)
{
	for (unsigned int i = 0; i < count; i++)
	{
		if (atom_data->dnt_pos[i] != (res + 1))
		{
			continue;
		}
		printf("%-4s %-4d %-4c %-4d ", atom_data->name[i], atom_data->index[i], atom_data->residue[i], atom_data->dnt_pos[i] + position);
		for (unsigned int j = 0; j < data_matrix->size2; j++)
		{
			printf("%-3.2f ", gsl_matrix_get(data_matrix, i, j));
		}
		putchar('\n');
	}
}

int RNAData::to_string_offset(int res, int position, char *s, int buffer_size, int string_index, uint32_t& idx_offset)
{
	for (unsigned int i = 0; i < count; i++)
	{
		if (atom_data->dnt_pos[i] != (res + 1))
		{
			continue;
		}
		idx_offset += 1;
		string_index += snprintf(&s[string_index], buffer_size - string_index, "%4s  %5d %-4s %3c  %4d    ", "ATOM", idx_offset, atom_data->name[i], atom_data->residue[i], atom_data->dnt_pos[i] + position);
		for (unsigned int j = 0; j < data_matrix->size2; j++)
		{
			string_index += snprintf(&s[string_index], buffer_size - string_index, "%8.3f", gsl_matrix_get(data_matrix, i, j));
		}
		string_index += snprintf(&s[string_index], buffer_size - string_index, "\n");
	}
	return string_index;
}

int RNAData::generate_string_prototype(int res, int position, char *s, int buffer_size, int string_index, uint32_t& idx_offset)
{
	for (unsigned int i = 0; i < count; i++)
	{
		if (atom_data->dnt_pos[i] != (res + 1))
		{
			continue;
		}
		idx_offset += 1;
		string_index += snprintf(&s[string_index], buffer_size - string_index, "%4s  %5d %-4s %3c  %4d    ", "ATOM", idx_offset, atom_data->name[i], atom_data->residue[i], atom_data->dnt_pos[i] + position);
		string_index += snprintf(&s[string_index], buffer_size - string_index, "%8.3f%8.3f%8.3f\n", 0.0F,0.0F,0.0F);
	}
	return string_index;
}

int RNAData::overwrite_string_prototype(int res, char *s, int buffer_size, int string_index, uint32_t& idx_offset)
{
	for(u_int32_t i = 0; i < count; i++)
	{
		if (atom_data->dnt_pos[i] != (res + 1))
		{
			continue;
		}
		idx_offset += 1;
		string_index += 30;
		string_index += snprintf(&s[string_index], buffer_size - string_index, "%8.3f%8.3f%8.3f", gsl_matrix_get(data_matrix,i,0),gsl_matrix_get(data_matrix,i,1),gsl_matrix_get(data_matrix,i,2));
		s[string_index] = '\n'; //snprintf will always attach a '\0' to the end of the string it creates. So we replace '\0' with '\n'
		string_index++;
	}
	return string_index;
}

void RNADataArray::initialize_ds(RunInfo& run_info, DimerLibArray& wc_library, size_t& arr_offset)
{
	const atom_id* wc_target1;
	const atom_id* wc_target2;

	ds_data.wc_rowsize[0] = get_target(wc_library[0]->name[0], &wc_target1);
	ds_data.wc_rowsize[1] = get_target(wc_library[0]->name[1], &wc_target2);
	ds_data.strand1_size = run_info.ds_strand1_n_frags;
	ds_data.iterator_max_Strand1 = run_info.ds_strand1_n_frags - 1;
	ds_data.wc_rows[0] = &ArrayMemBlock[arr_offset];
	arr_offset += ds_data.wc_rowsize[0];
	ds_data.wc_rows[1] = &ArrayMemBlock[arr_offset];
	for(int i = 0; i < wc_library[1]->atom_data->count; i++)
	{
		for(size_t j = 0; j < ds_data.wc_rowsize[0]; j++) 
		{
			if(wc_library[1]->atom_data->atom_ids[i] == wc_target1[j] && wc_library[1]->atom_data->dnt_pos[i] == 1)
			{
				ds_data.wc_rows[0][j] = (uint16_t)i;
			}
		}
		for(size_t j = 0; j < ds_data.wc_rowsize[1]; j++) 
		{
			if(wc_library[1]->atom_data->atom_ids[i] == wc_target2[j] && wc_library[1]->atom_data->dnt_pos[i] == 2)
			{
				ds_data.wc_rows[1][j] = (uint16_t)i;
			}
		}
	}
}

void RNADataArray::initialize(RunInfo& run_info, DimerLibArray& library, DimerLibArray& wc_library)
{
	int PLAC = library.PositiveAtomCount;
	int NLAC = library.NegativeAtomCount;
	iterator_max = run_info.n_fragments - 1;
	sequence = (RNAData *)malloc(sizeof(RNAData) * run_info.n_fragments);
	structure_type = run_info.structure_type;
	n_wc_pairs = run_info.n_wc_pairs;

	size_t matrix_memsize = 0;
	size_t array_memsize = 0;
	for(size_t i = 0; i < run_info.n_fragments; i++)
	{
		matrix_memsize   += calculate_matrix_memory_needed(library, i);
		array_memsize    += calculate_index_array_memory_needed(library, i);
	}
	matrix_memsize += (PLAC * NLAC); // Interaction Table Matrix
	array_memsize += sizeof(uint16_t) * run_info.n_fragments; // Only needed for Double Strand.
	MatrixMemBlock = gsl_block_alloc(matrix_memsize); //Allocate Block of memory to pass to individual RNAData's. This is to maintain contiguous memory.
	ArrayMemBlock = (uint16_t *)malloc(array_memsize); // Same as ^.

	size_t matrix_offset = 0;
	size_t array_offset = 0;

	for(size_t i = 0; i < run_info.n_fragments; i++)
	{
		sequence[i].initialize(library, i, 0, MatrixMemBlock, &matrix_offset, ArrayMemBlock, &array_offset);
	}

	if(structure_type == StrType::double_strand) { 
		initialize_ds(run_info, wc_library, array_offset);
	}

	*sequence[0]._flag = USED;
	iterator = 0;  //Start with first already DNT included
	count = 1;
	structure_energy = 0;
	COMS = gsl_matrix_alloc_from_block(MatrixMemBlock, matrix_offset, 2 * run_info.n_fragments, MATRIX_DIMENSION2, MATRIX_DIMENSION2);
	matrix_offset += (2 * run_info.n_fragments * MATRIX_DIMENSION2);
	InteractionTable = gsl_matrix_alloc_from_block(MatrixMemBlock, matrix_offset, PLAC, NLAC, NLAC);
	gsl_matrix_set_zero(InteractionTable); // Set all to 0
	PositiveAtomMap = library.PositiveAtomMap;
	NegativeAtomMap = library.NegativeAtomMap;
	TableRowCount = PLAC;
	TableColCount = NLAC;
	LAC = library.LargestAtomCount;

	uint32_t COMCheckCount = structure_type == StrType::double_strand ? run_info.n_fragments + 2 : run_info.n_fragments + 1;
	Radii = (double *)malloc(2 * run_info.n_fragments * sizeof(double));
	PassedCOMCheck = (bool*)malloc(COMCheckCount * sizeof(bool)); 
	if(n_wc_pairs != 0)
	{
		wc_rmsd = (float *)malloc(sizeof(float) * n_wc_pairs);
		wc_pairs = (IndexPair<size_t>*)malloc(sizeof(IndexPair<size_t>) * n_wc_pairs);
		wc_pair_res = (IndexPair<size_t>*)malloc(sizeof(IndexPair<size_t>) * n_wc_pairs);
		for(uint i = 0; i < n_wc_pairs; i++) {
			size_t idx1 = run_info.wc_pair_list[i].idx1;
			size_t idx2 = run_info.wc_pair_list[i].idx2;
			wc_pairs[i].idx1 = idx1 < 2 ? 0 : idx1 - 1;
			wc_pairs[i].idx2 = idx2 < 2 ? 0 : idx2 - 1;
			wc_pair_res[i].idx1 = idx1 == 0 ? 0 : 1;
			wc_pair_res[i].idx2 = idx2 == 0 ? 0 : 1;
			wc_rmsd[i] = 0.00;
		}
	}
	overwrite_initialize(0,run_info.frag_lib_bounds[0].idx1, library);
	generate_string_prototype(run_info);
}

void RNADataArray::overwrite(size_t LibIdx, size_t IdxInLib, DimerLibArray &library)
{
	sequence[LibIdx].overwrite(library, LibIdx, IdxInLib);
	Radii[LibIdx * 2] = library[LibIdx]->radii[0][IdxInLib];
	Radii[LibIdx * 2 + 1] = library[LibIdx]->radii[1][IdxInLib];
}

void RNADataArray::overwrite_initialize(size_t LibIdx, size_t IdxInLib, DimerLibArray &library)
{
	sequence[LibIdx].overwrite(library, LibIdx, IdxInLib);
	gsl_matrix_row_copy(COMS, LibIdx * 2, library[LibIdx]->data_matrices[IdxInLib], library[LibIdx]->atom_data->count + 0);
	gsl_matrix_row_copy(COMS, (LibIdx * 2) + 1, library[LibIdx]->data_matrices[IdxInLib], library[LibIdx]->atom_data->count + 1);
	Radii[LibIdx * 2] = library[LibIdx]->radii[0][IdxInLib];
	Radii[LibIdx * 2 + 1] = library[LibIdx]->radii[1][IdxInLib];
}

void RNADataArray::destroy()
{
	for(int i = 0; i < count; i++)
	{
		sequence[i].destroy();
	}
	free(sequence);
	gsl_matrix_free(COMS);
	free(Radii);
	free(PassedCOMCheck);
	if (string_initialized)
	{	
		free(string_out);
	}
	gsl_matrix_free(InteractionTable);
	free(PositiveAtomMap);
	free(NegativeAtomMap);
	gsl_block_free(MatrixMemBlock);
	free(ArrayMemBlock);
}

uint64_t RNADataArray::calculate_matrix_memory_needed(DimerLibArray& L, int idx)
{
	uint64_t memsize = 0;
	size_t data_mat_nrows = L[idx]->data_matrices[0]->size1;

	memsize += data_mat_nrows * MATRIX_DIMENSION2;
	memsize += 2 * MATRIX_DIMENSION2; //For COM Matrix
	//IDK why there was doubled up lines but I commented out, will test later...
	//memsize += data_mat_nrows * MATRIX_DIMENSION2;
	//memsize += 2 * MATRIX_DIMENSION2; //For COM Matrix

	return memsize;
}

uint64_t RNADataArray::calculate_index_array_memory_needed(DimerLibArray& L, int idx)
{
	uint64_t memsize = 0;
	const atom_id *dummy = NULL;

	/* Size for Submatrices Rows*/
	memsize += get_target(L[idx]->name[0], &dummy);
	memsize += get_target(L[idx]->name[1], &dummy);

	/* Size for WCSubmatrices Rows */
	memsize += get_WC_target(L[idx]->name[0], &dummy);
	memsize += get_WC_target(L[idx]->name[1], &dummy);

	/* Size for Steric Rows */
	for(int i = 0; i < L[idx]->atom_data->count; i++)
	{
		if( L[idx]->atom_data->atom_ids[i] == P   || L[idx]->atom_data->atom_ids[i] == OP1 || 
				L[idx]->atom_data->atom_ids[i] == OP2 || L[idx]->atom_data->atom_ids[i] == C5p ||
				L[idx]->atom_data->atom_ids[i] == O5p)
		{
			memsize++;
		}
	}

	/* Size for Charged Rows */
	for(int i = 0; i < L[idx]->atom_data->count; i++)
	{
		if( L[idx]->atom_data->charges[i] != atom_charge::NEUTRAL)
		{
			memsize++;
		}
	}
	memsize *= sizeof(uint16_t);

	return memsize;
}

RNAData *RNADataArray::operator[](int i)
{
	return &sequence[i];
}

/* Unused and should not be used for any reason (will remove at some point)*/
void RNADataArray::add_copy(RNAData *A)
{
	sequence[++iterator] = *A;
	count++;
}

void RNADataArray::add_move(RNAData *A)
{
	iterator++;
	gsl_matrix_row_copy(COMS, iterator * 2, A->data_matrix, A->get_residue_COM_index(0));
	gsl_matrix_row_copy(COMS, iterator * 2 + 1, A->data_matrix, A->get_residue_COM_index(1));
	Radii[iterator * 2] = A->COM_Radii[0];
	Radii[iterator * 2 + 1] = A->COM_Radii[1];
	sequence[iterator] = *A;
	count++;
	*A->_flag = USED;
}

RNAData *RNADataArray::current()
{
	return &sequence[iterator];
}

void RNADataArray::keep()
{
	iterator++;	
	gsl_matrix_row_copy(COMS, iterator * 2, sequence[iterator].data_matrix, sequence[iterator].get_residue_COM_index(0));
	gsl_matrix_row_copy(COMS, (iterator * 2) + 1, sequence[iterator].data_matrix, sequence[iterator].get_residue_COM_index(1));
	count++;
	*sequence[iterator]._flag = USED;
}

void RNADataArray::rollback()
{
	iterator--;
	count--;
}

void RNADataArray::safe_rollback() // unused
{
	iterator--;
	count--;
}

void RNADataArray::rollback_by(int amount)
{
	iterator -= amount;
	count -= amount;
}

bool RNADataArray::is_empty()
{
	return (count == 0);
}

bool RNADataArray::should_prepare_right(int working_pos) //return true if in left strand
{
	return (working_pos == ds_data.iterator_max_Strand1 + 1);
}

void RNADataArray::update_WC_rmsd(float rmsd_val, size_t idx)
{
	wc_rmsd[idx] = rmsd_val;
}

void RNADataArray::update_energy()
{
	float energy_ = 0.0;
	int Interactions = 0;

	for (int i = 0; i < count; i++)
	{
		energy_ += sequence[i].energy;
	}

	for(int i = count - 1; i > 1; i--)
	{
		Interactions += FindInteraction(0, 0, &sequence[0], 1, i, &sequence[i],
				InteractionTable, PositiveAtomMap, NegativeAtomMap, LAC);
		Interactions += FindInteraction(1, 0, &sequence[0], 1, i, &sequence[i],
				InteractionTable, PositiveAtomMap, NegativeAtomMap, LAC);
	}

	for (int i = 1; i <= count - 3; i++)
	{
		for (int j = count - 1; j > i + 1; j--)
		{
			Interactions += FindInteraction(1, i, &sequence[i], 1, j, &sequence[j],
					InteractionTable, PositiveAtomMap, NegativeAtomMap, LAC);
		}
	}

	if(Interactions > 1)
	{
		Interactions = HK_GetMaxMatching(InteractionTable);
	}
	energy_ -= Interactions;
	structure_energy = energy_;

	gsl_matrix_set_zero(InteractionTable);
}

int FindInteraction(int A_ResId, int A_Idx, RNAData *AData, int B_ResId, int B_Idx, RNAData *BData, gsl_matrix *AdjMatrix, int *PAdjMap, int*NAdjMap, size_t DIM2)
{
	gsl_matrix *A = AData->data_matrix;
	gsl_matrix *B = BData->data_matrix;
	uint16_t *A_P_Rows = AData->EnergyIndices[A_ResId * 2];
	uint16_t *A_N_Rows = AData->EnergyIndices[A_ResId * 2 + 1];
	uint16_t *B_P_Rows = BData->EnergyIndices[B_ResId * 2];
	uint16_t *B_N_Rows = BData->EnergyIndices[B_ResId * 2 + 1];
	size_t A_P_NRow = AData->count_per_Energy[A_ResId * 2];
	size_t A_N_NRow = AData->count_per_Energy[A_ResId * 2 + 1];
	size_t B_P_NRow = BData->count_per_Energy[B_ResId * 2];
	size_t B_N_NRow = BData->count_per_Energy[B_ResId * 2 + 1];
	size_t IDX1, IDX2;
	int Interactions = 0;
	for(size_t i = 0; i < A_P_NRow; i++)
	{
		for(size_t j = 0; j < B_N_NRow; j++)
		{
			if(distance_mat2mat(A, A_P_Rows[i], B, B_N_Rows[j]) < INTERACTION_DISTANCE)
			{
				//printf("Interaction between Resid: %d, atom: %s and Resid: %d, atom %s\n", AData->atom_data->dnt_pos[A_P_Rows[i]] - 1 + AData->position_in_lib[0], AData->atom_data->name[A_P_Rows[i]], BData->atom_data->dnt_pos[B_N_Rows[j]] - 1 + BData->position_in_lib[0], BData->atom_data->name[B_N_Rows[j]]);
				IDX1 = IDX_FLAT2D(A_Idx, A_P_Rows[i], DIM2);
				IDX2 = IDX_FLAT2D(B_Idx, B_N_Rows[j], DIM2);
				gsl_matrix_set(AdjMatrix, PAdjMap[IDX1], NAdjMap[IDX2], 1);
				Interactions++;
			}
		}
	}
	for(size_t i = 0; i < A_N_NRow; i++)
	{
		for(size_t j = 0; j < B_P_NRow; j++)
		{
			if(distance_mat2mat(A, A_N_Rows[i], B, B_P_Rows[j]) < INTERACTION_DISTANCE)
			{
				//printf("Interaction between Resid: %d, atom: %s and Resid: %d, atom %s\n", AData->atom_data->dnt_pos[A_N_Rows[i]] - 1 + AData->position_in_lib[0], AData->atom_data->name[A_N_Rows[i]], BData->atom_data->dnt_pos[B_P_Rows[j]] - 1 + BData->position_in_lib[0], BData->atom_data->name[B_P_Rows[j]]);
				IDX1 = IDX_FLAT2D(A_Idx, A_N_Rows[i], DIM2);
				IDX2 = IDX_FLAT2D(B_Idx, B_P_Rows[j], DIM2);
				gsl_matrix_set(AdjMatrix, PAdjMap[IDX2], NAdjMap[IDX1], 1);
				Interactions++;
			}
		}
	}
	return Interactions;
}

void RNADataArray::printall()
{
	sequence[0].print();
	for (int i = 1; i < count; i++)
	{
		sequence[i].print_offset(1, i);
	}
}

void RNADataArray::initialize_string()
{
	get_atom_sum();
	string_buffer = 55 * atom_sum + GLOBAL_STANDARD_STRING_LENGTH; //55 = length of ATOM line in PDB format + 2kb extra for REMARK lines and ENDMDLs
	string_out = (char *)malloc(sizeof(char) * string_buffer);
	string_index = 0;
	model_count = 0;
	string_initialized = true;
}

int RNADataArray::get_atom_sum()
{
	atom_sum = sequence[0].atom_data->count_per_res[0];
	atom_sum *= 2; //for Double Stranded
	for (int i = 0; i <= iterator_max; i++)
	{
		atom_sum += sequence[i].atom_data->count_per_res[1];
	}
	return atom_sum;
}

int RNADataArray::out_string_header_coord(RunInfo& run_info)
{
	char FORMAT_STRING[100];
	int jump_pos1, jump_pos2;
	string_index += snprintf(&string_out[string_index], string_buffer - string_index, "MODEL %-100d\n", ++model_count);
	string_index += snprintf(&string_out[string_index], string_buffer - string_index, "REMARK ");
	string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%s ", run_info.sequence);
	snprintf(&FORMAT_STRING[0], 100, "%s%d%s%1s", "%", (iterator_max + 1) * 4 - 1, "s", " ");
	jump_pos1 = string_index;
	string_index += snprintf(&string_out[string_index], string_buffer - string_index, FORMAT_STRING, " ");
	jump_pos2 = string_index;
	string_index = jump_pos1;
	for (int i = 0; i < count; i++)
	{
		string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%d", sequence[i].position_in_lib[1]  + 1);
		if (i != count - 1)
		{
			string_index += snprintf(&string_out[string_index], string_buffer - string_index, "-");
		}
	}
	string_out[string_index] = ' ';
	string_index = jump_pos2;
	string_index += snprintf(&string_out[string_index], string_buffer - string_index, " %6.3f ", structure_energy);
	if (n_wc_pairs != 0)
	{
		for(size_t i = 0; i < n_wc_pairs; i++) {
			string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%6.3f", wc_rmsd[i]);
		}
	}
	string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\n");
	return string_index;
}

int RNADataArray::out_string_header(RunInfo& run_info)
{
	if (run_info.run_options & RunOpts::write_coordinates)
	{
		string_index = out_string_header_coord(run_info);
	}
	else
	{
		char FORMAT_STRING[100];
		int jump_pos1, jump_pos2;
		string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%s ", run_info.sequence);
		string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#INDEX ");
		snprintf(&FORMAT_STRING[0], 100, "%s%d%s%1s", "%", (iterator_max + 1) * 4 - 1, "s", " ");
		jump_pos1 = string_index;
		string_index += snprintf(&string_out[string_index], string_buffer - string_index, FORMAT_STRING, " ");
		jump_pos2 = string_index;
		string_index = jump_pos1;
		for (int i = 0; i < count; i++)
		{
			string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%d", sequence[i].position_in_lib[1] + 1);
			if (i != count - 1)
			{
				string_index += snprintf(&string_out[string_index], string_buffer - string_index, "-");
			}
		}
		string_out[string_index] = ' ';
		string_index = jump_pos2;

		string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#ENERGY");
		string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%6.3f    ", structure_energy);
		if (n_wc_pairs != 0)
		{
			string_index += snprintf(&string_out[string_index], string_buffer - string_index, "#WC-RMSD");
			for(size_t i = 0; i < n_wc_pairs; i++) {
				string_index += snprintf(&string_out[string_index], string_buffer - string_index, "%6.3f", wc_rmsd[i]);
			}
		}
		string_index += snprintf(&string_out[string_index], string_buffer - string_index, "\n");
	}
	return string_index;
}

void RNADataArray::gen_string_prototype_ds(uint32_t& idx_offset)
{
	uint32_t pos = 0;
	for(int i = 1; i <= iterator_max; i++, pos++)
	{
		if(i == ds_data.iterator_max_Strand1 + 1)
		{
			string_index += snprintf(&string_out[string_index], string_buffer - string_index, "TER\n");
			pos++;
			string_index = sequence[i].generate_string_prototype(0, pos, string_out, string_buffer, string_index, idx_offset);
			string_index = sequence[i].generate_string_prototype(1, pos, string_out, string_buffer, string_index, idx_offset);
		}
		else
		{
			string_index = sequence[i].generate_string_prototype(1, pos, string_out, string_buffer, string_index, idx_offset);
		}
	}
	string_index = snprintf(&string_out[string_index], string_buffer - string_index, "ENDMDL\n");
}

void RNADataArray::gen_string_prototype_ss(uint32_t& idx_offset)
{
	for(int i = 1; i <= iterator_max; i++)
	{
		string_index = sequence[i].generate_string_prototype(1, i, string_out, string_buffer, string_index, idx_offset);
	}
	string_index = snprintf(&string_out[string_index], string_buffer - string_index, "ENDMDL\n");
}

void RNADataArray::generate_string_prototype(RunInfo& run_info)
{
	uint32_t idx_offset = 0;
	initialize_string();
	string_index = out_string_header(run_info);
	string_index = sequence[0].generate_string_prototype(0, 0, string_out, string_buffer, string_index, idx_offset);
	string_index = sequence[0].generate_string_prototype(1, 0, string_out, string_buffer, string_index, idx_offset);
	idx_offset = sequence[0].count;
	if(structure_type == StrType::single_strand) {
		gen_string_prototype_ss(idx_offset);
	}	else {
		gen_string_prototype_ds(idx_offset);
	}
	model_count--;
}

void RNADataArray::overwrite_string_prototype_ds(uint32_t& idx_offset)
{
	for(int i = 1; i <= iterator_max; i++)
	{   
		if(i == ds_data.iterator_max_Strand1 + 1)
		{
			string_index += 4;
			string_index = sequence[i].overwrite_string_prototype(0, string_out, string_buffer, string_index, idx_offset);
			string_index = sequence[i].overwrite_string_prototype(1, string_out, string_buffer, string_index, idx_offset);
		}
		else
		{
			string_index = sequence[i].overwrite_string_prototype(1, string_out, string_buffer, string_index, idx_offset);
		}
	}
}

void RNADataArray::overwrite_string_prototype_ss(uint32_t& idx_offset)
{
	for(int i = 1; i <= iterator_max; i++)
	{
		string_index = sequence[i].overwrite_string_prototype(1, string_out, string_buffer, string_index, idx_offset);
	}
	string_out[string_index] = 'E';
}

char* RNADataArray::overwrite_string_prototype(RunInfo& run_info)
{
	uint32_t idx_offset = sequence[0].count;
	string_index = 0;
	string_index = out_string_header(run_info);
	string_out[string_index] = 'A'; //snprintf will always attach a '\0' to the end of the string it creates. So we replace '\0' with 'A' b/c next line starts w/ "ATOM"
	string_index = sequence[0].overwrite_string_prototype(0, string_out, string_buffer, string_index, idx_offset);
	string_index = sequence[0].overwrite_string_prototype(1, string_out, string_buffer, string_index, idx_offset);
	if(structure_type == StrType::double_strand) {
		overwrite_string_prototype_ds(idx_offset);
	}	else {
		overwrite_string_prototype_ss(idx_offset);
	}
	return string_out;
}

int* RNADataArray::get_index()
{
	int *ar = (int *)malloc(sizeof(int) * count);
	for (int i = 0; i < count; i++)
	{
		ar[i] = sequence[i].position_in_lib[1];
	}
	return ar; 
}

void RNADataArray::print_index(int offset)
{
	for (int i = 0; i < count - 1 + offset; i++)
	{
		printf("%d-", sequence[i].position_in_lib[1]);
	}
	printf("%d\n", sequence[count - 1 + offset].position_in_lib[1]);
}

void RNADataArray::print_index()
{
	for (int i = 0; i < count - 1; i++)
	{
		printf("%d-", sequence[i].position_in_lib[1]);
	}
	printf("%d\n", sequence[count - 1].position_in_lib[1]);
}
