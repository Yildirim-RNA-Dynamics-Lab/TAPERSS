#include "steric_clash_check.hpp"

constexpr double SQR_VDW_RADIUS = square(VDW_RADIUS);

AttachStatus scc_full_check_all(gsl_matrix *A, size_t A_start, size_t A_end, gsl_matrix *B, size_t B_start, size_t B_end)
{
	for(size_t i = A_start; i < A_end; i++) {
		for(size_t j = B_start; j < B_end; j++) {
			if(sqr_distance_mat2mat(A, i, B, j) < SQR_VDW_RADIUS) {
				return AttachStatus::FAILED;
			}
		}
	}
	return AttachStatus::ATTACHED;
}

AttachStatus scc_check_last(gsl_matrix *A, uint16_t *A_Rows, size_t A_Count, gsl_matrix *B, size_t B_start, size_t B_end)
{
	for(size_t i = 0; i < A_Count; i++) {
		for(size_t j = B_start; j < B_end; j++) {
			if(sqr_distance_mat2mat(A, A_Rows[i], B, j) < SQR_VDW_RADIUS) {
				return AttachStatus::FAILED;
			}
		}
	}
	return AttachStatus::ATTACHED;
}

void cg_scc_check_com(gsl_matrix* M_COMS, double* M_Radii, int Count, gsl_matrix* A, size_t A_index, double A_Radius, bool* PassArray)
{    
	PassArray[0] = (sqr_distance_mat2mat(M_COMS,  0, A, A_index) > square((M_Radii[0] + A_Radius)));
	PassArray[1] = (sqr_distance_mat2mat(M_COMS, 1, A, A_index) > square((M_Radii[1] + A_Radius)));

	for(int i = 1; i < Count; i++) {
		PassArray[i + 1] = (sqr_distance_mat2mat(M_COMS, i * 2 + 1 , A, A_index) > square((M_Radii[i * 2 + 1] + A_Radius)));
	}
}

AttachStatus cg_scc(RNADataArray& model, RNAData* new_frag)
{  
	AttachStatus status;
	cg_scc_check_com(model.COMS, model.Radii, model.count, new_frag->data_matrix, new_frag->get_residue_COM_index(1), new_frag->COM_Radii[1], 
			model.PassedCOMCheck);

	if(model.PassedCOMCheck[0] == false) {
		status = scc_full_check_all(model[0]->data_matrix, model[0]->ResBoundaries[0], model[0]->ResBoundaries[1],
				new_frag->data_matrix, new_frag->ResBoundaries[2], new_frag->ResBoundaries[3]);
		if(status == AttachStatus::FAILED) {return status;}
	}

	if(model.PassedCOMCheck[model.count] == false) {
		status = scc_check_last(model[model.count - 1]->data_matrix, model[model.count - 1]->StericIndices[1], 
				model[model.count - 1]->count_per_Steric[1], new_frag->data_matrix, 
				new_frag->ResBoundaries[2], new_frag->ResBoundaries[3]);
		if(status == AttachStatus::FAILED) {return status;}
	}

	for(int i = 0; i < model.count - 1; i++) {
		if(model.PassedCOMCheck[i + 1] == false) {
			status = scc_full_check_all(model[i]->data_matrix, model[i]->ResBoundaries[2], model[i]->ResBoundaries[3],
					new_frag->data_matrix, new_frag->ResBoundaries[2], new_frag->ResBoundaries[3]);
			if(status == AttachStatus::FAILED) {return status;}
		}
	}
	return AttachStatus::ATTACHED;
}

template <bool IN2NDSTRAND>void cg_scc_ds_check_com(gsl_matrix* M_COMS, double* M_Radii, int Count, gsl_matrix* A, size_t A_index, double A_Radius, bool* PassArray, int size1)
{    
	PassArray[0] = (sqr_distance_mat2mat(M_COMS,  0, A, A_index) > square((M_Radii[0] + A_Radius)));
	PassArray[1] = (sqr_distance_mat2mat(M_COMS, 1, A, A_index) > square((M_Radii[1] + A_Radius)));
	if constexpr (IN2NDSTRAND == true) {
		for(int i = 1; i < size1; i++) {
			PassArray[i + 1] = (sqr_distance_mat2mat(M_COMS, i * 2 + 1 , A, A_index) > square((M_Radii[i * 2 + 1] + A_Radius)));
		}
		PassArray[size1] = (sqr_distance_mat2mat(M_COMS, size1 * 2, A, A_index) > square((M_Radii[size1] + A_Radius)));
		PassArray[size1 + 1] = (sqr_distance_mat2mat(M_COMS, size1 * 2 + 1, A, A_index) > square((M_Radii[size1 + 1] + A_Radius)));
		for(int i = size1 + 2; i < Count; i++) {
			PassArray[i + 1] = (sqr_distance_mat2mat(M_COMS, i * 2 + 1, A, A_index) > square((M_Radii[i * 2 + 1] + A_Radius)));
		}
		return;
	}   
	if constexpr (IN2NDSTRAND == false) {
		for(int i = 1; i < Count; i++) {
			PassArray[i + 1] = (sqr_distance_mat2mat(M_COMS, i * 2 + 1 , A, A_index) > square((M_Radii[i * 2 + 1] + A_Radius)));
		}
		return;
	}
}
template void cg_scc_ds_check_com<true>(gsl_matrix* M_COMS, double* M_Radii, int Count, gsl_matrix* A, size_t A_index, double A_Radius, bool* PassArray, int size1);
template void cg_scc_ds_check_com<false>(gsl_matrix* M_COMS, double* M_Radii, int Count, gsl_matrix* A, size_t A_index, double A_Radius, bool* PassArray, int size1);

/**
 * @brief Performs a steric clash check on the newly started right side strand after proper alignment. Checks both 5' and 3' residues of the new nucleotide
 * using COM check first, then a further check on any failed COM checks.
 * @param model - RNADataArray reference.
 * @param new_frag  - RNAData (new nucleotides) to be attached
 * @return AttachStatus - FAILED if steric clash is detected, ATTACHED otherwise
 */
AttachStatus cg_scc_ds_check_new_2nd(RNADataArray& model, RNAData* new_frag)
{
	cg_scc_ds_check_com<false>(model.COMS, model.Radii, model.count, new_frag->data_matrix, new_frag->get_residue_COM_index(0), new_frag->COM_Radii[0], 
			model.PassedCOMCheck, model.ds_data.strand1_size);
	if(model.PassedCOMCheck[0] == false) {
		if(scc_full_check_all(model[0]->data_matrix, model[0]->ResBoundaries[0], model[0]->ResBoundaries[1],
					new_frag->data_matrix, new_frag->ResBoundaries[0], new_frag->ResBoundaries[1]) == false) {
			return AttachStatus::FAILED;
		}
	}
	for(int i = 0; i < model.count; i++) {
		if(model.PassedCOMCheck[i + 1] == false) {
			if(scc_full_check_all(model[i]->data_matrix, model[i]->ResBoundaries[2], model[i]->ResBoundaries[3],
						new_frag->data_matrix, new_frag->ResBoundaries[0], new_frag->ResBoundaries[1]) == false) {
				return AttachStatus::FAILED;
			}
		}
	}
	cg_scc_ds_check_com<false>(model.COMS, model.Radii, model.count, new_frag->data_matrix, new_frag->get_residue_COM_index(1), new_frag->COM_Radii[1], 
			model.PassedCOMCheck, model.ds_data.strand1_size);
	if(model.PassedCOMCheck[0] == false) {
		if(scc_full_check_all(model[0]->data_matrix, model[0]->ResBoundaries[0], model[0]->ResBoundaries[1],
					new_frag->data_matrix, new_frag->ResBoundaries[2], new_frag->ResBoundaries[3]) == false) {
			return AttachStatus::FAILED;
		}
	}
	for(int i = 0; i < model.count; i++) {
		if(model.PassedCOMCheck[i + 1] == false) {
			if(scc_full_check_all(model[i]->data_matrix, model[i]->ResBoundaries[2], model[i]->ResBoundaries[3],
						new_frag->data_matrix, new_frag->ResBoundaries[2], new_frag->ResBoundaries[3]) == false) {
				return AttachStatus::FAILED;
			}
		}
	}
	return AttachStatus::ATTACHED;
}

/**
 * @brief Performs steric clash check for addition of a new (non 5') nucleotide to internal loop. First the COM distances are calculated, 
 * then any which fail COM check have a full atom-to-atom check. When performing a steric clash check on the previous nucleotide (i.e. 4th being added
 * residue is checked with 3rd) only the phosphate groups are checked with scc_check_last(...).
 * @tparam IN2NDSTRAND - template bool indicating whether or not the nucleotide being added is on the left side or right side. If true, nt is on the right
 * @param model - RNADataArray reference.
 * @param new_frag  - RNAData (new nucleotide) to be attached
 * @return AttachStatus - FAILED if steric clash is detected, ATTACHED otherwise
 */
template <bool IN2NDSTRAND>AttachStatus cg_scc_ds(RNADataArray& model, RNAData* new_frag)
{
	cg_scc_ds_check_com<IN2NDSTRAND>(model.COMS, model.Radii, model.count, new_frag->data_matrix, new_frag->get_residue_COM_index(1), new_frag->COM_Radii[1], 
			model.PassedCOMCheck, model.ds_data.strand1_size);
	if(model.PassedCOMCheck[0] == false) {
		if(scc_full_check_all(model[0]->data_matrix, model[0]->ResBoundaries[0], model[0]->ResBoundaries[1],
					new_frag->data_matrix, new_frag->ResBoundaries[2], new_frag->ResBoundaries[3]) == false) {
			return AttachStatus::FAILED;
		}
	}
	if(model.PassedCOMCheck[model.count] == false) {
		if(scc_check_last(model[model.count - 1]->data_matrix, model[model.count - 1]->StericIndices[1], 
					model[model.count - 1]->count_per_Steric[1], new_frag->data_matrix, 
					new_frag->ResBoundaries[2], new_frag->ResBoundaries[3]) == false) {
			return AttachStatus::FAILED;
		}
	}

	if constexpr(IN2NDSTRAND == true) {
		if(model.PassedCOMCheck[model.ds_data.strand1_size] == false) {
			if(scc_full_check_all(model[model.ds_data.strand1_size]->data_matrix, model[model.ds_data.strand1_size]->ResBoundaries[0], model[model.ds_data.strand1_size]->ResBoundaries[1],
						new_frag->data_matrix, new_frag->ResBoundaries[2], new_frag->ResBoundaries[3]) == false) {
				return AttachStatus::FAILED;
			}
		}
		for(int i = 0; i < model.count - 1; i++) {
			if(model.PassedCOMCheck[i + 1] == false) {
				if(scc_full_check_all(model[i]->data_matrix, model[i]->ResBoundaries[2], model[i]->ResBoundaries[3],
							new_frag->data_matrix, new_frag->ResBoundaries[2], new_frag->ResBoundaries[3]) == false) {
					return AttachStatus::FAILED;
				}
			}
		}
	} else if constexpr(IN2NDSTRAND == false) {
		for(int i = 0; i < model.count - 1; i++) {
			if(model.PassedCOMCheck[i + 1] == false) {
				if(scc_full_check_all(model[i]->data_matrix, model[i]->ResBoundaries[2], model[i]->ResBoundaries[3],
							new_frag->data_matrix, new_frag->ResBoundaries[2], new_frag->ResBoundaries[3]) == false) {
					return AttachStatus::FAILED;
				}
			}
		}
	}
	return AttachStatus::ATTACHED;
}
template AttachStatus cg_scc_ds<true>(RNADataArray& model, RNAData* new_frag);
template AttachStatus cg_scc_ds<false>(RNADataArray& model, RNAData* new_frag);
