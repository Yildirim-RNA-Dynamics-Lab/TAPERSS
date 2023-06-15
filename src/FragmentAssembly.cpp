#include "FragmentAssembly.hpp"

AttachStatus rotate(RNAData *reference, RNAData *rotated)
{
	double rmsd;
	double COMP[] = {0, 0, 0};
	double COMQ[] = {0, 0, 0};
	AttachStatus status = FAILED;

	gsl_matrix *MODEL = rotated->data_matrix;
	gsl_matrix *P = kabsch_prepare_matrix<KABSCH_MATRIX_P>(rotated->count_per_sub[0], MATRIX_DIMENSION2, rotated->submatrix_rows[0], rotated->data_matrix);
	gsl_matrix *Q = kabsch_prepare_matrix<KABSCH_MATRIX_Q>(reference->count_per_sub[1], MATRIX_DIMENSION2, reference->submatrix_rows[1], reference->data_matrix);
	gsl_matrix *P_WORK = kabsch_get_work_matrix(MATRIX_DIMENSION2, rotated->count_per_sub[0]);
	gsl_matrix *R;

	kabsch_calculate_rotation_matrix_Nx3fast(P, Q, P_WORK, COMP, COMQ);
	R = kabsch_get_rotation_matrix();
	rmsd = rmsd_generic(P, Q);

	if (rmsd <= GLOBAL_RMSD_LIMIT)
	{
		status = ATTACHED;
		translate_matrix(COMP, MODEL, -1.0F);
		apply_rotation_matrix(R, MODEL);
		translate_matrix(COMQ, MODEL, 1.0F);
	}
	return status;
}

template <uint32_t OPTS> AttachStatus prepare_right(RNAData* to_be_assembled, DimerLibArray &WC_Lib, RNADataArray &assembled) 
{
	RNAData *assembled_ref = assembled[assembled.ds_data.iterator_max_Strand1];
	gsl_matrix *R1, *R2;
	double COMP[] = {0, 0, 0}, COMQ[] = {0, 0, 0};

	gsl_matrix *WC_MODEL = WC_Lib[1]->data_matrices[0];
	gsl_matrix *MODEL = to_be_assembled->data_matrix;

	gsl_matrix *P1 = kabsch_prepare_matrix<KABSCH_MATRIX_P>(assembled.ds_data.wc_rowsize[0], MATRIX_DIMENSION2, assembled.ds_data.wc_rows[0], WC_MODEL);
	gsl_matrix *Q1 = kabsch_prepare_matrix<KABSCH_MATRIX_Q>(assembled_ref->count_per_sub[1], MATRIX_DIMENSION2, assembled_ref->submatrix_rows[1], assembled_ref->data_matrix);
	gsl_matrix *P1WORK = kabsch_get_work_matrix(MATRIX_DIMENSION2, assembled.ds_data.wc_rowsize[0]);

	kabsch_calculate_rotation_matrix_Nx3fast(P1, Q1, P1WORK, COMP, COMQ);
	R1 = kabsch_get_rotation_matrix();

	gsl_matrix *P2 = kabsch_prepare_matrix<KABSCH_MATRIX_P>(to_be_assembled->count_per_sub[0], MATRIX_DIMENSION2, to_be_assembled->submatrix_rows[0], MODEL);
	gsl_matrix *Q2 = kabsch_prepare_matrix<KABSCH_MATRIX_Q>(assembled.ds_data.wc_rowsize[1], MATRIX_DIMENSION2, assembled.ds_data.wc_rows[1], WC_MODEL);

	translate_matrix(COMP, Q2, -1.0F);
	apply_rotation_matrix(R1, Q2);
	translate_matrix(COMQ, Q2, 1.0F); 

	gsl_matrix *P2WORK = kabsch_get_work_matrix(MATRIX_DIMENSION2, to_be_assembled->count_per_sub[0]);

	COMP[0] = COMP[1] = COMP[2] = 0;
	COMQ[0] = COMQ[1] = COMQ[2] = 0;

	kabsch_calculate_rotation_matrix_Nx3fast(P2, Q2, P2WORK, COMP, COMQ);

	R2 = kabsch_get_rotation_matrix();

	translate_matrix(COMP, MODEL, -1.0F);
	apply_rotation_matrix(R2, MODEL);
	translate_matrix(COMQ, MODEL, 1.0F);

	COMP[0] = COMP[1] = COMP[2] = 0;
	COMQ[0] = COMQ[1] = COMQ[2] = 0;
	if(OPTS & RunOpts::str_filter_uses_ds_closing_bp) {
		WC_prepare_structure_matrix(0, assembled_ref, 1, to_be_assembled, 0);
		double RMSD = WC_check_pair(0);
		if(RMSD <= GLOBAL_WC_RMSD_LIMIT)
		{
			assembled.update_WC_rmsd(RMSD, 0);
			return cg_scc_ds_check_new_2nd(assembled, to_be_assembled);
		}
		return AttachStatus::FAILED;
	} else {
		return cg_scc_ds_check_new_2nd(assembled, to_be_assembled);
	}
}

AttachStatus fragment_assembly_ss(DimerLibArray &Lib, RNADataArray &assembled, CMB_Manager &manager)
{
	uint32_t dnmp_position = assembled.iterator + 1; // Position in sequence where new DNMP will be attached.
	uint32_t library_position = manager.attach_attempted[dnmp_position] + 1; //Continue from last attempted frag in lib + 1

	for (int i = library_position; i < manager.lib_bounds[dnmp_position].idx2; i++) {
		assembled.overwrite(dnmp_position, i, Lib);
		manager.attach_attempt(dnmp_position, i);
		if (rotate(assembled.current(), assembled[dnmp_position]) != AttachStatus::ATTACHED) {
			continue;
		}
		if (cg_scc(assembled, assembled[dnmp_position]) == AttachStatus::ATTACHED) {
			return AttachStatus::ATTACHED;
		}
	}
	return AttachStatus::FAILED;
}

template <uint32_t OPTS> AttachStatus fragment_assembly_ds(DimerLibArray &Lib, DimerLibArray &WC_Lib, RNADataArray &assembled, CMB_Manager &manager)
{
	uint32_t dnmp_position = assembled.iterator + 1; // Position in sequence where new DNMP will be attached.
	uint32_t library_position = manager.attach_attempted[dnmp_position] + 1;

	if (assembled.should_prepare_right(dnmp_position) == true) {
		for (int i = library_position; i < manager.lib_bounds[dnmp_position].idx2; i++) {
			assembled.overwrite(dnmp_position, i, Lib);
			manager.attach_attempt(dnmp_position, i);
			if (prepare_right<OPTS>(assembled[dnmp_position], WC_Lib, assembled) == AttachStatus::ATTACHED) {
				return AttachStatus::ATTACHED;
			}
			continue;
		}
		return AttachStatus::FAILED;
	}

	if(assembled.count > assembled.ds_data.strand1_size) {
		for (int i = library_position; i < manager.lib_bounds[dnmp_position].idx2; i++)
		{
			assembled.overwrite(dnmp_position, i, Lib);
			manager.attach_attempt(dnmp_position, i);
			if (rotate(assembled.current(), assembled[dnmp_position]) != AttachStatus::ATTACHED) {
				continue;
			}
			if(cg_scc_ds<true>(assembled, assembled[dnmp_position]) == AttachStatus::ATTACHED) {
				return AttachStatus::ATTACHED;
			}
			continue;
		}
		return AttachStatus::FAILED;
	} else {
		for (int i = library_position; i < manager.lib_bounds[dnmp_position].idx2; i++) {
			assembled.overwrite(dnmp_position, i, Lib);
			manager.attach_attempt(dnmp_position, i);
			if (rotate(assembled.current(), assembled[dnmp_position]) != AttachStatus::ATTACHED) {
				continue;
			}
			if(cg_scc_ds<false>(assembled, assembled[dnmp_position]) == AttachStatus::ATTACHED) {
				return AttachStatus::ATTACHED;
			}
			continue;
		}
		return AttachStatus::FAILED;
	}
}

template <uint32_t OPTS> AttachStatus fragment_assembly(DimerLibArray &lib, DimerLibArray &wc_lib, RNADataArray &assembled, CMB_Manager &manager)
{
	if(OPTS & RunOpts::strtype_ds) {
		if(OPTS & RunOpts::str_filter_uses_ds_closing_bp) {
			return fragment_assembly_ds<RunOpts::str_filter_uses_ds_closing_bp>(lib, wc_lib, assembled, manager);
		} else {
			return fragment_assembly_ds<0>(lib, wc_lib, assembled, manager);
		}
	} else {
		return fragment_assembly_ss(lib, assembled, manager);
	}
}
template AttachStatus fragment_assembly<RunOpts::strtype_ds | RunOpts::str_filter_uses_ds_closing_bp>(DimerLibArray &lib, DimerLibArray &wc_lib, RNADataArray &assembled, CMB_Manager &manager);
template AttachStatus fragment_assembly<RunOpts::strtype_ds>(DimerLibArray &lib, DimerLibArray &wc_lib, RNADataArray &assembled, CMB_Manager &manager);
template AttachStatus fragment_assembly<0>(DimerLibArray &lib, DimerLibArray &wc_lib, RNADataArray &assembled, CMB_Manager &manager);
