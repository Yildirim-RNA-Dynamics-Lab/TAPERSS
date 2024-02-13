#include "Combinatorial_Addition.hpp"
bool early_exit_handler(RunInfo& run_info)
{
	if(run_info.run_options & RunOpts::use_structure_filter) {
		return run_info.n_filter_structs_built == run_info.build_limit;
	} else {
		return run_info.n_total_structs_built == run_info.build_limit;
	}
}

void output_handler(RNADataArray &assembled, OutputString &o_string, RunInfo& run_info)
{
	run_info.n_total_structs_built++;
	if(run_info.run_options & RunOpts::use_structure_filter) {
		uint start = run_info.run_options & RunOpts::str_filter_uses_ds_closing_bp ? 1 : 0;
		if(WC_full_check<true>(assembled, start, run_info) != true) { return; }
	}
	assembled.update_energy();
	if(run_info.run_options & RunOpts::build_limit_by_energy) {
		int idx;
		if((idx = add_to_n_lowest_E(assembled.structure_energy, assembled)) != -1) {
			assembled.model_count = idx + 1;
			//o_string.insert_string(assembled.overwrite_string_prototype(run_info), idx);
		}
		
	} else {
		o_string.add_string(assembled.overwrite_string_prototype(run_info));
	}
	run_info.n_filter_structs_built++;
	return;
}

void fill_if_empty(CMB_Manager &manager, RNADataArray &assembled, DimerLibArray &Lib)
{
	uint32_t dnmp_pos = assembled.iterator + 1;
	uint32_t lib_pos = manager.attach_attempted[dnmp_pos] + 1;
	//DEBUG_PRINT("Filling Empty RNAData with index (%d:%d)\n", dnmp_pos, lib_pos);
	assembled.overwrite(dnmp_pos, lib_pos, Lib);
	assembled.keep();
	manager.attach_attempt(dnmp_pos, lib_pos);
}

template <AttachStatus status> bool rollback_handler(CMB_Manager &manager, RNADataArray &assembled, DimerLibArray &Lib)
{
	constexpr uint roll_back_offset = status == AttachStatus::ATTACHED ? 1 : 0;//+1 b/c we called keep(),we need to rollback an extra fragment to move on.
	if (manager.is_at_end()) {
		if (manager.check_lib_completion()) {
			return true;
		} else {
			assembled.rollback_by(manager.rollback_count + roll_back_offset); 
			if(assembled.is_empty()) { fill_if_empty(manager, assembled, Lib); }
		}
	} else {
		assembled.rollback();
	}  
	return false;
}
template bool rollback_handler<AttachStatus::ATTACHED>(CMB_Manager &manager, RNADataArray &assembled,  DimerLibArray &Lib);
template bool rollback_handler<AttachStatus::FAILED>(CMB_Manager &manager, RNADataArray &assembled,  DimerLibArray &Lib);

template <uint32_t OPTS> void run_combinatorial_internal(RNADataArray &assembled, DimerLibArray& frag_lib, 
		DimerLibArray& wc_lib, OutputString &o_string, RunInfo& run_info)
{
	CMB_Manager manager(frag_lib, run_info);
	bool run_done = false;
	while(run_done == false) {
		switch(fragment_assembly<OPTS>(frag_lib, wc_lib, assembled, manager, run_info)) {
			case AttachStatus::ATTACHED:
				assembled.keep();
				break;
			case AttachStatus::FAILED:
				run_done = rollback_handler<AttachStatus::FAILED>(manager, assembled, frag_lib);
				continue;
		}
		if (assembled.is_complete()) {
			output_handler(assembled, o_string, run_info);
			if (run_info.run_options & RunOpts::blind_build_limit){ if(early_exit_handler(run_info)) break; }
			run_done = rollback_handler<AttachStatus::ATTACHED>(manager, assembled, frag_lib);
		}
	}
}

void run_combinatorial(RNADataArray &assembled, DimerLibArray& frag_lib, DimerLibArray& wc_lib, OutputString &o_string, RunInfo& run_info)
{
	CMB_Manager manager(frag_lib, run_info);
	if(run_info.run_options & RunOpts::strtype_ds) {
		if(run_info.run_options & RunOpts::str_filter_uses_ds_closing_bp) {
			run_combinatorial_internal<RunOpts::strtype_ds | RunOpts::str_filter_uses_ds_closing_bp>(assembled, frag_lib, wc_lib, o_string, run_info);
		} else {
			run_combinatorial_internal<RunOpts::strtype_ds>(assembled, frag_lib, wc_lib, o_string, run_info);
		}
	} else {
		run_combinatorial_internal<0>(assembled, frag_lib, wc_lib, o_string, run_info);
	}
}
