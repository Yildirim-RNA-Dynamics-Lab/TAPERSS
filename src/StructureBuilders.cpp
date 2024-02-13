#include "StructureBuilders.hpp"

void build_structure_from_index(RNADataArray &assembled, DimerLibArray& frag_lib, OutputString &o_string, RunInfo& run_info)
{
	assembled.overwrite_initialize(0, run_info.index[0], frag_lib);
	for (int i = 1; i <= assembled.iterator_max; i++) {
		assembled.overwrite(i, run_info.index[i], frag_lib);
		attach_5p_to_3p(assembled.current(), assembled[i], run_info);
		assembled.keep();
	}
	if (run_info.run_options & RunOpts::use_structure_filter) {
		WC_full_check<false>(assembled, 0, run_info); 
	}
	assembled.update_energy();
	o_string.add_string(assembled.overwrite_string_prototype(run_info));
	run_info.n_total_structs_built++;
	run_info.n_filter_structs_built++;
}

void build_structure_from_index_ds(RNADataArray &assembled, DimerLibArray &frag_lib, DimerLibArray &wc_lib, OutputString &o_string, RunInfo& run_info)
{
	uint start = run_info.run_options & RunOpts::str_filter_uses_ds_closing_bp ? 1 : 0;
	assembled.overwrite_initialize(0, run_info.index[0], frag_lib);
	for (int i = 1; i <= assembled.iterator_max; i++) {
		assembled.overwrite(i, run_info.index[i], frag_lib);
		if (assembled.should_prepare_right(i) == true) {
			if(run_info.run_options & RunOpts::str_filter_uses_ds_closing_bp) {
				prepare_right<RunOpts::str_filter_uses_ds_closing_bp>(assembled[i], wc_lib, assembled, run_info);
			} else {
				prepare_right<0>(assembled[i], wc_lib, assembled, run_info);
			}
		} else {
			attach_5p_to_3p(assembled.current(), assembled[i], run_info);
		}
		assembled.keep();
	}
	if(run_info.run_options & RunOpts::use_structure_filter) { WC_full_check<false>(assembled, start, run_info); }
	o_string.add_string(assembled.overwrite_string_prototype(run_info));
	run_info.n_total_structs_built++;
	run_info.n_filter_structs_built++;
}

void build_structure_from_index_list(RNADataArray &assembled, DimerLibArray &frag_lib, DimerLibArray& wc_lib, OutputString &o_string, RunInfo& run_info)
{
	FILE* index_list = fopen(run_info.index_list_file, "r");
	if(index_list == NULL) {
		fprintf(stderr, "Error: Could not open index list file: %s\n", run_info.index_list_file);
		exit(2);
	}
	while(get_next_index_from_file(index_list, run_info)) {
		if(run_info.structure_type == StrType::double_strand) {
			build_structure_from_index_ds(assembled, frag_lib, wc_lib, o_string, run_info);
		} else {
			build_structure_from_index(assembled, frag_lib, o_string, run_info);
		}
		assembled.rollback_by(assembled.iterator);
	}
	fclose(index_list);
}
