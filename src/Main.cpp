//Started by Ivan Riveros on: **03/5/21**
#include "TAPERSS.hpp"
#include "Atom_info.hpp"
#include "DimerLib.hpp"
#include "RNAData.hpp"
#include "Combinatorial_Addition.hpp"
#include "OutputString.hpp"
#include "WatsonCrickPair.hpp"
#include "HBondDetector.hpp"
#include "InputHandler.hpp"
#include "StructureBuilders.hpp"
#include "RNADataArrayInternalLoop.hpp"
#include "n_lowest_energy.hpp"
#include <unistd.h>

using namespace std;

void initialize(RNADataArray& model, DimerLibArray& frag_library, DimerLibArray& wc_library, OutputString& output, RunInfo& run_info)
{
	run_info.init_timer.start_timer();

	load_libs(run_info, frag_library, false);
	frag_library.get_charged_atom_map();
	if(run_info.run_options & RunOpts::use_structure_filter || run_info.structure_type == StrType::double_strand) {
		load_libs(run_info, wc_library, true);
		WC_create(wc_library);
	}
	model.initialize(run_info, frag_library, wc_library);
	output.initialize(run_info, model.string_out);
	kabsch_create(frag_library.LargestAtomCount, MATRIX_DIMENSION2);
	HK_create(frag_library.PositiveAtomCount, frag_library.NegativeAtomCount);
	if(run_info.run_options & RunOpts::build_limit_by_energy) {
		create_n_lowest_E(run_info.build_limit, model);
	}

	run_info.init_timer.stop_timer();
	printf("Initialization Time: ");
	run_info.init_timer.print();
}

void run(RNADataArray& model, DimerLibArray& frag_library, DimerLibArray& wc_library, OutputString& output, RunInfo& run_info)
{
	printf("Starting Run...\n");
	run_info.run_timer.start_timer();
	switch(run_info.run_type)
	{
		case RunType::combinatorial:
			run_combinatorial(model, frag_library, wc_library, output, run_info);
			break;
		case RunType::build_from_index:
			if(run_info.structure_type == StrType::double_strand) {
				build_structure_from_index_ds(model, frag_library, wc_library, output, run_info);
			} else {
				build_structure_from_index(model, frag_library, output, run_info);
			}
			break;
		case RunType::build_from_index_list:
			build_structure_from_index_list(model, frag_library, wc_library, output, run_info);
			break;
		case RunType::runtype_undef:
			break;
	}

	if(run_info.run_options & RunOpts::build_limit_by_energy) {
		size_t old_models = run_info.n_total_structs_built;
		size_t ubound = run_info.build_limit;
		run_info.n_filter_structs_built = 0;
		model.model_count = 0;
		if(ubound > old_models) {
			ubound = old_models;
		}
		for(uint i = 0; i < ubound; i++) {
			model.rollback_by(model.iterator);
			update_run_info_index(run_info);
			build_structure_from_index(model, frag_library, output, run_info);
		}
		run_info.n_total_structs_built = old_models;
	}

	printf("Done.\n");
	run_info.run_timer.stop_timer();
	printf("Calculation Time: ");
	run_info.run_timer.print();
	run_info.print_stats();
}

void destroy_run_info(RunInfo& run_info)
{
	for(uint i = 0; i < run_info.n_fragments; i++) {
		free(run_info.fragment_lib_list[i]);
	}
	free(run_info.fragment_lib_list);
	free(run_info.lib_duplicate_record);
	free(run_info.frag_lib_bounds);
	if(run_info.run_options & RunOpts::use_structure_filter || run_info.structure_type == StrType::double_strand) {
		for(uint i = 0; i < run_info.n_wc_pairs; i++) {
			free(run_info.wc_lib_list[i]);
		}
		free(run_info.wc_lib_duplicate_record);
	}
	if(run_info.run_type == RunType::build_from_index) {
		free(run_info.index);
	}
}

void destroy(RNADataArray& model, DimerLibArray& frag_library, DimerLibArray& wc_library, OutputString& output, RunInfo& run_info)
{
	model.destroy();
	frag_library.destroy();
	output.destroy();
	if(run_info.run_options & RunOpts::use_structure_filter || run_info.structure_type == StrType::double_strand) {
		wc_library.destroy();
		WC_destroy();
	}
	kabsch_destroy();
	HK_destroy();
	if(run_info.run_options & RunOpts::build_limit_by_energy) {
		destroy_n_lowest_E();
	}
	destroy_run_info(run_info);
}

int main(int argc, char *ARGV[])
{
	RunInfo run_info;
	RNADataArray model;
	DimerLibArray frag_library;
	DimerLibArray wc_library;
	OutputString  output;

	input_handler(argc, ARGV, run_info);
	initialize(model, frag_library, wc_library, output, run_info);
	run(model, frag_library, wc_library, output, run_info);
	destroy(model, frag_library, wc_library, output, run_info);
}
