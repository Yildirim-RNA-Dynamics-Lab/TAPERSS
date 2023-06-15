#include "InputHandler.hpp"
#include <cstdint>

void trim_whitespace(char* s)
{
	const char whitespace[] = {' ', '\t', '\n'};
	size_t str_len = strlen(s);
	size_t match_loc = strcspn(s, whitespace);
	unsigned int idx = 0;

	if(match_loc == 0) {
		while(isspace(s[idx])) {
			idx++;
		}
		memmove(s, &s[idx], str_len + 1 - idx);
		str_len -= idx;
	}
	idx = str_len - 1;	
	while(isspace(s[idx])) {
		idx--;
	}
	s[idx + 1] = '\0';
}

void get_DNMP_lib_names(RunInfo& run_info)
{
	int name_position = strcspn(run_info.library_prototype, "X");
	char* seq = run_info.sequence;
	size_t lib_proto_len = strlen(run_info.library_prototype) + 1;
	run_info.fragment_lib_list = (char **)malloc(sizeof(char *) * run_info.n_fragments);
	if(run_info.structure_type == StrType::single_strand) {
		for (size_t i = 0; i < run_info.n_fragments; i++)
		{
			run_info.fragment_lib_list[i] = (char *)malloc(lib_proto_len);
			memcpy(run_info.fragment_lib_list[i], run_info.library_prototype, lib_proto_len);
			memcpy(&run_info.fragment_lib_list[i][name_position], seq++, sizeof(char) * 2);
		}
	}
	if(run_info.structure_type == StrType::double_strand) {
		for (size_t i = 0; i < run_info.ds_strand1_n_frags; i++)
		{
			run_info.fragment_lib_list[i] = (char *)malloc(lib_proto_len);
			memcpy(run_info.fragment_lib_list[i], run_info.library_prototype, lib_proto_len);
			memcpy(&run_info.fragment_lib_list[i][name_position], seq++, sizeof(char) * 2);
		}
		seq = &seq[2];
		for(size_t i = run_info.ds_strand1_n_frags; i < run_info.n_fragments; i++)
		{
			run_info.fragment_lib_list[i] = (char *)malloc(lib_proto_len);
			memcpy(run_info.fragment_lib_list[i], run_info.library_prototype, lib_proto_len);
			memcpy(&run_info.fragment_lib_list[i][name_position], seq++, sizeof(char) * 2);
		}
	}
}

uint32_t find_duplicates(RunInfo& run_info, bool for_WC)
{
	uint32_t num_duplicates = 0;
	uint32_t N = for_WC ? run_info.n_wc_pairs : run_info.n_fragments;
	char**   lib = for_WC ? run_info.wc_lib_list : run_info.fragment_lib_list;
	size_t*	 dup; 
	printf("N = %u\n", N);
	dup = (size_t *)calloc(N, sizeof(size_t));
	for (uint32_t i = N - 1; i > 0; i--)
	{
		dup[i] = i;
		for(uint32_t j = 0; j < i; j++)
		{
			if(!strcmp(lib[i], lib[j])) {
				dup[i] = j;
				num_duplicates++;
				break;
			}
		}
	}
	if(for_WC) { 
		run_info.wc_lib_duplicate_record = dup;
	} else {
		run_info.lib_duplicate_record = dup;
	}
	return num_duplicates;
}

uint32_t get_WC_partner_single(RunInfo& run_info, size_t idx1, size_t idx2, size_t lib_pos, size_t x_pos)
{
	char WC_pair[3];
	WC_pair[0] = run_info.sequence[idx1];
	WC_pair[1] = run_info.sequence[idx2];
	WC_pair[2] = '\0';

	if (WC_pair[0] == 'A') {
		if (WC_pair[1] != 'U') {
			printf("Error: Non-Watson Crick pair attempted: %s\n", WC_pair);
			return 1;
		}
	}
	else if (WC_pair[0] == 'U') {
		if (WC_pair[1] != 'A' && WC_pair[1] != 'G') {
			printf("Error: Non-Watson Crick pair attempted: %s\n", WC_pair);
			return 1;
		}
	}
	else if (WC_pair[0] == 'C') {
		if (WC_pair[1] != 'G') {
			printf("Error: Non-Watson Crick pair attempted: %s\n", WC_pair);
			return 1;
		}
	}
	else if (WC_pair[0] == 'G') {
		if (WC_pair[1] != 'C' && WC_pair[1] != 'U') {
			printf("Error: Non-Watson Crick pair attempted: %s\n", WC_pair);
			return 1;
		}
	}

	run_info.wc_lib_list[lib_pos] = (char *)malloc(strlen(run_info.wc_library_prototype) + 1);
	memcpy(run_info.wc_lib_list[lib_pos], run_info.wc_library_prototype, strlen(run_info.wc_library_prototype) + 1);
	memcpy(&run_info.wc_lib_list[lib_pos][x_pos], WC_pair, sizeof(char) * 2);
	return 0;
}

uint32_t get_WC_partner_full(RunInfo& run_info, size_t len, uint32_t* pair_record)
{
	uint32_t err_count = 0;
	uint32_t name_position = strcspn(run_info.wc_library_prototype, "X");
	uint32_t pairs_made = 0;
	uint32_t bypass_idx = run_info.n_fragments + 1;
	if(run_info.structure_type == StrType::double_strand)	{
		uint32_t x_loc = run_info.ds_strand1_n_frags + 1;
		bypass_idx = x_loc - 1;
		err_count += get_WC_partner_single(run_info, x_loc - 1, x_loc + 1, 0, name_position);
		if(pair_record[x_loc - 1] != x_loc - 1) {
			run_info.wc_pair_list[pairs_made] = {x_loc - 1, x_loc + 1};
			pairs_made++;
			run_info.run_options |= RunOpts::str_filter_uses_ds_closing_bp;
		}
	}
	for(uint32_t i = 0; i < len; i++) {
		if(pairs_made == run_info.n_wc_pairs) {break;}
		if(pair_record[i] != i && i != bypass_idx) {
			err_count += get_WC_partner_single(run_info, i, pair_record[i], pairs_made, name_position);
			run_info.wc_pair_list[pairs_made] = {i, pair_record[i]};
			pairs_made++;
		}
	}
	return err_count;
}

uint32_t get_index_int(char *index, uint32_t *indices)
{
	char buffer[5];
	uint32_t count = 0, buf_c = 0;
	for (char c = *index; c != '\0'; c = *(++index)) {
		if (c == '-') {
			buffer[buf_c] = '\0';
			indices[count] = atoi(buffer);
			buf_c = 0;
			count++;
			continue;
		}
		buffer[buf_c++] = *index;
	}
	buffer[buf_c] = '\0';
	indices[count] = atoi(buffer);
	return count;
}

template <typename T> uint32_t get_index_int_pair(char *index, IndexPair<T> *indices, uint8_t idx)
{
	char buffer[5];
	uint32_t count = 0, buf_c = 0;
	for (char c = *index; c != '\0'; c = *(++index)) {
		if (c == '-') {
			buffer[buf_c] = '\0';
			if(idx == 1)
				indices[count].idx1 = atoi(buffer);
			else
				indices[count].idx2 = atoi(buffer);
			buf_c = 0;
			count++;
			continue;
		}
		buffer[buf_c++] = *index;
	}
	buffer[buf_c] = '\0';
	if(idx == 1)
		indices[count].idx1 = atoi(buffer);
	else
		indices[count].idx2 = atoi(buffer);
	return count;
}

void parse_parallel_input(RunInfo& run_info) 
{
	run_info.frag_lib_bounds = (IndexPair<size_t>*)malloc(sizeof(IndexPair<size_t>) * run_info.n_fragments);
	if(run_info.parallel_lib_len[0] == '\0' && run_info.parallel_lib_idx[0] == '\0') {
		for(uint i = 0; i < run_info.n_fragments; i++) {
			run_info.frag_lib_bounds[i].idx1 = 0;
			run_info.frag_lib_bounds[i].idx2 = 0;
		}
	} else if (run_info.parallel_lib_len[0] == '\0' || run_info.parallel_lib_idx[0] == '\0') {
		fprintf(stderr, "Parallel run setup incomplete: missing either library segments lengths or indices!\n");
		exit(3);
	} else {
		uint32_t idx_len = get_index_int_pair(run_info.parallel_lib_idx, run_info.frag_lib_bounds, 1);
		uint32_t lib_len = get_index_int_pair(run_info.parallel_lib_len, run_info.frag_lib_bounds, 2);
		if(idx_len != lib_len || idx_len != run_info.n_fragments) {
			fprintf(stderr, "Parallel run setup error: incorrect number of library segments lengths or indices!\n");
			exit(3);
		}
	}
}

uint8_t parse_dot_backet(RunInfo& run_info, size_t len, uint32_t* pair_record)
{
	uint32_t open_brkt = 0;
	uint32_t clse_brkt = 0;

	for(size_t i = 0; i < len; i++) {
		if(run_info.dot_bracket[i] == '(') {
			open_brkt++;
		} else if (run_info.dot_bracket[i] == ')') {
			clse_brkt++;
		}
	}

	if(clse_brkt != open_brkt) {
		fprintf(stderr, "Error: Unpaired \"Bracket\" nucleotide found in input dot-bracket structure!\n");
		return 0;
	}

	for(size_t i = 0; i < len; i++) { pair_record[i] = i; }
	for(size_t i = 0; i < len; i++) {
		if(run_info.dot_bracket[i] == '(') {
			open_brkt = 0;
			clse_brkt = 0;
			for(size_t j = i + 1; j < len; j++) {
				if(run_info.dot_bracket[j] == '(') {
					open_brkt++;
				} if(run_info.dot_bracket[j] == ')') {
					if(open_brkt != clse_brkt) { clse_brkt++; }
					else { 
						run_info.n_wc_pairs++; 
						pair_record[i] = j;
					}
				}
			}
		}
	}
	return 1;
}

bool get_next_index_from_file(FILE* file, RunInfo& run_info)
{
	char line[GLOBAL_STANDARD_STRING_LENGTH];
	char indices[GLOBAL_STANDARD_STRING_LENGTH];
	if(fgets(line, sizeof(line), file) != NULL) {
		trim_whitespace(line);
		strcpy(indices, line);
	} else {
		return false;
	}
	get_index_int(indices, run_info.index);
	return true;
}

void read_input_file(char *file_name, RunInfo& run_info, char* idx_tmp)
{
	FILE *input = fopen(file_name, "r");
	if(input == nullptr) {
		printf("Could not open file: %s\n", file_name);
		exit(3);
	}
	GLOBAL_OUTPUT_FILE[0] = '\0';
	char line[GLOBAL_STANDARD_STRING_LENGTH];
	while (fgets(line, sizeof(line), input)) {
		char *header = strtok(line, "=");
		char *str1 = strtok(NULL, "=");
		trim_whitespace(header);
		trim_whitespace(str1);
		if (!strcasecmp(header, "sequence")) {
			strcpy(run_info.sequence, str1);
			continue;
		} else if (!strcasecmp(header, "OUTPUT-FILE")) {
			strcpy(run_info.output_file, str1);
			continue;
		} else if (!strcasecmp(header, "RUN-TYPE")) {
			if (!strcasecmp(str1, "cmb") || !strcasecmp(str1, "combinatorial")) {
				run_info.run_type = RunType::combinatorial;
			} else if (!strcasecmp(str1, "idx") || !strcasecmp(str1,"from-index")) {
				run_info.run_type = RunType::build_from_index;
			} else if (!strcasecmp(str1, "idx-list") || !strcasecmp(str1, "from-index-list")) {
				run_info.run_type = RunType::build_from_index_list;
			} else {
				printf("Unrecognized option: %s = %s\n", header, str1);
			}
			continue;
		} else if (!strcasecmp(header, "INDEX-LIST-FILE")) {
			strcpy(run_info.index_list_file, str1);
			continue;
		} else if (!strcasecmp(header, "structure-count-limit")) {
			run_info.build_limit = atoi(str1);
		} else if (!strcasecmp(header, "structure-count-limit-type")) {
			if (!strcasecmp(str1, "energy")) {
				run_info.run_options |= RunOpts::build_limit_by_energy;
			} else if (!strcasecmp(str1, "blind")) {
				run_info.run_options |= RunOpts::blind_build_limit;
			} else {
				printf("Unrecognized option: %s = %s\n", header, str1);
			}
			continue;
		} else if (!strcasecmp(header, "single-index-set")) {
			strcpy(idx_tmp, str1);
			continue;
		} else if (!strcasecmp(header, "overlap-rmsd-limit")) {
			run_info.rmsd_limit = atof(str1);
			continue;
		} else if (!strcasecmp(header, "watson-crick-rmsd-limit")) {
			run_info.wc_rmsd_limit = atof(str1);
			continue;
		} else if (!strcasecmp(header, "library-prototype")) {
			strcpy(run_info.library_prototype, str1);
			continue;
		} else if (!strcasecmp(header, "watson-crick-library-prototype")) {
			strcpy(run_info.wc_library_prototype, str1);
			continue;
		} else if (!strcasecmp(header, "WRITE-COORDINATES")) {
			if (!strcasecmp(str1, "TRUE")) {
				run_info.run_options |= RunOpts::write_coordinates;
			}
		} else if (!strcasecmp(header, "secondary-structure-filter")) {
			strcpy(run_info.dot_bracket, str1);
			run_info.run_options |= RunOpts::use_structure_filter;
			continue;
		} else {
			printf("Warning: Ignored unknown option in input file: %s\n", header);
		}
	}
	fclose(input);
}

uint32_t complete_run_info(RunInfo& run_info, char* idx_tmp, uint32_t err_count)
{
	if(strchr(run_info.sequence, 'x') != NULL) {
		uint32_t name_position = strcspn(run_info.wc_library_prototype, "X");
		uint32_t x_loc = strcspn(run_info.sequence, "x");
		run_info.structure_type = StrType::double_strand;	
		run_info.n_fragments = strlen(run_info.sequence) - 3;
		run_info.ds_strand1_n_frags = x_loc - 1;
		run_info.ds_strand2_n_frags = run_info.n_fragments - run_info.ds_strand1_n_frags;
		if(!(run_info.run_options & RunOpts::use_structure_filter)) {
			run_info.wc_lib_list = (char**)malloc(sizeof(char *));
			err_count += get_WC_partner_single(run_info, x_loc - 1, x_loc + 1, 0, name_position);
			run_info.n_wc_pairs = 1;
		}
	} else {
		run_info.structure_type = StrType::single_strand;
		run_info.n_fragments = strlen(run_info.sequence) - 1;
	}

	get_DNMP_lib_names(run_info);
	find_duplicates(run_info, false);
	parse_parallel_input(run_info);
	if(run_info.run_options & RunOpts::use_structure_filter) {
		size_t len_db = strlen(run_info.dot_bracket);
		if(len_db != strlen(run_info.sequence)) {
			fprintf(stderr, "Error: Number of nucleotides in secondary structure does not match sequence!\n");
			err_count++;
			return err_count;
		}
		uint32_t* tmp_pair_record = (uint32_t*)malloc(sizeof(uint32_t) * len_db);
		if(!parse_dot_backet(run_info, len_db, tmp_pair_record)) {err_count++; return err_count;}
		run_info.wc_lib_list = (char **)malloc(sizeof(char *) * run_info.n_wc_pairs);
		run_info.wc_pair_list = (IndexPair<size_t>*)malloc(sizeof(IndexPair<size_t>) * run_info.n_wc_pairs);
		err_count += get_WC_partner_full(run_info, len_db, tmp_pair_record);
		free(tmp_pair_record);
		find_duplicates(run_info, true);
	}

	if(run_info.run_type == RunType::build_from_index || run_info.run_type == RunType::build_from_index_list) {
		run_info.index = (uint32_t*)malloc(sizeof(uint32_t) * (strlen(run_info.sequence) - 1));
	}
	if(run_info.run_type == RunType::build_from_index) {
		get_index_int(idx_tmp, run_info.index);
	}
	return err_count;
}

void validate_run_info(RunInfo& run_info, char* tmp_idx)
{
	uint32_t err_count = 0;
	uint32_t warn_count = 0;
	if(run_info.sequence[0] == '\0') {
		fprintf(stderr, "Error: No sequence given in input file or as argument!\n");
		err_count++;
	}
	if(run_info.library_prototype[0] == '\0') {
		fprintf(stderr, "Error: No library file name prototype given in input file or as argument!\n");
		err_count++;
	}
	if(run_info.run_type == RunType::runtype_undef) {
		fprintf(stderr, "Error: Run type is not specified!\n");
		err_count++;
	}
	if((run_info.run_options & RunOpts::use_structure_filter) || (strchr(run_info.sequence, 'x') != NULL)) {
		if(run_info.wc_library_prototype[0] == '\0') {
			fprintf(stderr, "Error: No watson crick library file name prototype given in input file or as argument!\n");
			err_count++;
		}
	}
	if(run_info.run_type == RunType::build_from_index && tmp_idx[0] == '\0') {
		fprintf(stderr, "Error: Run is of type: \"Build from index\" but no index set is given!\n");
		err_count++;
	}
	if(run_info.run_type == RunType::build_from_index_list && run_info.index_list_file[0] == '\0') {
		fprintf(stderr, "Error: Run is of type: \"Build from index list\" but no index list file is given!\n");
		err_count++;
	}

	if(run_info.rmsd_limit == 0) {
		printf("Warning: Overlap RMSD Limit not given or set to 0. Using default value instead...\n");
		run_info.rmsd_limit = DEFAULT_RMSD_LIMIT;
		warn_count++;
	}
	if(run_info.wc_rmsd_limit == 0 && (run_info.run_options & RunOpts::use_structure_filter)) {
		printf("Warning: Watson-Crick RMSD Limit not given or set to 0. Using default value instead...\n");
		run_info.wc_rmsd_limit = DEFAULT_WC_RMSD_LIMIT;
		warn_count++;
	}
	if(run_info.build_limit != 0) {
		if(!((run_info.run_options & RunOpts::blind_build_limit) | (run_info.run_options & RunOpts::build_limit_by_energy)))
		{
			printf("Warning: Structure Build Limit specified, but limit type not given. Assuming blind limit...\n");
			run_info.run_options |= RunOpts::blind_build_limit;
			warn_count++;
		}
	}
	if((run_info.run_options & RunOpts::blind_build_limit) || (run_info.run_options & RunOpts::build_limit_by_energy)) {
		if(run_info.build_limit == 0) {
			fprintf(stderr, "Error: Build limit type specified but build limit is not set or set to 0!\n");
			err_count++;
		}
	}
	if(run_info.output_file[0] == '\0') {
		printf("Warning: No output file given. Using default instead...\n");
		strcpy(run_info.output_file, DEFAULT_OUTPUT_FILENAME);
		warn_count++;
	}
	if(run_info.memory_limit == 0) {
		printf("Warning: Memory limit not set or set to 0. Using default value instead...\n");
		run_info.memory_limit = DEFAULT_MEMORY_LIMIT;
		warn_count++;
	}

	if(err_count > 0) {
		fprintf(stderr, "Errors encountered in setup, aborting.\n");
		exit(2);
	}
	err_count += complete_run_info(run_info, tmp_idx, err_count);
	if(err_count > 0) {
		fprintf(stderr, "Errors encountered in setup, aborting.\n");
		exit(2);
	}

	if(warn_count > 0) {
		printf("Warnings issued in setup\n\n");
	}
}

void print_run_info(RunInfo& run_info)
{
	printf("Run Details:\n");
	printf("\tRun Type: ");
	switch(run_info.run_type) {
		case RunType::combinatorial:
			printf("Combinatorial run\n");
			break;
		case RunType::build_from_index:
			printf("Build from index set\n");
			break;
		case RunType::build_from_index_list:
			printf("Build multiple from index set file\n");
			break;
		default:
			printf("Undefined\n");
	}
	printf("\tOutput File: %s\n", run_info.output_file);
	printf("\tSequence: %s\n", run_info.sequence);
	switch(run_info.run_type) {
		case RunType::build_from_index:
			printf("Index Set: ");
			for(uint32_t i = 0; i < run_info.n_fragments; i++) {
				printf("%u", run_info.index[i]);
				if(i != run_info.n_fragments - 1) { printf("-"); }
			}
			break;
		case RunType::build_from_index_list:
			printf("Index Set File: %s\n", run_info.index_list_file);
			break;
		default:
			break;
	}
	printf("\tStructure Type: ");
	switch(run_info.structure_type) {
		case StrType::single_strand:
			printf("Single Stranded\n");
			printf("\tNumber of fragments per structure: %lu\n", run_info.n_fragments);
			break;
		case StrType::double_strand:
			printf("Double Stranded\n");
			printf("\tNumber of fragments per structure: Strand 1 = %lu, Strand 2 = %lu\n", 
					run_info.ds_strand1_n_frags, run_info.ds_strand2_n_frags);
			break;
		default:
			break;
	}
	printf("\tFragment Library Prototype: %s\n", run_info.library_prototype);
	printf("\tFragment Overlap RMSD Limit: %f\n", run_info.rmsd_limit);
	if(run_info.run_options & RunOpts::use_structure_filter) {
		printf("\tSecondary Structure: %s\n", run_info.dot_bracket);
		printf("\tNumber of Watson-Crick pairs: %lu\n", run_info.n_wc_pairs);
		printf("\tWatson-Crick Pair Library Prototype: %s\n", run_info.wc_library_prototype);
		printf("\tWatson-Crick Pair Overlap RMSD Limit: %f\n", run_info.wc_rmsd_limit);
	}
	printf("\tNumber of Structures: ");
	if(run_info.run_options & RunOpts::blind_build_limit) {
		printf("%lu. (Exits once %lu structures are built)\n", run_info.build_limit, run_info.build_limit);
	} else if (run_info.run_options & RunOpts::build_limit_by_energy) {
		printf("%lu. (Will build all structures, but print the lowest %lu)\n", run_info.build_limit, run_info.build_limit);
	} else {
		printf("No limit\n");
	}
	printf("\tBuffer Memory Limit (kB): %lu\n", run_info.memory_limit / 1000);
}

void parse_long_cmdline_opt(char* ARGV[], RunInfo& run_info, size_t idx, char* idx_tmp)
{
	char* opt = &ARGV[idx][2];
	char* field = ARGV[idx + 1];
	if(!strcasecmp(opt, "index")) {	
		strcpy(idx_tmp, field);
	} else if(!strcasecmp(opt, "index-list")) {
		strcpy(run_info.index_list_file, field);
	} else if(!strcasecmp(opt, "rmsd-lim")) {
		run_info.rmsd_limit = atof(field);
	} else if(!strcasecmp(opt, "wc-rmsd-lim")) {
		run_info.wc_rmsd_limit = atof(field);
	} else if (!strcasecmp(opt, "run-type")) {
		if (!strcasecmp(field, "cmb") || !strcasecmp(field, "combinatorial")) {
			run_info.run_type = RunType::combinatorial;
		} else if (!strcasecmp(field, "idx") || !strcasecmp(field,"from-index")) {
			run_info.run_type = RunType::build_from_index;
		}
		if (!strcasecmp(field, "idx-list") || !strcasecmp(field, "from-index-list")) {
			run_info.run_type = RunType::build_from_index_list;
		}
	} else if (!strcasecmp(opt, "str-count-limit")) {
		run_info.build_limit = atoi(field);
	} else if (!strcasecmp(opt, "str-count-limit-type")) {
		if (!strcasecmp(field, "energy"))	{
			run_info.run_options |= RunOpts::build_limit_by_energy;
		} else if (!strcasecmp(field, "none"))	{
			run_info.run_options |= RunOpts::blind_build_limit;
		}
	} else if (!strcasecmp(opt, "lib-prototype")) {
		strcpy(run_info.library_prototype, field);
	} else if (!strcasecmp(opt, "wc-lib-prototype")) {
		strcpy(run_info.wc_library_prototype, field);
	} else if (!strcasecmp(opt, "secondary-structure-filter")) {
		strcpy(run_info.dot_bracket, field);
	} else if (!strcasecmp(opt, "parallel-lib-lengths")) {
		strcpy(run_info.parallel_lib_len, field);
	} else if (!strcasecmp(opt, "parallel-lib-index")) {
		strcpy(run_info.parallel_lib_idx, field);
	} else {
		printf("Unknown argument flag: %s. Ignoring...\n", ARGV[idx]);
	}
}

void input_handler(int argc, char *ARGV[], RunInfo& run_info)
{
	char idx_tmp[GLOBAL_STANDARD_STRING_LENGTH] = {'\0'};
	if (argc == 1) {
		printf("No Input\n");
		exit(1);
	}
	for (int i = 0; i < argc; i++)
	{
		if (ARGV[i][0] == '-') {
			switch (ARGV[i][1])
			{
				case 'i':
					read_input_file(ARGV[i + 1], run_info, idx_tmp);
					break;
				case 'o':	
					strcpy(run_info.output_file, ARGV[i + 1]);
					break;
				case 's':
					strcpy(run_info.sequence, ARGV[i + 1]);
					break;
				case 't':
					if (!strcasecmp(ARGV[i + 1], "ss")) {
						run_info.structure_type = StrType::single_strand;
					}
					else if (!strcasecmp(ARGV[i + 1], "ds")) {
						run_info.structure_type = StrType::double_strand;
					}
					break;
				case 'c':
					run_info.run_options |= RunOpts::write_coordinates;
					break;
				case '-':
					parse_long_cmdline_opt(ARGV, run_info, i, idx_tmp);
					break;
				default:
					printf("Unknown argument flag: %s. Ignoring...\n", ARGV[i]);
					break;
			}
		}
	}
	validate_run_info(run_info, idx_tmp);
	print_run_info(run_info);
}

void free_libs(char **A, int s)
{
	for (int i = 0; i < s; i++)
		free(A[i]);
	free(A);
}
