#include "InputHandler.hpp"

void get_diNt_names(char *sequence, char** rtn, int N_diNts)
{
	int name_position = 0;
	//*N_diNts = strlen(sequence) - 1;
	//char **rtn = (char **)malloc(sizeof(char *) * *N_diNts);
	while (LIBRARY_FILENAME_PROTOTYPE[name_position] != 'X')
	{
		name_position++;
	};
	//name_position -= 1;

	for (int i = 0; i < N_diNts; i++)
	{
		rtn[i] = (char *)malloc(sizeof(LIBRARY_FILENAME_PROTOTYPE));
		memcpy(rtn[i], LIBRARY_FILENAME_PROTOTYPE, sizeof(LIBRARY_FILENAME_PROTOTYPE));
		memcpy(&rtn[i][name_position], sequence++, sizeof(char) * 2);
	}
}

int find_duplicates(char** rtn, int *duplicates, int N_diNts)
{
	uint num_duplicates = 0;
	for (int i = N_diNts - 1; i > -1; i--)
	{
		for(int j = 0; j < i; j++)
		{
			if(!strcmp(rtn[i], rtn[j]))
			{
				duplicates[i] = j;
				num_duplicates++;
			}
		}
	}
	//Comment
	return num_duplicates;
}

void get_WC_partner(char *sequence, char** rtn, int Nt1, int Nt2)
{
	char WC_pair[3];

	int name_position = 0;
	while (WATSON_CRICK_LIBRARY_PROTOTYPE[name_position] != 'X')
	{
		name_position++;
	}

	WC_pair[0] = sequence[Nt1];
	WC_pair[1] = sequence[Nt2];
	WC_pair[2] = '\0';

	if (WC_pair[0] == 'A')
	{
		if (WC_pair[1] != 'C' && WC_pair[1] != 'U')
		{
			printf("Not potential hairpin, exiting... ( this is a temporary fix while program is dumb )\n");
			printf("Pair Attempted: %s\n", WC_pair);
			exit(2);
		}
	}
	else if (WC_pair[0] == 'U')
	{
		if (WC_pair[1] != 'A' && WC_pair[1] != 'G')
		{
			printf("Not potential hairpin, exiting... ( this is a temporary fix while program is dumb )\n");
			printf("Pair Attempted: %s\n", WC_pair);
			exit(2);
		}
	}
	else if (WC_pair[0] == 'C')
	{
		if (WC_pair[1] != 'A' && WC_pair[1] != 'G')
		{
			printf("Not potential hairpin, exiting... ( this is a temporary fix while program is dumb )\n");
			printf("Pair Attempted: %s\n", WC_pair);
			exit(2);
		}
	}
	else if (WC_pair[0] == 'G')
	{
		if (WC_pair[1] != 'C' && WC_pair[1] != 'U')
		{
			printf("Not potential hairpin, exiting... ( this is a temporary fix while program is dumb. Instead it should look through whole sequence to find h-bond pairs and then check if hairpin is possible)\n");
			printf("Pair Attempted: %s\n", WC_pair);
			exit(2);
		}
	}

	*rtn = (char *)malloc(sizeof(WATSON_CRICK_LIBRARY_PROTOTYPE));
	memcpy(*rtn, WATSON_CRICK_LIBRARY_PROTOTYPE, sizeof(WATSON_CRICK_LIBRARY_PROTOTYPE));
	memcpy(&rtn[0][name_position], WC_pair, sizeof(char) * 2);
}

void get_index_int(char *index, uint32_t *indices)
{
	char buffer[5];
	int count = 0, buf_c = 0;
	for (char c = *index; c != '\0'; c = *(++index))
	{
		if (c == '-')
		{
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
}

int read_input_index_file(char *file_name, int n_DiNts)
{
	FILE *input = fopen(file_name, "r");
	char line[GLOBAL_STANDARD_STRING_LENGTH];
	char indices[GLOBAL_STANDARD_STRING_LENGTH];
	int str_count = 0;

	while (fgets(line, sizeof(line), input)) // Get line count
	{
		char *header;
		header = strtok(line, " ");
		if (!strcmp(header, GLOBAL_INPUT_SEQUENCE))
			str_count++;
	}
	rewind(input);

	GLOBAL_INPUT_INDICES_LIST = (uint32_t **)malloc(str_count * sizeof(int *));

	for (int i = 0; i < str_count; i++)
	{
		GLOBAL_INPUT_INDICES_LIST[i] = (uint32_t *)malloc(n_DiNts * sizeof(int));
		if(fgets(line, sizeof(line), input) != NULL)
		{
			char *str1;
			char *header = strtok(line, " ");
			if (!strcmp(header, GLOBAL_INPUT_SEQUENCE))
			{
				str1 = strtok(NULL, " ");
				str1 = strtok(NULL, " ");
				str1[strcspn(str1, "\t")] = '\0'; // equivalent to chomp() from perl
				strcpy(indices, str1);
			}
			get_index_int(indices, GLOBAL_INPUT_INDICES_LIST[i]);
		}
	}
	fclose(input);
	return str_count;
}

void trim_whitespace(char* s)
{
	const char whitespace[] = {' ', '\t', '\n'};
	size_t str_len = strlen(s);
	size_t match_loc = strcspn(s, whitespace);
	unsigned int idx = 0;

	if(match_loc == 0)
	{
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

void read_input_file(char *file_name, RunInfo& run_info, char* idx_tmp)
{
	FILE *input = fopen(file_name, "r");
	if(input == nullptr)
	{
		printf("Could not open file: %s\n", file_name);
		exit(3);
	}
	GLOBAL_OUTPUT_FILE[0] = '\0';
	char line[GLOBAL_STANDARD_STRING_LENGTH];
	while (fgets(line, sizeof(line), input))
	{
		char *header = strtok(line, "=");
		char *str1 = strtok(NULL, "=");
		trim_whitespace(header);
		trim_whitespace(str1);
		if (!strcasecmp(header, "sequence"))
		{
			strcpy(run_info.sequence, str1);
			continue;
		}
		else if (!strcasecmp(header, "OUTPUT-FILE"))
		{
			strcpy(run_info.output_file, str1);
			continue;
		}
		else if (!strcasecmp(header, "RUN-TYPE"))
		{
			if (!strcasecmp(str1, "cmb") || !strcasecmp(str1, "combinatorial"))
			{
				run_info.run_type = RunType::combinatorial;
			}
			else if (!strcasecmp(str1, "idx") || !strcasecmp(str1,"from-index"))
			{
				run_info.run_type = RunType::build_from_index;
			}
			if (!strcasecmp(str1, "idx-list") || !strcasecmp(str1, "from-index-list"))
			{
				run_info.run_type = RunType::build_from_index_list;
			}
			continue;
		}
		else if (!strcasecmp(header, "RNA-TYPE"))
		{
			if (!strcasecmp(str1, "singe-strand"))
			{
				run_info.structure_type = StrType::single_strand;
			}
			if (!strcasecmp(str1, "double-strand"))	{
				run_info.structure_type = StrType::double_strand;
			}
		}
		else if (!strcasecmp(header, "INDEX-LIST-FILE"))
		{
			strcpy(run_info.index_list_file, str1);
			continue;
		}
		else if (!strcasecmp(header, "structure-count-limit"))
		{
			run_info.build_limit = atoi(str1);
		}
		else if (!strcasecmp(header, "structure-count-limit-type"))
		{
			if (!strcasecmp(str1, "energy"))	{
				run_info.run_options |= RunOpts::build_limit_by_energy;
			}
			else if (!strcasecmp(str1, "none"))	{
				run_info.run_options |= RunOpts::blind_build_limit;
			}
			continue;
		}
		else if (!strcasecmp(header, "single-index-set"))
		{
			strcpy(idx_tmp, str1);
			continue;
		}
		else if (!strcasecmp(header, "overlap-rmsd-limit"))
		{
			run_info.rmsd_limit = atof(str1);
			continue;
		}
		else if (!strcasecmp(header, "watson-crick-rmsd-limit"))
		{
			run_info.wc_rmsd_limit = atof(str1);
			continue;
		}
		else if (!strcasecmp(header, "library-prototype"))
		{
			strcpy(run_info.library_prototype, str1);
			continue;
		}
		else if (!strcasecmp(header, "watson-crick-library-prototype"))
		{
			strcpy(run_info.wc_library_prototype, str1);
			continue;
		}
		else if (!strcasecmp(header, "WRITE-COORDINATES"))
		{
			if (!strcasecmp(str1, "TRUE"))
			{
				run_info.run_options |= RunOpts::write_coordinates;
			}
		}
		else
		{
			printf("Warning: Ignored unknown option in input file: %s\n", header);
		}
	}
	
	fclose(input);
}

void validate_run_info(RunInfo& run_info, char* tmp_idx)
{
	size_t err_count = 0;
	size_t warn_count = 0;
	if(run_info.sequence[0] == '\0')
	{
		fprintf(stderr, "Error: No sequence given in input file or as argument!\n");
		err_count++;
	}
	if(run_info.library_prototype[0] == '\0')
	{
		fprintf(stderr, "Error: No library file name prototype given in input file or as argument!\n");
		err_count++;
	}
	if((run_info.run_options & RunOpts::use_h_bond_filter) || (strchr(run_info.sequence, 'x') != NULL))
	{
		if(run_info.wc_library_prototype[0] == '\0')
		{
			fprintf(stderr, "Error: No watson crick library file name prototype given in input file or as argument!\n");
			err_count++;
		}
	}
	if(run_info.run_type == RunType::build_from_index && tmp_idx[0] == '\0')
	{
		fprintf(stderr, "Error: Run is of type: \"Build from index\" but no index set is given!\n");
		err_count++;
	}
	if(run_info.run_type == RunType::build_from_index_list && run_info.index_list_file[0] == '\0')
	{
		fprintf(stderr, "Error: Run is of type: \"Build from index list\" but no index list file is given!\n");
		err_count++;
	}

	if(run_info.rmsd_limit == 0)
	{
		printf("Warning: Overlap RMSD Limit not given or set to 0. Using default value instead...\n");
		run_info.rmsd_limit = DEFAULT_RMSD_LIMIT;
		warn_count++;
	}
	if(run_info.wc_rmsd_limit == 0)
	{
		printf("Warning: Watson-Crick RMSD Limit not given or set to 0. Using default value instead...\n");
		run_info.wc_rmsd_limit = DEFAULT_WC_RMSD_LIMIT;
		warn_count++;
	}
	if(run_info.build_limit != 0)
	{
		if(!(run_info.run_options & RunOpts::blind_build_limit) || !(run_info.run_options & RunOpts::build_limit_by_energy))
		{
			printf("Warning: Structure Build Limit specified, but limit type not given. Assuming blind limit...\n");
			run_info.run_options |= RunOpts::blind_build_limit;
			warn_count++;
		}
	}
	if((run_info.run_options & RunOpts::blind_build_limit) || (run_info.run_options & RunOpts::build_limit_by_energy))
	{
		if(run_info.build_limit == 0)
		{
			fprintf(stderr, "Error: Build limit type specified but build limit is not set or set to 0!\n");
			err_count++;
		}
	}
	if(run_info.output_file[0] == '\0')
	{
		printf("Warning: No output file given. Using default instead...\n");
		warn_count++;
	}
}

void complete_run_info(RunInfo& run_info)
{

}

void parse_long_cmdline_opt(char* ARGV[], RunInfo& run_info, size_t idx, char* idx_tmp)
{
	char* opt = &ARGV[idx][2];
	char* field = ARGV[idx + 1];
	if(!strcasecmp(opt, "index")) 
	{	
		strcpy(idx_tmp, field);
	}
	else if(!strcasecmp(opt, "index-list"))
	{
		strcpy(run_info.index_list_file, field);
	}
	else if(!strcasecmp(opt, "rmsd-lim"))
	{
		run_info.rmsd_limit = atof(field);
	}
	else if(!strcasecmp(opt, "wc-rmsd-lim"))
	{
		run_info.wc_rmsd_limit = atof(field);
	}
	else if (!strcasecmp(opt, "run-type"))
	{
		if (!strcasecmp(field, "cmb") || !strcasecmp(field, "combinatorial"))
		{
			run_info.run_type = RunType::combinatorial;
		}
		else if (!strcasecmp(field, "idx") || !strcasecmp(field,"from-index"))
		{
			run_info.run_type = RunType::build_from_index;
		}
		if (!strcasecmp(field, "idx-list") || !strcasecmp(field, "from-index-list"))
		{
			run_info.run_type = RunType::build_from_index_list;
		}
	}
	else if (!strcasecmp(opt, "str-count-limit"))
	{
		run_info.build_limit = atoi(field);
	}
	else if (!strcasecmp(opt, "str-count-limit-type"))
	{
		if (!strcasecmp(field, "energy"))	{
			run_info.run_options |= RunOpts::build_limit_by_energy;
		}
		else if (!strcasecmp(field, "none"))	{
			run_info.run_options |= RunOpts::blind_build_limit;
		}
	}
	else if (!strcasecmp(opt, "lib-prototype"))
	{
		strcpy(run_info.library_prototype, field);
	}
	else if (!strcasecmp(opt, "wc-lib-prototype"))
	{
		strcpy(run_info.wc_library_prototype, field);
	}
}

void input_handler(int argc, char *ARGV[], RunInfo& run_info)
{
	char idx_tmp[GLOBAL_STANDARD_STRING_LENGTH] = {'\0'};
	if (argc == 1)
	{
		printf("No Input\n");
		exit(1);
	}
	for (int i = 0; i < argc; i++)
	{
		if (ARGV[i][0] == '-')
		{
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
					if (!strcasecmp(ARGV[i + 1], "ss"))
					{
						run_info.structure_type = StrType::single_strand;
					}
					else if (!strcasecmp(ARGV[i + 1], "ds"))
					{
						run_info.structure_type = StrType::double_strand;
					}
					break;
				case 'c':
					run_info.run_options |= RunOpts::write_coordinates;
					break;
				case '-':
					parse_long_cmdline_opt(ARGV, run_info, i, idx_tmp);
					break;
			}
		}
	}

	if(idx_tmp[0] != '\0' && run_info.sequence[0] != '\0' && run_info.run_type == RunType::build_from_index)
	{
		run_info.index = (uint32_t*)malloc(sizeof(uint32_t) * (strlen(run_info.sequence) - 1));
		get_index_int(idx_tmp, run_info.index);
	}
}

void free_libs(char **A, int s)
{
	for (int i = 0; i < s; i++)
		free(A[i]);
	free(A);
}
