#include "n_lowest_energy.hpp"

float* NLowestE;
size_t* NLowestIdx;
size_t N;
size_t count;
size_t output_iterator;
void create_n_lowest_E(size_t n, RNADataArray& model)
{
	N = n;
	count = model.iterator_max + 1;
	output_iterator = 0;
	NLowestE = (float *)malloc(sizeof(float) * N);
	NLowestIdx = (size_t *)malloc(sizeof(size_t) * N * (model.iterator_max + 1));
	for(uint i = 0; i < n; i++) {
		NLowestE[i] = (float)INT64_MAX;
		for(uint j = 0; j < count; j++) {
			NLowestIdx[j + (i * count)] = INT64_MAX;
		}
	}
}

void destroy_n_lowest_E()
{
	free(NLowestE);
	free(NLowestIdx);
}

void shift_following(float E, int i, RNADataArray& model)
{
	float prev = NLowestE[i];
	float next;
	memmove(&NLowestIdx[(i + 1) * (count)], &NLowestIdx[i * (count)], sizeof(size_t) * ((N - (1 + i)) * (count)));
	for(uint j = 0; j < count; j++) {
		uint idx = j  + (i * count);
		NLowestIdx[idx]	= model[j]->position_in_lib[1];
	}
	NLowestE[i] = E;
	for(uint j = i + 1; j < N; j++) {
		next = NLowestE[j];
		NLowestE[j] = prev;
		prev = next;
	}
}

void print_lowest()
{
	DEBUG_PRINT("Lowest List:\n");
	for(uint i = 0; i < N; i++) {
		printf("[%d] = %f: (", i, NLowestE[i]);
		printf("%lu", NLowestIdx[i * count]);
		for(uint j = 1; j < count; j++) {
			printf("-%lu", NLowestIdx[j + (i * count)]);
		}
		printf(")\n");
	}
}

void update_run_info_index(RunInfo& run_info) {
	for(uint i = 0; i < count; i++) {
		run_info.index[i] = NLowestIdx[i + (output_iterator * count)];
	}
	output_iterator++;
}

int add_to_n_lowest_E(float E, RNADataArray& model)
{
	uint i = 0;
	if(E < NLowestE[0]) {
		//DEBUG_PRINT("Replacing [%d] = %f with: %f\n", i, NLowestE[i], E);
		shift_following(E, 0, model);
		//print_lowest();
		return 0;
	}
	for(i = N-1; i >= 1; i--) {
		if(E > NLowestE[i]) {
			break;	
		}
	}
	i++;
	if(i <= N-1) {
		//DEBUG_PRINT("Replacing [%d] = %f with: %f\n", i, NLowestE[i], E);
		shift_following(E, i, model);
		//print_lowest();
		return i;
	}
	return -1;
}
