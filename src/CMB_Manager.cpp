#include "CMB_Manager.hpp"
/**
 * @brief Constructs a new CMB_Manager based on a previously constructed DimerLibArray.
 * 
 * @param &LA reference to DimerLibArray used for fragment libraries (Not WatsonCrick Pairs!!)
**/
CMB_Manager::CMB_Manager(DimerLibArray &LA, RunInfo& run_info)
{
	count = LA.count;

	lib_bounds = (IndexPair<int32_t> *)malloc(sizeof(IndexPair<int32_t>) * count);
	attach_attempted = (int32_t *)malloc(count * sizeof(int32_t));
	libs_completed = (bool *)calloc(count, sizeof(bool));

	memset(attach_attempted, -1, count * sizeof(int32_t));
	for (uint32_t i = 0; i < count; i++) {
		lib_bounds[i].idx1 = run_info.frag_lib_bounds[i].idx1;
		lib_bounds[i].idx2 = run_info.frag_lib_bounds[i].idx2;
		attach_attempted[i] = lib_bounds[i].idx1 - 1;
	}

	attach_attempted[0] = lib_bounds[0].idx1;
	last_attempted[0] = 0;
	last_attempted[1] = lib_bounds[0].idx1;
}

CMB_Manager::~CMB_Manager()
{
	free(attach_attempted);
	free(lib_bounds);
	free(libs_completed);
}

void CMB_Manager::attach_attempt(int32_t i, int32_t j)
{
	attach_attempted[i] = j;
	last_attempted[0] = i;
	last_attempted[1] = j;
}

bool CMB_Manager::check_lib_completion()
{
	rollback_count = 0;
	memset(libs_completed, false, count * sizeof(bool));

	for (int i = last_attempted[0]; i >= 0; i--) {
		if (attach_attempted[i] == lib_bounds[i].idx2 - 1) {// Checks which libraries have had all conformations tested
			//DEBUG_PRINT("Lib %d is at end (%d vs %d); removing and resesting to (%d)\n", i, attach_attempted[i], lib_bounds[i].idx2 - 1, lib_bounds[i].idx1 - 1);
			libs_completed[i] = true;
			attach_attempted[i] = lib_bounds[i].idx1 - 1;                   // Reset attempts for completed lib
			rollback_count++;                           // Track how many libraries have been fully tested.
		} else {                                           // First lib that has not been completed is position to continue from 
			break;
		}
	}
	return rollback_count == last_attempted[0] + 1; //Check if the all libraries have been fully tested (aka all viable combinations have been tested)
}

/* Deprecated */
void CMB_Manager::clear_attempts()
{
	for (uint64_t i = 0; i < count; i++)
	{
		if (libs_completed[i] == true)
		{
			attach_attempted[i] = -1;
		}
	}
}

/* Deprecated */
uint64_t CMB_Manager::get_reset_count()
{
	uint64_t n_reset = 0;
	for (int64_t i = 0; i < last_attempted[0] + 1; i++)
	{
			if (libs_completed[i] == true)
			{
				n_reset++;
			}
	}
	//printf("n_reset = %d\n", n_reset);
	return n_reset;
}

