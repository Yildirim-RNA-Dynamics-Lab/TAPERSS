#ifndef CMB_MANAGER_HPP
#define CMB_MANAGER_HPP

#include "RNACMB.hpp"
#include "Atom_info.hpp"
#include "DimerLib.hpp"
#include "RNAData.hpp"
#include "OutputString.hpp"
#include "WatsonCrickPair.hpp"

struct CMB_Manager
{
	IndexPair<int32_t> *lib_bounds;         //Number of structures in each library
	int32_t *attach_attempted;      //Tracks if all structures in each respective library has been tested
	int32_t last_attempted[2];
	uint32_t count;
	int32_t rollback_count;
	bool *libs_completed;

	CMB_Manager(DimerLibArray& LA, RunInfo& run_info);
	~CMB_Manager();
	void attach_attempt(int i, int j);
	inline bool is_at_end();
	bool check_lib_completion();
	void clear_attempts();
	uint64_t get_reset_count();
};


inline bool CMB_Manager::is_at_end()
{
	return last_attempted[1] == lib_bounds[last_attempted[0]].idx2 - 1;
}

#endif
