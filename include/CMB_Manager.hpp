#ifndef CMB_MANAGER_HPP
#define CMB_MANAGER_HPP

#include "RNACMB.hpp"
#include "Atom_info.hpp"
#include "DimerLib.hpp"
#include "RNAData.hpp"
#include "output_string.hpp"
#include "WatsonCrickPair.hpp"

struct CMB_Manager
{
    int32_t *count_per_lib;         //Number of structures in each library
    int32_t *attach_attempted;      //Tracks if all structures in each respective library has been tested
    int32_t last_attempted[2];
    uint32_t count;
    int32_t rollback_count;
    bool *libs_completed;
    uint64_t strs_built;
    uint64_t hairpins_built;
    uint64_t internal_loops_built;                  

    CMB_Manager(DimerLibArray& LA);
    ~CMB_Manager();
    void attach_attempt(int i, int j);
    inline bool is_at_end();
    bool check_lib_completion();
    void clear_attempts();
    uint64_t get_reset_count();
    void successful_construction();
};


inline bool CMB_Manager::is_at_end()
{
    return last_attempted[1] == count_per_lib[last_attempted[0]] - 1;
}

#endif