#ifndef CMB_MANAGER_HPP
#define CMB_MANAGER_HPP

#include "RNACMB.hpp"
#include "Atom_info.hpp"
#include "DimerLib.hpp"
#include "RNAData.hpp"
#include "output_string.hpp"
#include "Hairpin.hpp"

struct CMB_Manager
{
    uint64_t *count_per_lib;         //Number of structures in each library
    bool **attach_attempted;    //Tracks if all structures in each respective library has been tested
    uint64_t last_attempted[2];
    bool *libs_completed;
    uint_fast64_t strs_built;
    uint_fast64_t hairpins_built;
    uint_fast64_t internal_loops_built;
    uint64_t count;                  //For deallocation

    CMB_Manager(DimerLibArray& LA);
    ~CMB_Manager();
    void attach_attempt(int i, int j);
    bool is_at_end();
    void check_lib_completion();
    void clear_attempts();
    uint64_t get_reset_count();
    void successful_construction();
};

#endif