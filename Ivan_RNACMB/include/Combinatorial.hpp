#ifndef COMBINATORIAL_HPP
#define COMBINATORIAL_HPP

#include "RNACMB.hpp"
#include "Atom_info.hpp"
#include "DimerLib.hpp"
#include "RNA_data.hpp"
#include "output_string.hpp"
#include "Hairpin.hpp"

enum attach_status{FAILED, ATTACHED, NOT_CHECKED};

struct CMB_Manager
{
    int *count_per_lib;         //Number of structures in each library
    bool **attach_attempted;    //Tracks if all structures in each respective library has been tested
    int last_attempted[2];
    bool *libs_completed;
    int strs_built;
    int hairpins_built;
    int count;                  //For deallocation

    CMB_Manager(DimerLibArray& LA);
    ~CMB_Manager();
    void attach_attempt(int i, int j);
    bool is_at_end();
    void check_lib_completion();
    void clear_attempts();
    int get_reset_count();
    void successful_construction();
};

void update_energy(RNA_data_array& sequence);
attach_status rotate(RNA_data *reference, RNA_data *rotated);
bool overlap_check(RNA_data_array& sequence, RNA_data *attach);
attach_status check_attachment(RNA_data_array& sequence, RNA_data *attach);
bool combinatorial_addition(DimerLibArray& Lib, RNA_data_array& assembled, CMB_Manager& manager, output_string& o_string, DimerLibArray& WC_Lib);

#endif