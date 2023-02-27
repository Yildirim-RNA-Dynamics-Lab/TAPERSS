#include "Combinatorial_Addition.hpp"
template <STRUCTFILTER_TYPE FILTER> bool early_exit_handler(CMB_Manager &manager)
{
    if constexpr (FILTER == HAIRPIN)
    {
        if (manager.hairpins_built == GLOBAL_STRUCTURE_LIMIT_COUNT)
        {
            return true;
        }
        return false;
    }
    if constexpr (FILTER == INTERNAL_LOOP)
    {
        if (manager.internal_loops_built == GLOBAL_STRUCTURE_LIMIT_COUNT)
        {
            return true;
        }
        return false;
    }
    if constexpr (FILTER == NONE)
    {
        if (manager.strs_built == GLOBAL_STRUCTURE_LIMIT_COUNT)
        {
            return true;
        }
        return false;
    }
}
template bool early_exit_handler<STRUCTFILTER_TYPE::HAIRPIN>(CMB_Manager &manager);
template bool early_exit_handler<STRUCTFILTER_TYPE::INTERNAL_LOOP>(CMB_Manager &manager);
template bool early_exit_handler<STRUCTFILTER_TYPE::NONE>(CMB_Manager &manager);

template <STRUCTFILTER_TYPE FILTER> void output_handler(CMB_Manager &manager, RNADataArray &assembled, output_string &o_string)
{
    if constexpr (FILTER == HAIRPIN)
    {
        double RMSD;
        WC_prepare_structure_matrix(0, assembled[0], 0, assembled[assembled.iterator_max], 1);
        RMSD = WC_check_pair(0);
        if (RMSD <= GLOBAL_WC_RMSD_LIMIT)
        {
            assembled.update_WC_rmsd(RMSD);
            assembled.update_energy();
            o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
            manager.hairpins_built++;
        }
        manager.strs_built++;
        return;
    }
    if constexpr (FILTER == INTERNAL_LOOP)
    {
        double RMSD;
        WC_prepare_structure_matrix(0, assembled[0], 0, assembled[assembled.iterator_max], 1);
        RMSD = WC_check_pair(0);
        if (RMSD <= GLOBAL_WC_RMSD_LIMIT)
        {
            assembled.update_WC_rmsd(RMSD);
            assembled.update_energy();
            o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
            manager.internal_loops_built++;
        }
        manager.strs_built++;
        return;
    }
    if constexpr (FILTER == NONE)
    {
        assembled.update_energy();
        o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
        manager.strs_built++;
    }
}
template void output_handler<STRUCTFILTER_TYPE::HAIRPIN>(CMB_Manager &manager, RNADataArray &assembled, output_string &o_string);
template void output_handler<STRUCTFILTER_TYPE::INTERNAL_LOOP>(CMB_Manager &manager, RNADataArray &assembled, output_string &o_string);
template void output_handler<STRUCTFILTER_TYPE::NONE>(CMB_Manager &manager, RNADataArray &assembled, output_string &o_string);

void fill_if_empty(CMB_Manager &manager, RNADataArray &assembled, DimerLibArray &Lib)
{
    if (assembled.is_empty())
    {
        uint32_t dnmp_pos = assembled.iterator + 1;
        uint32_t lib_pos = manager.attach_attempted[dnmp_pos] + 1;
        assembled.overwrite(dnmp_pos, lib_pos, Lib);
        assembled.keep();
        manager.attach_attempt(dnmp_pos, lib_pos);
    }
}

template <attach_status status> bool rollback_handler(CMB_Manager &manager, RNADataArray &assembled, DimerLibArray &Lib)
{
    if (manager.is_at_end())
    {
        if (manager.check_lib_completion())
        {
            return true;
        }
        if constexpr (status == attach_status::ATTACHED) 
        { 
            assembled.rollback_by(manager.rollback_count + 1); //+1 b/c we called keep(), so we need to rollback an extra dnmp to move on.
            fill_if_empty(manager, assembled, Lib);
            return false; 
        } 
        if constexpr (status == attach_status::FAILED)   
        { 
            assembled.rollback_by(manager.rollback_count); 
            fill_if_empty(manager, assembled, Lib);
            return false; 
        }
    }
    else
    {
        assembled.rollback();
        return false;
    }  
}
template bool rollback_handler<attach_status::ATTACHED>(CMB_Manager &manager, RNADataArray &assembled,  DimerLibArray &Lib);
template bool rollback_handler<attach_status::FAILED>(CMB_Manager &manager, RNADataArray &assembled,  DimerLibArray &Lib);

template <STRUCTFILTER_TYPE FILTER> bool combinatorial_addition(DimerLibArray &Lib, RNADataArray &assembled, CMB_Manager &manager, output_string &o_string)
{
    switch(fragment_assembly(Lib, assembled, manager))
    {
    case attach_status::ATTACHED:
        assembled.keep();
        break;
    case attach_status::FAILED:
        return rollback_handler<attach_status::FAILED>(manager, assembled, Lib);
    default:      //This condition will never be reached..only here to stop compiler warnings
        break;
    }

    if (assembled.is_complete())
    {
        output_handler<FILTER>(manager, assembled, o_string);
        if constexpr (STRUCTURE_BUILD_LIMIT){ if(early_exit_handler<FILTER>(manager)){return true;} }
        return rollback_handler<attach_status::ATTACHED>(manager, assembled, Lib);
    }
    return false;
}
template bool combinatorial_addition<STRUCTFILTER_TYPE::HAIRPIN>(DimerLibArray &Lib, RNADataArray &assembled, CMB_Manager &manager, output_string &o_string);
template bool combinatorial_addition<STRUCTFILTER_TYPE::NONE>(DimerLibArray &Lib, RNADataArray &assembled, CMB_Manager &manager, output_string &o_string);

bool combinatorial_addition_IL(DimerLibArray &Lib, RNADataArrayInternalLoop &assembled, CMB_Manager &manager, output_string &o_string, DimerLibArray &WC_Lib)
{
    
    switch(fragment_assembly_IL(Lib, WC_Lib, assembled, manager))
    {
    case attach_status::ATTACHED:
        assembled.keep();
        break;
    case attach_status::FAILED:
        return rollback_handler<attach_status::FAILED>(manager, assembled, Lib);
    default:      //This condition will never be reached..only here to stop compiler warnings
        break;
    }

    if (assembled.is_complete())
    {
        output_handler<STRUCTFILTER_TYPE::INTERNAL_LOOP>(manager, assembled, o_string);
        if constexpr (STRUCTURE_BUILD_LIMIT){ if(early_exit_handler<STRUCTFILTER_TYPE::INTERNAL_LOOP>(manager)){return true;} }
        return rollback_handler<attach_status::ATTACHED>(manager, assembled, Lib);
    }
    return false;
}
