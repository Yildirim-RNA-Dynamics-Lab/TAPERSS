#include "StructureBuilders.hpp"

void create_custom_structure(DimerLibArray &Lib, DimerLibArray &WC_Lib, RNADataArray &assembled, output_string &o_string, int *indices)
{
    RNAData *custom_attach;
    RNAData *base;
    attach_status status;
    assembled.rollback();
    assembled.add_move(nullptr);//new RNAData(Lib, 0, indices[0]));
    for (int i = 1; i <= assembled.iterator_max; i++)
    {
        base = assembled.current();
        custom_attach = nullptr;//new RNAData(Lib, i, indices[i]);
        
        rotate(base, custom_attach);
        custom_attach->update_submatrices();
        status = check_attachment(assembled, custom_attach);
        assembled.add_move(custom_attach);
        if (status == ATTACHED)
            printf("Successful attachment on: %d\n", i);
        else
            printf("Unsuccessful attachment on: %d\n", i);
    }
    if (GLOBAL_PERFORM_STRUCTCHECK == HAIRPIN)
    {
        is_WC_pair(assembled, WC_Lib, 0, assembled.iterator_max, 0);
    }
    assembled.update_energy();
    o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
}

void create_custom_structure_IL(DimerLibArray &Lib, DimerLibArray &WC_Lib, RNADataArrayInternalLoop &assembled, output_string &o_string, int *indices)
{
    RNAData *custom_attach;
    RNAData *base;
    attach_status status;
    assembled.rollback();
    assembled.add_move(nullptr);//new RNAData(Lib, 0, indices[0]));
    for (int i = 1; i <= assembled.iterator_max; i++)
    {
        base = assembled.current();
        custom_attach = nullptr;//new RNAData(Lib, i, indices[i]);
        printf("i = %d\n", i);
        if(assembled.inLeft_or_inRight(i) == true) 
        {
            printf("i = %d\n", i);
            assembled.prepare_right(custom_attach, WC_Lib); 
            assembled.add_move(custom_attach);
            continue;
            /*
            o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
            o_string.add_string((char *)"!!!!!!Unexpected steric clash!!!!!!\n", sizeof("!!!!!!Unexpected steric clash!!!!!!\n"));
            */
        }
        rotate(base, custom_attach);
        custom_attach->update_submatrices();
        status = check_attachment(assembled, custom_attach);
        assembled.add_move(custom_attach);
        if (status == ATTACHED)
            printf("Successful attachment on: %d\n", i);
        else
            printf("Unsuccessful attachment on: %d\n", i);
    }
    is_WC_pair(assembled, WC_Lib, 0, assembled.iterator_max, 0);
    assembled.update_energy();
    o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
}

template <bool PerformChecks, STRUCTCHECK_TYPE StructCheck> void create_custom_structure_list(DimerLibArray &Lib, DimerLibArray &WC_Lib, RNADataArray &assembled, output_string &o_string, int num_strs)
{
    RNAData *custom_attach;
    RNAData *base;
    attach_status status = ATTACHED;
    RNAData *WC_pair = nullptr;//new RNAData(WC_Lib, 0, 0, true);
    gsl_matrix *WC_model_matrix = make_WC_submatrix(WC_pair, WC_pair);
    gsl_matrix *Sequence_matrix = make_WC_submatrix(WC_pair, WC_pair); //Preallocation of Sequence matrix only... This will be improved later

    if(PerformChecks == true)
    {
        for (int i = 0; i < num_strs; i++)
        {
            assembled.rollback_by(assembled.iterator);
            assembled.add_move(nullptr);//new RNAData(Lib, 0, GLOBAL_INPUT_INDICES_LIST[i][0]));
            for (int j = 1; j <= assembled.iterator_max; j++)
            {
                //printf("j = %d\n", j);
                //printf("Iterator = %d\n", assembled.iterator);
                base = assembled.current();
                custom_attach = nullptr;//new RNAData(Lib, j, GLOBAL_INPUT_INDICES_LIST[i][j]);
                rotate(base, custom_attach);
                custom_attach->update_submatrices();
                status = check_attachment(assembled, custom_attach);
                if(status != ATTACHED)
                {
                    //assembled.rollback_by(j + 1);
                    //printf("Failed Attach\n");
                    break;
                }
                else
                {
                    //printf("Attached\n");
                    assembled.add_move(custom_attach);
                }            
            }
            if(status == ATTACHED)
            {
                if constexpr (StructCheck == HAIRPIN)
                {
                    assembled[0]->make_WC_submatrices();
                    assembled[1]->make_WC_submatrices();
                    overwrite_WC_submatrix_gsl(assembled[0]->get_WC_target_matrix(0), assembled[assembled.iterator_max]->get_WC_target_matrix(1), Sequence_matrix);
                    RMSD_WC_pair_gsl(WC_model_matrix, Sequence_matrix);
                }
                assembled.update_energy();                  
                o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
            }
        }
    }
    else
    {
        for (int i = 0; i < num_strs; i++)
        {
            assembled.rollback_by(assembled.iterator);
            assembled.add_move(nullptr);//new RNAData(Lib, 0, GLOBAL_INPUT_INDICES_LIST[i][0]));
            for (int j = 1; j <= assembled.iterator_max; j++)
            {                    
                base = assembled.current();
                custom_attach = nullptr;//new RNAData(Lib, j, GLOBAL_INPUT_INDICES_LIST[i][j]);
                rotate(base, custom_attach);
                custom_attach->update_submatrices();
                assembled.add_move(custom_attach);           
            }
            if constexpr (StructCheck == HAIRPIN)
            {
                assembled[0]->make_WC_submatrices();
                assembled[assembled.iterator_max]->make_WC_submatrices();
                overwrite_WC_submatrix_gsl(assembled[0]->get_WC_target_matrix(0), assembled[assembled.iterator_max]->get_WC_target_matrix(1), Sequence_matrix);
                assembled.update_WC_rmsd(RMSD_WC_pair_gsl(WC_model_matrix, Sequence_matrix));
            }
            assembled.update_energy();                  
            o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
        }
    }
    gsl_matrix_free(WC_model_matrix);
    gsl_matrix_free(Sequence_matrix);
    delete(WC_pair);
}
template void create_custom_structure_list<PERFORM_CHECKS_ON_CUSTOM_BUILD, HAIRPIN>(DimerLibArray& Lib, DimerLibArray& WC_Lib,RNADataArray& assembled, output_string& o_string, int num_strs);
template void create_custom_structure_list<PERFORM_CHECKS_ON_CUSTOM_BUILD, NONE>(DimerLibArray& Lib, DimerLibArray& WC_Lib,RNADataArray& assembled, output_string& o_string, int num_strs);

void create_custom_structure_list_testing(DimerLibArray &Lib, RNADataArray &assembled, output_string &o_string, int num_strs)
{
    RNAData *custom_attach;
    RNAData *base;
    for (int i = 0; i < num_strs; i++)
    {
        assembled.rollback_by(assembled.iterator);
        assembled.add_move(nullptr);//new RNAData(Lib, 0, GLOBAL_INPUT_INDICES_LIST[i][0]));
        for (int j = 1; j <= assembled.iterator_max; j++)
        {
            base = assembled.current();
            custom_attach = nullptr;//new RNAData(Lib, j, GLOBAL_INPUT_INDICES_LIST[i][j]);
            rotate(base, custom_attach);
            custom_attach->update_submatrices();
            SCC_record_COM_distance(assembled, custom_attach);
            assembled.add_move(custom_attach);
        }
        //exit(1);
        o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
    }
}