#include "StructureBuilders.hpp"

void create_custom_structure(DimerLibArray &Lib, DimerLibArray &WC_Lib, RNADataArray &assembled, output_string &o_string, int *indices)
{
    RNAData *custom_attach;
    RNAData *base;
    attach_status status;
    assembled.rollback();
    assembled.add_move(new RNAData(Lib, 0, indices[0]));
    for (int i = 1; i <= assembled.iterator_max; i++)
    {
        base = assembled.current();
        custom_attach = new RNAData(Lib, i, indices[i]);
        
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
    assembled.add_move(new RNAData(Lib, 0, indices[0]));
    for (int i = 1; i <= assembled.iterator_max; i++)
    {
        base = assembled.current();
        custom_attach = new RNAData(Lib, i, indices[i]);
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

void create_custom_structure_list(DimerLibArray &Lib, DimerLibArray &WC_Lib, RNADataArray &assembled, output_string &o_string, int num_strs)
{
    RNAData *custom_attach;
    RNAData *base;
    for (int i = 0; i < num_strs; i++)
    {
        assembled.rollback_by(assembled.iterator);
        assembled.add_move(new RNAData(Lib, 0, GLOBAL_INPUT_INDICES_LIST[i][0]));
        for (int j = 1; j <= assembled.iterator_max; j++)
        {
            base = assembled.current();
            custom_attach = new RNAData(Lib, j, GLOBAL_INPUT_INDICES_LIST[i][j]);
            rotate(base, custom_attach);
            custom_attach->update_submatrices();
            assembled.add_move(custom_attach);
        }
        if (GLOBAL_PERFORM_STRUCTCHECK == HAIRPIN)
        {
            is_WC_pair(assembled, WC_Lib, 0, assembled.iterator_max, 0);
        }
        assembled.update_energy();
        o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
    }
}

void create_custom_structure_list_testing(DimerLibArray &Lib, RNADataArray &assembled, output_string &o_string, int num_strs)
{
    RNAData *custom_attach;
    RNAData *base;
    for (int i = 0; i < num_strs; i++)
    {
        assembled.rollback_by(assembled.iterator);
        assembled.add_move(new RNAData(Lib, 0, GLOBAL_INPUT_INDICES_LIST[i][0]));
        for (int j = 1; j <= assembled.iterator_max; j++)
        {
            base = assembled.current();
            custom_attach = new RNAData(Lib, j, GLOBAL_INPUT_INDICES_LIST[i][j]);
            rotate(base, custom_attach);
            custom_attach->update_submatrices();
            SCC_record_COM_distance(assembled, custom_attach);
            assembled.add_move(custom_attach);
        }
        //exit(1);
        o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
    }
}