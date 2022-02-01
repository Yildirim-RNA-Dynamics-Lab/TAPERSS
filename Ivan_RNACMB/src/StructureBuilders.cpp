#include "StructureBuilders.hpp"

void create_custom_structure(DimerLibArray &Lib, DimerLibArray &WC_Lib, RNA_data_array &assembled, output_string &o_string, int *indices)
{
    RNA_data *custom_attach;
    RNA_data *base;
    attach_status status;
    assembled.rollback();
    assembled.add_move(new RNA_data(Lib, 0, indices[0]));
    for (int i = 1; i <= assembled.iterator_max; i++)
    {
        base = assembled.current();
        custom_attach = new RNA_data(Lib, i, indices[i]);
        rotate(base, custom_attach);
        custom_attach->update_submatrices();
        status = check_attachment(assembled, custom_attach);
        assembled.add_move(custom_attach);
        if (status == ATTACHED)
            printf("Successful attachment on: %d\n", i);
        else
            printf("Unsuccessful attachment on: %d\n", i);
    }
    if (GLOBAL_PERFORM_HAIRPIN_CHECK)
    {
        is_hairpin(assembled, WC_Lib);
    }
    assembled.update_energy();
    o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
}

void create_custom_structure_list(DimerLibArray &Lib, DimerLibArray &WC_Lib, RNA_data_array &assembled, output_string &o_string, int num_strs)
{
    RNA_data *custom_attach;
    RNA_data *base;
    for (int i = 0; i < num_strs; i++)
    {
        assembled.rollback_by(assembled.iterator);
        assembled.add_move(new RNA_data(Lib, 0, GLOBAL_INPUT_INDICES_LIST[i][0]));
        for (int j = 1; j <= assembled.iterator_max; j++)
        {
            base = assembled.current();
            custom_attach = new RNA_data(Lib, j, GLOBAL_INPUT_INDICES_LIST[i][j]);
            rotate(base, custom_attach);
            custom_attach->update_submatrices();
            assembled.add_move(custom_attach);
        }
        if (GLOBAL_PERFORM_HAIRPIN_CHECK)
        {
            is_hairpin(assembled, WC_Lib);
        }
        assembled.update_energy();
        o_string.add_string(assembled.to_string(), assembled.get_atom_sum());
    }
}