#include "CMB_Manager.hpp"

CMB_Manager::CMB_Manager(DimerLibArray &LA)
{
    count = LA.count;

    count_per_lib = (int *)malloc(sizeof(int) * LA.count);
    attach_attempted = (bool **)malloc(sizeof(bool *) * LA.count);

    for (int i = 0; i < LA.count; i++)
    {
        count_per_lib[i] = LA[i]->count;
        attach_attempted[i] = (bool *)calloc(LA[i]->count, sizeof(bool));
    }

    last_attempted[0] = 0;
    last_attempted[1] = 0;

    strs_built = 0;
    hairpins_built = 0;
    libs_completed = (bool *)calloc(count, sizeof(bool));
}

CMB_Manager::~CMB_Manager()
{
    for (int i = 0; i < count; i++)
    {
        free(attach_attempted[i]);
    }
    free(attach_attempted);
    free(count_per_lib);
    free(libs_completed);
}

void CMB_Manager::attach_attempt(int i, int j)
{
    attach_attempted[i][j] = true;
    last_attempted[0] = i;
    last_attempted[1] = j;
}

bool CMB_Manager::is_at_end()
{
    //printf("Library = %d, Model = %d, Total structures in library = %d\n", last_attempted[0], last_attempted[1], count_per_lib[last_attempted[0]] - 1);
    if (last_attempted[1] == count_per_lib[last_attempted[0]] - 1)
    {
        //printf("Library = %d, Model = %d, Total structures in library = %d\n", last_attempted[0], last_attempted[1], count_per_lib[last_attempted[0]] - 1);
        return true;
    }
    return false;
}

void CMB_Manager::check_lib_completion()
{
    int counter = 0;

    int first_completed = last_attempted[0];

    for (int i = last_attempted[0]; i >= 0; i--)
    {
        libs_completed[i] = false;
        if (attach_attempted[i][count_per_lib[i] - 1] == true)
        {
            libs_completed[i] = true;
            counter++;
            if (i < first_completed)
                first_completed = i;
            // printf("Lib %d complete! Max = %d\n", i, count_per_lib[i] - 1);
        }
        /*else
        {
            printf("Lib %d not complete! Max = %d\n", i, count_per_lib[i] - 1);
        }*/
    }
    if (counter != (count - first_completed))
    {
        bool swap = true;
        for (int i = last_attempted[0]; i >= first_completed; i--)
        {
            if (libs_completed[i] == false)
            {
                swap = false;
            }
            libs_completed[i] = swap;
        }
    }

    if (counter == last_attempted[0] + 1)
        ;
    else
        libs_completed[0] = false;
    return;
}

void CMB_Manager::clear_attempts()
{
    for (int i = 0; i < count; i++)
    {
        if (libs_completed[i] == true)
        {
            // printf("reseting attempts for %d\n", i);
            for (int j = 0; j < count_per_lib[i]; j++)
            {
                attach_attempted[i][j] = false;
            }
        }
    }
}

int CMB_Manager::get_reset_count()
{
    int n_reset = 0;
    for (int i = 0; i < last_attempted[0] + 1; i++)
    {
        if (libs_completed[i] == true)
        {
            n_reset++;
        }
    }
    //printf("n_reset = %d\n", n_reset);
    return n_reset;
}

void CMB_Manager::successful_construction()
{
    strs_built++;
}

