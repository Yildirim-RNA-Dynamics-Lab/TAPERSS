#include "Atom_info.hpp"

/* Gets atom ID */
atom_id get_atom_id(char *name)
{
    atom_id atom = NOTATOM;
    switch (name[0])
    {
    case 'N':
        switch (name[1])
        {
        case '1':
            atom = N1;
            break;
        case '2':
            atom = N2;
            break;
        case '3':
            atom = N3;
            break;
        case '4':
            atom = N4;
            break;
        case '6':
            atom = N6;
            break;
        case '7':
            atom = N7;
            break;
        case '9':
            atom = N9;
            break;
        }
        break;
    case 'O':
        switch (name[1])
        {
        case 'P':
            if (name[2] == '1')
                atom = OP1;
            else if (name[2] == '2')
                atom = OP2;
            break;
        case '1':
            atom = O1P;
            break;
        case '2':
            if (name[2] == '\'')
                atom = O2p;
            else if (name[2] == 'P')
                atom = O2P;
            else
                atom = O2;
            break;
        case '3':
            atom = O3p;
            break;
        case '4':
            if (name[2] == '\'')
                atom = O4p;
            else
                atom = O4;
            break;
        case '5':
            atom = O5p;
            break;
        case '6':
            atom = O6;
            break;
        }
        break;
    case 'C':
        switch (name[1])
        {
        case '1':
            atom = C1p;
            break;
        case '2':
            if (name[2] == '\'')
                atom = C2p;
            else
                atom = C2;
            break;
        case '3':
            atom = C3p;
            break;
        case '4':
            if (name[2] == '\'')
                atom = C4p;
            else
                atom = C4;
            break;
        case '5':
            if (name[2] == '\'')
                atom = C5p;
            else
                atom = C5;
            break;
        case '6':
            atom = C6;
            break;
        case '8':
            atom = C8;
            break;
        }
        break;
    case 'P':
        atom = P;
        break;
    default:
        printf("Unknown atom: %s\n", name);
        exit(2);
    }
    return atom;
}

atom_charge get_atom_charge(atom_id name, char res)
{
    atom_charge charge = NEUTRAL;
    switch (name)
    {
    case N1:
        switch (res)
        {
        case 'A':
            charge = NEGATIVE;
            break;
        case 'G':
            charge = POSITIVE;
            break;
        }
        break;
    case N2:
        charge = POSITIVE;
        break;
    case N3:
        switch (res)
        {
        case 'U':
            charge = POSITIVE;
            break;
        default:
            charge = NEGATIVE;
            break;
        }
        break;
    case N4:
        charge = POSITIVE;
        break;
    case N6:
        charge = POSITIVE;
        break;
    case N7:
        charge = NEGATIVE;
        break;
    case O2:
        charge = NEGATIVE;
        break;
    case O4:
        charge = NEGATIVE;
        break;
    case O6:
        charge = NEGATIVE;
        break;
    default:
        charge = NEUTRAL;
        break;
    }
    return charge;
}

atom_info::atom_info()
{
    iterator = 0;
}

atom_info::atom_info(int n)
{
    name = (char **)malloc(sizeof(char *) * n);
    index = (int *)malloc(sizeof(int) * n);
    residue = (char *)malloc(sizeof(char) * n);
    dnt_pos = (int *)malloc(sizeof(int) * n);
    atom_ids = (atom_id *)malloc(sizeof(atom_id) * n);
    charges = (atom_charge *)malloc(sizeof(atom_charge) * n);
    count = n;
    for (int i = 0; i < n; i++)
    {
        name[i] = (char *)malloc(sizeof(char) * 5);
    }
    iterator = 0;
}

atom_info::atom_info(char **n, int *in, char *r, int *d, int c)
{
    name = (char **)malloc(sizeof(char *) * c);
    for (int i = 0; i < c; i++)
    {
        name[i] = (char *)malloc(sizeof(char) * 5);
        strcpy(name[i], n[i]);
    }

    index = (int *)malloc(sizeof(int) * c);
    residue = (char *)malloc(sizeof(char) * c);
    dnt_pos = (int *)malloc(sizeof(int) * c);

    memcpy(index, in, sizeof(int) * c);
    memcpy(residue, r, sizeof(char) * c);
    memcpy(dnt_pos, d, sizeof(int) * c);

    count = c;
}

void atom_info::initialize_atom_info(int n)
{
    name = (char **)malloc(sizeof(char *) * n);
    index = (int *)malloc(sizeof(int) * n);
    residue = (char *)malloc(sizeof(char) * n);
    dnt_pos = (int *)malloc(sizeof(int) * n);
    count = n;
    for (int i = 0; i < n; i++)
    {
        name[i] = (char *)malloc(sizeof(char) * 5);
    }
    iterator = 0;
}

atom_info::~atom_info()
{
    // printf("Destructor Called!\n");
    for (int i = 0; i < count; i++)
    {
        free(name[i]);
    }
    free(index);
    free(residue);
    free(dnt_pos);
    free(atom_ids);
    free(charges);
    free(name);
}

void atom_info::add_atom(char *N, int i, char r, int p)
{
    strcpy(name[iterator], N);
    index[iterator] = i;
    residue[iterator] = r;
    dnt_pos[iterator] = p;
    atom_ids[iterator] = get_atom_id(N);
    charges[iterator] = get_atom_charge(atom_ids[iterator], r);
    iterator++;
}

void atom_info::clear()
{
    iterator = 0;
}

atom_info::atom_info(const atom_info &A)
{
    name = (char **)malloc(sizeof(char *) * A.count);
    for (int i = 0; i < A.count; i++)
    {
        name[i] = (char *)malloc(sizeof(char) * 5);
        strcpy(name[i], A.name[i]);
    }

    index = (int *)malloc(sizeof(int) * A.count);
    residue = (char *)malloc(sizeof(char) * A.count);
    dnt_pos = (int *)malloc(sizeof(int) * A.count);

    memcpy(index, A.index, sizeof(int) * A.count);
    memcpy(residue, A.residue, sizeof(char) * A.count);
    memcpy(dnt_pos, A.dnt_pos, sizeof(int) * A.count);

    count = A.count;
}

void atom_info::move(atom_info &A)
{
    name = A.name;
    index = A.index;
    residue = A.residue;
    dnt_pos = A.dnt_pos;
    count = A.count;
}

void atom_info::print_at(int line)
{
    printf("%-6s%5d %4s %3c  %4d    ", "ATOM", index[line], name[line], residue[line], dnt_pos[line]);
}

