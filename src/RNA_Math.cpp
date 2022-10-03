#include "RNA_Math.hpp"

double rmsd_generic(gsl_matrix *A, gsl_matrix *B)
{
    double tempx, tempy, tempz;
    double temprmsd = 0.0;
    int atom_count = 0;

    /*if (A->size1 != B->size1)
    {
        exit(1);
    }*/

    for (unsigned int i = 0; i < A->size1; i++)
    {
        tempx = square((gsl_matrix_get(A, i, 0) - gsl_matrix_get(B, i, 0)));
        tempy = square((gsl_matrix_get(A, i, 1) - gsl_matrix_get(B, i, 1)));
        tempz = square((gsl_matrix_get(A, i, 2) - gsl_matrix_get(B, i, 2)));
        temprmsd += tempx + tempy + tempz;
        atom_count++;
    }
    return sqrt(temprmsd / (double)atom_count);
}

void get_matrix_COM(gsl_matrix *M, double *COM)
{
    for (unsigned int i = 0; i < M->size1; i++)
    {
        COM[0] += gsl_matrix_get(M, i, 0);
        COM[1] += gsl_matrix_get(M, i, 1);
        COM[2] += gsl_matrix_get(M, i, 2);
    }

    COM[0] /= M->size1;
    COM[1] /= M->size1;
    COM[2] /= M->size1;
}

void center_matrix(gsl_matrix *M, double *COM)
{
    for (unsigned int i = 0; i < M->size1; i++)
    {
        gsl_matrix_set(M, i, 0, (gsl_matrix_get(M, i, 0) - COM[0]));
        gsl_matrix_set(M, i, 1, (gsl_matrix_get(M, i, 1) - COM[1]));
        gsl_matrix_set(M, i, 2, (gsl_matrix_get(M, i, 2) - COM[2]));
    }
}

void translate_matrix(double *__restrict__ tV, gsl_matrix *M, double scalar)
{
    for (unsigned int i = 0; i < M->size1; i++)
    {
        gsl_matrix_set(M, i, 0, (gsl_matrix_get(M, i, 0) + (tV[0] * scalar)));
        gsl_matrix_set(M, i, 1, (gsl_matrix_get(M, i, 1) + (tV[1] * scalar)));
        gsl_matrix_set(M, i, 2, (gsl_matrix_get(M, i, 2) + (tV[2] * scalar)));
    }
}

void apply_rotation_matrix(gsl_matrix *R, gsl_matrix *M)
{
    gsl_matrix *M_TEMP = kabsch_get_work_matrix(M->size2, M->size1);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, R, M, 0.0, M_TEMP);
    gsl_matrix_transpose_memcpy(M, M_TEMP);
}

void gsl_matrix_print(gsl_matrix *M, const char* name)
{
    printf("---------------Matrix name: %s\n", name);
    for (unsigned int i = 0; i < M->size1; i++)
    {
        for (unsigned int j = 0; j < M->size2; j++)
            printf("%3.3f ", gsl_matrix_get(M, i, j));
        putchar('\n');
    }
}

void gsl_matrix_print_row(gsl_matrix *M, size_t i)
{
    for (unsigned int j = 0; j < M->size2; j++)
        printf("%3.3f ", gsl_matrix_get(M, i, j));
    putchar('\n');

}

int array_min_idx_for_energy(int *Arr, int N)
{
    int MinVal = N;
    int MinIdx = -1;
    for(int i = 0; i < N; i++)
    {
        if(Arr[i] < MinVal && Arr[i] != 0)
        {
            MinVal = Arr[i];
            MinIdx = i;
        }
    }
    return MinIdx;
}

void gsl_vector_print(gsl_vector* V)
{
    printf("%f %f %f\n", gsl_vector_get(V, 0), gsl_vector_get(V, 1), gsl_vector_get(V, 2));
}