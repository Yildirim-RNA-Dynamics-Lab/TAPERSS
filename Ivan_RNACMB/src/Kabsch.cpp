#include "Kabsch.hpp"

double get_determinant(gsl_matrix *A, bool inPlace) 
{
    double det;
    int signum = 0;
    gsl_permutation *p = gsl_permutation_alloc(A->size1);
    gsl_matrix *tmpA = nullptr;

    if (inPlace)
       tmpA = A;
    else {
      gsl_matrix *tmpA = gsl_matrix_alloc(A->size1, A->size2);
      gsl_matrix_memcpy(tmpA , A);
    }

    gsl_linalg_LU_decomp(tmpA , p , &signum);
    det = gsl_linalg_LU_det(tmpA , signum);
    gsl_permutation_free(p);
    if (!inPlace)
    {
       gsl_matrix_free(tmpA);
    }
    
    return det;
}

gsl_matrix *kabsch_get_rotation_matrix_generic(gsl_matrix *P, gsl_matrix *Q, double * __restrict__ COMP, double * __restrict__ COMQ)
{
    gsl_matrix *P_TEMP = gsl_matrix_alloc(P->size2, P->size1);
    const int dimension2 = P->size2;

    gsl_matrix *H = gsl_matrix_alloc(dimension2, dimension2);
    gsl_matrix *V = gsl_matrix_alloc(dimension2, dimension2);

    gsl_vector *S = gsl_vector_alloc(dimension2);     
    gsl_vector *work = gsl_vector_alloc(dimension2);     

    gsl_matrix *VUt = gsl_matrix_alloc(dimension2, dimension2);
    gsl_matrix *DIA = gsl_matrix_alloc(dimension2, dimension2);

    gsl_matrix *TEMP = gsl_matrix_alloc(dimension2, dimension2);
    gsl_matrix *R = gsl_matrix_alloc(dimension2, dimension2);

    center_matrix(P, COMP); 
    center_matrix(Q, COMQ);

    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, P, Q, 0.0, H);    

    gsl_linalg_SV_decomp(H, V, S, work);
    //gsl_linalg_SV_decomp_jacobi(H, V, S);

    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, V, H, 0.0, VUt);

    double det = get_determinant(VUt, true);
    
    gsl_matrix_set_identity( DIA );
    gsl_matrix_set (DIA, dimension2 - 1, dimension2 - 1, det);

    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, DIA, H, 0.0, TEMP);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V, TEMP, 0.0, R);    
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, R, P, 0.0, P_TEMP);            
    
    gsl_matrix_transpose_memcpy(P, P_TEMP);

    gsl_matrix_free(H);
    gsl_matrix_free(V);

    gsl_vector_free(S);     
    gsl_vector_free(work);     

    gsl_matrix_free(VUt);
    gsl_matrix_free(DIA);

    gsl_matrix_free(TEMP);
    gsl_matrix_free(P_TEMP);
    return R;
}