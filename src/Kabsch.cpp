#include "Kabsch.hpp"


gsl_block  *KABSCH_MEMBLOCK;
gsl_matrix *H_KABSCH;   
gsl_matrix *V_KABSCH; 
gsl_matrix *VUt_KABSCH; 
gsl_matrix *DIA_KABSCH;
gsl_matrix *TEMP_KABSCH;
gsl_matrix *R_KABSCH; //Rotation Matrix, values are updated after kabsch_get_rotation_matrix_generic_fast(...)
gsl_matrix *P_KABSCH;
gsl_matrix *PWORK_KABSCH;
gsl_matrix *Q_KABSCH;
gsl_vector *S_KABSCH;
gsl_vector *WORK_KABSCH;
uint32_t   P_OFFSET_KABSCH;
uint32_t   PWORK_OFFSET_KABSCH;
uint32_t   Q_OFFSET_KABSCH;


double get_determinant_3x3fast(gsl_matrix *A)
{
  double det1,det2,det3; //Using 3 variables to try to use automatic parallelization (vectorization) by CPU (using -march=native flag in GCC)
  det1 = gsl_matrix_get(A, 0, 0) * ((gsl_matrix_get(A, 1, 1) * gsl_matrix_get(A, 2, 2)) - (gsl_matrix_get(A, 1, 2) * gsl_matrix_get(A, 2, 1)));
  det2 = gsl_matrix_get(A, 0, 1) * ((gsl_matrix_get(A, 1, 0) * gsl_matrix_get(A, 2, 2)) - (gsl_matrix_get(A, 1, 2) * gsl_matrix_get(A, 2, 0)));
  det3 = gsl_matrix_get(A, 0, 2) * ((gsl_matrix_get(A, 1, 0) * gsl_matrix_get(A, 2, 1)) - (gsl_matrix_get(A, 1, 1) * gsl_matrix_get(A, 2, 0)));
  return det1 - det2 + det3;
}

double get_determinant(gsl_matrix *A, bool inPlace)
{
  double det;
  int signum = 0;
  gsl_permutation *p = gsl_permutation_alloc(A->size1);
  gsl_matrix *tmpA = nullptr;

  if (inPlace)
    tmpA = A;
  else
  {
    gsl_matrix *tmpA = gsl_matrix_alloc(A->size1, A->size2);
    gsl_matrix_memcpy(tmpA, A);
  }

  gsl_linalg_LU_decomp(tmpA, p, &signum);
  det = gsl_linalg_LU_det(tmpA, signum);
  gsl_permutation_free(p);
  if (!inPlace)
  {
    gsl_matrix_free(tmpA);
  }

  return det;
}

void kabsch_create(size_t M, size_t N)
{
  size_t real_mem_size = (6 * MATRIX_DIMENSION2 * MATRIX_DIMENSION2) + (2 * MATRIX_DIMENSION2) + (N * M * 3);//6 3x3 Matrices + 2 3-long vectors + additional needed
  constexpr size_t tda = MATRIX_DIMENSION2;
  KABSCH_MEMBLOCK = gsl_block_alloc(real_mem_size);
  H_KABSCH    = gsl_matrix_alloc_from_block(KABSCH_MEMBLOCK,  0, MATRIX_DIMENSION2, MATRIX_DIMENSION2, tda);   
  V_KABSCH    = gsl_matrix_alloc_from_block(KABSCH_MEMBLOCK,  9, MATRIX_DIMENSION2, MATRIX_DIMENSION2, tda); 
  VUt_KABSCH  = gsl_matrix_alloc_from_block(KABSCH_MEMBLOCK, 18, MATRIX_DIMENSION2, MATRIX_DIMENSION2, tda); 
  DIA_KABSCH  = gsl_matrix_alloc_from_block(KABSCH_MEMBLOCK, 27, MATRIX_DIMENSION2, MATRIX_DIMENSION2, tda);
  TEMP_KABSCH = gsl_matrix_alloc_from_block(KABSCH_MEMBLOCK, 36, MATRIX_DIMENSION2, MATRIX_DIMENSION2, tda);
  R_KABSCH    = gsl_matrix_alloc_from_block(KABSCH_MEMBLOCK, 45, MATRIX_DIMENSION2, MATRIX_DIMENSION2, tda);
  S_KABSCH    = gsl_vector_alloc_from_block(KABSCH_MEMBLOCK, 54, MATRIX_DIMENSION2, 1);
  WORK_KABSCH = gsl_vector_alloc_from_block(KABSCH_MEMBLOCK, 57, MATRIX_DIMENSION2, 1);
  P_OFFSET_KABSCH = 60;
  PWORK_OFFSET_KABSCH = (P_OFFSET_KABSCH + (N * M));
  Q_OFFSET_KABSCH = (PWORK_OFFSET_KABSCH + (N * M));
  P_KABSCH = gsl_matrix_alloc_from_block(KABSCH_MEMBLOCK,  P_OFFSET_KABSCH, M, N, tda);
  PWORK_KABSCH = gsl_matrix_alloc_from_block(KABSCH_MEMBLOCK,  PWORK_OFFSET_KABSCH, N, M, M);
  Q_KABSCH = gsl_matrix_alloc_from_block(KABSCH_MEMBLOCK,  Q_OFFSET_KABSCH, M, N, tda);
  gsl_matrix_set_identity(DIA_KABSCH);
}

void kabsch_prepare_matrix_internal(size_t M, size_t N, uint16_t *Rows, gsl_matrix *source, gsl_matrix *dest)
{
    dest->size1 = M;
    dest->size2 = N;
    dest->tda   = N;
    for(size_t i = 0; i < M; i++)
    {
      gsl_matrix_row_copy(dest, i, source, Rows[i]);
    }
}

template <kabsch_matrix M_REQUEST> gsl_matrix* kabsch_prepare_matrix(size_t M, size_t N, uint16_t *Rows, gsl_matrix *source)
{
  if constexpr (M_REQUEST == kabsch_matrix::KABSCH_MATRIX_P)
  {
    kabsch_prepare_matrix_internal(M, N, Rows, source, P_KABSCH);
    return P_KABSCH;
  }
  else if constexpr (M_REQUEST == kabsch_matrix::KABSCH_MATRIX_Q)
  {
    kabsch_prepare_matrix_internal(M, N, Rows, source, Q_KABSCH);
    return Q_KABSCH;
  }
}
template gsl_matrix* kabsch_prepare_matrix<KABSCH_MATRIX_P>(size_t M, size_t N, uint16_t *Rows, gsl_matrix *source);
template gsl_matrix* kabsch_prepare_matrix<KABSCH_MATRIX_Q>(size_t M, size_t N, uint16_t *Rows, gsl_matrix *source);

gsl_matrix* kabsch_get_work_matrix(size_t M, size_t N)
{
  PWORK_KABSCH->size1 = M;
  PWORK_KABSCH->size2 = N;
  PWORK_KABSCH->tda   = N;
  return PWORK_KABSCH;
}

void kabsch_destroy()
{
    gsl_matrix_free(H_KABSCH);
    gsl_matrix_free(V_KABSCH);
    gsl_matrix_free(VUt_KABSCH);
    gsl_matrix_free(DIA_KABSCH);
    gsl_matrix_free(TEMP_KABSCH);
    gsl_matrix_free(R_KABSCH);
    gsl_matrix_free(P_KABSCH);
    gsl_matrix_free(PWORK_KABSCH);
    gsl_matrix_free(Q_KABSCH);
    gsl_vector_free(S_KABSCH);
    gsl_vector_free(WORK_KABSCH);
    gsl_block_free(KABSCH_MEMBLOCK);
}

gsl_matrix* kabsch_allocate_work_matrix(gsl_matrix *P)
{
  return gsl_matrix_alloc(P->size2, P->size1);
}

gsl_matrix* kabsch_get_rotation_matrix()
{
  return R_KABSCH;
}

void kabsch_calculate_rotation_matrix_Nx3fast(gsl_matrix *P, gsl_matrix *Q, gsl_matrix *P_WORK, double *__restrict__ COMP, double *__restrict__ COMQ)
{
  double Determinant;
  get_matrix_COM(P, COMP);
  center_matrix(P, COMP);
  get_matrix_COM(Q, COMQ);
  center_matrix(Q, COMQ);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, P, Q, 0.0, H_KABSCH);
  gsl_linalg_SV_decomp(H_KABSCH, V_KABSCH, S_KABSCH, WORK_KABSCH);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, V_KABSCH, H_KABSCH, 0.0, VUt_KABSCH);
  Determinant = get_determinant_3x3fast(VUt_KABSCH);
  gsl_matrix_set(DIA_KABSCH, MATRIX_DIMENSION2 - 1, MATRIX_DIMENSION2 - 1, Determinant);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, DIA_KABSCH, H_KABSCH, 0.0, TEMP_KABSCH);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, V_KABSCH, TEMP_KABSCH, 0.0, R_KABSCH);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, R_KABSCH, P, 0.0, P_WORK);
  gsl_matrix_transpose_memcpy(P, P_WORK);
}


/* To be deprecrated. Use kabsch_calculate_rotation_matrix_Nx3fast instead*/
gsl_matrix *kabsch_get_rotation_matrix_generic(gsl_matrix *P, gsl_matrix *Q, double *__restrict__ COMP, double *__restrict__ COMQ)
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

  get_matrix_COM(P, COMP);
  center_matrix(P, COMP);
  get_matrix_COM(Q, COMQ);
  center_matrix(Q, COMQ);

  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, P, Q, 0.0, H);

  gsl_linalg_SV_decomp(H, V, S, work);
  // gsl_linalg_SV_decomp_jacobi(H, V, S);

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, V, H, 0.0, VUt);

  double det = get_determinant(VUt, true);

  gsl_matrix_set_identity(DIA);
  gsl_matrix_set(DIA, dimension2 - 1, dimension2 - 1, det);

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