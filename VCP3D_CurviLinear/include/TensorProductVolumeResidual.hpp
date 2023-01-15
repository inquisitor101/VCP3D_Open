// Template function to carry out the actual tensor product to accumulate the
// volume residual. See TensorProductVolumeResidual for the details.
// The template construction makes sure that the values of K and M are known
// at compile time, such that the compiler can optimize much better.
template<int K, int M>
void TPAVR(const int N,
           su2double *ATx,
           su2double *ATy,
           su2double *ATz,
           su2double **B,
           su2double **C)
{
  // Determine the corresponding padded value of K.
  const int KP = ((K+vecLen1D-1)/vecLen1D)*vecLen1D;

  // Cast the one dimensional input arrays to multi-dimensional arrays.
  // Note that C++ stores multi-dimensional arrays in row major order,
  // hence the indices are reversed compared to the column major order
  // definition explanation in TensorProductVolumeResidual.
  const su2double (*atx)[KP] = (const su2double (*)[KP]) ATx;
  const su2double (*aty)[KP] = (const su2double (*)[KP]) ATy;
  const su2double (*atz)[KP] = (const su2double (*)[KP]) ATz;

  // Outer loop over N.
  for(int l=0; l<N; ++l)
  {
    // Cast the current value of the input arrays to multi-dimensional
    // arrays.
    const su2double (*b)[M][M] = (const su2double (*)[M][M]) B[l];
    su2double (*c)[K][K]       = (su2double (*)[K][K]) C[l];

    // Define the variables to store the intermediate results.
    su2double tmpK[M][M][KP];
    su2double tmpJ[K][M][KP];
    su2double tmpI[K][K][KP];

    // Tensor product in k-direction to obtain the
    // solution in the K points in k-direction.
    for(int i=0; i<M; ++i)
    {
      for(int j=0; j<M; ++j)
      {
#pragma omp simd
        for(int k=0; k<KP; ++k) tmpK[i][j][k] = zero;
        for(int kk=0; kk<M; ++kk)
        {
#pragma omp simd
          for(int k=0; k<KP; ++k)
            tmpK[i][j][k] += atz[kk][k] * b[kk][j][i];
        }
      }
    }

    // Tensor product in j-direction to obtain the
    // solution in the K-points in j-direction.
    for(int k=0; k<K; ++k)
    {
      for(int i=0; i<M; ++i)
      {
#pragma omp simd
        for(int j=0; j<KP; ++j) tmpJ[k][i][j] = zero;
        for(int jj=0; jj<M; ++jj)
        {
#pragma omp simd
          for(int j=0; j<KP; ++j)
            tmpJ[k][i][j] += aty[jj][j] * tmpK[i][jj][k];
        }
      }
    }

    // Tensor product in i-direction to obtain the
    // solution in the K-points in i-direction.
    for(int k=0; k<K; ++k)
    {
      for(int j=0; j<K; ++j)
      {
#pragma omp simd
        for(int i=0; i<KP; ++i) tmpI[k][j][i] = zero;
        for(int ii=0; ii<M; ++ii)
        {
#pragma omp simd
          for(int i=0; i<KP; ++i)
            tmpI[k][j][i] += atx[ii][i] * tmpJ[k][ii][j];
        }
      }
    }

    // Add tmpI to c.
    for(int k=0; k<K; ++k)
      for(int j=0; j<K; ++j)
#pragma omp simd safelen(K)
        for(int i=0; i<K; ++i)
          c[k][j][i] += tmpI[k][j][i];

  } // End of the loop over N.
}
