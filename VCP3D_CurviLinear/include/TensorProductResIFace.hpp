// Template function to carry out the actual tensor product to accumulate the
// residual with a contribution from a face in i-direction. See
// TensorProductResIFace for the details. The template construction makes
// sure that the values of K and M are known at compile time, such that the
// compiler can optimize much better.
template<int K, int M>
void TPRIF(const int N,
           su2double *AFacex,
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
  // definition explanation in TensorProductResIFace.
  const su2double (*aty)[KP] = (const su2double (*)[KP]) ATy;
  const su2double (*atz)[KP] = (const su2double (*)[KP]) ATz;

  // Define the variables to store the intermediate results.
  su2double tmpJ[K][KP], tmpK[M][KP];

  // Outer loop over N.
  for(int l=0; l<N; ++l)
  {
    // Cast the current value of B and C to multi-dimensional arrays.
    const su2double (*b)[M] = (const su2double (*)[M]) B[l];
    su2double (*c)[K][K]    = (su2double (*)[K][K]) C[l];

    // Tensor product in k-direction to obtain the
    // solution in the K-points in k-direction.
    for(int j=0; j<M; ++j)
    {
#pragma omp simd
      for(int k=0; k<KP; ++k) tmpK[j][k] = zero;
      for(int kk=0; kk<M; ++kk)
      {
#pragma omp simd
        for(int k=0; k<KP; ++k)
          tmpK[j][k] += atz[kk][k] * b[kk][j];
      }
    }

    // Tensor product in j-direction to obtain the solution in the K-points
    // in j-direction.
    for(int k=0; k<K; ++k)
    {
#pragma omp simd
      for(int j=0; j<KP; ++j) tmpJ[k][j] = zero;
      for(int jj=0; jj<M; ++jj)
      {
#pragma omp simd
        for(int j=0; j<KP; ++j)
          tmpJ[k][j] += aty[jj][j] * tmpK[jj][k];
      }
    }

    // Tensor product in i-direction to obtain the values in
    // the K-points in i-direction. This is the final result that
    // must be added to c.
    for(int k=0; k<K; ++k)
      for(int j=0; j<K; ++j)
#pragma omp simd safelen(K)
        for(int i=0; i<K; ++i)
          c[k][j][i] += AFacex[i]*tmpJ[k][j];

  } // End of the loop over N.
}
