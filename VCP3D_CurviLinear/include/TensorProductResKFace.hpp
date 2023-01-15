// Template function to carry out the actual tensor product to accumulate the
// residual with a contribution from a face in k-direction. See
// TensorProductResKFace for the details. The template construction makes
// sure that the values of K and M are known at compile time, such that the
// compiler can optimize much better.
template<int K, int M>
void TPRKF(const int N,
           su2double *ATx,
           su2double *ATy,
           su2double *AFacez,
           su2double **B,
           su2double **C)
{
  // Determine the corresponding padded value of K.
  const int KP = ((K+vecLen1D-1)/vecLen1D)*vecLen1D;

  // Cast the one dimensional input arrays to multi-dimensional arrays.
  // Note that C++ stores multi-dimensional arrays in row major order,
  // hence the indices are reversed compared to the column major order
  // definition explanation in TensorProductResIFace.
  const su2double (*atx)[KP] = (const su2double (*)[KP]) ATx;
  const su2double (*aty)[KP] = (const su2double (*)[KP]) ATy;

  // Define the variables to store the intermediate results.
  su2double tmpI[K][KP];
  su2double tmpJ[M][KP];

  // Outer loop over N.
  for(int l=0; l<N; ++l)
  {
    // Cast the current value of B and C to multi-dimensional arrays.
    const su2double (*b)[M] = (const su2double (*)[M]) B[l];
    su2double (*c)[K][K]    = (su2double (*)[K][K]) C[l];

    // Tensor product in j-direction to obtain the
    // solution in the K-points in j-direction.
    for(int i=0; i<M; ++i)
    {
#pragma omp simd
      for(int j=0; j<KP; ++j) tmpJ[i][j] = zero;
      for(int jj=0; jj<M; ++jj)
      {
#pragma omp simd
        for(int j=0; j<KP; ++j)
          tmpJ[i][j] += aty[jj][j] * b[jj][i];
      }
    }

    // Tensor product in i-direction to obtain the solution in the K-points
    // in i-direction.
    for(int j=0; j<K; ++j)
    {
#pragma omp simd
      for(int i=0; i<KP; ++i) tmpI[j][i] = zero;
      for(int ii=0; ii<M; ++ii)
      {
#pragma omp simd
        for(int i=0; i<KP; ++i)
          tmpI[j][i] += atx[ii][i] * tmpJ[ii][j];
      }
    }

    // Tensor product in k-direction to obtain the values in
    // the K-points in k-direction. This is the final result that
    // must be added to c.
    for(int k=0; k<K; ++k)
      for(int j=0; j<K; ++j)
#pragma omp simd safelen(K)
        for(int i=0; i<K; ++i)
          c[k][j][i] += AFacez[k]*tmpI[j][i];

  } // End of the loop over N.
}
