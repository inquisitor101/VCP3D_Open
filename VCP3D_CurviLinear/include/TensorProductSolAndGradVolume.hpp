// Template function to carry out the actual tensor product to compute the
// solution and its gradients in the volume. See TensorProductSolAndGradVolume
// for the details. The template construction makes sure that the values of K and
// M are known at compile time, such that the compiler can optimize much better.
template<int K, int M>
void TPSAGV(const int N,
            su2double *A,
            su2double *ADer,
            su2double **B,
            su2double **C,
            su2double **CDerX,
            su2double **CDerY,
            su2double **CDerZ)
{
  // Determine the corresponding padded value of M.
  const int MP = ((M+vecLen1D-1)/vecLen1D)*vecLen1D;

  // Cast the one dimensional input arrays to multi-dimensional arrays.
  // Note that C++ stores multi-dimensional arrays in row major order,
  // hence the indices are reversed compared to the column major order
  // definition explanation in TensorProductSolAndGradVolume.
  const su2double (*a)[MP]    = (const su2double (*)[MP]) A;
  const su2double (*aDer)[MP] = (const su2double (*)[MP]) ADer;

  // Outer loop over N.
  for(int l=0; l<N; ++l)
  {
    // Cast the current value of B to a multi-dimensional array.
    const su2double (*b)[K][K] = (const su2double (*)[K][K]) B[l];

    // Define the variables to store the intermediate results.
    su2double tmpK[K][K][MP];
    su2double tmpJ[M][K][MP];
    su2double tmpI[M][M][MP];

    // Tensor product in k-direction to obtain the
    // solution in the M points in k-direction.
    for(int i=0; i<K; ++i)
    {
      for(int j=0; j<K; ++j)
      {
#pragma omp simd
        for(int k=0; k<MP; ++k) tmpK[i][j][k] = zero;
        for(int kk=0; kk<K; ++kk)
        {
#pragma omp simd
          for(int k=0; k<MP; ++k)
            tmpK[i][j][k] += a[kk][k] * b[kk][j][i];
        }
      }
    }

    // Tensor product in j-direction to obtain the
    // solution in the M-points in j-direction.
    for(int k=0; k<M; ++k)
    {
      for(int i=0; i<K; ++i)
      {
#pragma omp simd
        for(int j=0; j<MP; ++j) tmpJ[k][i][j] = zero;
        for(int jj=0; jj<K; ++jj)
        {
#pragma omp simd
          for(int j=0; j<MP; ++j)
            tmpJ[k][i][j] += a[jj][j] * tmpK[i][jj][k];
        }
      }
    }

    // Check if C must be computed.
    if( C )
    {
      // Cast the current value of C to a multi-dimensional array.
      su2double (*c)[M][M] = (su2double (*)[M][M]) C[l];

      // Tensor product in i-direction to obtain the
      // solution in the M-points in i-direction.
      // This is the final result for the solution.
      for(int k=0; k<M; ++k)
      {
        for(int j=0; j<M; ++j)
        {
#pragma omp simd
          for(int i=0; i<MP; ++i) tmpI[k][j][i] = zero;
          for(int ii=0; ii<K; ++ii)
          {
#pragma omp simd
            for(int i=0; i<MP; ++i)
              tmpI[k][j][i] += a[ii][i] * tmpJ[k][ii][j];
          }
        }
      }

      // Copy the values to the appropriate location in c.
      for(int k=0; k<M; ++k)
        for(int j=0; j<M; ++j)
          for(int i=0; i<M; ++i)
            c[k][j][i] = tmpI[k][j][i];
    }

    // Check if CDerX must be computed.
    if( CDerX )
    {
      // Cast the current value of CDerX to a multi-dimensional array.
      su2double (*cDerX)[M][M] = (su2double (*)[M][M]) CDerX[l];

      // Tensor product in i-direction to obtain the x-derivative
      // of the solution in the M-points in i-direction. This is
      // the final result for the z-derivative. Note that tmpJ can
      // be reused, which saves two tensor products compared to a
      // full re-evaluation of this derivative.
      for(int k=0; k<M; ++k)
      {
        for(int j=0; j<M; ++j)
        {
#pragma omp simd
          for(int i=0; i<MP; ++i) tmpI[k][j][i] = zero;
          for(int ii=0; ii<K; ++ii)
          {
#pragma omp simd
            for(int i=0; i<MP; ++i)
              tmpI[k][j][i] += aDer[ii][i] * tmpJ[k][ii][j];
          }
        }
      }

      // Copy the values to the appropriate location in cDerX.
      for(int k=0; k<M; ++k)
        for(int j=0; j<M; ++j)
          for(int i=0; i<M; ++i)
            cDerX[k][j][i] = tmpI[k][j][i];
    }

    // Check if CDerY must be computed.
    if( CDerY )
    {
      // Cast the current value of CDerY to a multi-dimensional array.
      su2double (*cDerY)[M][M] = (su2double (*)[M][M]) CDerY[l];

      // Tensor product in j-direction to obtain the y-derivative
      // of the solution in the M-points in j-direction. Note that
      // tmpK can be reused, which saves one tensor product
      // compared to a full re-evaluation of this derivative.
      for(int k=0; k<M; ++k)
      {
        for(int i=0; i<K; ++i)
        {
#pragma omp simd
          for(int j=0; j<MP; ++j) tmpJ[k][i][j] = zero;
          for(int jj=0; jj<K; ++jj)
          {
#pragma omp simd
            for(int j=0; j<MP; ++j)
              tmpJ[k][i][j] += aDer[jj][j] * tmpK[i][jj][k];
          }
        }
      }

      // Tensor product in i-direction to obtain the
      // y-derivative of the solution in the M-points
      // in i-direction. This is the final result.
      for(int k=0; k<M; ++k)
      {
        for(int j=0; j<M; ++j)
        {
#pragma omp simd
          for(int i=0; i<MP; ++i) tmpI[k][j][i] = zero;
          for(int ii=0; ii<K; ++ii)
          {
#pragma omp simd
            for(int i=0; i<MP; ++i)
              tmpI[k][j][i] += a[ii][i] * tmpJ[k][ii][j];
          }
        }
      }

      // Copy the values to the appropriate location in cDerY.
      for(int k=0; k<M; ++k)
        for(int j=0; j<M; ++j)
          for(int i=0; i<M; ++i)
            cDerY[k][j][i] = tmpI[k][j][i];
    }

    // Check if CDerZ must be computed.
    if( CDerZ )
    {
      // Cast the current value of CDerZ to a multi-dimensional array.
      su2double (*cDerZ)[M][M] = (su2double (*)[M][M]) CDerZ[l];

      // Tensor product in k-direction to obtain the z-derivative
      // of the solution in the M-points in k-direction.
     for(int i=0; i<K; ++i)
      {
        for(int j=0; j<K; ++j)
        {
#pragma omp simd
          for(int k=0; k<MP; ++k) tmpK[i][j][k] = zero;
          for(int kk=0; kk<K; ++kk)
          {
#pragma omp simd
            for(int k=0; k<MP; ++k)
              tmpK[i][j][k] += aDer[kk][k] * b[kk][j][i];
          }
        }
      }

      // Tensor product in j-direction to obtain the z-derivative
      // of the solution in the M-points in j-direction.
      for(int k=0; k<M; ++k)
      {
        for(int i=0; i<K; ++i)
        {
#pragma omp simd
          for(int j=0; j<MP; ++j) tmpJ[k][i][j] = zero;
          for(int jj=0; jj<K; ++jj)
          {
#pragma omp simd
            for(int j=0; j<MP; ++j)
              tmpJ[k][i][j] += a[jj][j] * tmpK[i][jj][k];
          }
        }
      }

      // Tensor product in i-direction to obtain the
      // z-derivative of the solution in the M-points
      // in i-direction. This is the final result.
      for(int k=0; k<M; ++k)
      {
        for(int j=0; j<M; ++j)
        {
#pragma omp simd
          for(int i=0; i<MP; ++i) tmpI[k][j][i] = zero;
          for(int ii=0; ii<K; ++ii)
          {
#pragma omp simd
            for(int i=0; i<MP; ++i)
              tmpI[k][j][i] += a[ii][i] * tmpJ[k][ii][j];
          }
        }
      }

      // Copy the values to the appropriate location in cDerZ.
      for(int k=0; k<M; ++k)
        for(int j=0; j<M; ++j)
          for(int i=0; i<M; ++i)
            cDerZ[k][j][i] = tmpI[k][j][i];
    }
  } // End of the loop over N.
}
