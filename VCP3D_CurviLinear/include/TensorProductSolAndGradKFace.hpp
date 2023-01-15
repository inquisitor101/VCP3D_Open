// Template function to carry out the actual tensor product to compute the
// solution and its gradients on a face in k-direction. See
// TensorProductSolAndGradKFace for the details. The template construction makes
// sure that the values of K and M are known at compile time, such that the
// compiler can optimize much better.
template<int K, int M>
void TPSAGK(const int N,
            su2double *A,
            su2double *ADer,
            su2double *AFace,
            su2double *ADerFace,
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

  // Define the variables to store the intermediate results.
  su2double tmpI[K][MP], tmpJ[M][MP];
  su2double bk[K][K], bkDer[K][K];

  // Outer loop over N.
  for(int l=0; l<N; ++l)
  {
    // Cast the current value of B, C and CDer to multi-dimensional arrays.
    const su2double (*b)[K][K] = (const su2double (*)[K][K]) B[l];

		// Compute the tensor product in k-direction for both AFace
    // and ADerFace. The tensor products in the k-direction are
    // different, because both AFace and ADerFace are rank 1 tensors.
    for(int j=0; j<K; ++j)
    {
#pragma omp simd safelen(K)
      for(int i=0; i<K; ++i)
      {
        bk[j][i]    = zero;
        bkDer[j][i] = zero;
      }
    }

    for(int k=0; k<K; ++k)
    {
      for(int j=0; j<K; ++j)
      {
#pragma omp simd safelen(K)
        for(int i=0; i<K; ++i)
        {
          bk[j][i]    += AFace[k]   *b[k][j][i];
          bkDer[j][i] += ADerFace[k]*b[k][j][i];
        }
      }
    }

    // Tensor product in i-direction to obtain the
    // solution in the M-points in i-direction.
    for(int j=0; j<K; ++j)
    {
#pragma omp simd
      for(int i=0; i<MP; ++i) tmpI[j][i] = zero;
      for(int ii=0; ii<K; ++ii)
      {
#pragma omp simd
        for(int i=0; i<MP; ++i)
          tmpI[j][i] += a[ii][i] * bk[j][ii];
      }
    }

		// Check if C must be computed.
		if( C )
		{
			// Cast the current value of C to a multi-dimensional array.
			su2double (*c)[M] = (su2double (*)[M]) C[l];

    	// Tensor product in j-direction to obtain the
    	// solution in the M-points in j-direction.
    	// This is the final result for the solution.
    	for(int i=0; i<M; ++i)
    	{
#pragma omp simd
      	for(int j=0; j<MP; ++j) tmpJ[i][j] = zero;
      	for(int jj=0; jj<K; ++jj)
      	{
#pragma omp simd
    	    for(int j=0; j<MP; ++j)
    	      tmpJ[i][j] += a[jj][j] * tmpI[jj][i];
    	  }
    	}

    	// Copy the values to the appropriate location in c.
    	for(int j=0; j<M; ++j)
    	  for(int i=0; i<M; ++i)
    	    c[j][i] = tmpJ[i][j];
		}

		// Check if CDerY must be computed.
		if( CDerY )
		{
			// Cast the current value of CDerY to a multi-dimensional array.
			su2double (*cDerY)[M] = (su2double (*)[M]) CDerY[l];

    	// The values of tmpI can be reused to compute the derivatives of
    	// the solution in j-direction. Carry out the tensor product to
    	// compute these values.
    	for(int i=0; i<M; ++i)
    	{
#pragma omp simd
      	for(int j=0; j<MP; ++j) tmpJ[i][j] = zero;
      	for(int jj=0; jj<K; ++jj)
      	{
#pragma omp simd
    	    for(int j=0; j<MP; ++j)
    	      tmpJ[i][j] += aDer[jj][j] * tmpI[jj][i];
    	  }
    	}

    	// Copy the values to the appropriate location in cDerY.
    	for(int j=0; j<M; ++j)
    	  for(int i=0; i<M; ++i)
    	    cDerY[j][i] = tmpJ[i][j];
		}

		// Check if CDerX must be computed.
		if( CDerX )
		{
			// Cast the current value of CDerX to a multi-dimensional array.
			su2double (*cDerX)[M] = (su2double (*)[M]) CDerX[l];

    	// Tensor product in i-direction to obtain the x-derivatives
    	// of the solution in the M-points in i-direction.
    	for(int j=0; j<K; ++j)
    	{
#pragma omp simd
      	for(int i=0; i<MP; ++i) tmpI[j][i] = zero;
      	for(int ii=0; ii<K; ++ii)
      	{
#pragma omp simd
    	    for(int i=0; i<MP; ++i)
    	      tmpI[j][i] += aDer[ii][i] * bk[j][ii];
    	  }
    	}

    	// Tensor product in j-direction to obtain the x-derivatives
    	// of the solution in the M-points in j-direction.
    	// This is the final result for the x-derivatives of the solution.
    	for(int i=0; i<M; ++i)
    	{
#pragma omp simd
      	for(int j=0; j<MP; ++j) tmpJ[i][j] = zero;
      	for(int jj=0; jj<K; ++jj)
      	{
#pragma omp simd
    	    for(int j=0; j<MP; ++j)
    	      tmpJ[i][j] += a[jj][j] * tmpI[jj][i];
    	  }
    	}

    	// Copy the values to the appropriate location in cDerX.
    	for(int j=0; j<M; ++j)
    	  for(int i=0; i<M; ++i)
    	    cDerX[j][i] = tmpJ[i][j];
		}

		// Check if CDerZ must be computed.
		if( CDerZ )
		{
			// Cast the current value of CDerZ to a multi-dimensional array.
			su2double (*cDerZ)[M] = (su2double (*)[M]) CDerZ[l];

    	// Tensor product in i-direction to obtain the z-derivatives
    	// of the solution in the M-points in i-direction.
    	for(int j=0; j<K; ++j)
    	{
#pragma omp simd
      	for(int i=0; i<MP; ++i) tmpI[j][i] = zero;
      	for(int ii=0; ii<K; ++ii)
      	{
#pragma omp simd
    	    for(int i=0; i<MP; ++i)
    	      tmpI[j][i] += a[ii][i] * bkDer[j][ii];
    	  }
    	}

    	// Tensor product in j-direction to obtain the z-derivatives
    	// of the solution in the M-points in j-direction.
    	// This is the final result for the z-derivatives of the solution.
    	for(int i=0; i<M; ++i)
    	{
#pragma omp simd
      	for(int j=0; j<MP; ++j) tmpJ[i][j] = zero;
      	for(int jj=0; jj<K; ++jj)
      	{
#pragma omp simd
    	    for(int j=0; j<MP; ++j)
    	      tmpJ[i][j] += a[jj][j] * tmpI[jj][i];
    	  }
    	}

    	// Copy the values to the appropriate location in cDerZ.
    	for(int j=0; j<M; ++j)
    	  for(int i=0; i<M; ++i)
    	    cDerZ[j][i] = tmpJ[i][j];
		}
  }  // End of the loop over N.
}
