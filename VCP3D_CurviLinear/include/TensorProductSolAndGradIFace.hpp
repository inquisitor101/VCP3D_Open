// Template function to carry out the actual tensor product to compute the
// solution and its gradients on a face in i-direction. See
// TensorProductSolAndGradIFace for the details. The template construction makes
// sure that the values of K and M are known at compile time, such that the
// compiler can optimize much better.
template<int K, int M>
void TPSAGI(const int N,
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
  su2double tmpJ[K][MP], tmpK[M][MP];
  su2double bi[K][K], biDer[K][K];

  // Outer loop over N.
  for(int l=0; l<N; ++l)
  {
    // Cast the current value of B, C and CDer to multi-dimensional arrays.
    const su2double (*b)[K][K] = (const su2double (*)[K][K]) B[l];

		// Compute the tensor product in i-direction for both AFace
    // and ADerFace. The tensor products in the i-direction are
    // different, because both AFace and ADerFace are rank 1 tensors.
    for(int k=0; k<K; ++k)
    {
#pragma omp simd safelen(K)
      for(int j=0; j<K; ++j)
      {
        bi[k][j]    = zero;
        biDer[k][j] = zero;
      }
    }

    for(int k=0; k<K; ++k)
    {
      for(int i=0; i<K; ++i)
      {
#pragma omp simd safelen(K)
        for(int j=0; j<K; ++j)
        {
          bi[k][j]    += AFace[i]   *b[k][j][i];
          biDer[k][j] += ADerFace[i]*b[k][j][i];
        }
      }
    }

    // Tensor product in j-direction to obtain the
    // solution in the M-points in j-direction.
    for(int k=0; k<K; ++k)
    {
#pragma omp simd
      for(int j=0; j<MP; ++j) tmpJ[k][j] = zero;
      for(int jj=0; jj<K; ++jj)
      {
#pragma omp simd
        for(int j=0; j<MP; ++j)
          tmpJ[k][j] += a[jj][j] * bi[k][jj];
      }
    }


		// Check if C must be computed.
		if( C )
		{
			// Cast the current value of C to a multi-dimensional array.
    	su2double (*c)[M] = (su2double (*)[M]) C[l];

    	// Tensor product in k-direction to obtain the
    	// solution in the M-points in k-direction.
    	// This is the final result for the solution.
    	for(int j=0; j<M; ++j)
    	{
#pragma omp simd
      	for(int k=0; k<MP; ++k) tmpK[j][k] = zero;
      	for(int kk=0; kk<K; ++kk)
      	{
#pragma omp simd
    	    for(int k=0; k<MP; ++k)
    	      tmpK[j][k] += a[kk][k] * tmpJ[kk][j];
    	  }
    	}

    	// Copy the values to the appropriate location in c.
    	for(int k=0; k<M; ++k)
    	  for(int j=0; j<M; ++j)
    	    c[k][j] = tmpK[j][k];
		}

		// Check if CDerZ must be computed.
		if( CDerZ )
		{
			// Cast the current value of CDerZ to a multi-dimensional array.
    	su2double (*cDerZ)[M] = (su2double (*)[M]) CDerZ[l];
    	
			// The values of tmpJ can be reused to compute the derivatives of
    	// the solution in k-direction. Carry out the tensor product to
    	// compute these values.
    	for(int j=0; j<M; ++j)
    	{
#pragma omp simd
      	for(int k=0; k<MP; ++k) tmpK[j][k] = zero;
      	for(int kk=0; kk<K; ++kk)
      	{
#pragma omp simd
    	    for(int k=0; k<MP; ++k)
    	      tmpK[j][k] += aDer[kk][k] * tmpJ[kk][j];
    	  }
    	}

    	// Copy the values to the appropriate location in cDerZ.
    	for(int k=0; k<M; ++k)
    	  for(int j=0; j<M; ++j)
    	    cDerZ[k][j] = tmpK[j][k];
		}

		// Check if CDerY must be computed.
		if( CDerY )
		{
			// Cast the current value of CDerY to a multi-dimensional array.
			su2double (*cDerY)[M] = (su2double (*)[M]) CDerY[l];
    	
			// Tensor product in j-direction to obtain the y-derivatives
    	// of the solution in the M-points in j-direction.
    	for(int k=0; k<K; ++k)
    	{
#pragma omp simd
      	for(int j=0; j<MP; ++j) tmpJ[k][j] = zero;
      	for(int jj=0; jj<K; ++jj)
      	{
#pragma omp simd
    	    for(int j=0; j<MP; ++j)
    	      tmpJ[k][j] += aDer[jj][j] * bi[k][jj];
    	  }
    	}

    	// Tensor product in k-direction to obtain the y-derivatives
    	// of the solution in the M-points in k-direction.
    	// This is the final result for the y-derivatives of the solution.
    	for(int j=0; j<M; ++j)
    	{
#pragma omp simd
      	for(int k=0; k<MP; ++k) tmpK[j][k] = zero;
      	for(int kk=0; kk<K; ++kk)
      	{
#pragma omp simd
    	    for(int k=0; k<MP; ++k)
    	      tmpK[j][k] += a[kk][k] * tmpJ[kk][j];
    	  }
    	}

    	// Copy the values to the appropriate location in cDerY.
    	for(int k=0; k<M; ++k)
    	  for(int j=0; j<M; ++j)
    	    cDerY[k][j] = tmpK[j][k];
		}

		// Check if CDerX must be computed.
		if( CDerX )
		{
			// Cast the current value of CDerX to a multi-dimensional array.
			su2double (*cDerX)[M] = (su2double (*)[M]) CDerX[l];

    	// Tensor product in j-direction to obtain the x-derivatives
    	// of the solution in the M-points in j-direction.
    	for(int k=0; k<K; ++k)
    	{
#pragma omp simd
      	for(int j=0; j<MP; ++j) tmpJ[k][j] = zero;
      	for(int jj=0; jj<K; ++jj)
      	{
#pragma omp simd
    	    for(int j=0; j<MP; ++j)
    	      tmpJ[k][j] += a[jj][j] * biDer[k][jj];
    	  }
    	}

    	// Tensor product in k-direction to obtain the x-derivatives
    	// of the solution in the M-points in k-direction.
    	// This is the final result for the x-derivatives of the solution.
    	for(int j=0; j<M; ++j)
    	{
#pragma omp simd
      	for(int k=0; k<MP; ++k) tmpK[j][k] = zero;
      	for(int kk=0; kk<K; ++kk)
      	{
#pragma omp simd
    	    for(int k=0; k<MP; ++k)
    	      tmpK[j][k] += a[kk][k] * tmpJ[kk][j];
    	  }
    	}

    	// Copy the values to the appropriate location in cDerX.
    	for(int k=0; k<M; ++k)
    	  for(int j=0; j<M; ++j)
    	    cDerX[k][j] = tmpK[j][k];
		}
  }  // End of the loop over N.
}
