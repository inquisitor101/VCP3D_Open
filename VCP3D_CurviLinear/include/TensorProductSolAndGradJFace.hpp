// Template function to carry out the actual tensor product to compute the
// solution and its gradients on a face in j-direction. See
// TensorProductSolAndGradJFace for the details. The template construction makes
// sure that the values of K and M are known at compile time, such that the
// compiler can optimize much better.
template<int K, int M>
void TPSAGJ(const int N,
            su2double *A,
            su2double *ADer,
            su2double *AFace,
            su2double *ADerFace,
            su2double **B,
            su2double **C,
            su2double **CDerX,
            su2double **CDerY,
            su2double  **CDerZ)
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
  su2double tmpI[K][MP], tmpK[M][MP];
  su2double bj[K][K], bjDer[K][K];

  // Outer loop over N.
  for(int l=0; l<N; ++l)
  {
    // Cast the current value of B, C and CDer to multi-dimensional arrays.
    const su2double (*b)[K][K] = (const su2double (*)[K][K]) B[l];

		// Compute the tensor product in j-direction for both AFace
    // and ADerFace. The tensor products in the j-direction are
    // different, because both AFace and ADerFace are rank 1 tensors.
    for(int k=0; k<K; ++k)
    {
#pragma omp simd safelen(K)
      for(int i=0; i<K; ++i)
      {
        bj[k][i]    = zero;
        bjDer[k][i] = zero;
      }
    }

    for(int k=0; k<K; ++k)
    {
      for(int j=0; j<K; ++j)
      {
#pragma omp simd safelen(K)
        for(int i=0; i<K; ++i)
        {
          bj[k][i]    += AFace[j]   *b[k][j][i];
          bjDer[k][i] += ADerFace[j]*b[k][j][i];
        }
      }
    }

    // Tensor product in i-direction to obtain the
    // solution in the M-points in i-direction.
    for(int k=0; k<K; ++k)
    {
#pragma omp simd
      for(int i=0; i<MP; ++i) tmpI[k][i] = zero;
      for(int ii=0; ii<K; ++ii)
      {
#pragma omp simd
        for(int i=0; i<MP; ++i)
          tmpI[k][i] += a[ii][i] * bj[k][ii];
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
    	for(int i=0; i<M; ++i)
    	{
#pragma omp simd
      	for(int k=0; k<MP; ++k) tmpK[i][k] = zero;
      	for(int kk=0; kk<K; ++kk)
      	{
#pragma omp simd
    	    for(int k=0; k<MP; ++k)
    	      tmpK[i][k] += a[kk][k] * tmpI[kk][i];
    	  }
    	}

    	// Copy the values to the appropriate location in c.
    	for(int k=0; k<M; ++k)
    	  for(int i=0; i<M; ++i)
    	    c[k][i] = tmpK[i][k];
		}

		// Check if CDerZ must be computed.
		if( CDerZ )
		{
			// Cast the current value of CDerZ to a multi-dimensional array.
			su2double (*cDerZ)[M] = (su2double (*)[M]) CDerZ[l];
    
			// The values of tmpI can be reused to compute the derivatives of
    	// the solution in k-direction. Carry out the tensor product to
    	// compute these values.
    	for(int i=0; i<M; ++i)
    	{
#pragma omp simd
      	for(int k=0; k<MP; ++k) tmpK[i][k] = zero;
      	for(int kk=0; kk<K; ++kk)
      	{
#pragma omp simd
    	    for(int k=0; k<MP; ++k)
    	      tmpK[i][k] += aDer[kk][k] * tmpI[kk][i];
    	  }
    	}

    	// Copy the values to the appropriate location in cDerZ.
    	for(int k=0; k<M; ++k)
    	  for(int i=0; i<M; ++i)
    	    cDerZ[k][i] = tmpK[i][k];
		}

		// Check if CDerX must be computed.
		if( CDerX )
		{
			// Cast the current value of CDerX to a multi-dimensional array.
			su2double (*cDerX)[M] = (su2double (*)[M]) CDerX[l];

    	// Tensor product in i-direction to obtain the x-derivatives
    	// of the solution in the M-points in i-direction.
    	for(int k=0; k<K; ++k)
    	{
#pragma omp simd
      	for(int i=0; i<MP; ++i) tmpI[k][i] = zero;
      	for(int ii=0; ii<K; ++ii)
      	{
#pragma omp simd
    	    for(int i=0; i<MP; ++i)
    	      tmpI[k][i] += aDer[ii][i] * bj[k][ii];
    	  }
    	}

    	// Tensor product in k-direction to obtain the x-derivatives
    	// of the solution in the M-points in k-direction.
    	// This is the final result for the x-derivatives of the solution.
    	for(int i=0; i<M; ++i)
    	{
#pragma omp simd
      	for(int k=0; k<MP; ++k) tmpK[i][k] = zero;
      	for(int kk=0; kk<K; ++kk)
      	{
#pragma omp simd
    	    for(int k=0; k<MP; ++k)
    	      tmpK[i][k] += a[kk][k] * tmpI[kk][i];
    	  }
    	}

    	// Copy the values to the appropriate location in cDerX.
    	for(int k=0; k<M; ++k)
    	  for(int i=0; i<M; ++i)
    	    cDerX[k][i] = tmpK[i][k];
		}

		// Check if CDerY must be computed.
		if( CDerY )
		{
			// Cast the current value of CDerY to a multi-dimensional array.
			su2double (*cDerY)[M] = (su2double (*)[M]) CDerY[l];

    	// Tensor product in i-direction to obtain the y-derivatives
    	// of the solution in the M-points in i-direction.
    	for(int k=0; k<K; ++k)
    	{
#pragma omp simd
      	for(int i=0; i<MP; ++i) tmpI[k][i] = zero;
      	for(int ii=0; ii<K; ++ii)
      	{
#pragma omp simd
    	    for(int i=0; i<MP; ++i)
    	      tmpI[k][i] += a[ii][i] * bjDer[k][ii];
    	  }
    	}

    	// Tensor product in k-direction to obtain the y-derivatives
    	// of the solution in the M-points in k-direction.
    	// This is the final result for the y-derivatives of the solution.
    	for(int i=0; i<M; ++i)
    	{
#pragma omp simd
      	for(int k=0; k<MP; ++k) tmpK[i][k] = zero;
      	for(int kk=0; kk<K; ++kk)
      	{
#pragma omp simd
    	    for(int k=0; k<MP; ++k)
    	      tmpK[i][k] += a[kk][k] * tmpI[kk][i];
    	  }
    	}

    	// Copy the values to the appropriate location in cDerY.
    	for(int k=0; k<M; ++k)
    	  for(int i=0; i<M; ++i)
    	    cDerY[k][i] = tmpK[i][k];
		}
  }  // End of the loop over N.
}
