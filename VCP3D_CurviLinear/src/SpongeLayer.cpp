//------------------------------------------------------------------------------
// File, which contains the implementation of the member functions
// of element and solver classes, which correspond to sponge layer methods.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"


//-----------------------------------------------------------------------------

// Function, which initializes the sponge layer elements.
void SolverClass::InitSpongeLayer(void)
{
	// Error flag, default false.
	bool error      = false;
	// Extract the number of volume integration points.
	const int nInt  = mStandardHex.mNIntegration;
	// Explicitly extract the flags for readability.
	const bool imin = mInputParam->mSpongeLayerSpecified[0];
	const bool imax = mInputParam->mSpongeLayerSpecified[1];
	const bool jmin = mInputParam->mSpongeLayerSpecified[2];
	const bool jmax = mInputParam->mSpongeLayerSpecified[3];
	const bool kmin = mInputParam->mSpongeLayerSpecified[4];
	const bool kmax = mInputParam->mSpongeLayerSpecified[5];

	// Loop over the elements. This loop can be parallelized with OpenMP.
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(3), schedule(static)
#endif
  for(int k=1; k<=mNElemPerRankK; ++k)
    for(int j=1; j<=mNElemPerRankJ; ++j)
      for(int i=1; i<=mNElemPerRankI; ++i)
				mElements[i][j][k]->ConfigureSpongeLayer(mInputParam, &mStandardHex);

	// Decide which block side(s) this sponge layer belongs to and initialize it.
	
	// First, the x-direction.
	if( imin ) error = SetUpSpongeLayerIMIN();
	if( imax ) error = SetUpSpongeLayerIMAX();

	// Then, the y-direction.
	if( jmin ) error = SetUpSpongeLayerJMIN();
	if( jmax ) error = SetUpSpongeLayerJMAX();

	// Then, the z-direction.
	if( kmin ) error = SetUpSpongeLayerKMIN();
	if( kmax ) error = SetUpSpongeLayerKMAX();

	// Terminate if an error was found.
	if( error ) 
		TerminateAll("SolverClass::InitSpongeLayer", __FILE__, __LINE__,
				         "Inconsistent sponge layer location.");

	// Zero all non-used damping function dimensions, since they temporarily
	// store the coordinates. The padded vaues are already zeroed.
	
	// Loop over the elements and zero all non-used damping function dimensions. 
	// This is needed, since they temporarily store the coordinates. Note, the 
	// padded values are already zeroed. 
	// This loop can be parallelized with OpenMP.
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(3), schedule(static)
#endif
  for(int k=1; k<=mNElemPerRankK; ++k)
	{
    for(int j=1; j<=mNElemPerRankJ; ++j)
		{
      for(int i=1; i<=mNElemPerRankI; ++i)
			{
				if( mElements[i][j][k]->mSpongeLayerElement )
				{
					if( !mElements[i][j][k]->mSpongeXMIN && !mElements[i][j][k]->mSpongeXMAX ) 
						for(int l=0; l<nInt; ++l) mElements[i][j][k]->mDampingFunction[0][l] = zero;
					if( !mElements[i][j][k]->mSpongeYMIN && !mElements[i][j][k]->mSpongeYMAX ) 
						for(int l=0; l<nInt; ++l) mElements[i][j][k]->mDampingFunction[1][l] = zero;
					if( !mElements[i][j][k]->mSpongeZMIN && !mElements[i][j][k]->mSpongeZMAX ) 
						for(int l=0; l<nInt; ++l) mElements[i][j][k]->mDampingFunction[2][l] = zero;
				}
			}
		}
	}
	
	// Read the averaged solution from file. This is used as the target damping state.
	ReadAveragedSolutionSpongeLayer();
}

//-----------------------------------------------------------------------------

// Function, which sets up a sponge layer in the IMIN direction.
bool SolverClass::SetUpSpongeLayerIMIN(void)
{
	// Initialize error flag to false.
	bool error         = false;
	// Extract the damping constant and exponent in this layer.
	const su2double ss = mInputParam->mDampingConstant[0];
	const su2double aa = mInputParam->mDampingExponent[0];

	// Determine the index of the min and max nodal coordinate per element.
	const int I0 = 0;
	const int I1 = mStandardHex.mNPoly;


	// Initialize the min and max value of the relevant coordinate.
	su2double xmin = valLarge, xmax = -valLarge;

	// Loop over all elements and extract the min and max threshold values in the DOFs.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(3), reduction(min: xmin), reduction(max: xmax)
#endif
	for(int k=1; k<=mNElemPerRankK; ++k)
	{
		for(int j=1; j<=mNElemPerRankJ; ++j)
		{
			for(int i=1; i<=mNElemPerRankI; ++i)
			{
				if( mElements[i][j][k]->mSpongeXMIN )
				{
					xmin = std::min(xmin, mElements[i][j][k]->mCoorNodalGridDOFs[0][I0]);
					xmax = std::max(xmax, mElements[i][j][k]->mCoorNodalGridDOFs[0][I1]);
				}
			}
		}
	}


#ifdef HAVE_MPI
	{
		// Temporary buffer.
		su2double locBuf_xmin = xmin;
		su2double locBuf_xmax = xmax;

		// Send the local max and min on each rank and compute the global
		// min and max values and distribute them to all ranks.
		MPI_Allreduce(&locBuf_xmin, &xmin, 1, MPI_SU2DOUBLE, MPI_MIN, MPI_COMM_WORLD); 
		MPI_Allreduce(&locBuf_xmax, &xmax, 1, MPI_SU2DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
#endif

	// Deduce the inverse of the layer width.
	const su2double ovdx = one/(xmax - xmin);
	// Deduce the interface coordinate.
	const su2double x0   = xmax;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(3)
#endif
	for(int k=1; k<=mNElemPerRankK; ++k)
	{
		for(int j=1; j<=mNElemPerRankJ; ++j)
		{
			for(int i=1; i<=mNElemPerRankI; ++i)
			{
				if( mElements[i][j][k]->mSpongeXMIN )
				{
					// Initialize the damping function on the volume 
					// integration points in the x-direction.
#pragma omp simd
					for(int l=0; l<mStandardHex.mNIntegration; l++)
					{
						// Extract the temporarily stored x-coordinate.
						const su2double x = mElements[i][j][k]->mDampingFunction[0][l];
						// Construct the spatial damping function.
						mElements[i][j][k]->mDampingFunction[0][l] = ss*POW( ovdx*FABS(x - x0), aa ); 
					
						// Consistency check: make sure the element is within the interface threshold.
						if( x > x0 ) error = true;
					}

					// If a characteristic matching function is specified, deduce the scaling factor. 
					// This scaling is used to normalize the damping function, such that its range is
					// transformed from [0, sigma] into [0, 1].
					if( mInputParam->mCharacteristicMatchingLayer ) 
						mElements[i][j][k]->mMatchingScaleXMIN = ( FABS(ss) < 1.0e-4 ) ? one : one/ss;
				}
			}
		}
	}

	// Return the error flag.
	return error;
}

//-----------------------------------------------------------------------------

// Function, which sets up a sponge layer in the IMAX direction.
bool SolverClass::SetUpSpongeLayerIMAX(void)
{
	// Initialize error flag to false.
	bool error         = false;
	// Extract the damping constant and exponent in this layer.
	const su2double ss = mInputParam->mDampingConstant[1];
	const su2double aa = mInputParam->mDampingExponent[1];

	// Determine the index of the min and max nodal coordinate per element.
	const int I0 = 0;
	const int I1 = mStandardHex.mNPoly;


	// Initialize the min and max value of the relevant coordinate.
	su2double xmin = valLarge, xmax = -valLarge;

	// Loop over all elements and extract the min and max threshold values in the DOFs.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(3), reduction(min: xmin), reduction(max: xmax)
#endif
	for(int k=1; k<=mNElemPerRankK; ++k)
	{
		for(int j=1; j<=mNElemPerRankJ; ++j)
		{
			for(int i=1; i<=mNElemPerRankI; ++i)
			{
				if( mElements[i][j][k]->mSpongeXMAX )
				{
					xmin = std::min(xmin, mElements[i][j][k]->mCoorNodalGridDOFs[0][I0]);
					xmax = std::max(xmax, mElements[i][j][k]->mCoorNodalGridDOFs[0][I1]);
				}
			}
		}
	}


#ifdef HAVE_MPI
	{
		// Temporary buffer.
		su2double locBuf_xmin = xmin;
		su2double locBuf_xmax = xmax;

		// Send the local max and min on each rank and compute the global
		// min and max values and distribute them to all ranks.
		MPI_Allreduce(&locBuf_xmin, &xmin, 1, MPI_SU2DOUBLE, MPI_MIN, MPI_COMM_WORLD); 
		MPI_Allreduce(&locBuf_xmax, &xmax, 1, MPI_SU2DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
#endif

	// Deduce the inverse of the layer width.
	const su2double ovdx = one/(xmax - xmin);
	// Deduce the interface coordinate.
	const su2double x0   = xmin;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(3)
#endif
	for(int k=1; k<=mNElemPerRankK; ++k)
	{
		for(int j=1; j<=mNElemPerRankJ; ++j)
		{
			for(int i=1; i<=mNElemPerRankI; ++i)
			{
				if( mElements[i][j][k]->mSpongeXMAX )
				{
					// Initialize the damping function on the volume 
					// integration points in the x-direction.
#pragma omp simd
					for(int l=0; l<mStandardHex.mNIntegration; l++)
					{
						// Extract the temporarily stored x-coordinate.
						const su2double x = mElements[i][j][k]->mDampingFunction[0][l];
						// Construct the spatial damping function.
						mElements[i][j][k]->mDampingFunction[0][l] = ss*POW( ovdx*FABS(x - x0), aa ); 
					
						// Consistency check: make sure the element is within the interface threshold.
						if( x < x0 ) error = true;
					}
					
					// If a characteristic matching function is specified, deduce the scaling factor. 
					// This scaling is used to normalize the damping function, such that its range is
					// transformed from [0, sigma] into [0, 1].
					if( mInputParam->mCharacteristicMatchingLayer ) 
						mElements[i][j][k]->mMatchingScaleXMAX = ( FABS(ss) < 1.0e-4 ) ? one : one/ss;
				}
			}
		}
	}

	// Return the error flag.
	return error;
}

//-----------------------------------------------------------------------------

// Function, which sets up a sponge layer in the JMIN direction.
bool SolverClass::SetUpSpongeLayerJMIN(void)
{
	// Initialize error flag to false.
	bool error         = false;
	// Extract the damping constant and exponent in this layer.
	const su2double ss = mInputParam->mDampingConstant[2];
	const su2double aa = mInputParam->mDampingExponent[2];

	// Determine the index of the min and max nodal coordinate per element.
	const int J0 = 0;
	const int J1 = mStandardHex.mNPoly*mStandardHex.mNDOFs1D;


	// Initialize the min and max value of the relevant coordinate.
	su2double ymin = valLarge, ymax = -valLarge;

	// Loop over all elements and extract the min and max threshold values in the DOFs.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(3), reduction(min: ymin), reduction(max: ymax)
#endif
	for(int k=1; k<=mNElemPerRankK; ++k)
	{
		for(int j=1; j<=mNElemPerRankJ; ++j)
		{
			for(int i=1; i<=mNElemPerRankI; ++i)
			{
				if( mElements[i][j][k]->mSpongeYMIN )
				{
					ymin = std::min(ymin, mElements[i][j][k]->mCoorNodalGridDOFs[1][J0]);
					ymax = std::max(ymax, mElements[i][j][k]->mCoorNodalGridDOFs[1][J1]);
				}
			}
		}
	}


#ifdef HAVE_MPI
	{
		// Temporary buffer.
		su2double locBuf_ymin = ymin;
		su2double locBuf_ymax = ymax;

		// Send the local max and min on each rank and compute the global
		// min and max values and distribute them to all ranks.
		MPI_Allreduce(&locBuf_ymin, &ymin, 1, MPI_SU2DOUBLE, MPI_MIN, MPI_COMM_WORLD); 
		MPI_Allreduce(&locBuf_ymax, &ymax, 1, MPI_SU2DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
#endif

	// Deduce the inverse of the layer width.
	const su2double ovdy = one/(ymax - ymin);
	// Deduce the interface coordinate.
	const su2double y0   = ymax;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(3)
#endif
	for(int k=1; k<=mNElemPerRankK; ++k)
	{
		for(int j=1; j<=mNElemPerRankJ; ++j)
		{
			for(int i=1; i<=mNElemPerRankI; ++i)
			{
				if( mElements[i][j][k]->mSpongeYMIN )
				{
					// Initialize the damping function on the volume 
					// integration points in the y-direction.
#pragma omp simd
					for(int l=0; l<mStandardHex.mNIntegration; l++)
					{
						// Extract the temporarily stored y-coordinate.
						const su2double y = mElements[i][j][k]->mDampingFunction[1][l];
						// Construct the spatial damping function.
						mElements[i][j][k]->mDampingFunction[1][l] = ss*POW( ovdy*FABS(y - y0), aa ); 
					
						// Consistency check: make sure the element is within the interface threshold.
						if( y > y0 ) error = true;
					}
	
					// If a characteristic matching function is specified, deduce the scaling factor. 
					// This scaling is used to normalize the damping function, such that its range is
					// transformed from [0, sigma] into [0, 1].
					if( mInputParam->mCharacteristicMatchingLayer ) 
						mElements[i][j][k]->mMatchingScaleYMIN = ( FABS(ss) < 1.0e-4 ) ? one : one/ss;
				}
			}
		}
	}

	// Return the error flag.
	return error;
}

//-----------------------------------------------------------------------------

// Function, which sets up a sponge layer in the JMAX direction.
bool SolverClass::SetUpSpongeLayerJMAX(void)
{
	// Initialize error flag to false.
	bool error         = false;
	// Extract the damping constant and exponent in this layer.
	const su2double ss = mInputParam->mDampingConstant[3];
	const su2double aa = mInputParam->mDampingExponent[3];

	// Determine the index of the min and max nodal coordinate per element.
	const int J0 = 0;
	const int J1 = mStandardHex.mNPoly*mStandardHex.mNDOFs1D;


	// Initialize the min and max value of the relevant coordinate.
	su2double ymin = valLarge, ymax = -valLarge;

	// Loop over all elements and extract the min and max threshold values in the DOFs.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(3), reduction(min: ymin), reduction(max: ymax)
#endif
	for(int k=1; k<=mNElemPerRankK; ++k)
	{
		for(int j=1; j<=mNElemPerRankJ; ++j)
		{
			for(int i=1; i<=mNElemPerRankI; ++i)
			{
				if( mElements[i][j][k]->mSpongeYMAX )
				{
					ymin = std::min(ymin, mElements[i][j][k]->mCoorNodalGridDOFs[1][J0]);
					ymax = std::max(ymax, mElements[i][j][k]->mCoorNodalGridDOFs[1][J1]);
				}
			}
		}
	}


#ifdef HAVE_MPI
	{
		// Temporary buffer.
		su2double locBuf_ymin = ymin;
		su2double locBuf_ymax = ymax;

		// Send the local max and min on each rank and compute the global
		// min and max values and distribute them to all ranks.
		MPI_Allreduce(&locBuf_ymin, &ymin, 1, MPI_SU2DOUBLE, MPI_MIN, MPI_COMM_WORLD); 
		MPI_Allreduce(&locBuf_ymax, &ymax, 1, MPI_SU2DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
#endif

	// Deduce the inverse of the layer width.
	const su2double ovdy = one/(ymax - ymin);
	// Deduce the interface coordinate.
	const su2double y0   = ymin;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(3)
#endif
	for(int k=1; k<=mNElemPerRankK; ++k)
	{
		for(int j=1; j<=mNElemPerRankJ; ++j)
		{
			for(int i=1; i<=mNElemPerRankI; ++i)
			{
				if( mElements[i][j][k]->mSpongeYMAX )
				{
					// Initialize the damping function on the volume 
					// integration points in the y-direction.
#pragma omp simd
					for(int l=0; l<mStandardHex.mNIntegration; l++)
					{
						// Extract the temporarily stored y-coordinate.
						const su2double y = mElements[i][j][k]->mDampingFunction[1][l];
						// Construct the spatial damping function.
						mElements[i][j][k]->mDampingFunction[1][l] = ss*POW( ovdy*FABS(y - y0), aa ); 
					
						// Consistency check: make sure the element is within the interface threshold.
						if( y < y0 ) error = true;
					}

					// If a characteristic matching function is specified, deduce the scaling factor. 
					// This scaling is used to normalize the damping function, such that its range is
					// transformed from [0, sigma] into [0, 1].
					if( mInputParam->mCharacteristicMatchingLayer ) 
						mElements[i][j][k]->mMatchingScaleYMAX = ( FABS(ss) < 1.0e-4 ) ? one : one/ss;
				}
			}
		}
	}

	// Return the error flag.
	return error;
}

//-----------------------------------------------------------------------------

// Function, which sets up a sponge layer in the KMIN direction.
bool SolverClass::SetUpSpongeLayerKMIN(void)
{
	// Initialize error flag to false.
	bool error         = false;
	// Extract the damping constant and exponent in this layer.
	const su2double ss = mInputParam->mDampingConstant[4];
	const su2double aa = mInputParam->mDampingExponent[4];

	// Determine the index of the min and max nodal coordinate per element.
	const int K0 = 0;
	const int K1 = mStandardHex.mNPoly*mStandardHex.mNDOFs1D*mStandardHex.mNDOFs1D;


	// Initialize the min and max value of the relevant coordinate.
	su2double zmin = valLarge, zmax = -valLarge;

	// Loop over all elements and extract the min and max threshold values in the DOFs.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(3), reduction(min: kmin), reduction(max: kmax)
#endif
	for(int k=1; k<=mNElemPerRankK; ++k)
	{
		for(int j=1; j<=mNElemPerRankJ; ++j)
		{
			for(int i=1; i<=mNElemPerRankI; ++i)
			{
				if( mElements[i][j][k]->mSpongeZMIN )
				{
					zmin = std::min(zmin, mElements[i][j][k]->mCoorNodalGridDOFs[2][K0]);
					zmax = std::max(zmax, mElements[i][j][k]->mCoorNodalGridDOFs[2][K1]);
				}
			}
		}
	}


#ifdef HAVE_MPI
	{
		// Temporary buffer.
		su2double locBuf_zmin = zmin;
		su2double locBuf_zmax = zmax;

		// Send the local max and min on each rank and compute the global
		// min and max values and distribute them to all ranks.
		MPI_Allreduce(&locBuf_zmin, &zmin, 1, MPI_SU2DOUBLE, MPI_MIN, MPI_COMM_WORLD); 
		MPI_Allreduce(&locBuf_zmax, &zmax, 1, MPI_SU2DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
#endif

	// Deduce the inverse of the layer width.
	const su2double ovdz = one/(zmax - zmin);
	// Deduce the interface coordinate.
	const su2double z0   = zmax;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(3)
#endif
	for(int k=1; k<=mNElemPerRankK; ++k)
	{
		for(int j=1; j<=mNElemPerRankJ; ++j)
		{
			for(int i=1; i<=mNElemPerRankI; ++i)
			{
				if( mElements[i][j][k]->mSpongeZMIN )
				{
					// Initialize the damping function on the volume 
					// integration points in the z-direction.
#pragma omp simd
					for(int l=0; l<mStandardHex.mNIntegration; l++)
					{
						// Extract the temporarily stored z-coordinate.
						const su2double z = mElements[i][j][k]->mDampingFunction[2][l];
						// Construct the spatial damping function.
						mElements[i][j][k]->mDampingFunction[2][l] = ss*POW( ovdz*FABS(z - z0), aa ); 
					
						// Consistency check: make sure the element is within the interface threshold.
						if( z > z0 ) error = true;
					}
					
					// If a characteristic matching function is specified, deduce the scaling factor. 
					// This scaling is used to normalize the damping function, such that its range is
					// transformed from [0, sigma] into [0, 1].
					if( mInputParam->mCharacteristicMatchingLayer ) 
						mElements[i][j][k]->mMatchingScaleZMIN = ( FABS(ss) < 1.0e-4 ) ? one : one/ss;
				}
			}
		}
	}

	// Return the error flag.
	return error;
}

//-----------------------------------------------------------------------------

// Function, which sets up a sponge layer in the KMAX direction.
bool SolverClass::SetUpSpongeLayerKMAX(void)
{
	// Initialize error flag to false.
	bool error         = false;
	// Extract the damping constant and exponent in this layer.
	const su2double ss = mInputParam->mDampingConstant[5];
	const su2double aa = mInputParam->mDampingExponent[5];

	// Determine the index of the min and max nodal coordinate per element.
	const int K0 = 0;
	const int K1 = mStandardHex.mNPoly*mStandardHex.mNDOFs1D*mStandardHex.mNDOFs1D;


	// Initialize the min and max value of the relevant coordinate.
	su2double zmin = valLarge, zmax = -valLarge;

	// Loop over all elements and extract the min and max threshold values in the DOFs.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(3), reduction(min: zmin), reduction(max: zmax)
#endif
	for(int k=1; k<=mNElemPerRankK; ++k)
	{
		for(int j=1; j<=mNElemPerRankJ; ++j)
		{
			for(int i=1; i<=mNElemPerRankI; ++i)
			{
				if( mElements[i][j][k]->mSpongeZMAX )
				{
					zmin = std::min(zmin, mElements[i][j][k]->mCoorNodalGridDOFs[2][K0]);
					zmax = std::max(zmax, mElements[i][j][k]->mCoorNodalGridDOFs[2][K1]);
				}
			}
		}
	}


#ifdef HAVE_MPI
	{
		// Temporary buffer.
		su2double locBuf_zmin = zmin;
		su2double locBuf_zmax = zmax;

		// Send the local max and min on each rank and compute the global
		// min and max values and distribute them to all ranks.
		MPI_Allreduce(&locBuf_zmin, &zmin, 1, MPI_SU2DOUBLE, MPI_MIN, MPI_COMM_WORLD); 
		MPI_Allreduce(&locBuf_zmax, &zmax, 1, MPI_SU2DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}
#endif

	// Deduce the inverse of the layer width.
	const su2double ovdz = one/(zmax - zmin);
	// Deduce the interface coordinate.
	const su2double z0   = zmin;

#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(3)
#endif
	for(int k=1; k<=mNElemPerRankK; ++k)
	{
		for(int j=1; j<=mNElemPerRankJ; ++j)
		{
			for(int i=1; i<=mNElemPerRankI; ++i)
			{
				if( mElements[i][j][k]->mSpongeZMAX )
				{
					// Initialize the damping function on the volume 
					// integration points in the z-direction.
#pragma omp simd
					for(int l=0; l<mStandardHex.mNIntegration; l++)
					{
						// Extract the temporarily stored z-coordinate.
						const su2double z = mElements[i][j][k]->mDampingFunction[2][l];
						// Construct the spatial damping function.
						mElements[i][j][k]->mDampingFunction[2][l] = ss*POW( ovdz*FABS(z - z0), aa ); 
					
						// Consistency check: make sure the element is within the interface threshold.
						if( z < z0 ) error = true;
					}
					
					// If a characteristic matching function is specified, deduce the scaling factor. 
					// This scaling is used to normalize the damping function, such that its range is
					// transformed from [0, sigma] into [0, 1].
					if( mInputParam->mCharacteristicMatchingLayer ) 
						mElements[i][j][k]->mMatchingScaleZMAX = ( FABS(ss) < 1.0e-4 ) ? one : one/ss;
				}
			}
		}
	}

	// Return the error flag.
	return error;
}

//-----------------------------------------------------------------------------

// Function, which reads the averaged solution in the sponge layer.
void SolverClass::ReadAveragedSolutionSpongeLayer(void)
{
  // Easier storage of the total number of elements in the three index directions.
  const int nElemI = mInputParam->mSubfaces[1][0]->mElemIBeg;
  const int nElemJ = mInputParam->mSubfaces[3][0]->mElemJBeg;
  const int nElemK = mInputParam->mSubfaces[5][0]->mElemKBeg;

  // Determine the total number of DOFs in the simulation.
  const int nDOFs    = mStandardHex.mNDOFs;
  const int nDOFsTot = nElemI*nElemJ*nElemK*nDOFs;

  // Determine the precision of the restart file.
  short precisionRestartFile = 0;
  if(rank == 0)
    precisionRestartFile = DeterminePrecisionRestartFile(mInputParam->mAveragedSolutionFile.c_str());

#ifdef HAVE_MPI
  // Broadcast the precision to all ranks.
  MPI_Bcast(&precisionRestartFile, 1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

  // Open the restart file for reading and check if it went OK.
  MPI_File fh;
  if(MPI_File_open(MPI_COMM_WORLD, mInputParam->mAveragedSolutionFile.c_str(),
                   MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) != MPI_SUCCESS)
    TerminateAll("SolverClass::ReadAveragedSolutionSpongeLayer", __FILE__, __LINE__,
                 "Averaged solution file could not be opened for reading");

  // Read the first 5 integer variables of the solution file. All processors
  // read them, because a collective read is nonsense here.
  int header[5];
  MPI_File_read(fh, header, 5, MPI_INT, MPI_STATUS_IGNORE);

  // Check if byte swapping must be applied.
  bool byteSwap = false;
  if(header[0] != SU2_MagicNumber)
  {
    byteSwap = true;
    SwapBytes(header, sizeof(int), 5);
  }

  // Check if this is indeed an SU2 restart file.
  if(header[0] != SU2_MagicNumber)
    TerminateAll("SolverClass::ReadAveragedSolutionSpongeLayer", __FILE__, __LINE__,
                 "Averaged solution file is not an SU2 restart file");

  // Easier storage of the number of variables and DOFS in the restart file.
  const int nVarRestart  = header[1];
  const int nDOFsRestart = header[2];

  // Check if the number of DOFs is correct.
  if(nDOFsRestart != nDOFsTot)
    TerminateAll("SolverClass::ReadAveragedSolutionSpongeLayer", __FILE__, __LINE__,
                 "Wrong number of DOFs in the averaged solution file");

  // Determine the byte offset where the data for the solution
  // variables starts.
  MPI_Offset offset = 5*sizeof(int) + nVarRestart*CGNS_STRING_SIZE;

  // Make a distinction between the precision and call ReadSolutionFromRestartFile
  // with the appropriate template parameter to read the restart file accordingly.
  switch( precisionRestartFile )
  {
    case 4:
    {
      ReadSpongeLayerSolutionFromRestartFile<float>(fh, nVarRestart, byteSwap, offset);
      break;
    }

    case 8:
    {
      ReadSpongeLayerSolutionFromRestartFile<double>(fh, nVarRestart, byteSwap, offset);
      break;
    }

    default:
      TerminateAll("SolverClass::ReadAveragedSolutionSpongeLayer", __FILE__, __LINE__,
                   "Unknown precision encountered in the averaged solution file");
  }

  // Close the file.
  MPI_File_close(&fh);

#else
  // No MPI is used, so the reading of the file can be done sequentially.
  // Open the averaged solution file for binary reading.
  FILE *restartFile = std::fopen(mInputParam->mAveragedSolutionFile.c_str(), "rb");
  if( !restartFile )
    Terminate("SolverClass::ReadAveragedSolutionSpongeLayer", __FILE__, __LINE__,
              "Averaged solution file could not be opened for binary reading");

  // Read the first 5 integer variables of the solution file.
  int header[5];
  if(std::fread(header, sizeof(int), 5, restartFile) != 5)
    Terminate("SolverClass::ReadAveragedSolutionSpongeLayer", __FILE__, __LINE__,
              "File header could not be read");

  // Check if byte swapping must be applied.
  bool byteSwap = false;
  if(header[0] != SU2_MagicNumber)
  {
    byteSwap = true;
    SwapBytes(header, sizeof(int), 5);
  }

  // Check if this is indeed an SU2 restart file.
  if(header[0] != SU2_MagicNumber)
    Terminate("SolverClass::ReadAveragedSolutionSpongeLayer", __FILE__, __LINE__,
              "Averaged solution file is not an SU2 restart file");

  // Easier storage of the number of variables and DOFS in the restart file.
  const int nVarRestart  = header[1];
  const int nDOFsRestart = header[2];

  // Check if the number of DOFs is correct.
  if(nDOFsRestart != nDOFsTot)
    Terminate("SolverClass::ReadAveragedSolutionSpongeLayer", __FILE__, __LINE__,
              "Wrong number of DOFs in the averaged solution file");

  // Read the names of the solution variables.
  for(int i=0; i<nVarRestart; ++i)
  {
    char varName[CGNS_STRING_SIZE];
    if(std::fread(varName, sizeof(char), CGNS_STRING_SIZE, restartFile) != CGNS_STRING_SIZE)
      Terminate("SolverClass::ReadAveragedSolutionSpongeLayer", __FILE__, __LINE__,
                "Variable name could not be read");
  }

  // Make a distinction between the precision and call ReadSpongeLayerSolutionFromRestartFile
  // with the appropriate template parameter to read the restart file accordingly.
  switch( precisionRestartFile )
  {
    case 4:
    {
      ReadSpongeLayerSolutionFromRestartFile<float>(restartFile, nVarRestart, byteSwap);
      break;
    }

    case 8:
    {
      ReadSpongeLayerSolutionFromRestartFile<double>(restartFile, nVarRestart, byteSwap);
      break;
    }

    default:
      Terminate("SolverClass::ReadAveragedSolutionSpongeLayer", __FILE__, __LINE__,
                "Unknown precision encountered in the averaged solution file");
  }

  // Close the file again.
  std::fclose(restartFile);
#endif

	// So far, the nodal averaged solution is read and stored in mSolOld. Now,
	// we need to interpolate it onto the volume integration nodes and store it 
	// inside mDampingState.
	// This loop can be parallelized with OpenMP.
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(3), schedule(static)
#endif
  for(int k=1; k<=mNElemPerRankK; ++k)
	{
    for(int j=1; j<=mNElemPerRankJ; ++j)
		{
      for(int i=1; i<=mNElemPerRankI; ++i)
			{
				if( mElements[i][j][k]->mSpongeLayerElement )
				{
					mElements[i][j][k]->InitializeDampingState(mInputParam, &mStandardHex);
				}
			}
		}
	}
}

//------------------------------------------------------------------------------

// Function, which configures a sponge layer element.
void ElementClass::ConfigureSpongeLayer(const InputParamClass      *inputParam,
		                                    const StandardElementClass *standardHex)
{
	// Extract global element indices.
	const int I = mGlobalInd[0];
	const int J = mGlobalInd[1];
	const int K = mGlobalInd[2];

	// Extract total number of elements per dimension.
	const int nxElem = inputParam->mSubfaces[1][0]->mElemIBeg;
	const int nyElem = inputParam->mSubfaces[3][0]->mElemJBeg;
	const int nzElem = inputParam->mSubfaces[5][0]->mElemKBeg;

	// Extract total number of sponge layer elements, per x-dimension.
	const int nbIMin = inputParam->mNSpongeLayerElements[0];
	const int nbIMax = inputParam->mNSpongeLayerElements[1];

	// Extract total number of sponge layer elements, per y-dimension.
	const int nbJMin = inputParam->mNSpongeLayerElements[2];
	const int nbJMax = inputParam->mNSpongeLayerElements[3];

	// Extract total number of sponge layer elements, per z-dimension.
	const int nbKMin = inputParam->mNSpongeLayerElements[4];
	const int nbKMax = inputParam->mNSpongeLayerElements[5];

	// Check if current element belongs in the x-dimension.
	mSpongeXMIN = ( I <   nbIMin         ) ? true : false;  
  mSpongeXMAX = ( I >= (nxElem-nbIMax) ) ? true : false;

	// Check if current element belongs in the y-dimension.
	mSpongeYMIN = ( J <   nbJMin         ) ? true : false;
	mSpongeYMAX = ( J >= (nyElem-nbJMax) ) ? true : false;

	// Check if current element belongs in the z-dimension.
	mSpongeZMIN = ( K <   nbKMin         ) ? true : false;
	mSpongeZMAX = ( K >= (nzElem-nbKMax) ) ? true : false;

	// In case this element is not a sponge layer, return from function.
	// Otherwise, flag the current element as a sponge layer element and proceed.
	if( !mSpongeXMIN && !mSpongeXMAX && 
			!mSpongeYMIN && !mSpongeYMAX && 
			!mSpongeZMIN && !mSpongeZMAX ) return;
	else mSpongeLayerElement = true;

	// Determine number of (padded) integration points.
	const int nIntPad = standardHex->mNIntegrationPad;
	// Determine number of          integration points.
	const int nInt    = standardHex->mNIntegration;
	// Determine number of       1D integration points.
	const int nInt1D  = standardHex->mNIntegration1D;
	// Determine number of       1D nodal       points.
	const int nDOFs1D = standardHex->mNDOFs1D;

	// Determine the number of bytes that must be allocated by
	// each of the member variables.
	const size_t nBytes = nIntPad * sizeof(su2double);

	// Allocate memory for the damping function for the arrays needed. 
	// Note, the damping is general in all 3-dimensions for now.
	mDampingFunction.resize(3, NULL);
	
	for(int i=0; i<3; i++)
	{
		mDampingFunction[i] = (su2double *) AllocateMemory(nBytes);

		if( !mDampingFunction[i] )
			Terminate("ElementClass::ConfigureSpongeLayer", __FILE__, __LINE__,
								"Memory allocation for mDampingFunction failed.");
	}

	// Determine the coordinates in the volume integration points. 
	// Note, these are temporarily stored in mDampingFunction.
  TensorProductSolAndGradVolume(nInt1D, 3, nDOFs1D, standardHex->mLagrangeInt1D,
                                NULL, mCoorNodalGridDOFs.data(),
                                mDampingFunction.data(), NULL, NULL, NULL);

	// Initialize the padded values to avoid problems later on.
	for(int l=0; l<3; ++l) 
		for(int m=nInt; m<nIntPad; ++m) 
			mDampingFunction[l][m] = zero;


	// Allocate memory for the target-state being damped against.
	mDampingState.resize(nVar, NULL);

	for(int i=0; i<nVar; i++)
	{
		mDampingState[i] = (su2double *) AllocateMemory(nBytes);

		if( !mDampingState[i] )
			Terminate("ElementClass::ConfigureSpongeLayer", __FILE__, __LINE__,
								"Memory allocation for mDampingState failed.");
	}

}

//------------------------------------------------------------------------------

// Function, which initializes the damping state in a sponge layer.
void ElementClass::InitializeDampingState(const InputParamClass      *inputParam,
		                                      const StandardElementClass *standardHex)
{
	// Determine number of (padded) integration points.
	const int nIntPad = standardHex->mNIntegrationPad;
	// Determine number of          integration points.
	const int nInt    = standardHex->mNIntegration;
	// Determine number of       1D integration points.
	const int nInt1D  = standardHex->mNIntegration1D;
	// Determine number of       1D nodal       points.
	const int nDOFs1D = standardHex->mNDOFs1D;

	// Invert the reference values.
  const su2double rhoRefInv = one/rhoRef;
  const su2double uRefInv   = one/uRef;
  const su2double pRefInv   = one/pRef;


  // Call TensorProductSolAndGradVolume to carry out the actual tensor product
  // multiplication to compute the solution in the integration points. 
	// Note, the averaged nodal solution was temporarily stored in mSolOld.
  TensorProductSolAndGradVolume(nInt1D, nVar, nDOFs1D, standardHex->mLagrangeInt1D,
                                NULL, mSolOld.data(),
                                mDampingState.data(), NULL, NULL, NULL);

	// Initialize the padded values to avoid problems.
	for(int l=0; l<nVar; ++l)
		for(int m=nInt; m<nIntPad; ++m) 
			mDampingState[l][m] = mDampingState[l][0];


	// Convert the primitive variables into conservative variables.
	// Easier storage of some expressions involving gamma.
	const su2double gm1   = GamConstant - one;
	const su2double ovgm1 = one/gm1;

	// Loop over the (padded) number of integration points.
#pragma omp simd
 	for(int l=0; l<nIntPad; ++l)
 	{
		// Extract the averaged primitive state and non-dimensionalize it.
		const su2double rho   = mDampingState[0][l]*rhoRefInv;
		const su2double u     = mDampingState[1][l]*uRefInv;
		const su2double v     = mDampingState[2][l]*uRefInv;
		const su2double w     = mDampingState[3][l]*uRefInv;
		const su2double p     = mDampingState[4][l]*pRefInv;

	  // Assemble the target-state variables in conservative form.
		mDampingState[0][l]   = rho;
		mDampingState[1][l]   = rho*u;
		mDampingState[2][l]   = rho*v;
		mDampingState[3][l]   = rho*w;
		mDampingState[4][l]   = p*ovgm1 + half*rho*(u*u + v*v + w*w);
	}
}



