//------------------------------------------------------------------------------
// File, which contains the implementation of the member functions of the
// the class StandardElementClass.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"

//------------------------------------------------------------------------------

// Constructor. Initialize the pointer variables.
StandardElementClass::StandardElementClass()
{
  mIntWeights                        = NULL;
  mLegendreDOFs1D                    = NULL;
  mDerLegendreDOFs1D                 = NULL;
  mLegendreInt1D                     = NULL;
  mLegendreInt1DTranspose            = NULL;
  mDerLegendreInt1D                  = NULL;
  mDerLegendreInt1DTranspose         = NULL;
  mIntWeightsFace                    = NULL;
  mLegendreMinFace1D                 = NULL;
  mLegendreMaxFace1D                 = NULL;
  mDerLegendreMinFace1D              = NULL;
  mDerLegendreMaxFace1D              = NULL;
  mVandermonde1D                     = NULL;
  mVandermonde1DInverse              = NULL;
  mVandermonde1DTranspose            = NULL;
  mBasisCenter                       = NULL;
	mLagrangeInt1D                     = NULL;
	mLagrangeInt1DTranspose            = NULL;
	mDerLagrangeInt1D                  = NULL;
	mDerLagrangeInt1DTranspose         = NULL;
	mDerLagrangeDOFs1D                 = NULL;
	mDerLagrangeDOFs1DTranspose        = NULL;
	mPolynomialCorrectionMatrix        = NULL;
	mDerLagrangeInt1D_BaseInt          = NULL;
	mDerLagrangeInt1D_BaseIntTranspose = NULL;
}

//------------------------------------------------------------------------------

// Destructor. Release the allocate memory.
StandardElementClass::~StandardElementClass()
{
  FreeMemory((void **) &mIntWeights);
  FreeMemory((void **) &mLegendreDOFs1D);
  FreeMemory((void **) &mDerLegendreDOFs1D);
  FreeMemory((void **) &mLegendreInt1D);
  FreeMemory((void **) &mLegendreInt1DTranspose);
  FreeMemory((void **) &mDerLegendreInt1D);
  FreeMemory((void **) &mDerLegendreInt1DTranspose);
  FreeMemory((void **) &mIntWeightsFace);
  FreeMemory((void **) &mLegendreMinFace1D);
  FreeMemory((void **) &mLegendreMaxFace1D);
  FreeMemory((void **) &mDerLegendreMinFace1D);
  FreeMemory((void **) &mDerLegendreMaxFace1D);
  FreeMemory((void **) &mVandermonde1D);
  FreeMemory((void **) &mVandermonde1DInverse);
  FreeMemory((void **) &mVandermonde1DTranspose);
  FreeMemory((void **) &mBasisCenter);
	FreeMemory((void **) &mLagrangeInt1D);
	FreeMemory((void **) &mLagrangeInt1DTranspose);
	FreeMemory((void **) &mDerLagrangeInt1D);
	FreeMemory((void **) &mDerLagrangeInt1DTranspose);
	FreeMemory((void **) &mDerLagrangeDOFs1D);
	FreeMemory((void **) &mDerLagrangeDOFs1DTranspose);
	FreeMemory((void **) &mPolynomialCorrectionMatrix);
	FreeMemory((void **) &mDerLagrangeInt1D_BaseInt);
	FreeMemory((void **) &mDerLagrangeInt1D_BaseIntTranspose);
}

//------------------------------------------------------------------------------

// Function, which computes the values of the basis functions in the
// center of the element.
void StandardElementClass::BasisFunctionsCenter(void)
{
  // Allocate the memory for mBasisCenter.
  mBasisCenter = (su2double *) AllocateMemory(mNDOFs*sizeof(su2double));
  if( !mBasisCenter )
    Terminate("StandardElementClass::BasisFunctionsCenter", __FILE__, __LINE__,
              "Memory allocation failure for mBasisCenter.");

  // Determine the values of the Legendre polynomials in the center
  // of the element.
  std::vector<su2double> rCenter(1);
  rCenter[0] = zero;

  std::vector<su2double> V(mNDOFs1D);
  Vandermonde1D(mNDOFs1D, rCenter, V);

  // Loop over the three dimensions and compute the 3D basis functions
  // in the center of the element.
  int ii = 0;
  for(int k=0; k<mNDOFs1D; ++k)
    for(int j=0; j<mNDOFs1D; ++j)
      for(int i=0; i<mNDOFs1D; ++i, ++ii)
        mBasisCenter[ii] = V[i]*V[j]*V[k];
}

//------------------------------------------------------------------------------

// Function, which determines the data for this standard element.
void StandardElementClass::DataStandardElement(const InputParamClass *inputParam)
{
  // Set the polynomial degree for this standard element and the number of DOFs.
  mNPoly   = inputParam->mNPolySolDOFs;
  mNDOFs1D = mNPoly+1;
  mNDOFs   = mNDOFs1D*mNDOFs1D*mNDOFs1D;

  // Determine the padded values for the number of DOFs, both in 1D and 3D.
  mNDOFs1DPad = ((mNDOFs1D+vecLen1D-1)/vecLen1D)*vecLen1D;
  mNDOFsPad   = ((mNDOFs  +vecLen3D-1)/vecLen3D)*vecLen3D;

  // Determine the 1D location of the DOFs. The 3D location can be obtained
  // via a tensor product.
  LocationDOFs1D(mNDOFs1D, inputParam->mDOFLocation, mRDOFs1D);

  // Determine the 1D location of the integration points and their weights.
  // The 3D location can be obtained via a tensor product.
  LocationIntegration1D(inputParam->mNPolyIntRule);

  // Allocate the memory for the integration weights in 3D and determine them.
  // Note that the memory is allocation for the padded value of the number
  // of 3D integration points. The dummy values are initialized to zero
  // to avoid problems.
  mNIntegration    = mNIntegration1D*mNIntegration1D*mNIntegration1D;
  mNIntegrationPad = ((mNIntegration+vecLen3D-1)/vecLen3D)*vecLen3D;

  mIntWeights = (su2double *) AllocateMemory(mNIntegrationPad*sizeof(su2double));
  if( !mIntWeights )
    Terminate("StandardElementClass::DataStandardElement", __FILE__, __LINE__,
              "Memory allocation failure for mIntWeights.");

  int ii = 0;
  for(int k=0; k<mNIntegration1D; ++k)
    for(int j=0; j<mNIntegration1D; ++j)
      for(int i=0; i<mNIntegration1D; ++i, ++ii)
        mIntWeights[ii] = mIntegrationWeights1D[i]*mIntegrationWeights1D[j]
                        * mIntegrationWeights1D[k];

  for(int i=mNIntegration; i<mNIntegrationPad; ++i)
    mIntWeights[i] = zero;

  // Determine the, possibly padded, 1D Vandermonde matrix in the DOFs.
  VandermondeDOFs1D();

  // Determine the values of the 1D Legendre basis functions and its
  // derivatives in the 1D DOFs.
  ValuesLegendreBasisDOFs1D();

  // Determine the values of the 1D Legendre basis functions and its
  // derivatives in the 1D integration points.
  ValuesLegendreBasisInt1D();

  // Determine the weights of the 1D integration rule based on the DOFs.
  Weights1DIntegrationRuleDOFs();

  // Determine the values of the 3D basis functions in the center of
  // the element.
  BasisFunctionsCenter();

	// Determine the face data.
  DataStandardFace();

	// Determine the values of the 1D Lagrange basis and its derivatives 
	// in the integration points and 1D DOFs of the element.
	ValuesLagrangeBasis1D();

	// If an NSCBC is specified, compute the least-squares polynomial-correction
	// matrix and the two on the max boundary.
	if( inputParam->mNSCBC_Specified )
		AssemblePolynomialCorrectionOperators();
}

//------------------------------------------------------------------------------

// Function, which determines the values of the polynomial-correction
// operators, in case an NSCBC is specified.
void StandardElementClass::AssemblePolynomialCorrectionOperators(void)
{
	// Specify the number of points needed in this function.
	const int nInt1D     = mNIntegration1D;
	const int nInt1DPad  = mNIntegration1DPad;
	const int nDOFs1D    = mNDOFs1D;
	const int nDOFs1DPad = mNDOFs1DPad;
	const int nDOFs2D    = mNDOFs1D*mNDOFs1D;
	const int nInt2D     = mNIntegration2D;
	const int nInt2DPad  = mNIntegration2DPad;

	// Allocate the memory for mPolynomialCorrectionMatrix. Note that for 
	// efficiency reasons the number of integration points is padded.
	mPolynomialCorrectionMatrix = (su2double *) AllocateMemory(nInt2DPad*nInt2D*sizeof(su2double));

	// Check if memory allocation failed.
	if( !mPolynomialCorrectionMatrix )
		Terminate("StandardElementClass::AssemblePolynomialCorrectionOperators", __FILE__, __LINE__,
				      "Memory allocation failure for mPolynomialCorrectionMatrix");

	// Initialize these variables to zero, such that the padded values get a value.
	for(int i=0; i<nInt2D*nInt2DPad; ++i)
		mPolynomialCorrectionMatrix[i] = zero;

	// Form a temporary 2D array to store the 2D Lagrange interpolation matrix that 
	// converts nodal to integration points.
	std::vector< std::vector<su2double> > L2D_p2q;
	
	// Number of rows and columns used in L2D_p2q.
	const int nRow = nInt2D;
	const int nCol = nDOFs2D;

	// Initialize the 2D Lagrange interpolation function, nodal to integration points.
	L2D_p2q.resize(nRow);
	for(int i=0; i<nRow; i++)
		L2D_p2q[i].resize(nCol, zero);

	// Note, below uses Edwin's transposed definition which is actually the 
	// standard definition used in MATLAB notation. More specifically,
	// mLagrangeInt1DTranspose[iRow*nDOFs1DPad+iCol] == L1D_p2q[nRow][nCol], 
	//   where nRow: nInt1D and nCol: nDOFs1D.
	int iRow = 0;
	for(int i1=0; i1<nInt1D; i1++){
		for(int i2=0; i2<nInt1D; i2++){
			int iCol = 0;
			for(int j1=0; j1<nDOFs1D; j1++){
				for(int j2=0; j2<nDOFs1D; j2++){
					L2D_p2q[iRow][iCol++] = mLagrangeInt1DTranspose[i1*nDOFs1DPad+j1]
						                     *mLagrangeInt1DTranspose[i2*nDOFs1DPad+j2];
				}
			}
			iRow++;
		}
	}

	// External-coefficient differentiation matrix.
	std::vector< std::vector<su2double> > dL2D_p2q = L2D_p2q;
	// Factor in the MAX-face derivative coefficient. 
	// Note, the dr_max = -dr_min.
	for(int i=0; i<nRow; i++) 
		for(int j=0; j<nCol; j++) 
			dL2D_p2q[i][j] *= -mDerLagrangeDOFs1D[0];

	// Initialize least-squares matrix as a 1D vector, since we need to 
	// invert it using the InverseMatrix function. Note, make sure it is
	// initialized to zero.
	std::vector<su2double> LS(nDOFs2D*nDOFs2D, zero);

	// Compute the least-squares matrix.
	for(int i=0; i<nDOFs2D; i++){
		for(int j=0; j<nDOFs2D; j++){
			for(int k=0; k<nInt2D; k++){
				LS[i*nDOFs2D+j] += dL2D_p2q[k][i]*dL2D_p2q[k][j];
			}
		}
	}
	
	// Invert the least-squares matrix.
	InverseMatrix(nDOFs2D, LS);

	// Compute the product of LS*dL2D_p2q^T. Note, here transpose uses the 
	// standard definition also in MATLAB, unlike Edwin's transpose.
	for(int i=0; i<nDOFs2D; i++){
		for(int j=0; j<nInt2D; j++){
			su2double tmp = zero;
			for(int k=0; k<nDOFs2D; k++){
				tmp += LS[i*nDOFs2D+k]*dL2D_p2q[j][k];
			} 
			mPolynomialCorrectionMatrix[i*nInt2DPad+j] = tmp;
		}
	}

	// Copy the data from mPolynomialCorrectionMatrix to dL2D_p2q since it 
	// is no longer needed. Note, these are stored as [j][i] because of the
	// dimensions in dL2D_p2q.
	for(int i=0; i<nDOFs2D; i++)
		for(int j=0; j<nInt2D; j++)
			dL2D_p2q[j][i] = mPolynomialCorrectionMatrix[i*nInt2DPad+j];

	// Finalize the PC-matrix by multiplying from the left by the Lagrange
	// interpolation matrix from nodal to integration points in 2D.
	for(int i=0; i<nInt2D; i++){
		for(int j=0; j<nInt2D; j++){
			su2double tmp = zero;
			for(int k=0; k<nDOFs2D; k++){
				tmp += L2D_p2q[i][k]*dL2D_p2q[j][k]; 
			}
			mPolynomialCorrectionMatrix[i*nInt2DPad+j] = tmp;
		}
	}

	// Allocate the memory for mDerLagrangeInt1D_BaseInt and mDerLagrangeInt1D_BaseIntTranspose. 
	// Note that for efficiency reasons the number of integration points is padded.
	mDerLagrangeInt1D_BaseInt          = (su2double *) AllocateMemory(nInt1DPad*nInt1D*sizeof(su2double));
	mDerLagrangeInt1D_BaseIntTranspose = (su2double *) AllocateMemory(nInt1DPad*nInt1D*sizeof(su2double));

	// Check if memory allocation failed.
	if( !mDerLagrangeInt1D_BaseInt || !mDerLagrangeInt1D_BaseIntTranspose )
		Terminate("StandardElementClass::AssemblePolynomialCorrectionOperators", __FILE__, __LINE__,
				      "Memory allocation failure for mDerLagrangeInt1D_BaseInt and mDerLagrangeInt1D_BaseIntTranspose");

	// Initialize these variables to zero, such that the padded values get a value.
	for(int i=0; i<nInt1D*nInt1DPad; ++i)
		mDerLagrangeInt1D_BaseInt[i] = mDerLagrangeInt1D_BaseIntTranspose[i] = zero;

	// Compute the derivative of the Lagrange basis function on the 1D integration points using 
	// the integration points as the basis.
	LagrangianBasisFunctions(mRIntegration1D, mRIntegration1D, 
			                     NULL, mDerLagrangeInt1D_BaseInt,
													 NULL, NULL);

	// Set the relevant values of mDerLagrangeInt1D_BaseIntTranspose by taking the 
	// transpose of mDerLagrangeInt1D_BaseInt.
  for(int j=0; j<nInt1D; ++j)
  {
    for(int i=0; i<nInt1D; ++i)
    {
      const int ii  = i*nInt1DPad + j;
      const int iiT = j*nInt1DPad + i;

      mDerLagrangeInt1D_BaseIntTranspose[iiT] = mDerLagrangeInt1D_BaseInt[ii];
    }
  }
}

//------------------------------------------------------------------------------

// Function, which determines the values of the 1D Lagrange basis functions
// and its derivatives in the 1D interpolation points and 1D DOFs.
void StandardElementClass::ValuesLagrangeBasis1D(void)
{
  // Allocate the memory for mLagrangeInt1D and mDerLagrangeInt1D. Note that
  // for efficiency reasons the number of integration points is padded for these arrays.
  mLagrangeInt1D    = (su2double *) AllocateMemory(mNIntegration1DPad*mNDOFs1D*sizeof(su2double));
  mDerLagrangeInt1D = (su2double *) AllocateMemory(mNIntegration1DPad*mNDOFs1D*sizeof(su2double));

	// Check if memory allocation failed.
	if( !mLagrangeInt1D || !mDerLagrangeInt1D )
		Terminate("StandardElementClass::ValuesLagrangeBasis1D", __FILE__, __LINE__,
				      "Memory allocation failure for mLagrangeInt1D and mDerLagrangeInt1D");

	// Initialize these variables to zero, such that the padded values get a value.
	for(int i=0; i<mNIntegration1DPad*mNDOFs1D; ++i)
		mLagrangeInt1D[i] = mDerLagrangeInt1D[i] = zero;


	// Compute the Lagrange basis function and its derivate for the 1D integration points.
	LagrangianBasisFunctions(mRDOFs1D, mRIntegration1D, 
			                     mLagrangeInt1D, mDerLagrangeInt1D,
													 NULL, NULL);


  // Allocate the memory for mLagrangeInt1DTranspose and mDerLagrangeInt1DTranspose. 
	// Note that for efficiency reasons the number of DOFs is padded for these arrays.
  mLagrangeInt1DTranspose    = (su2double *) AllocateMemory(mNIntegration1D*mNDOFs1DPad*sizeof(su2double));
  mDerLagrangeInt1DTranspose = (su2double *) AllocateMemory(mNIntegration1D*mNDOFs1DPad*sizeof(su2double));

	// Check if memory allocation failed.
	if( !mLagrangeInt1DTranspose || !mDerLagrangeInt1DTranspose )
		Terminate("StandardElementClass::ValuesLagrangeBasis1D", __FILE__, __LINE__,
				      "Memory allocation failure for mLagrangeInt1DTranspose and mDerLagrangeInt1DTranspose");

	// Initialize these variables to zero, such that the padded values get a value.
	for(int i=0; i<mNIntegration1D*mNDOFs1DPad; ++i)
		mLagrangeInt1DTranspose[i] = mDerLagrangeInt1DTranspose[i] = zero;

	// Set the relevant values of mLagrangeInt1DTranspose and
  // mDerLagrangeInt1DTranspose by taking the transpose of mLagrangeInt1D
  // and mDerLagrangeInt1D, respectively.
  for(int j=0; j<mNIntegration1D; ++j)
  {
    for(int i=0; i<mNDOFs1D; ++i)
    {
      const int ii  = i*mNIntegration1DPad + j;
      const int iiT = j*mNDOFs1DPad + i;

      mLagrangeInt1DTranspose[iiT]    = mLagrangeInt1D[ii];
      mDerLagrangeInt1DTranspose[iiT] = mDerLagrangeInt1D[ii];
    }
  }

  // Allocate the memory for mDerLagrangeDOFs1D and its transpose mDerLagrangeDOFs1DTranspose. 
	// Note that for efficiency reasons the number of DOFs is padded for this array.
	mDerLagrangeDOFs1D          = (su2double *) AllocateMemory(mNDOFs1DPad*mNDOFs1D*sizeof(su2double));
	mDerLagrangeDOFs1DTranspose = (su2double *) AllocateMemory(mNDOFs1DPad*mNDOFs1D*sizeof(su2double));

	// Check if memory allocation failed.
	if( !mDerLagrangeDOFs1D || !mDerLagrangeDOFs1DTranspose )
		Terminate("StandardElementClass::ValuesLagrangeBasis1D", __FILE__, __LINE__,
				      "Memory allocation failure for mDerLagrangeDOFs1D and mDerLagrangeDOFs1DTranspose");

	// Initialize these variables to zero, such that the padded values get a value.
	for(int i=0; i<mNDOFs1DPad*mNDOFs1D; ++i)
		mDerLagrangeDOFs1D[i] = mDerLagrangeDOFs1DTranspose[i] = zero;

	// Compute the Lagrange basis function and its derivate for the 1D DOFs.
	LagrangianBasisFunctions(mRDOFs1D, mRDOFs1D, 
			                     NULL, mDerLagrangeDOFs1D,
													 NULL, NULL);

	// Set the relevant values of mDerLagrangeDOFs1DTranspose by taking the 
	// transpose of mDerLagrangeDOFs1D.
  for(int j=0; j<mNDOFs1D; ++j)
  {
    for(int i=0; i<mNDOFs1D; ++i)
    {
      const int ii  = i*mNDOFs1DPad + j;
      const int iiT = j*mNDOFs1DPad + i;

      mDerLagrangeDOFs1DTranspose[iiT] = mDerLagrangeDOFs1D[ii];
    }
  }
}

//------------------------------------------------------------------------------

// Function, which determines the standard element data needed for the
// computation of the face integrals.
void StandardElementClass::DataStandardFace(void)
{
  // Allocate the memory for the integration weights of the face.
  // Note that this array is padded to allow for vectorization.
  mNIntegration2D    = mNIntegration1D * mNIntegration1D;
  mNIntegration2DPad = ((mNIntegration2D+vecLen2D-1)/vecLen2D)*vecLen2D;

  mIntWeightsFace = (su2double *) AllocateMemory(mNIntegration2DPad*sizeof(su2double));
  if( !mIntWeightsFace )
    Terminate("StandardElementClass::DataStandardFace", __FILE__, __LINE__,
              "Memory allocation failure for mIntWeightsFace.");

  // Determine the integration weights of the face.
  // The dummy values are initialized to zero.
  int ii = 0;
  for(int j=0; j<mNIntegration1D; ++j)
    for(int i=0; i<mNIntegration1D; ++i, ++ii)
      mIntWeightsFace[ii] = mIntegrationWeights1D[i] * mIntegrationWeights1D[j];

  for(int i=mNIntegration2D; i<mNIntegration2DPad; ++i)
    mIntWeightsFace[i] = zero;

  // Allocate the memory for the value of the Legendre basis functions and its
  // derivatives on the min and max faces. Initialize the values to zero.
  mLegendreMinFace1D    = (su2double *) AllocateMemory(mNDOFs1DPad*sizeof(su2double));
  mLegendreMaxFace1D    = (su2double *) AllocateMemory(mNDOFs1DPad*sizeof(su2double));
  mDerLegendreMinFace1D = (su2double *) AllocateMemory(mNDOFs1DPad*sizeof(su2double));
  mDerLegendreMaxFace1D = (su2double *) AllocateMemory(mNDOFs1DPad*sizeof(su2double));

  if(!mLegendreMinFace1D || !mLegendreMaxFace1D)
    Terminate("StandardElementClass::DataStandardFace", __FILE__, __LINE__,
              "Memory allocation failure for mLegendreMinFace1D and mLegendreMaxFace1D");

  if(!mDerLegendreMinFace1D || !mDerLegendreMaxFace1D)
    Terminate("StandardElementClass::DataStandardFace", __FILE__, __LINE__,
              "Memory allocation failure for mDerLegendreMinFace1D and mDerLegendreMaxFace1D");

  for(int i=0; i<mNDOFs1DPad; ++i)
    mLegendreMinFace1D[i]    = mLegendreMaxFace1D[i]    =
    mDerLegendreMinFace1D[i] = mDerLegendreMaxFace1D[i] = zero;

  // Store the parametric coordinate of the min face in a vector.
  // This is needed for the calls to Vandermonde1D and GradVandermonde1D.
  std::vector<su2double> rFace(1);
  rFace[0] = -one;

  // Determine the value of the Legendre basis functions on the min boundary.
  std::vector<su2double> V(mNDOFs1D);
  Vandermonde1D(mNDOFs1D, rFace, V);

  for(int k=0; k<mNDOFs1D; ++k)
    mLegendreMinFace1D[k] = V[k];

  // Determine the gradients of the Legendre basis functions on the min boundary.
  GradVandermonde1D(mNDOFs1D, rFace, V);

  for(int k=0; k<mNDOFs1D; ++k)
    mDerLegendreMinFace1D[k] = V[k];

  // Set the parametric coordinate of the max face 
  rFace[0] = one;

  // Determine the value of the Legendre basis functions on the max boundary.
  Vandermonde1D(mNDOFs1D, rFace, V);

  for(int k=0; k<mNDOFs1D; ++k)
    mLegendreMaxFace1D[k] = V[k];

  // Determine the gradients of the Legendre basis functions on the max boundary.
  GradVandermonde1D(mNDOFs1D, rFace, V);

  for(int k=0; k<mNDOFs1D; ++k)
    mDerLegendreMaxFace1D[k] = V[k];
}

//------------------------------------------------------------------------------

// Function, which determines the 1D Lagrangian basis functions for the
// given 1D location of the DOFs and points.
void StandardElementClass::LagrangianBasisFunctions(const std::vector<su2double> &rDOFs1D,
                                                    const std::vector<su2double> &rPoints1D,
                                                    su2double                    *lagrangePoints1D,
                                                    su2double                    *derLagrangePoints1D,
                                                    su2double                    *derLagrangeMinFace1D,
                                                    su2double                    *derLagrangeMaxFace1D)
{
  // Determine the number of DOFs associated with rDOFs1D.
  const int nDOFs1D = (int) rDOFs1D.size();

  // Determine the inverse of the 1D Vandermonde matrix in the 1D DOFs.
  std::vector<su2double> VInv(nDOFs1D*nDOFs1D);

  Vandermonde1D(nDOFs1D, rDOFs1D, VInv);
  InverseMatrix(nDOFs1D, VInv);

  // Determine the number of points associated with rPoints1D and its padded value.
  const int nPoints1D    = (int) rPoints1D.size();
  const int nPoints1DPad = ((nPoints1D+vecLen1D-1)/vecLen1D)*vecLen1D;

  // Determine the 1D Vandermonde matrix in the 1D points.
  std::vector<su2double> V(nDOFs1D*nPoints1D);
  Vandermonde1D(nDOFs1D, rPoints1D, V);

  // Check if the Lagrange polynomials in the points must be computed.
  if( lagrangePoints1D )
  {
    // The actual values for the Lagrangian basis functions in the
    // points are obtained by carrying out the matrix product V*VInv.
    for(int j=0; j<nDOFs1D; ++j)
    {
      for(int i=0; i<nPoints1D; ++i)
      {
        const int ii = j*nPoints1DPad + i;

        for(int k=0; k<nDOFs1D; ++k)
        {
          const int indV    = k*nPoints1D + i;
          const int indVInv = j*nDOFs1D + k;

          lagrangePoints1D[ii] += V[indV]*VInv[indVInv];
        }
      }
    }

    // Check if the Lagrangian basis functions are computed correctly.
    // The check is that for every point the basis functions should sum up to 1.
    for(int j=0; j<nPoints1D; ++j)
    {
      su2double val = zero;
      for(int i=0; i<nDOFs1D; ++i)
      {
        const int ii = i*nPoints1DPad + j;
        val += lagrangePoints1D[ii];
      }

      // Check if the Lagrangian basis functions sum up to 1.
      if(FABS(val-one) > epsThreshold)
        Terminate("StandardElementClass::LagrangianBasisFunctions", __FILE__,
                  __LINE__, "Lagrangian basis function do not sum up to 1.");
    }
  }

  // Check if the derivatives of the Lagrange polynomials in the points
  // must be computed.
  if( derLagrangePoints1D )
  {
    // Determine the gradients of the 1D Vandermonde matrix in the 1D points.
    GradVandermonde1D(nDOFs1D, rPoints1D, V);

    // The actual values for the derivatives of the Lagrangian basis functions in
    // the points are obtained by carrying out the matrix product V*VInv.
    for(int j=0; j<nDOFs1D; ++j)
    {
      for(int i=0; i<nPoints1D; ++i)
      {
        const int ii = j*nPoints1DPad + i;

        for(int k=0; k<nDOFs1D; ++k)
        {
          const int indV    = k*nPoints1D + i;
          const int indVInv = j*nDOFs1D + k;

          derLagrangePoints1D[ii] += V[indV]*VInv[indVInv];
        }
      }
    }

    // Check if the Lagrangian basis functions are computed correctly.
    // The check is that for every point the derivatives should sum up to 0.
    for(int j=0; j<nPoints1D; ++j)
    {
      su2double valDer = zero;
      for(int i=0; i<nDOFs1D; ++i)
      {
        const int ii = i*nPoints1DPad + j;
        valDer += derLagrangePoints1D[ii];
      }

      // Check if the derivatives of the Lagrangian basis functions sum up to 0.
      if(FABS(valDer) > epsThreshold)
        Terminate("StandardElementClass::LagrangianBasisFunctions", __FILE__, __LINE__,
                  "Derivatives of the Lagrangian basis function do not sum up to 0.");
    }
  }

  // Check if the derivatives of the Lagrangian basis functions must
  // be computed at the min and max faces.
  if(derLagrangeMinFace1D && derLagrangeMaxFace1D)
  {
    // Store the parametric coordinate of the min face in a vector.
    // This is needed for the call to GradVandermonde1D.
    std::vector<su2double> rFace(1);
    rFace[0] = -one;

    // Determine the gradients of the 1D Vandermonde matrix in rFace.
    V.resize(nDOFs1D);
    GradVandermonde1D(nDOFs1D, rFace, V);

    // The actual values for the derivatives of the Lagrangian basis functions at
    // the face are obtained by carrying out the matrix product V*VInv.
    for(int j=0; j<nDOFs1D; ++j)
    {
      for(int k=0; k<nDOFs1D; ++k)
      {
        const int indVInv = j*nDOFs1D + k;
        derLagrangeMinFace1D[j] += V[k]*VInv[indVInv];
      }
    }

    // Also set the parametric coordinate of the max face and determine
    // the gradients of the Vandermonde matrix.
    rFace[0] = one;
    GradVandermonde1D(nDOFs1D, rFace, V);

    // Compute the derivatives of the Lagrangian basis functions at the face
    // by carrying out the matrix product V*VInv.
    for(int j=0; j<nDOFs1D; ++j)
    {
      for(int k=0; k<nDOFs1D; ++k)
      {
        const int indVInv = j*nDOFs1D + k;
        derLagrangeMaxFace1D[j] += V[k]*VInv[indVInv];
      }
    }

    // Check if the derivatives of the Lagrangian basis functions are computed
    // correctly.  The check is that the derivatives should sum up to 0.
    su2double valDerMin = zero, valDerMax = zero;
    for(int i=0; i<nDOFs1D; ++i)
    {
      valDerMin += derLagrangeMinFace1D[i];
      valDerMax += derLagrangeMaxFace1D[i];
    }

    if((FABS(valDerMin) > epsThreshold) || (FABS(valDerMax) > epsThreshold))
      Terminate("StandardElementClass::LagrangianBasisFunctions", __FILE__, __LINE__,
                "Derivatives of the Lagrangian basis function at the face do not sum up to 0.");
  }
}

//------------------------------------------------------------------------------

// Function, which determines the 1D location of the integration points.
void StandardElementClass::LocationIntegration1D(const int nPolyIntRule)
{
  // Determine the number of 1D integration points and its padded value.
  mNIntegration1D    = nPolyIntRule/2 + 1;
  mNIntegration1DPad = ((mNIntegration1D+vecLen1D-1)/vecLen1D)*vecLen1D;

  // Allocate the memory for the 1D integration points and weights.
  // This is not done with mRIntegration1D and mIntegrationWeights1D,
  // because su2double may not be equal to double.
  std::vector<double> rIntegration1D(mNIntegration1D);
  std::vector<double> integrationWeights1D(mNIntegration1D);

  // Determine the location of the integration points and their weights.
  CGaussJacobiQuadrature GaussJacobi;
  GaussJacobi.GetQuadraturePoints(0.0, 0.0, -1.0, 1.0, rIntegration1D,
                                  integrationWeights1D);

  // Allocate the memory for the 1D integration points and weights and
  // copy the data from rIntegration1D and integrationWeights1D.
  mRIntegration1D.resize(mNIntegration1D);
  mIntegrationWeights1D.resize(mNIntegration1D);

  for(int i=0; i<mNIntegration1D; ++i)
  {
    mRIntegration1D[i]       = (su2double) rIntegration1D[i];
    mIntegrationWeights1D[i] = (su2double) integrationWeights1D[i];
  }
}

//------------------------------------------------------------------------------

// Function, which determines the values of the 1D Legendre basis functions
// and its derivatives in the 1D DOFs.
void StandardElementClass::ValuesLegendreBasisDOFs1D(void)
{
  // Allocate the memory for mLegendreDOFs1D and mDerLegendreDOFs1D. Note that
  // for efficiency reasons the number of DOFs is padded for this array.
  // Initialize the values of these arrays to zero.
  mLegendreDOFs1D    = (su2double *) AllocateMemory(mNDOFs1DPad*mNDOFs1D*sizeof(su2double));
  mDerLegendreDOFs1D = (su2double *) AllocateMemory(mNDOFs1DPad*mNDOFs1D*sizeof(su2double));

  for(int i=0; i<(mNDOFs1DPad*mNDOFs1D); ++i)
    mLegendreDOFs1D[i] = mDerLegendreDOFs1D[i] = zero;

  // Determine the 1D Vandermonde matrix and its derivatives in the 1D DOFs.
  std::vector<su2double> V(mNDOFs1D*mNDOFs1D), VDer(mNDOFs1D*mNDOFs1D);

  Vandermonde1D(mNDOFs1D, mRDOFs1D, V);
  GradVandermonde1D(mNDOFs1D, mRDOFs1D, VDer);

  // Copy the data from V and VDer into mLegendreDOFs1D and mDerLegendreDOFs1D.
  for(int j=0; j<mNDOFs1D; ++j)
  {
    for(int i=0; i<mNDOFs1D; ++i)
    {
      const int ii   = j*mNDOFs1DPad + i;
      const int indV = j*mNDOFs1D + i;

      mLegendreDOFs1D[ii]    = V[indV];
      mDerLegendreDOFs1D[ii] = VDer[indV];
    }
  }
}

//------------------------------------------------------------------------------

// Function, which determines the values of the 1D Legendre basis functions
// and its derivatives in the 1D integration points.
void StandardElementClass::ValuesLegendreBasisInt1D(void)
{
  // Allocate the memory for mLegendreInt1D and mDerLegendreInt1D. Note that for
  // efficiency reasons the number of integration points is padded for these
  // arrays. Initialize their values to zero, such that also the padded values
  // are initialized.
  mLegendreInt1D    = (su2double *) AllocateMemory(mNIntegration1DPad*mNDOFs1D*sizeof(su2double));
  mDerLegendreInt1D = (su2double *) AllocateMemory(mNIntegration1DPad*mNDOFs1D*sizeof(su2double));

  for(int i=0; i<(mNIntegration1DPad*mNDOFs1D); ++i)
    mLegendreInt1D[i] = mDerLegendreInt1D[i] = zero;

  // Determine the 1D Vandermonde matrix in the 1D integration points.
  std::vector<su2double> V(mNDOFs1D*mNIntegration1D);
  Vandermonde1D(mNDOFs1D, mRIntegration1D, V);

  // Copy the data from V into mLegendreInt1D.
  for(int j=0; j<mNDOFs1D; ++j)
  {
    for(int i=0; i<mNIntegration1D; ++i)
    {
      const int ii   = j*mNIntegration1DPad + i;
      const int indV = j*mNIntegration1D    + i;

      mLegendreInt1D[ii] = V[indV];
    }
  }

  // Determine the gradients of the 1D Vandermonde matrix in the
  // 1D integration points.
  GradVandermonde1D(mNDOFs1D, mRIntegration1D, V);

  // Copy the data from V into mDerLegendreInt1D.
  for(int j=0; j<mNDOFs1D; ++j)
  {
    for(int i=0; i<mNIntegration1D; ++i)
    {
      const int ii   = j*mNIntegration1DPad + i;
      const int indV = j*mNIntegration1D    + i;

      mDerLegendreInt1D[ii] = V[indV];
    }
  }

  // Allocate the memory for mLegendreInt1DTranspose and mDerLegendreInt1DTranspose.
  // Note that for efficiency reasons the number of DOFs is padded for these arrays.
  // Initialize the values of these arrays to zero, such that also the padded
  // values are initialized.
  mLegendreInt1DTranspose    = (su2double *) AllocateMemory(mNIntegration1D*mNDOFs1DPad*sizeof(su2double));
  mDerLegendreInt1DTranspose = (su2double *) AllocateMemory(mNIntegration1D*mNDOFs1DPad*sizeof(su2double));

  for(int i=0; i<(mNIntegration1D*mNDOFs1DPad); ++i)
    mLegendreInt1DTranspose[i] = mDerLegendreInt1DTranspose[i] = zero;

  // Set the relevant values of mLegendreInt1DTranspose and
  // mDerLegendreInt1DTranspose by taking the transpose of mLegendreInt1D
  // and mDerLegendreInt1D, respectively.
  for(int j=0; j<mNIntegration1D; ++j)
  {
    for(int i=0; i<mNDOFs1D; ++i)
    {
      const int ii  = i*mNIntegration1DPad + j;
      const int iiT = j*mNDOFs1DPad + i;

      mLegendreInt1DTranspose[iiT]    = mLegendreInt1D[ii];
      mDerLegendreInt1DTranspose[iiT] = mDerLegendreInt1D[ii];
    }
  }
}

//------------------------------------------------------------------------------

// Function, which determines the values of the, possibly padded, 1D
// Vandermonde matrix in the DOFs.
void StandardElementClass::VandermondeDOFs1D(void)
{
  // Determine the 1D Vandermonde matrix in the 1D DOFs.
  std::vector<su2double> V(mNDOFs1D*mNDOFs1D);

  Vandermonde1D(mNDOFs1D, mRDOFs1D, V);

  // Allocate the memory to store the 1D Vandermonde matrix and, its inverse and its
  // transpose. Take padding into account for the first dimension. Initialize the
  // data to zero, such that also the padded values are initialized.
  mVandermonde1D          = (su2double *) AllocateMemory(mNDOFs1DPad*mNDOFs1D*sizeof(su2double));
  mVandermonde1DInverse   = (su2double *) AllocateMemory(mNDOFs1DPad*mNDOFs1D*sizeof(su2double));
  mVandermonde1DTranspose = (su2double *) AllocateMemory(mNDOFs1DPad*mNDOFs1D*sizeof(su2double));

  if(!mVandermonde1D || !mVandermonde1DInverse || !mVandermonde1DTranspose)
    Terminate("StandardElementClass::VandermondeDOFs1D", __FILE__, __LINE__,
              "Memory allocation failure for the Vandermonde matrices");

  for(int i=0; i<(mNDOFs1DPad*mNDOFs1D); ++i)
    mVandermonde1D[i] = mVandermonde1DInverse[i] = mVandermonde1DTranspose[i] = zero;

  // Copy the values from V into mVandermonde1D and mVandermonde1DTranspose.
  // Take the padding into account.
  for(int j=0; j<mNDOFs1D; ++j)
  {
    for(int i=0; i<mNDOFs1D; ++i)
    {
      const int ii  = j*mNDOFs1DPad + i;
      const int iiT = i*mNDOFs1DPad + j;
      const int jj  = j*mNDOFs1D    + i;

      mVandermonde1D[ii]           = V[jj];
      mVandermonde1DTranspose[iiT] = V[jj];
    }
  }

  // Compute the inverse of the Vandermonde matrix.
  InverseMatrix(mNDOFs1D, V);

  // Copy the values of V into mVandermonde1DInverse.
  // Take the padding into account.
  for(int j=0; j<mNDOFs1D; ++j)
  {
    for(int i=0; i<mNDOFs1D; ++i)
    {
      const int ii = j*mNDOFs1DPad + i;
      const int jj = j*mNDOFs1D    + i;

      mVandermonde1DInverse[ii] = V[jj];
    }
  }
}

//------------------------------------------------------------------------------

// Function, which determine the weights of the 1D integration rule
// based on the DOFs.
void StandardElementClass::Weights1DIntegrationRuleDOFs(void)
{
  // Determine the inverse of the 1D Vandermonde matrix in the 1D DOFs.
  std::vector<su2double> VInv(mNDOFs1D*mNDOFs1D);

  Vandermonde1D(mNDOFs1D, mRDOFs1D, VInv);
  InverseMatrix(mNDOFs1D, VInv);

  // Determine the 1D Vandermonde matrix in the 1D integration points.
  std::vector<su2double> V(mNDOFs1D*mNIntegration1D);
  Vandermonde1D(mNDOFs1D, mRIntegration1D, V);

  // Allocate the memory for the 1D Lagrangian basis functions in the
  // integration points. For consistency reasons the number of integration
  // points is padded.
  std::vector<su2double> lagrangeInt1D(mNIntegration1DPad*mNDOFs1D, zero);

  // The actual values for the Lagrangian basis functions in the integration
  // points is obtained by carrying out the matrix product V*VInv.
  for(int j=0; j<mNDOFs1D; ++j)
  {
    for(int i=0; i<mNIntegration1D; ++i)
    {
      const int ii = j*mNIntegration1DPad + i;

      for(int k=0; k<mNDOFs1D; ++k)
      {
        const int indV    = k*mNIntegration1D + i;
        const int indVInv = j*mNDOFs1D + k;

        lagrangeInt1D[ii] += V[indV]*VInv[indVInv];
      }
    }
  }

  // Initialize the integration weights for the 1D integration rule based
  // on the DOFs to zero.
  mIntegrationWeights1DDOFs.resize(mNDOFs1D, zero);

  // Loop over the 1D DOFs and compute the integration weights of the
  // integration rule based on the DOFs.
  for(int j=0; j<mNDOFs1D; ++j)
  {
    for(int i=0; i<mNIntegration1D; ++i)
    {
      const int ii = j*mNIntegration1DPad + i;
      mIntegrationWeights1DDOFs[j] += mIntegrationWeights1D[i]
                                    * lagrangeInt1D[ii];
    }
  }
}
