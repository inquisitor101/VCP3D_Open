//------------------------------------------------------------------------------
// File, which contains the implementation of the member functions of the
// class ElementClass.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"

//-----------------------------------------------------------------------------

// Local class, which stores the information of a donor element for the
// exchange points of the wall model.
class CExchangeInfoWMClass
{
public:
  //--------------------------------------------------
  // Constructor and destructor.
  //--------------------------------------------------

  // Constructor, nothing to be done.
  CExchangeInfoWMClass(){}

  // Destructor, nothing to be done.
  ~CExchangeInfoWMClass(){}

  // Copy constructor.
  CExchangeInfoWMClass(const CExchangeInfoWMClass &other){Copy(other);}

  //--------------------------------------------------
  // Operators.
  //--------------------------------------------------

  // Assignment operator.
  CExchangeInfoWMClass& operator=(const CExchangeInfoWMClass &other)
  {
    Copy(other);
    return (*this);
  }

  // Less than operator, needed for the sorting.
  bool operator<(const CExchangeInfoWMClass &other) const
  {
    if(mK     != other.mK)     return (mK     < other.mK);
    if(mJ     != other.mJ)     return (mJ     < other.mJ);
    if(mI     != other.mI)     return (mI     < other.mI);
    if(mIntID != other.mIntID) return (mIntID < other.mIntID);

    return false; // Objects are considered identical.
  }

  //--------------------------------------
  // Public member variables.
  //--------------------------------------

  // Local element indices of the donor element.
  int mK, mJ, mI;

  // ID of the corresponding integration point.
  int mIntID;

  // Parametric coordinates of the exchange point in the donor element.
  su2double mParCoor[3];

private:
  //--------------------------------------
  // Private member functions.
  //--------------------------------------

  // Copy function.
  void Copy(const CExchangeInfoWMClass &other)
  {
    mK     = other.mK;
    mJ     = other.mJ;
    mI     = other.mI;
    mIntID = other.mIntID;

    mParCoor[0] = other.mParCoor[0];
    mParCoor[1] = other.mParCoor[1];
    mParCoor[2] = other.mParCoor[2];
  }
};

//------------------------------------------------------------------------------

// Overloaded constructor. Initialize the member variables.
ElementClass::ElementClass(const InputParamClass      *inputParam,
                           const StandardElementClass *standardHex,
                           const int                  *globInd,
                           const ENUM_ELEM_TYPE       elemType)
{
  // Initialize the local indices to 0. These will be overwritten later
  // for owned elements.
  mLocalInd[0] = mLocalInd[1] = mLocalInd[2] = 0;

  // Copy the element type and the global indices.
  mElemType = elemType;

  mGlobalInd[0] = globInd[0];
  mGlobalInd[1] = globInd[1];
  mGlobalInd[2] = globInd[2];

  // Initialize the periodic transformations to zero.
  mTransIMin[0] = mTransIMin[1] = mTransIMin[2] = zero;
  mTransIMax[0] = mTransIMax[1] = mTransIMax[2] = zero;
  mTransJMin[0] = mTransJMin[1] = mTransJMin[2] = zero;
  mTransJMax[0] = mTransJMax[1] = mTransJMax[2] = zero;
  mTransKMin[0] = mTransKMin[1] = mTransKMin[2] = zero;
  mTransKMax[0] = mTransKMax[1] = mTransKMax[2] = zero;

  // Initialize the pointers to NULL.
  mBCIMin = NULL;
  mBCIMax = NULL;
  mBCJMin = NULL;
  mBCJMax = NULL;
  mBCKMin = NULL;
  mBCKMax = NULL;

  mAveEddyVis = NULL;

  // Check if this is an internal element.
  if(elemType == INTERNAL_ELEMENT)
  {
    // Determine the number of bytes that must be allocated by
    // each of the member variables.
    const size_t nBytes     = standardHex->mNDOFsPad * sizeof(su2double);
    const size_t nBytesDUdV = standardHex->mNIntegrationPad * sizeof(su2double);

    // Allocate the memory for the arrays that must always be allocated.
    mSol.resize(nVar, NULL);
    mSolOld.resize(nVar, NULL);
    mRes.resize(nVar, NULL);

    for(int i=0; i<nVar; ++i)
    {
      mSol[i]     = (su2double *) AllocateMemory(nBytes);
      mSolOld[i]  = (su2double *) AllocateMemory(nBytes);
      mRes[i]     = (su2double *) AllocateMemory(nBytes);

      if(!mSol[i] || !mSolOld[i] || !mRes[i])
        Terminate("ElementClass::ElementClass", __FILE__, __LINE__,
                  "Memory allocation failure");
    }

    // Test if the average solution must be computed. If so, allocate
    // the appropriate memory, which is initialized to zero.
    if( inputParam->mComputeAverageSolution )
    {
      mAvePrim.resize(nVar, NULL);
      mAveVelProd.resize(6, NULL);

      for(int i=0; i<nVar; ++i)
      {
        mAvePrim[i] = (su2double *) AllocateMemory(nBytes);
        if( !mAvePrim[i] )
          Terminate("ElementClass::ElementClass", __FILE__, __LINE__,
                    "Memory allocation failure for mAvePrim");
#pragma omp simd
        for(int j=0; j<standardHex->mNDOFsPad; ++j)
          mAvePrim[i][j] = zero;
      }

      for(int i=0; i<6; ++i)
      {
        mAveVelProd[i] = (su2double *) AllocateMemory(nBytes);
        if( !mAveVelProd[i] )
          Terminate("ElementClass::ElementClass", __FILE__, __LINE__,
                    "Memory allocation failure for mAveVelProd");
#pragma omp simd
        for(int j=0; j<standardHex->mNDOFsPad; ++j)
          mAveVelProd[i][j] = zero;
      }

      mAveEddyVis = (su2double *) AllocateMemory(nBytes);
      if( !mAveEddyVis )
        Terminate("ElementClass::ElementClass", __FILE__, __LINE__,
                  "Memory allocation failure for mAveEddyVis");

#pragma omp simd
      for(int j=0; j<standardHex->mNDOFsPad; ++j)
        mAveEddyVis[j] = zero;
    }

    // Test if the memory for the transformation matrices dUdV in the
    // integration points must be allocated.
    if(inputParam->mFEMVariables == ENTROPY_VARIABLES)
    {
      const int sizeDUdV = nVar*(nVar+1)/2;
      mDUdVInt.resize(sizeDUdV);

      for(int i=0; i<sizeDUdV; ++i)
      {
        mDUdVInt[i] = (su2double *) AllocateMemory(nBytesDUdV);
        if( !mDUdVInt[i] )
          Terminate("ElementClass::ElementClass", __FILE__, __LINE__,
                    "Memory allocation failure for mDUdVInt");
      }
    }
  }
}

//------------------------------------------------------------------------------

// Destructor. Release the memory.
ElementClass::~ElementClass()
{
  for(unsigned long i=0; i<mCoorNodalGridDOFs.size(); ++i)
    FreeMemory((void **) &mCoorNodalGridDOFs[i]);

  for(unsigned long i=0; i<mCoorNodalSolDOFs.size(); ++i)
    FreeMemory((void **) &mCoorNodalSolDOFs[i]);

  for(unsigned long i=0; i<mSurfMetricIntIMin.size(); ++i)
    FreeMemory((void **) &mSurfMetricIntIMin[i]);

  for(unsigned long i=0; i<mSurfMetricIntIMax.size(); ++i)
    FreeMemory((void **) &mSurfMetricIntIMax[i]);

  for(unsigned long i=0; i<mSurfMetricIntJMin.size(); ++i)
    FreeMemory((void **) &mSurfMetricIntJMin[i]);

  for(unsigned long i=0; i<mSurfMetricIntJMax.size(); ++i)
    FreeMemory((void **) &mSurfMetricIntJMax[i]);

  for(unsigned long i=0; i<mSurfMetricIntKMin.size(); ++i)
    FreeMemory((void **) &mSurfMetricIntKMin[i]);

  for(unsigned long i=0; i<mSurfMetricIntKMax.size(); ++i)
    FreeMemory((void **) &mSurfMetricIntKMax[i]);

  for(unsigned long i=0; i<mVolMetricInt.size(); ++i)
    FreeMemory((void **) &mVolMetricInt[i]);

  for(unsigned long i=0; i<mVolMetricSolDOFs.size(); ++i)
    FreeMemory((void **) &mVolMetricSolDOFs[i]);

  for(unsigned long i=0; i<mSol.size(); ++i)
    FreeMemory((void **) &mSol[i]);

  for(unsigned long i=0; i<mAvePrim.size(); ++i)
    FreeMemory((void **) &mAvePrim[i]);

  for(unsigned long i=0; i<mAveVelProd.size(); ++i)
    FreeMemory((void **) &mAveVelProd[i]);

  FreeMemory((void **) &mAveEddyVis);

  for(unsigned long i=0; i<mSolOld.size(); ++i)
    FreeMemory((void **) &mSolOld[i]);

  for(unsigned long i=0; i<mRes.size(); ++i)
    FreeMemory((void **) &mRes[i]);

  for(unsigned long i=0; i<mResIMin.size(); ++i)
    FreeMemory((void **) &mResIMin[i]);

  for(unsigned long i=0; i<mResJMin.size(); ++i)
    FreeMemory((void **) &mResJMin[i]);

  for(unsigned long i=0; i<mResKMin.size(); ++i)
    FreeMemory((void **) &mResKMin[i]);

  for(unsigned long i=0; i<mDUdVInt.size(); ++i)
    FreeMemory((void **) &mDUdVInt[i]);

  for(unsigned long i=0; i<mPrescribedDataIMin.size(); ++i)
    FreeMemory((void **) &mPrescribedDataIMin[i]);

  for(unsigned long i=0; i<mPrescribedDataIMax.size(); ++i)
    FreeMemory((void **) &mPrescribedDataIMax[i]);

  for(unsigned long i=0; i<mPrescribedDataJMin.size(); ++i)
    FreeMemory((void **) &mPrescribedDataJMin[i]);

  for(unsigned long i=0; i<mPrescribedDataJMax.size(); ++i)
    FreeMemory((void **) &mPrescribedDataJMax[i]);

  for(unsigned long i=0; i<mPrescribedDataKMin.size(); ++i)
    FreeMemory((void **) &mPrescribedDataKMin[i]);

  for(unsigned long i=0; i<mPrescribedDataKMax.size(); ++i)
    FreeMemory((void **) &mPrescribedDataKMax[i]);
}

//------------------------------------------------------------------------------

// Function, which computes the surface area of each (physical) boundary
void ElementClass::ComputeSurfaceAreaBoundary(const InputParamClass      *inputParam,
																							const StandardElementClass *standardHex,
																							su2double                  *surfArea)
{
	// Surface area per boundary. 
	for(unsigned short iBoundary=0; iBoundary<6; iBoundary++)
	{
		// Boundary and metric pointer.
		SubfaceBaseClass *BC = NULL;
		su2double **metric   = NULL;
		switch(iBoundary){
			case(0): BC = mBCIMin; metric = mSurfMetricIntIMin.data() ; break;
			case(1): BC = mBCIMax; metric = mSurfMetricIntIMax.data() ; break;
			case(2): BC = mBCJMin; metric = mSurfMetricIntJMin.data() ; break;
			case(3): BC = mBCJMax; metric = mSurfMetricIntJMax.data() ; break;
			case(4): BC = mBCKMin; metric = mSurfMetricIntKMin.data() ; break;
			case(5): BC = mBCKMax; metric = mSurfMetricIntKMax.data() ; break;
			default: break;
		}
		
		if(BC != NULL && BC->mBoundaryID == iBoundary)
			for(int l=0; l<standardHex->mNIntegration2D; l++)
				surfArea[iBoundary] += standardHex->mIntWeightsFace[l]*metric[3][l];
	}
}

//------------------------------------------------------------------------------

// Function, which sets the tuning parameters required for a NSCBC.
void ElementClass::TuningParamNSCBC(const InputParamClass 		 *inputParam,
																		const StandardElementClass *standardHex)
{
	// Operate on a per-face basis.
	for(unsigned short iBoundary=0; iBoundary<6; iBoundary++)
	{
		// Process required pointers.
		SubfaceBaseClass *BC = NULL;	
		su2double **metric   = NULL;
		su2double dirSign    = one;
		switch(iBoundary){
			case(0): BC = mBCIMin; metric = mSurfMetricIntIMin.data(); dirSign = -one; break;
			case(1): BC = mBCIMax; metric = mSurfMetricIntIMax.data(); dirSign =  one; break;
			case(2): BC = mBCJMin; metric = mSurfMetricIntJMin.data(); dirSign = -one; break;
			case(3): BC = mBCJMax; metric = mSurfMetricIntJMax.data(); dirSign =  one; break;
			case(4): BC = mBCKMin; metric = mSurfMetricIntKMin.data(); dirSign = -one; break;
			case(5): BC = mBCKMax; metric = mSurfMetricIntKMax.data(); dirSign =  one; break;
			default: break;
		}

		// Assign values, if NSCBC iBoundary exists.
		if(  BC != NULL && BC->mBoundaryID == iBoundary) 
			if( BC->GetTypeBoundaryPrescribed() == BC_OUTFLOW_CHARACTERISTIC 
				  ||
				  BC->GetTypeBoundaryPrescribed() == BC_INFLOW_CHARACTERISTIC )
		{
			// Configure boundary.
			BC->ConfigureParamNSCBC(inputParam, standardHex, metric, dirSign);
		}
	}
}

//------------------------------------------------------------------------------

// Function, which configures the required averaged Mach on each NSCBC boundary.
void ElementClass::AverageMachNSCBC(const InputParamClass       *inputParam,
																	  const StandardElementClass  *standardHex,
																	  su2double									 **workArray,
																	  su2double 									*LocalMachWeighted)
{
	// Easier storage of the number of DOFs and integration points.
  const int nDOFs1D  = standardHex->mNDOFs1D;
  const int nInt1D   = standardHex->mNIntegration1D;
  const int nInt     = standardHex->mNIntegration2D;

  // Set the pointers for the left and right solution 
  // in the integration points of the face.
  su2double **sol = workArray;

	// Function pointer to relevant tensor-product.
	void (*TensorProduct)(const int M, const int N, const int K,
                        su2double *A,
                        su2double *ADer,
                        su2double *AFace,
                        su2double *ADerFace,
                        su2double **B,
                        su2double **C,
                        su2double **CDerX,
                        su2double **CDerY,
                        su2double **CDerZ);

	// Class pointer to relevant boundary face.
	SubfaceBaseClass *BC              = NULL;
	su2double       **metric          = NULL;
	su2double        *legendreFace    = NULL;
	su2double        *derLegendreFace = NULL;
	su2double         dirSign         = one;

	for(unsigned short iBoundary=0; iBoundary<6; iBoundary++)
	{
		// Process required pointers.
		switch(iBoundary){
			case(0):
			{ BC = mBCIMin; dirSign = -one;
				metric          = mSurfMetricIntIMin.data();
				legendreFace    = standardHex->mLegendreMinFace1D;
				derLegendreFace = standardHex->mDerLegendreMinFace1D;
				TensorProduct   = TensorProductSolAndGradIFace; 
				break;
			}
			case(1): 
			{
				BC = mBCIMax; dirSign =  one;
				metric 					= mSurfMetricIntIMax.data(); 
				legendreFace    = standardHex->mLegendreMaxFace1D;
				derLegendreFace = standardHex->mDerLegendreMaxFace1D;
				TensorProduct 	= TensorProductSolAndGradIFace; 
				break;
			}
			case(2): 
			{
				BC = mBCJMin; dirSign = -one;
				metric 					= mSurfMetricIntJMin.data();
				legendreFace    = standardHex->mLegendreMinFace1D;
				derLegendreFace = standardHex->mDerLegendreMinFace1D;
				TensorProduct 	= TensorProductSolAndGradJFace; 
				break;
			}
			case(3): 
			{
				BC = mBCJMax; dirSign =  one;
				metric 					= mSurfMetricIntJMax.data();
				legendreFace    = standardHex->mLegendreMaxFace1D;
				derLegendreFace = standardHex->mDerLegendreMaxFace1D;
				TensorProduct 	= TensorProductSolAndGradJFace; 
				break;
			}
			case(4): 
			{
				BC = mBCKMin; dirSign = -one;
				metric 					= mSurfMetricIntKMin.data();
				legendreFace    = standardHex->mLegendreMinFace1D;
				derLegendreFace = standardHex->mDerLegendreMinFace1D;
				TensorProduct 	= TensorProductSolAndGradKFace; 
				break;
			}
			case(5): 
			{
				BC = mBCKMax; dirSign =  one;
				metric 					= mSurfMetricIntKMax.data();
				legendreFace    = standardHex->mLegendreMaxFace1D;
				derLegendreFace = standardHex->mDerLegendreMaxFace1D;
				TensorProduct 	= TensorProductSolAndGradKFace; 
				break;
			}
			default: break;
		}

		// Execute averaging process, if NSCBC iBoundary exists.
		if( BC != NULL && BC->mBoundaryID == iBoundary )
			if( BC->GetTypeBoundaryPrescribed() == BC_OUTFLOW_CHARACTERISTIC
			    || 
					BC->GetTypeBoundaryPrescribed() == BC_INFLOW_CHARACTERISTIC )
		{
	  	// Compute the solution in the integration points of the face. 
			// The solution corresponds to the BC face of this element, 
			// hence the Legendre basis functions and its derivatives must be
	  	// taken at the relevant BC boundary.
			TensorProduct(nInt1D, nVar, nDOFs1D,
	  	              standardHex->mLegendreInt1D,
	  	              standardHex->mDerLegendreInt1D,
	  	              legendreFace,
	  	              derLegendreFace,
	  	              mSol.data(), sol, NULL, NULL, NULL);
	
			// Compute Average data required.
			LocalMachWeighted[iBoundary] += BC->WeightedMachElement(inputParam, nInt,
																		                         sol, dirSign, metric,
																		                         standardHex->mIntWeightsFace);
		}
	}
}

//------------------------------------------------------------------------------

// Function, which sets the averaged Mach number over the NSCBC boundary.
void ElementClass::SetAverageBoundaryMachNSCBC(const su2double *value)
{
	// Class pointer to relevant boundary face.
	SubfaceBaseClass *BC = NULL;

	for(unsigned short iBoundary=0; iBoundary<6; iBoundary++)
	{
		// Process required pointers.
		switch(iBoundary){
			case(0): BC = mBCIMin; break;
			case(1): BC = mBCIMax; break;
			case(2): BC = mBCJMin; break;
			case(3): BC = mBCJMax; break;
			case(4): BC = mBCKMin; break;
			case(5): BC = mBCKMax; break;
			default: break;
		}

		// Assign values, if NSCBC iBoundary exists.
		if( BC != NULL && BC->mBoundaryID == iBoundary )
			if( BC->GetTypeBoundaryPrescribed() == BC_OUTFLOW_CHARACTERISTIC
			    || 
					BC->GetTypeBoundaryPrescribed() == BC_INFLOW_CHARACTERISTIC )
		{
			// Assign data required.
			BC->SetMachAverageBoundary(value[iBoundary]);
		}
	}
}

//------------------------------------------------------------------------------

// Function, which adds fluctuations to the dimensional primitive variables.
void ElementClass::AddFluctuationsPrimVar(const InputParamClass      *inputParam,
                                          const StandardElementClass *standardHex,
                                          const unsigned long        timeStep,
                                          su2double                  **primVar)
{
  // Easier storage of the number of DOFs of the element.
  const int nDOFs = standardHex->mNDOFs;

  // Make a distinction between the cases.
  switch( inputParam->mFluctuationsInitSol )
  {
    case SYNTHETIC:
    case SINE_AND_COSINE:
    {
      TerminateAll("ElementClass::AddFluctuationsPrimVar", __FILE__, __LINE__,
                   "Option not implemented yet");
      break;
    }

    //--------------------------------------------------------------------------

    case RANDOM:
    {
      // Random fluctuations should be added. Determine the inverse of the
      // length in the three directions of the bounding box.
      const su2double lenXInv = one/(inputParam->mUpperCoorBoxFluctuations[0]
                              -      inputParam->mLowerCoorBoxFluctuations[0]);
      const su2double lenYInv = one/(inputParam->mUpperCoorBoxFluctuations[1]
                              -      inputParam->mLowerCoorBoxFluctuations[1]);
      const su2double lenZInv = one/(inputParam->mUpperCoorBoxFluctuations[2]
                              -      inputParam->mLowerCoorBoxFluctuations[2]);
      // Loop over the DOFs.
      for(int l=0; l<nDOFs; ++l)
      {
        // Easier storage of the coordinates.
        const su2double x = mCoorNodalSolDOFs[0][l];
        const su2double y = mCoorNodalSolDOFs[1][l];
        const su2double z = mCoorNodalSolDOFs[2][l];

        // Test if the coordinates are within the bounding box for which
        // the fluctuations must be applied.
        if((x >= inputParam->mLowerCoorBoxFluctuations[0]) &&
           (y >= inputParam->mLowerCoorBoxFluctuations[1]) &&
           (z >= inputParam->mLowerCoorBoxFluctuations[2]) &&
           (x <= inputParam->mUpperCoorBoxFluctuations[0]) &&
           (y <= inputParam->mUpperCoorBoxFluctuations[1]) &&
           (z <= inputParam->mUpperCoorBoxFluctuations[2]))
        {
          // Easier storage of the velocities.
          su2double u = primVar[1][l];
          su2double v = primVar[2][l];
          su2double w = primVar[3][l];

          // Scale the coordinates, such that they are between 0 and 1
          // inside the specified box.
          const su2double xx = lenXInv*(x - inputParam->mLowerCoorBoxFluctuations[0]);
          const su2double yy = lenYInv*(y - inputParam->mLowerCoorBoxFluctuations[1]);
          const su2double zz = lenZInv*(z - inputParam->mLowerCoorBoxFluctuations[2]);

          // Determine the seed of the random generator based on the coordinates.
          // This means that the random generator will always produce the same
          // number for the same coordinates, such that the element boundaries
          // are not visible after adding the perturbations.
          const int seed = (int) (36*xx + 821*yy + 18955*zz);
          std::default_random_engine rand_gen(seed);
          std::uniform_real_distribution<su2double> u01(-0.1,0.1);

          // Compute the perturbations for the three velocities.
          const su2double vMag = SQRT(u*u + v*v + w*w);
          const su2double utrp =  vMag * u01(rand_gen);
          const su2double vtrp =  vMag * u01(rand_gen);
          const su2double wtrp =  vMag * u01(rand_gen);

          // Store the perturbed velocities.
          primVar[1][l] = u + utrp;
          primVar[2][l] = v + vtrp;
          primVar[3][l] = w + wtrp;
        }
      }

      break;
    }

    //--------------------------------------------------------------------------

		case PRESSURE_PULSE:
		{
			// Extract pulse center.
			const su2double x0 = inputParam->mPulseCenter[0];
			const su2double y0 = inputParam->mPulseCenter[1];
			const su2double z0 = inputParam->mPulseCenter[2];

			// Extract pulse strength.
			const su2double d0 = inputParam->mPulseStrength;
			// Extract pulse width.
			const su2double w0 = inputParam->mPulseWidth;
	
			// Characteristic coefficient in Gaussian exponential.
			const su2double cc = half/(w0*w0); 

			// Loop over the DOFs.
			for(int l=0; l<nDOFs; ++l)
			{
				// Extract coordinates, with respect to pulse center.
				const su2double dx = mCoorNodalSolDOFs[0][l] - x0;
				const su2double dy = mCoorNodalSolDOFs[1][l] - y0;
				const su2double dz = mCoorNodalSolDOFs[2][l] - z0;

				// Compute radial distance squared.
				const su2double r2 = dx*dx + dy*dy + dz*dz;

				// Extract the background pressure.
				const su2double pinf = primVar[4][l];

				// Store the new perturbed pressure.
				primVar[4][l] = pinf*( one + d0*EXP(-r2*cc) );
			}

			break;
		}

    //--------------------------------------------------------------------------

    default:  // Just to avoid a compiler warning.
      break;
  }
}

//------------------------------------------------------------------------------

// Function, which allocates the residual and halo solution arrays.
void ElementClass::AllocateResidualAndHaloSolutionArrays(const StandardElementClass *standardHex)
{
  // Abbreviate the number of bytes that must be allocated.
  const size_t nBytes = standardHex->mNDOFsPad * sizeof(su2double);

  // Make a distinction between the element types.
  switch( mElemType )
  {
    case INTERNAL_ELEMENT:
    {
      // Internal element. If the iMin boundary is an internal boundary,
      // allocate the memory for mResIMin.
      if( !mBCIMin )
      {
        mResIMin.resize(nVar, NULL);
        for(int i=0; i<nVar; ++i)
        {
          mResIMin[i] = (su2double *) AllocateMemory(nBytes);
          if( !mResIMin[i] )
            Terminate("ElementClass::AllocateResidualAndHaloSolutionArrays",
                      __FILE__, __LINE__, "Memory allocation failure for mResIMin");
        }
      }

      // If the jMin boundary is an internal boundary, allocate the memory for mResJMin.
      if( !mBCJMin )
      {
        mResJMin.resize(nVar, NULL);
        for(int i=0; i<nVar; ++i)
        {
          mResJMin[i] = (su2double *) AllocateMemory(nBytes);
          if( !mResJMin[i] )
            Terminate("ElementClass::AllocateResidualAndHaloSolutionArrays",
                      __FILE__, __LINE__, "Memory allocation failure for mResJMin");
        }
      }

      // If the kMin boundary is an internal boundary, allocate the memory for mResKMin.
      if( !mBCKMin )
      {
        mResKMin.resize(nVar, NULL);
        for(int i=0; i<nVar; ++i)
        {
          mResKMin[i] = (su2double *) AllocateMemory(nBytes);
          if( !mResKMin[i] )
            Terminate("ElementClass::AllocateResidualAndHaloSolutionArrays",
                      __FILE__, __LINE__, "Memory allocation failure for mResKMin");
        }
      }

      break;
    }

    //--------------------------------------------------------------------------

    case HALO_ELEMENT_IMIN:
    case HALO_ELEMENT_JMIN:
    case HALO_ELEMENT_KMIN:
    {
      // Halo element on a min boundary. The solution for these elements
      // must be communicated. Just allocate the memory for mSol.
      mSol.resize(nVar, NULL);
      for(int i=0; i<nVar; ++i)
      {
        mSol[i] = (su2double *) AllocateMemory(nBytes);
        if( !mSol[i] )
          Terminate("ElementClass::AllocateResidualAndHaloSolutionArrays",
                    __FILE__, __LINE__, "Memory allocation failure for mSol");
      }

      break;
    }

    //--------------------------------------------------------------------------

    case HALO_ELEMENT_IMAX:
    {
      // Halo element on a iMax boundary. This halo is used to store the
      // residual of the iMin boundary of the element. Hence only allocate
      // the memory for mResIMin.
      mResIMin.resize(nVar);
      for(int i=0; i<nVar; ++i)
      {
        mResIMin[i] = (su2double *) AllocateMemory(nBytes);
        if( !mResIMin[i] )
          Terminate("ElementClass::AllocateResidualAndHaloSolutionArrays",
                    __FILE__, __LINE__, "Memory allocation failure for mResIMin");
      }

      break;
    }

    //--------------------------------------------------------------------------

    case HALO_ELEMENT_JMAX:
    {
      // Halo element on a jMax boundary. This halo is used to store the
      // residual of the jMin boundary of the element. Hence only allocate
      // the memory for mResJMin.
      mResJMin.resize(nVar);
      for(int i=0; i<nVar; ++i)
      {
        mResJMin[i] = (su2double *) AllocateMemory(nBytes);
        if( !mResJMin[i] )
          Terminate("ElementClass::AllocateResidualAndHaloSolutionArrays",
                    __FILE__, __LINE__, "Memory allocation failure for mResJMin");
      }

      break;
    }

    //--------------------------------------------------------------------------

    case HALO_ELEMENT_KMAX:
    {
      // Halo element on a kMax boundary. This halo is used to store the
      // residual of the kMin boundary of the element. Hence only allocate
      // the memory for mResKMin.
      mResKMin.resize(nVar);
      for(int i=0; i<nVar; ++i)
      {
        mResKMin[i] = (su2double *) AllocateMemory(nBytes);
        if( !mResKMin[i] )
          Terminate("ElementClass::AllocateResidualAndHaloSolutionArrays",
                    __FILE__, __LINE__, "Memory allocation failure for mResKMin");
      }

      break;
    }
  }
}

//------------------------------------------------------------------------------

// Function, which computes the coordinates of the nodal solution DOFs
// and the metric terms in these DOFs.
void ElementClass::ComputeCoorNodalSolDOFs(const int nDOFs1DGrid,
                                           const int nDOFs1DSol,
                                           su2double *lagrangeDOFs1D,
                                           su2double *derLagrangeDOFs1D)
{
  // Determine the number of solution DOFs in 3D and its padded value.
  const int nDOFs3D    = nDOFs1DSol*nDOFs1DSol*nDOFs1DSol;
  const int nDOFs3DPad = ((nDOFs3D+vecLen3D-1)/vecLen3D)*vecLen3D;

  // Allocate the memory for the coordinates of the nodal solution DOFs.
  mCoorNodalSolDOFs.resize(3, NULL);
  for(unsigned int m=0; m<mCoorNodalSolDOFs.size(); ++m)
  {
    mCoorNodalSolDOFs[m] = (su2double *) AllocateMemory(nDOFs3DPad*sizeof(su2double));
    if( !mCoorNodalSolDOFs[m] )
      Terminate("ElementClass::ComputeCoorNodalSolDOFs", __FILE__, __LINE__, 
                "Memory allocation failure for mCoorNodalSolDOFs");
  }

  // Determine the coordinates of the nodal solution DOFs.
  TensorProductSolAndGradVolume(nDOFs1DSol, 3, nDOFs1DGrid, lagrangeDOFs1D,
                                derLagrangeDOFs1D, mCoorNodalGridDOFs.data(),
                                mCoorNodalSolDOFs.data(), NULL, NULL, NULL);

  // Initialize the padded values to avoid problems.
  for(int m=0; m<3; ++m)
    for(int l=nDOFs3D; l<nDOFs3DPad; ++l)
      mCoorNodalSolDOFs[m][l] = mCoorNodalSolDOFs[m][0];

  // Allocate the memory for the metric terms of the nodal solution DOFs.
  mVolMetricSolDOFs.resize(10, NULL);
  for(unsigned int m=0; m<mVolMetricSolDOFs.size(); ++m)
  {
    mVolMetricSolDOFs[m] = (su2double *) AllocateMemory(nDOFs3DPad*sizeof(su2double));
    if( !mVolMetricSolDOFs[m] )
      Terminate("ElementClass::ComputeCoorNodalSolDOFs", __FILE__, __LINE__, 
                "Memory allocation failure for mVolMetricSolDOFs");
  }

  // Use mVolMetricSolDOFs as a buffer to store the derivatives of the coordinates
  // in the solution DOFs. Set the pointers accordingly.
  su2double **dcoordr = mVolMetricSolDOFs.data();
  su2double **dcoords = dcoordr + 3;
  su2double **dcoordt = dcoords + 3;

  // Determine the derivatives of the coordinates in the volume solution DOFs.
  TensorProductSolAndGradVolume(nDOFs1DSol, 3, nDOFs1DGrid, lagrangeDOFs1D,
                                derLagrangeDOFs1D, mCoorNodalGridDOFs.data(),
                                NULL, dcoordr, dcoords, dcoordt);

  // Initialize the padded values to avoid problems.
  for(int m=0; m<3; ++m)
  {
    for(int l=nDOFs3D; l<nDOFs3DPad; ++l)
    {
      dcoordr[m][l] = dcoordr[m][0];
      dcoords[m][l] = dcoords[m][0];
      dcoordt[m][l] = dcoordt[m][0];
    }
  }

  // Loop over the solution DOFs and compute the metric terms.
#pragma omp simd
  for(int l=0; l<nDOFs3DPad; ++l)
  {
    // Copy the data from dcoordr, dcoords and dcoordt.
    const su2double dxdr = dcoordr[0][l];
    const su2double dydr = dcoordr[1][l];
    const su2double dzdr = dcoordr[2][l];

    const su2double dxds = dcoords[0][l];
    const su2double dyds = dcoords[1][l];
    const su2double dzds = dcoords[2][l];

    const su2double dxdt = dcoordt[0][l];
    const su2double dydt = dcoordt[1][l];
    const su2double dzdt = dcoordt[2][l];

    // Compute the Jacobian and its inverse.
    mVolMetricSolDOFs[0][l] = dxdr*(dyds*dzdt - dzds*dydt)
                            - dxds*(dydr*dzdt - dzdr*dydt)
                            + dxdt*(dydr*dzds - dzdr*dyds);
    const su2double Jinv = one/mVolMetricSolDOFs[0][l];

    // Compute and store the other matrix terms.
    mVolMetricSolDOFs[1][l] = Jinv*(dyds*dzdt - dzds*dydt);  // drdx
    mVolMetricSolDOFs[2][l] = Jinv*(dzds*dxdt - dxds*dzdt);  // drdy
    mVolMetricSolDOFs[3][l] = Jinv*(dxds*dydt - dyds*dxdt);  // drdz

    mVolMetricSolDOFs[4][l] = Jinv*(dzdr*dydt - dydr*dzdt);  // dsdx
    mVolMetricSolDOFs[5][l] = Jinv*(dxdr*dzdt - dzdr*dxdt);  // dsdy
    mVolMetricSolDOFs[6][l] = Jinv*(dydr*dxdt - dxdr*dydt);  // dsdz

    mVolMetricSolDOFs[7][l] = Jinv*(dydr*dzds - dzdr*dyds);  // dtdx
    mVolMetricSolDOFs[8][l] = Jinv*(dzdr*dxds - dxdr*dzds);  // dtdy
    mVolMetricSolDOFs[9][l] = Jinv*(dxdr*dyds - dydr*dxds);  // dtdz
  }

  // Check if negative Jacobians are present.
  for(int l=0; l<nDOFs3D; ++l)
  {
    if(mVolMetricSolDOFs[0][l] <= zero)
      Terminate("ElementClass::ComputeCoorNodalSolDOFs", __FILE__, __LINE__,
                "Negative or zero Jacobians encountered");
  }
}

//------------------------------------------------------------------------------

// Function, which computes the different length scales of the element.
void ElementClass::ComputeLengthScales(const InputParamClass      *inputParam,
                                       const StandardElementClass *standardHex)
{
  // Determine the volume of the element.
  mVolElem = zero;
  for(int l=0; l<standardHex->mNIntegration; ++l)
    mVolElem += standardHex->mIntWeights[l]*mVolMetricInt[0][l];

  // Determine the length scale of the element based on the volume.
  // From this value the LES length scale can be computed, which takes
  // the polynomial degree of the solution into account.
  const su2double third = one/three;
  mLenScaleVolume = POW(mVolElem, third);
  mLenScaleLES    = mLenScaleVolume/inputParam->mNPolySolDOFs;

  // Compute the sum of the areas of the iMin and iMax faces
  // as well as the normals on these faces.
  su2double sumAreasFaces = zero;
  mAveNormIDir[0] = mAveNormIDir[1] = mAveNormIDir[2] = zero;

  for(int l=0; l<standardHex->mNIntegration2D; ++l)
  {
    const su2double w = standardHex->mIntWeightsFace[l];

    mAveNormIDir[0] += w*(mSurfMetricIntIMin[0][l]*mSurfMetricIntIMin[3][l]
                     +    mSurfMetricIntIMax[0][l]*mSurfMetricIntIMax[3][l]);
    mAveNormIDir[1] += w*(mSurfMetricIntIMin[1][l]*mSurfMetricIntIMin[3][l]
                     +    mSurfMetricIntIMax[1][l]*mSurfMetricIntIMax[3][l]);
    mAveNormIDir[2] += w*(mSurfMetricIntIMin[2][l]*mSurfMetricIntIMin[3][l]
                     +    mSurfMetricIntIMax[2][l]*mSurfMetricIntIMax[3][l]);

    sumAreasFaces   += w*(mSurfMetricIntIMin[3][l] + mSurfMetricIntIMax[3][l]);
  }

  // Determine the length scale in i-direction of the element.
  mLenScaleIDir = two*mVolElem/sumAreasFaces;

  // Determine the averaged unit normal in i-direction.
  mAveNormIDir[0] /= sumAreasFaces;
  mAveNormIDir[1] /= sumAreasFaces;
  mAveNormIDir[2] /= sumAreasFaces;

  // Compute the sum of the areas of the jMin and jMax faces
  // as well as the normals on these faces.
  sumAreasFaces = zero;
  mAveNormJDir[0] = mAveNormJDir[1] = mAveNormJDir[2] = zero;

  for(int l=0; l<standardHex->mNIntegration2D; ++l)
  {
    const su2double w = standardHex->mIntWeightsFace[l];

    mAveNormJDir[0] += w*(mSurfMetricIntJMin[0][l]*mSurfMetricIntJMin[3][l]
                     +    mSurfMetricIntJMax[0][l]*mSurfMetricIntJMax[3][l]);
    mAveNormJDir[1] += w*(mSurfMetricIntJMin[1][l]*mSurfMetricIntJMin[3][l]
                     +    mSurfMetricIntJMax[1][l]*mSurfMetricIntJMax[3][l]);
    mAveNormJDir[2] += w*(mSurfMetricIntJMin[2][l]*mSurfMetricIntJMin[3][l]
                     +    mSurfMetricIntJMax[2][l]*mSurfMetricIntJMax[3][l]);

    sumAreasFaces += w*(mSurfMetricIntJMin[3][l] + mSurfMetricIntJMax[3][l]);
  }

  // Determine the length scale in j-direction of the element.
  mLenScaleJDir = two*mVolElem/sumAreasFaces;

  // Determine the averaged unit normal in j-direction.
  mAveNormJDir[0] /= sumAreasFaces;
  mAveNormJDir[1] /= sumAreasFaces;
  mAveNormJDir[2] /= sumAreasFaces;

  // Compute the sum of the areas of the kMin and kMax faces
  // as well as the normals on these faces.
  sumAreasFaces = zero;
  mAveNormKDir[0] = mAveNormKDir[1] = mAveNormKDir[2] = zero;

  for(int l=0; l<standardHex->mNIntegration2D; ++l)
  {
    const su2double w = standardHex->mIntWeightsFace[l];

    mAveNormKDir[0] += w*(mSurfMetricIntKMin[0][l]*mSurfMetricIntKMin[3][l]
                     +    mSurfMetricIntKMax[0][l]*mSurfMetricIntKMax[3][l]);
    mAveNormKDir[1] += w*(mSurfMetricIntKMin[1][l]*mSurfMetricIntKMin[3][l]
                     +    mSurfMetricIntKMax[1][l]*mSurfMetricIntKMax[3][l]);
    mAveNormKDir[2] += w*(mSurfMetricIntKMin[2][l]*mSurfMetricIntKMin[3][l]
                     +    mSurfMetricIntKMax[2][l]*mSurfMetricIntKMax[3][l]);

    sumAreasFaces += w*(mSurfMetricIntKMin[3][l] + mSurfMetricIntKMax[3][l]);
  }

  // Determine the length scale in k-direction of the element.
  mLenScaleKDir = two*mVolElem/sumAreasFaces;

  // Determine the averaged unit normal in k-direction.
  mAveNormKDir[0] /= sumAreasFaces;
  mAveNormKDir[1] /= sumAreasFaces;
  mAveNormKDir[2] /= sumAreasFaces;
}

//------------------------------------------------------------------------------

// Function, which computes the metric terms.
void ElementClass::ComputeMetricTerms(const int nDOFs1DGrid,
                                      const int nInt1D,
                                      su2double *lagrangeInt1D,
                                      su2double *derLagrangeInt1D,
                                      su2double *lagrangeMinFace1D,
                                      su2double *lagrangeMaxFace1D,
                                      su2double *derLagrangeMinFace1D,
                                      su2double *derLagrangeMaxFace1D)
{
  // Determine the number of integration points in 2D and 3D and
  // its padded values.
  const int nInt2D = nInt1D*nInt1D;
  const int nInt3D = nInt2D*nInt1D;

  const int nInt2DPad = ((nInt2D+vecLen2D-1)/vecLen2D)*vecLen2D;
  const int nInt3DPad = ((nInt3D+vecLen3D-1)/vecLen3D)*vecLen3D;

  //------------------------------------------------------
  //  The volume metric terms in the integration points.
  //------------------------------------------------------

  // Allocate the memory for the metric terms in the volume integration points.
  mVolMetricInt.resize(13, NULL);
  for(unsigned int m=0; m<mVolMetricInt.size(); ++m)
  {
    mVolMetricInt[m] = (su2double *) AllocateMemory(nInt3DPad*sizeof(su2double));
    if( !mVolMetricInt[m] )
      Terminate("ElementClass::ComputeMetricTerms", __FILE__, __LINE__, 
                "Memory allocation failure for mVolMetricInt");
  }

  // Use mVolMetricInt as a buffer to store the interpolated coordinates
  // and the derivatives of the coordinates in the integration points.
  // Set the pointers accordingly. Note that the coordinates of the
  // integration points are stored in the last positions of mVolMetricInt.
  su2double **coor    = mVolMetricInt.data() + 10;
  su2double **dcoordr = mVolMetricInt.data();
  su2double **dcoords = dcoordr + 3;
  su2double **dcoordt = dcoords + 3;

  // Determine the coordinates and the derivatives of the coordinates in
  // the volume integration points.
  TensorProductSolAndGradVolume(nInt1D, 3, nDOFs1DGrid, lagrangeInt1D,
                                derLagrangeInt1D, mCoorNodalGridDOFs.data(),
                                coor, dcoordr, dcoords, dcoordt);

  // Initialize the padded values to avoid problems.
  for(int m=0; m<3; ++m)
  {
    for(int l=nInt3D; l<nInt3DPad; ++l)
    {
      coor[m][l]    = coor[m][0];
      dcoordr[m][l] = dcoordr[m][0];
      dcoords[m][l] = dcoords[m][0];
      dcoordt[m][l] = dcoordt[m][0];
    }
  }

  // Loop over the volume integration points to compute the required
  // metric terms.
#pragma omp simd
  for(int l=0; l<nInt3DPad; ++l)
  {
    // Copy the data from dcoordr, dcoords and dcoordt.
    const su2double dxdr = dcoordr[0][l];
    const su2double dydr = dcoordr[1][l];
    const su2double dzdr = dcoordr[2][l];

    const su2double dxds = dcoords[0][l];
    const su2double dyds = dcoords[1][l];
    const su2double dzds = dcoords[2][l];

    const su2double dxdt = dcoordt[0][l];
    const su2double dydt = dcoordt[1][l];
    const su2double dzdt = dcoordt[2][l];

    // Compute and store the Jacobian of the transformation.
    mVolMetricInt[0][l] = dxdr*(dyds*dzdt - dzds*dydt)
                        - dxds*(dydr*dzdt - dzdr*dydt)
                        + dxdt*(dydr*dzds - dzdr*dyds); // J

    // Compute and store the other matrix terms.
    const su2double Jinv = one/mVolMetricInt[0][l];

    mVolMetricInt[1][l] = Jinv*(dyds*dzdt - dzds*dydt);  // drdx
    mVolMetricInt[2][l] = Jinv*(dzds*dxdt - dxds*dzdt);  // drdy
    mVolMetricInt[3][l] = Jinv*(dxds*dydt - dyds*dxdt);  // drdz

    mVolMetricInt[4][l] = Jinv*(dzdr*dydt - dydr*dzdt);  // dsdx
    mVolMetricInt[5][l] = Jinv*(dxdr*dzdt - dzdr*dxdt);  // dsdy
    mVolMetricInt[6][l] = Jinv*(dydr*dxdt - dxdr*dydt);  // dsdz

    mVolMetricInt[7][l] = Jinv*(dydr*dzds - dzdr*dyds);  // dtdx
    mVolMetricInt[8][l] = Jinv*(dzdr*dxds - dxdr*dzds);  // dtdy
    mVolMetricInt[9][l] = Jinv*(dxdr*dyds - dydr*dxds);  // dtdz
  }

  // Check if negative Jacobians are present.
  for(int l=0; l<nInt3D; ++l)
  {
    if(mVolMetricInt[0][l] <= zero)
      Terminate("ElementClass::ComputeMetricTerms", __FILE__, __LINE__,
                "Negative or zero Jacobians encountered");
  }

  //------------------------------------------------------
  //    The surface metric terms of the iMin boundary.
  //------------------------------------------------------

  // Allocate the memory for the metric terms in the surface integration
  // points of the iMin boundary.
  mSurfMetricIntIMin.resize(16, NULL);
  for(unsigned int m=0; m<mSurfMetricIntIMin.size(); ++m)
  {
    mSurfMetricIntIMin[m] = (su2double *) AllocateMemory(nInt2DPad*sizeof(su2double));
    if( !mSurfMetricIntIMin[m] )
      Terminate("ElementClass::ComputeMetricTerms", __FILE__, __LINE__, 
                "Memory allocation failure for mSurfMetricIntIMin");
  }

  // Use mSurfMetricIntIMin as a buffer to store the derivatives of the
  // coordinates in the integration points. Set the pointers accordingly.
  // Store the coordinates of the integration points in the 3 entries.
  coor    = mSurfMetricIntIMin.data() + 13;
  dcoordr = mSurfMetricIntIMin.data();
  dcoords = dcoordr + 3;
  dcoordt = dcoords + 3;

  // Determine the coordinates and its derivatives in the surface
  // integration points of the iMin boundary.
  TensorProductSolAndGradIFace(nInt1D, 3, nDOFs1DGrid, lagrangeInt1D,
                               derLagrangeInt1D, lagrangeMinFace1D,
                               derLagrangeMinFace1D,
                               mCoorNodalGridDOFs.data(), coor,
                               dcoordr, dcoords, dcoordt);

  // Initialize the padded values to avoid problems.
  for(int m=0; m<3; ++m)
  {
    for(int l=nInt2D; l<nInt2DPad; ++l)
    {
      coor[m][l]    = coor[m][0];
      dcoordr[m][l] = dcoordr[m][0];
      dcoords[m][l] = dcoords[m][0];
      dcoordt[m][l] = dcoordt[m][0];
    }
  }

  // Compute the actual metric data on the iMin boundary.
  FinalMetricIFace(nInt2DPad, mSurfMetricIntIMin);

  //------------------------------------------------------
  //    The surface metric terms of the iMax boundary.
  //------------------------------------------------------

  // Allocate the memory for the metric terms in the surface integration
  // points of the iMax boundary.
  mSurfMetricIntIMax.resize(16, NULL);
  for(unsigned int m=0; m<mSurfMetricIntIMax.size(); ++m)
  {
    mSurfMetricIntIMax[m] = (su2double *) AllocateMemory(nInt2DPad*sizeof(su2double));
    if( !mSurfMetricIntIMax[m] )
      Terminate("ElementClass::ComputeMetricTerms", __FILE__, __LINE__, 
                "Memory allocation failure for mSurfMetricIntIMax");
  }

  // Use mSurfMetricIntIMax as a buffer to store the derivatives of the
  // coordinates in the integration points. Set the pointers accordingly.
  // Store the coordinates of the integration points in the 3 entries.
  coor    = mSurfMetricIntIMax.data() + 13;
  dcoordr = mSurfMetricIntIMax.data();
  dcoords = dcoordr + 3;
  dcoordt = dcoords + 3;

  // Determine the derivatives of the coordinates in the surface
  // integration points of the iMax boundary.
  TensorProductSolAndGradIFace(nInt1D, 3, nDOFs1DGrid, lagrangeInt1D,
                               derLagrangeInt1D, lagrangeMaxFace1D,
                               derLagrangeMaxFace1D,
                               mCoorNodalGridDOFs.data(), coor,
                               dcoordr, dcoords, dcoordt);

  // Initialize the padded values to avoid problems.
  for(int m=0; m<3; ++m)
  {
    for(int l=nInt2D; l<nInt2DPad; ++l)
    {
      coor[m][l]    = coor[m][0];
      dcoordr[m][l] = dcoordr[m][0];
      dcoords[m][l] = dcoords[m][0];
      dcoordt[m][l] = dcoordt[m][0];
    }
  }

  // Compute the actual metric data on the iMax boundary.
  FinalMetricIFace(nInt2DPad, mSurfMetricIntIMax);

  //------------------------------------------------------
  //    The surface metric terms of the jMin boundary.
  //------------------------------------------------------

  // Allocate the memory for the metric terms in the surface integration
  // points of the jMin boundary.
  mSurfMetricIntJMin.resize(16, NULL);
  for(unsigned int m=0; m<mSurfMetricIntJMin.size(); ++m)
  {
    mSurfMetricIntJMin[m] = (su2double *) AllocateMemory(nInt2DPad*sizeof(su2double));
    if( !mSurfMetricIntJMin[m] )
      Terminate("ElementClass::ComputeMetricTerms", __FILE__, __LINE__, 
                "Memory allocation failure for mSurfMetricIntJMin");
  }

  // Use mSurfMetricIntJMin as a buffer to store the derivatives of the
  // coordinates in the integration points. Set the pointers accordingly.
  // Store the coordinates of the integration points in the 3 entries.
  coor    = mSurfMetricIntJMin.data() + 13;
  dcoordr = mSurfMetricIntJMin.data();
  dcoords = dcoordr + 3;
  dcoordt = dcoords + 3;

  // Determine the derivatives of the coordinates in the surface
  // integration points of the jMin boundary.
  TensorProductSolAndGradJFace(nInt1D, 3, nDOFs1DGrid, lagrangeInt1D,
                               derLagrangeInt1D, lagrangeMinFace1D,
                               derLagrangeMinFace1D,
                               mCoorNodalGridDOFs.data(), coor,
                               dcoordr, dcoords, dcoordt);

  // Initialize the padded values to avoid problems.
  for(int m=0; m<3; ++m)
  {
    for(int l=nInt2D; l<nInt2DPad; ++l)
    {
      coor[m][l]    = coor[m][0];
      dcoordr[m][l] = dcoordr[m][0];
      dcoords[m][l] = dcoords[m][0];
      dcoordt[m][l] = dcoordt[m][0];
    }
  }

  // Compute the actual metric data on the jMin boundary.
  FinalMetricJFace(nInt2DPad, mSurfMetricIntJMin);

  //------------------------------------------------------
  //    The surface metric terms of the jMax boundary.
  //------------------------------------------------------

  // Allocate the memory for the metric terms in the surface integration
  // points of the jMax boundary.
  mSurfMetricIntJMax.resize(16, NULL);
  for(unsigned int m=0; m<mSurfMetricIntJMax.size(); ++m)
  {
    mSurfMetricIntJMax[m] = (su2double *) AllocateMemory(nInt2DPad*sizeof(su2double));
    if( !mSurfMetricIntJMax[m] )
      Terminate("ElementClass::ComputeMetricTerms", __FILE__, __LINE__, 
                "Memory allocation failure for mSurfMetricIntJMax");
  }

  // Use mSurfMetricIntJMax as a buffer to store the derivatives of the
  // coordinates in the integration points. Set the pointers accordingly.
  // Store the coordinates of the integration points in the 3 entries.
  coor    = mSurfMetricIntJMax.data() + 13;
  dcoordr = mSurfMetricIntJMax.data();
  dcoords = dcoordr + 3;
  dcoordt = dcoords + 3;

  // Determine the derivatives of the coordinates in the surface
  // integration points of the jMax boundary.
  TensorProductSolAndGradJFace(nInt1D, 3, nDOFs1DGrid, lagrangeInt1D,
                               derLagrangeInt1D, lagrangeMaxFace1D,
                               derLagrangeMaxFace1D,
                               mCoorNodalGridDOFs.data(), coor,
                               dcoordr, dcoords, dcoordt);

  // Initialize the padded values to avoid problems.
  for(int m=0; m<3; ++m)
  {
    for(int l=nInt2D; l<nInt2DPad; ++l)
    {
      coor[m][l]    = coor[m][0];
      dcoordr[m][l] = dcoordr[m][0];
      dcoords[m][l] = dcoords[m][0];
      dcoordt[m][l] = dcoordt[m][0];
    }
  }

  // Compute the actual metric data on the jMax boundary.
  FinalMetricJFace(nInt2DPad, mSurfMetricIntJMax);

  //------------------------------------------------------
  //    The surface metric terms of the kMin boundary.
  //------------------------------------------------------

  // Allocate the memory for the metric terms in the surface integration
  // points of the kMin boundary.
  mSurfMetricIntKMin.resize(16, NULL);
  for(unsigned int m=0; m<mSurfMetricIntKMin.size(); ++m)
  {
    mSurfMetricIntKMin[m] = (su2double *) AllocateMemory(nInt2DPad*sizeof(su2double));
    if( !mSurfMetricIntKMin[m] )
      Terminate("ElementClass::ComputeMetricTerms", __FILE__, __LINE__, 
                "Memory allocation failure for mSurfMetricIntKMin");
  }

  // Use mSurfMetricIntKMin as a buffer to store the derivatives of the
  // coordinates in the integration points. Set the pointers accordingly.
  // Store the coordinates of the integration points in the 3 entries.
  coor    = mSurfMetricIntKMin.data() + 13;
  dcoordr = mSurfMetricIntKMin.data();
  dcoords = dcoordr + 3;
  dcoordt = dcoords + 3;

  // Determine the derivatives of the coordinates in the surface
  // integration points of the kMin boundary.
  TensorProductSolAndGradKFace(nInt1D, 3, nDOFs1DGrid, lagrangeInt1D,
                               derLagrangeInt1D, lagrangeMinFace1D,
                               derLagrangeMinFace1D,
                               mCoorNodalGridDOFs.data(), coor,
                               dcoordr, dcoords, dcoordt);

  // Initialize the padded values to avoid problems.
  for(int m=0; m<3; ++m)
  {
    for(int l=nInt2D; l<nInt2DPad; ++l)
    {
      coor[m][l]    = coor[m][0];
      dcoordr[m][l] = dcoordr[m][0];
      dcoords[m][l] = dcoords[m][0];
      dcoordt[m][l] = dcoordt[m][0];
    }
  }

  // Compute the actual metric data on the kMin boundary.
  FinalMetricKFace(nInt2DPad, mSurfMetricIntKMin);

  //------------------------------------------------------
  //    The surface metric terms of the kMax boundary.
  //------------------------------------------------------

  // Allocate the memory for the metric terms in the surface integration
  // points of the kMax boundary.
  mSurfMetricIntKMax.resize(16, NULL);
  for(unsigned int m=0; m<mSurfMetricIntKMax.size(); ++m)
  {
    mSurfMetricIntKMax[m] = (su2double *) AllocateMemory(nInt2DPad*sizeof(su2double));
    if( !mSurfMetricIntKMax[m] )
      Terminate("ElementClass::ComputeMetricTerms", __FILE__, __LINE__, 
                "Memory allocation failure for mSurfMetricIntKMax");
  }

  // Use mSurfMetricIntKMax as a buffer to store the derivatives of the
  // coordinates in the integration points. Set the pointers accordingly.
  // Store the coordinates of the integration points in the 3 entries.
  coor    = mSurfMetricIntKMax.data() + 13;
  dcoordr = mSurfMetricIntKMax.data();
  dcoords = dcoordr + 3;
  dcoordt = dcoords + 3;

  // Determine the derivatives of the coordinates in the surface
  // integration points of the kMax boundary.
  TensorProductSolAndGradKFace(nInt1D, 3, nDOFs1DGrid, lagrangeInt1D,
                               derLagrangeInt1D, lagrangeMaxFace1D,
                               derLagrangeMaxFace1D,
                               mCoorNodalGridDOFs.data(), coor,
                               dcoordr, dcoords, dcoordt);

  // Initialize the padded values to avoid problems.
  for(int m=0; m<3; ++m)
  {
    for(int l=nInt2D; l<nInt2DPad; ++l)
    {
      coor[m][l]    = coor[m][0];
      dcoordr[m][l] = dcoordr[m][0];
      dcoords[m][l] = dcoords[m][0];
      dcoordt[m][l] = dcoordt[m][0];
    }
  }

  // Compute the actual metric data on the kMax boundary.
  FinalMetricKFace(nInt2DPad, mSurfMetricIntKMax);
}

//------------------------------------------------------------------------------

// Function, which computes the dimensional primitive variables
// for the given number of items.
void ElementClass::ComputePrimitiveVariables(const InputParamClass *inputParam,
                                             const int             nItems,
                                             su2double             **primVar)
{
  // Make a distinction between the working variables.
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      // The conservative variables are used as working variables.
      // Easier storage of gamma-1.
      const su2double gm1 = GamConstant - one;

      // Loop over the number of items.
#pragma omp simd
      for(int l=0; l<nItems; ++l)
      {
        // Compute the primitive variables from the conservative ones,
        // which are currently stored in primVar.
        const su2double rho    = primVar[0][l];
        const su2double rhoInv = one/rho;
        const su2double u      = rhoInv*primVar[1][l];
        const su2double v      = rhoInv*primVar[2][l];
        const su2double w      = rhoInv*primVar[3][l];
        const su2double p      = gm1*(primVar[4][l] - half*rho*(u*u + v*v + w*w));

        // Store the dimensional primitive variables.
        primVar[0][l] = rhoRef*rho;
        primVar[1][l] = uRef  *u;
        primVar[2][l] = uRef  *v;
        primVar[3][l] = uRef  *w;
        primVar[4][l] = pRef  *p;
      }

      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      // The entropy variables are used as working variables.
      // Easier storage of some expressions involving gamma.
      const su2double gm1   =  GamConstant - one;
      const su2double ov1mg = -one/gm1;

      // Loop over the number of items.
#pragma omp simd
      for(int l=0; l<nItems; ++l)
      {
        // Compute the primitive variables from the entropy ones,
        // which are currently stored in primVar.
        const su2double V4Inv =  one/primVar[4][l];
        const su2double u     = -V4Inv*primVar[1][l];
        const su2double v     = -V4Inv*primVar[2][l];
        const su2double w     = -V4Inv*primVar[3][l];
        const su2double eKin  =  half*(u*u + v*v + w*w);
        const su2double s     =  GamConstant - gm1*(primVar[0][l] - primVar[4][l]*eKin);
        const su2double tmp   = -primVar[4][l]*EXP(s);
        const su2double rho   =  POW(tmp, ov1mg);
        const su2double p     = -rho*V4Inv;

        // Store the dimensional primitive variables.
        primVar[0][l] = rhoRef*rho;
        primVar[1][l] = uRef  *u;
        primVar[2][l] = uRef  *v;
        primVar[3][l] = uRef  *w;
        primVar[4][l] = pRef  *p;
      }

      break;
    }

    //--------------------------------------------------------------------------

    default:
    {
       // This is just to avoid a compiler warning.
       break;
    }
  }
}

//------------------------------------------------------------------------------

// Function, which computes the dimensional primitive variables in
// the nodal DOFs.
void ElementClass::ComputePrimitiveVariablesNodalDOFs(const InputParamClass      *inputParam,
                                                      const StandardElementClass *standardHex,
                                                      su2double                  **primVar)
{
  // Interpolate the working variables to the nodal DOFs.
  const int nDOFs1D = standardHex->mNDOFs1D;
  TensorProductSolAndGradVolume(nDOFs1D, nVar, nDOFs1D, standardHex->mLegendreDOFs1D,
                                standardHex->mDerLegendreDOFs1D, mSol.data(), 
                                primVar, NULL, NULL, NULL);

  // Call ComputePrimitiveVariables to carry out the conversion to primitive variables.
  ComputePrimitiveVariables(inputParam, standardHex->mNDOFs, primVar);

  // Initialize the padded values to the first value.
  for(int m=0; m<nVar; ++m)
    for(int l=standardHex->mNDOFs; l<standardHex->mNDOFsPad; ++l)
      primVar[m][l] = primVar[m][0];
}

//------------------------------------------------------------------------------

// Function, which computes the allowable time step for the element.
su2double ElementClass::ComputeTimeStep(const InputParamClass      *inputParam,
                                        const StandardElementClass *standardHex,
                                        const su2double            *f1,
                                        const su2double            *f2,
                                        su2double                  **solNodal,
                                        su2double                  **gradSolNodal)
{
  // Compute the dimensional viscosity. Note that the reference length is taken as one.
  const su2double muDim = mu*pRef/uRef;

  // Easier storage of the inverse of the length scale in the three directions
  // as well as the square of these values. Sum these values.
  const su2double dliInv = one/mLenScaleIDir;
  const su2double dljInv = one/mLenScaleJDir;
  const su2double dlkInv = one/mLenScaleKDir;

  const su2double dli2Inv = dliInv*dliInv, dlj2Inv = dljInv*dljInv, dlk2Inv = dlkInv*dlkInv;
  const su2double dLen2Inv = dli2Inv + dlj2Inv + dlk2Inv;

  // Determine the viscous spectral radius over nu for both the laminar and the
  // turbulent part.
  const su2double radVisOverNuLam  = std::max(std::max(one, two+lambdaOverMu), factHeatFlux_Lam)*dLen2Inv;
  const su2double radVisOverNuTurb = std::max(std::max(one, two+lambdaOverMu), factHeatFlux_Turb)*dLen2Inv;

  // Easier storage of the total number of DOFs and its padded value.
  const int nDOFs    = standardHex->mNDOFs;
  const int nDOFsPad = standardHex->mNDOFsPad;

  // Compute the primitive variables in the nodal DOFs.
  ComputePrimitiveVariablesNodalDOFs(inputParam, standardHex, solNodal);

  // Initialize the padded values to avoid problems.
  for(int m=0; m<nVar; ++m)
    for(int l=nDOFs; l<nDOFsPad; ++l)
      solNodal[m][l] = solNodal[m][0];

  // Set the pointer for the eddy viscosity. Note that it is assumed that the first
  // index of solNodal is allocated with dimension nVar+1.
  su2double *eddyVis = solNodal[nVar];

  // Determine whether or not a subgrid scale model is used.
  if(inputParam->mSGSModelType != NO_SGS_MODEL)
  {
    // Subgrid scale model is used. First compute the
    // gradients of the velocities in the nodal DOFs.
    ComputeVelocityGradientsNodalDOFs(inputParam, standardHex, solNodal,
                                      gradSolNodal);

    // Compute the eddy viscosities in the nodal DOFs.
    inputParam->mSGSModel->EddyViscosity(nDOFsPad, mLenScaleLES, solNodal,
                                         gradSolNodal, eddyVis);
  }
  else
  {
    // No subgrid scale model used. Set the eddy viscosity to zero.
#pragma omp simd
    for(int l=0; l<nDOFsPad; ++l)
      eddyVis[l] = zero;
  }

  // Loop over the padded number of nodal DOFs and compute the inviscid
  // and viscous spectral radius and the speed of sound squared.
  // These values are stored in the first three indices of solNodal.
#pragma omp simd
  for(int l=0; l<nDOFsPad; ++l)
  {
    // Easier storage of the primitive variables.
    const su2double rho = solNodal[0][l];
    const su2double u   = solNodal[1][l];
    const su2double v   = solNodal[2][l];
    const su2double w   = solNodal[3][l];
    const su2double p   = solNodal[4][l];

    // Compute the speed of sound squared.
    const su2double rhoInv = one/rho;
    const su2double a2 = GamConstant*p*rhoInv;

    // Compute the velocities in i-, j- and k-direction.
    const su2double velI = u*mAveNormIDir[0] + v*mAveNormIDir[1] + w*mAveNormIDir[2];
    const su2double velJ = u*mAveNormJDir[0] + v*mAveNormJDir[1] + w*mAveNormJDir[2];
    const su2double velK = u*mAveNormKDir[0] + v*mAveNormKDir[1] + w*mAveNormKDir[2];

    // Compute the inviscid and viscous spectral radius for this DOF.
    const su2double a = SQRT(FABS(a2));
    const su2double radInv = (FABS(velI)+a)*dliInv + (FABS(velJ)+a)*dljInv
                           + (FABS(velK)+a)*dlkInv;
    const su2double radVis = rhoInv*(muDim*radVisOverNuLam
                           +         eddyVis[l]*radVisOverNuTurb);

    // Store the the spectral radii and a2 in solNodal.
    solNodal[0][l] = radInv;
    solNodal[1][l] = radVis;
    solNodal[2][l] = a2;
  }

  // Initialization of the minimum value of the speed of sound squared
  // and the maximim values for the spectral radii.
  su2double a2Min = (su2double) 1.e+10;
  su2double radInvMax = zero, radVisMax = zero;

  // Loop over the DOFs and compute the data for a2Min, radInvMax and radVisMax.
  for(int l=0; l<nDOFs; ++l)
  {
    radInvMax = std::max(radInvMax, solNodal[0][l]);
    radVisMax = std::max(radVisMax, solNodal[1][l]);
    a2Min     = std::min(a2Min,     solNodal[2][l]);
  }

  // If the minimum of the speed of sound squared is negative, terminate.
  if(a2Min < zero)
    Terminate("ElementClass::ComputeTimeStep",
              __FILE__, __LINE__, "Negative square speed of sound");

  // Compute the inverse of the local time step per unit CFL number.
  const int nPoly = standardHex->mNPoly;
  const su2double dtInv = f1[nPoly]*radInvMax + f2[nPoly]*radVisMax;

  // Compute the time step for this element, which is returned.
  // Note that it must be non-dimensionalized, where the reference
  // length is one.
  const su2double dtElem = inputParam->mCFL*uRef/dtInv;
  return dtElem;
}

//------------------------------------------------------------------------------

// Function, which computes the velocity gradients in the nodal DOFs.
void ElementClass::ComputeVelocityGradientsNodalDOFs(const InputParamClass      *inputParam,
                                                     const StandardElementClass *standardHex,
                                                     su2double                  **primVarNodal,
                                                     su2double                  **gradSolNodal)
{
  // Store the inverse of the reference value of the density and velocity.
  const su2double rhoRefInv = one/rhoRef;
  const su2double uRefInv   = one/uRef;

  // Easier storage of the number of DOFs and its padded values.
  const int nDOFs    = standardHex->mNDOFs;
  const int nDOFsPad = standardHex->mNDOFsPad;

  // Set the pointers for the gradients in xi-, eta- and zeta-direction.
  su2double **dSolDr = gradSolNodal;
  su2double **dSolDs = dSolDr + nVar;
  su2double **dSolDt = dSolDs + nVar;

  // Compute the gradients of the working variables in the nodal DOFs.
  const int nDOFs1D = standardHex->mNDOFs1D;
  TensorProductSolAndGradVolume(nDOFs1D, nVar, nDOFs1D, standardHex->mLegendreDOFs1D,
                                standardHex->mDerLegendreDOFs1D, mSol.data(),
                                NULL, dSolDr, dSolDs, dSolDt);

  // Initialize the padded values to zero.
  for(int m=0; m<nVar; ++m)
    for(int l=nDOFs; l<nDOFsPad; ++l)
      dSolDr[m][l] = dSolDs[m][l] = dSolDt[m][l] = zero;

  // Make a distinction between the working variables.
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      // Conservative variables are used. Loop over the padded number of DOFs.
#pragma omp simd
      for(int l=0; l<nDOFsPad; ++l)
      {
        // Compute the non-dimensional density and velocities.
        const su2double rho = primVarNodal[0][l]*rhoRefInv;
        const su2double u   = primVarNodal[1][l]*uRefInv;
        const su2double v   = primVarNodal[2][l]*uRefInv;
        const su2double w   = primVarNodal[3][l]*uRefInv;

        // Compute the non-dimensional velocity gradients in computational space.
        const su2double rhoInv = one/rho;

        const su2double dudr = rhoInv*(dSolDr[1][l] - u*dSolDr[0][l]);
        const su2double duds = rhoInv*(dSolDs[1][l] - u*dSolDs[0][l]);
        const su2double dudt = rhoInv*(dSolDt[1][l] - u*dSolDt[0][l]);

        const su2double dvdr = rhoInv*(dSolDr[2][l] - v*dSolDr[0][l]);
        const su2double dvds = rhoInv*(dSolDs[2][l] - v*dSolDs[0][l]);
        const su2double dvdt = rhoInv*(dSolDt[2][l] - v*dSolDt[0][l]);

        const su2double dwdr = rhoInv*(dSolDr[3][l] - w*dSolDr[0][l]);
        const su2double dwds = rhoInv*(dSolDs[3][l] - w*dSolDs[0][l]);
        const su2double dwdt = rhoInv*(dSolDt[3][l] - w*dSolDt[0][l]);

        // Easier storage of the metric terms.
        const su2double drdx = mVolMetricSolDOFs[1][l];
        const su2double drdy = mVolMetricSolDOFs[2][l];
        const su2double drdz = mVolMetricSolDOFs[3][l];

        const su2double dsdx = mVolMetricSolDOFs[4][l];
        const su2double dsdy = mVolMetricSolDOFs[5][l];
        const su2double dsdz = mVolMetricSolDOFs[6][l];

        const su2double dtdx = mVolMetricSolDOFs[7][l];
        const su2double dtdy = mVolMetricSolDOFs[8][l];
        const su2double dtdz = mVolMetricSolDOFs[9][l];

        // Store the dimensional Cartesian velocity gradients in gradSolNodal.
        gradSolNodal[0][l] = uRef*(dudr*drdx + duds*dsdx + dudt*dtdx);
        gradSolNodal[1][l] = uRef*(dvdr*drdx + dvds*dsdx + dvdt*dtdx);
        gradSolNodal[2][l] = uRef*(dwdr*drdx + dwds*dsdx + dwdt*dtdx);

        gradSolNodal[3][l] = uRef*(dudr*drdy + duds*dsdy + dudt*dtdy);
        gradSolNodal[4][l] = uRef*(dvdr*drdy + dvds*dsdy + dvdt*dtdy);
        gradSolNodal[5][l] = uRef*(dwdr*drdy + dwds*dsdy + dwdt*dtdy);

        gradSolNodal[6][l] = uRef*(dudr*drdz + duds*dsdz + dudt*dtdz);
        gradSolNodal[7][l] = uRef*(dvdr*drdz + dvds*dsdz + dvdt*dtdz);
        gradSolNodal[8][l] = uRef*(dwdr*drdz + dwds*dsdz + dwdt*dtdz);
      }

      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      // Entropy variables are used. Loop over the padded number of DOFs.
#pragma omp simd
      for(int l=0; l<nDOFsPad; ++l)
      {
        // Compute the non-dimensional velocities and V4Inv.
        const su2double u     =  primVarNodal[1][l]*uRefInv;
        const su2double v     =  primVarNodal[2][l]*uRefInv;
        const su2double w     =  primVarNodal[3][l]*uRefInv;
        const su2double V4Inv = -primVarNodal[4][l]*uRefInv*uRefInv/primVarNodal[0][l];

        // Compute the non-dimensional velocity gradients in computational space.
        const su2double dudr = -V4Inv*(dSolDr[1][l] + u*dSolDr[4][l]);
        const su2double dvdr = -V4Inv*(dSolDr[2][l] + v*dSolDr[4][l]);
        const su2double dwdr = -V4Inv*(dSolDr[3][l] + w*dSolDr[4][l]);

        const su2double duds = -V4Inv*(dSolDs[1][l] + u*dSolDs[4][l]);
        const su2double dvds = -V4Inv*(dSolDs[2][l] + v*dSolDs[4][l]);
        const su2double dwds = -V4Inv*(dSolDs[3][l] + w*dSolDs[4][l]);

        const su2double dudt = -V4Inv*(dSolDt[1][l] + u*dSolDt[4][l]);
        const su2double dvdt = -V4Inv*(dSolDt[2][l] + v*dSolDt[4][l]);
        const su2double dwdt = -V4Inv*(dSolDt[3][l] + w*dSolDt[4][l]);

        // Easier storage of the metric terms.
        const su2double drdx = mVolMetricSolDOFs[1][l];
        const su2double drdy = mVolMetricSolDOFs[2][l];
        const su2double drdz = mVolMetricSolDOFs[3][l];

        const su2double dsdx = mVolMetricSolDOFs[4][l];
        const su2double dsdy = mVolMetricSolDOFs[5][l];
        const su2double dsdz = mVolMetricSolDOFs[6][l];

        const su2double dtdx = mVolMetricSolDOFs[7][l];
        const su2double dtdy = mVolMetricSolDOFs[8][l];
        const su2double dtdz = mVolMetricSolDOFs[9][l];

        // Store the dimensional Cartesian velocity gradients in gradSolNodal.
        gradSolNodal[0][l] = uRef*(dudr*drdx + duds*dsdx + dudt*dtdx);
        gradSolNodal[1][l] = uRef*(dvdr*drdx + dvds*dsdx + dvdt*dtdx);
        gradSolNodal[2][l] = uRef*(dwdr*drdx + dwds*dsdx + dwdt*dtdx);

        gradSolNodal[3][l] = uRef*(dudr*drdy + duds*dsdy + dudt*dtdy);
        gradSolNodal[4][l] = uRef*(dvdr*drdy + dvds*dsdy + dvdt*dtdy);
        gradSolNodal[5][l] = uRef*(dwdr*drdy + dwds*dsdy + dwdt*dtdy);

        gradSolNodal[6][l] = uRef*(dudr*drdz + duds*dsdz + dudt*dtdz);
        gradSolNodal[7][l] = uRef*(dvdr*drdz + dvds*dsdz + dvdt*dtdz);
        gradSolNodal[8][l] = uRef*(dwdr*drdz + dwds*dsdz + dwdt*dtdz);
      }

      break;
    }

    //--------------------------------------------------------------------------

    default:
    {
      // This is just to avoid a compiler warning.
      break;
    }
  }
}
    
//------------------------------------------------------------------------------

// Function, which computes the working variables from the primitive variables.
void ElementClass::ComputeWorkingVariables(const InputParamClass      *inputParam,
                                           const StandardElementClass *standardHex,
                                           su2double                  **primVar)
{
  // Easier storage of the number of DOFs of the element and compute
  // the inverse of the reference values.
  const int nDOFs = standardHex->mNDOFs;

  const su2double rhoRefInv = one/rhoRef;
  const su2double uRefInv   = one/uRef;
  const su2double pRefInv   = one/pRef;

  // Convert the dimensional primitive variables in the nodal DOFs to the working
  // variables in the nodal DOFs. This depends on the working variables used.
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      // The conservative variables are used as working variables.
      // Compute 1/(gam-1).
      const su2double ovgm1 = GamConstant - one;

      // Loop over the DOFs.
#pragma omp simd
      for(int l=0; l<nDOFs; ++l)
      {
        // Compute the non-dimensional primitive variables.
        const su2double rho = primVar[0][l]*rhoRefInv;
        const su2double u   = primVar[1][l]*uRefInv;
        const su2double v   = primVar[2][l]*uRefInv;
        const su2double w   = primVar[3][l]*uRefInv;
        const su2double p   = primVar[4][l]*pRefInv;

        // Compute the conservative variables, which are
        // stored again primVar.
        primVar[0][l] = rho;
        primVar[1][l] = rho*u;
        primVar[2][l] = rho*v;
        primVar[3][l] = rho*w;
        primVar[4][l] = p*ovgm1 + half*rho*(u*u + v*v + w*w);
      }

      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      // The entropy variables are used as working variables.
      // Easier storage of some expressions involving gamma.
      const su2double gm1   = GamConstant - one;
      const su2double ovgm1 = one/gm1;

      // Loop over the DOFs.
#pragma omp simd
      for(int l=0; l<nDOFs; ++l)
      {
        // Compute the non-dimensional primitive variables.
        const su2double rho = primVar[0][l]*rhoRefInv;
        const su2double u   = primVar[1][l]*uRefInv;
        const su2double v   = primVar[2][l]*uRefInv;
        const su2double w   = primVar[3][l]*uRefInv;
        const su2double p   = primVar[4][l]*pRefInv;

        // Compute the entropy variables, which are
        // stored again in primVar.
        const su2double s    = LOG(p/POW(rho,GamConstant));
        const su2double pInv = one/p;

        primVar[0][l] =  (GamConstant-s)*ovgm1 - half*rho*(u*u + v*v + w*w)*pInv;
        primVar[1][l] =  rho*u*pInv;
        primVar[2][l] =  rho*v*pInv;
        primVar[3][l] =  rho*w*pInv;
        primVar[4][l] = -rho*pInv;
      }

      break;
    }

    //--------------------------------------------------------------------------

    default:
    {
       // This is just to avoid a compiler warning.
       break;
    }
  }

  // Convert the values in the nodal DOFs to the modal form.
  const int nDOFs1D = standardHex->mNDOFs1D;
  TensorProductSolAndGradVolume(nDOFs1D, nVar, nDOFs1D, standardHex->mVandermonde1DInverse,
                                standardHex->mVandermonde1DInverse, primVar, mSol.data(), 
                                NULL, NULL, NULL);
}

//------------------------------------------------------------------------------

// Function, which creates the final metric terms for a face in i-direction.
void ElementClass::FinalMetricIFace(const int                nItems,
                                    std::vector<su2double *> &surfMetricInt)
{
  // Loop over the number of items for which the final metric terms
  // must be computed.
#pragma omp simd
  for(int l=0; l<nItems; ++l)
  {
    // On entry the derivatives of the physical coordinates w.r.t. the
    // parametric coordinates are stored in surfMetricInt. Copy these values.
    const su2double dxdr = surfMetricInt[0][l];
    const su2double dydr = surfMetricInt[1][l];
    const su2double dzdr = surfMetricInt[2][l];

    const su2double dxds = surfMetricInt[3][l];
    const su2double dyds = surfMetricInt[4][l];
    const su2double dzds = surfMetricInt[5][l];

    const su2double dxdt = surfMetricInt[6][l];
    const su2double dydt = surfMetricInt[7][l];
    const su2double dzdt = surfMetricInt[8][l];

    // Carry out the cross product dxds X dxdt, where x stands for the
    // vector (x,y,z). This gives the normal of the face pointing in
    // positive r-direction.
    const su2double nx = dyds*dzdt - dydt*dzds;
    const su2double ny = dxdt*dzds - dxds*dzdt;
    const su2double nz = dxds*dydt - dxdt*dyds;

    // Determine the length and the inverse of the length of this normal.
    const su2double lenNorm    = SQRT(nx*nx + ny*ny + nz*nz);
    const su2double invLenNorm = one/lenNorm;

    // Store the unit normal and the length in the first 4 entries of surfMetricInt.
    surfMetricInt[0][l] = nx*invLenNorm;
    surfMetricInt[1][l] = ny*invLenNorm;
    surfMetricInt[2][l] = nz*invLenNorm;
    surfMetricInt[3][l] = lenNorm;

    // Compute the derivatives of the parametric coordinates w.r.t. the
    // Cartesian coordinates, i.e. drdx, drdy, etc.
    const su2double Jinv = one/(dxdr*(dyds*dzdt - dzds*dydt)
                         -      dxds*(dydr*dzdt - dzdr*dydt)
                         +      dxdt*(dydr*dzds - dzdr*dyds));

    surfMetricInt[4][l] = (dyds*dzdt - dzds*dydt)*Jinv;  // drdx
    surfMetricInt[5][l] = (dzds*dxdt - dxds*dzdt)*Jinv;  // drdy
    surfMetricInt[6][l] = (dxds*dydt - dyds*dxdt)*Jinv;  // drdz

    surfMetricInt[7][l] = (dzdr*dydt - dydr*dzdt)*Jinv;  // dsdx
    surfMetricInt[8][l] = (dxdr*dzdt - dzdr*dxdt)*Jinv;  // dsdy
    surfMetricInt[9][l] = (dydr*dxdt - dxdr*dydt)*Jinv;  // dsdz

    surfMetricInt[10][l] = (dydr*dzds - dzdr*dyds)*Jinv;  // dtdx
    surfMetricInt[11][l] = (dzdr*dxds - dxdr*dzds)*Jinv;  // dtdy
    surfMetricInt[12][l] = (dxdr*dyds - dydr*dxds)*Jinv;  // dtdz
  }
}

//------------------------------------------------------------------------------

// Function, which creates the final metric terms for a face in j-direction.
void ElementClass::FinalMetricJFace(const int                nItems,
                                    std::vector<su2double *> &surfMetricInt)
{
  // Loop over the number of items for which the final metric terms
  // must be computed.
#pragma omp simd
  for(int l=0; l<nItems; ++l)
  {
    // On entry the derivatives of the physical coordinates w.r.t. the
    // parametric coordinates are stored in surfMetricInt. Copy these values.
    const su2double dxdr = surfMetricInt[0][l];
    const su2double dydr = surfMetricInt[1][l];
    const su2double dzdr = surfMetricInt[2][l];

    const su2double dxds = surfMetricInt[3][l];
    const su2double dyds = surfMetricInt[4][l];
    const su2double dzds = surfMetricInt[5][l];

    const su2double dxdt = surfMetricInt[6][l];
    const su2double dydt = surfMetricInt[7][l];
    const su2double dzdt = surfMetricInt[8][l];

    // Carry out the cross product dxdt X dxdr, where x stands for the
    // vector (x,y,z). This gives the normal of the face pointing in
    // positive s-direction.
    const su2double nx = dydt*dzdr - dydr*dzdt;
    const su2double ny = dxdr*dzdt - dxdt*dzdr;
    const su2double nz = dxdt*dydr - dxdr*dydt;

    // Determine the length and the inverse of the length of this normal.
    const su2double lenNorm    = SQRT(nx*nx + ny*ny + nz*nz);
    const su2double invLenNorm = one/lenNorm;

    // Store the unit normal and the length in the first 4 entries of surfMetricInt.
    surfMetricInt[0][l] = nx*invLenNorm;
    surfMetricInt[1][l] = ny*invLenNorm;
    surfMetricInt[2][l] = nz*invLenNorm;
    surfMetricInt[3][l] = lenNorm;

    // Compute the derivatives of the parametric coordinates w.r.t. the
    // Cartesian coordinates, i.e. drdx, drdy, etc.
    const su2double Jinv = one/(dxdr*(dyds*dzdt - dzds*dydt)
                         -      dxds*(dydr*dzdt - dzdr*dydt)
                         +      dxdt*(dydr*dzds - dzdr*dyds));

    surfMetricInt[4][l] = (dyds*dzdt - dzds*dydt)*Jinv;  // drdx
    surfMetricInt[5][l] = (dzds*dxdt - dxds*dzdt)*Jinv;  // drdy
    surfMetricInt[6][l] = (dxds*dydt - dyds*dxdt)*Jinv;  // drdz

    surfMetricInt[7][l] = (dzdr*dydt - dydr*dzdt)*Jinv;  // dsdx
    surfMetricInt[8][l] = (dxdr*dzdt - dzdr*dxdt)*Jinv;  // dsdy
    surfMetricInt[9][l] = (dydr*dxdt - dxdr*dydt)*Jinv;  // dsdz

    surfMetricInt[10][l] = (dydr*dzds - dzdr*dyds)*Jinv;  // dtdx
    surfMetricInt[11][l] = (dzdr*dxds - dxdr*dzds)*Jinv;  // dtdy
    surfMetricInt[12][l] = (dxdr*dyds - dydr*dxds)*Jinv;  // dtdz
  }
}

//------------------------------------------------------------------------------

// Function, which creates the final metric terms for a face in k-direction.
void ElementClass::FinalMetricKFace(const int                nItems,
                                    std::vector<su2double *> &surfMetricInt)
{
  // Loop over the number of items for which the final metric terms
  // must be computed.
#pragma omp simd
  for(int l=0; l<nItems; ++l)
  {
    // On entry the derivatives of the physical coordinates w.r.t. the
    // parametric coordinates are stored in surfMetricInt. Copy these values.
    const su2double dxdr = surfMetricInt[0][l];
    const su2double dydr = surfMetricInt[1][l];
    const su2double dzdr = surfMetricInt[2][l];

    const su2double dxds = surfMetricInt[3][l];
    const su2double dyds = surfMetricInt[4][l];
    const su2double dzds = surfMetricInt[5][l];

    const su2double dxdt = surfMetricInt[6][l];
    const su2double dydt = surfMetricInt[7][l];
    const su2double dzdt = surfMetricInt[8][l];

    // Carry out the cross product dxdr X dxds, where x stands for the
    // vector (x,y,z). This gives the normal of the face pointing in
    // positive t-direction.
    const su2double nx = dydr*dzds - dyds*dzdr;
    const su2double ny = dxds*dzdr - dxdr*dzds;
    const su2double nz = dxdr*dyds - dxds*dydr;

    // Determine the length and the inverse of the length of this normal.
    const su2double lenNorm    = SQRT(nx*nx + ny*ny + nz*nz);
    const su2double invLenNorm = one/lenNorm;

    // Store the unit normal and the length in the first 4 entries of surfMetricInt.
    surfMetricInt[0][l] = nx*invLenNorm;
    surfMetricInt[1][l] = ny*invLenNorm;
    surfMetricInt[2][l] = nz*invLenNorm;
    surfMetricInt[3][l] = lenNorm;

    // Compute the derivatives of the parametric coordinates w.r.t. the
    // Cartesian coordinates, i.e. drdx, drdy, etc.
    const su2double Jinv = one/(dxdr*(dyds*dzdt - dzds*dydt)
                         -      dxds*(dydr*dzdt - dzdr*dydt)
                         +      dxdt*(dydr*dzds - dzdr*dyds));

    surfMetricInt[4][l] = (dyds*dzdt - dzds*dydt)*Jinv;  // drdx
    surfMetricInt[5][l] = (dzds*dxdt - dxds*dzdt)*Jinv;  // drdy
    surfMetricInt[6][l] = (dxds*dydt - dyds*dxdt)*Jinv;  // drdz

    surfMetricInt[7][l] = (dzdr*dydt - dydr*dzdt)*Jinv;  // dsdx
    surfMetricInt[8][l] = (dxdr*dzdt - dzdr*dxdt)*Jinv;  // dsdy
    surfMetricInt[9][l] = (dydr*dxdt - dxdr*dydt)*Jinv;  // dsdz

    surfMetricInt[10][l] = (dydr*dzds - dzdr*dyds)*Jinv;  // dtdx
    surfMetricInt[11][l] = (dzdr*dxds - dxdr*dzds)*Jinv;  // dtdy
    surfMetricInt[12][l] = (dxdr*dyds - dydr*dxds)*Jinv;  // dtdz
  }
}

//------------------------------------------------------------------------------

// Function, which carries out a containment search in a high order element,
// where an initial guess can be obtained from linear sub-elements.
void ElementClass::HighOrderContainmentSearch(const su2double              *coorExchange,
                                              const unsigned short         subElem,
                                              const StandardElementClass   *standardHex,
                                              const su2double              *weightsInterpol,
                                              const std::vector<su2double> &VInv,
                                              const std::vector<su2double> &rSample1D,
                                              su2double                    *parCoor)
{
  // Definition of the maximum number of iterations in the Newton solver
  // and the tolerance level.
  const int       maxIt     = 50;
  const su2double tolNewton = (su2double) 1.e-10;

  // Store the number of 1D grid DOFs a bit easier.
  const int nDOFs1DGrid = (int) standardHex->mRDOFsGrid1D.size();

  // Determine the i-, j- and k-index of the subelement.
  const int nElem1D = (int) (rSample1D.size()-1);
  const int nElem2D = nElem1D*nElem1D;

  const int k   = subElem/nElem2D;
  const int rem = subElem - k*nElem2D;
  const int j   = rem/nElem1D;
  const int i   = rem - j*nElem1D;

  // Easier storage of the 1D parametric coordinates of the sample points
  // and interpolation weights.
  const su2double *r = rSample1D.data();
  const su2double *w = weightsInterpol;

  // Determine the initial guess of the parametric coordinates.
  parCoor[0] = (w[0]+w[3]+w[4]+w[7])*r[i] + (w[1]+w[2]+w[5]+w[6])*r[i+1];
  parCoor[1] = (w[0]+w[1]+w[4]+w[5])*r[j] + (w[2]+w[3]+w[6]+w[7])*r[j+1];
  parCoor[2] = (w[0]+w[1]+w[2]+w[3])*r[k] + (w[4]+w[5]+w[6]+w[7])*r[k+1];

  // Define the vectors to either store or compute the
  // Lagrangian basis functions and its derivatives.
  std::vector<su2double> rPoints(1), V(nDOFs1DGrid), dV(nDOFs1DGrid);
  std::vector<su2double> lagR(nDOFs1DGrid), dLagR(nDOFs1DGrid);
  std::vector<su2double> lagS(nDOFs1DGrid), dLagS(nDOFs1DGrid);
  std::vector<su2double> lagT(nDOFs1DGrid), dLagT(nDOFs1DGrid);

  // Loop over the maximum number of iterations.
  int itCount;
  for(itCount=0; itCount<maxIt; ++itCount)
  {
    // Determine the Lagrangian basis function and its derivatives
    // in r-direction.
    rPoints[0] = parCoor[0];
    Vandermonde1D(nDOFs1DGrid, rPoints, V);
    GradVandermonde1D(nDOFs1DGrid, rPoints, dV);

    for(int j=0; j<nDOFs1DGrid; ++j)
    {
      lagR[j] = dLagR[j] = zero;
      for(int k=0; k<nDOFs1DGrid; ++k)
      {
        const int indVInv = j*nDOFs1DGrid + k;
        lagR[j]  +=  V[k]*VInv[indVInv];
        dLagR[j] += dV[k]*VInv[indVInv];
      }
    }

    // Determine the Lagrangian basis function and its derivatives
    // in s-direction.
    rPoints[0] = parCoor[1];
    Vandermonde1D(nDOFs1DGrid, rPoints, V);
    GradVandermonde1D(nDOFs1DGrid, rPoints, dV);

    for(int j=0; j<nDOFs1DGrid; ++j)
    {
      lagS[j] = dLagS[j] = zero;
      for(int k=0; k<nDOFs1DGrid; ++k)
      {
        const int indVInv = j*nDOFs1DGrid + k;
        lagS[j]  +=  V[k]*VInv[indVInv];
        dLagS[j] += dV[k]*VInv[indVInv];
      }
    }

    // Determine the Lagrangian basis function and its derivatives
    // in t-direction.
    rPoints[0] = parCoor[2];
    Vandermonde1D(nDOFs1DGrid, rPoints, V);
    GradVandermonde1D(nDOFs1DGrid, rPoints, dV);

    for(int j=0; j<nDOFs1DGrid; ++j)
    {
      lagT[j] = dLagT[j] = zero;
      for(int k=0; k<nDOFs1DGrid; ++k)
      {
        const int indVInv = j*nDOFs1DGrid + k;
        lagT[j]  +=  V[k]*VInv[indVInv];
        dLagT[j] += dV[k]*VInv[indVInv];
      }
    }

    // Initialize the functional to be computed and its Jacobian.
    su2double f0 = coorExchange[0], f1 = coorExchange[1], f2 = coorExchange[2];
    su2double a00 = zero, a01 = zero, a02 = zero;
    su2double a10 = zero, a11 = zero, a12 = zero;
    su2double a20 = zero, a21 = zero, a22 = zero;

    // Loop over the grid DOFs to compute the functional and minus
    // the Jacobian.
    int ind = 0;
    for(int k=0; k<nDOFs1DGrid; ++k)
    {
      for(int j=0; j<nDOFs1DGrid; ++j)
      {
        const su2double  lslt =  lagS[j]* lagT[k];
        const su2double dlslt = dLagS[j]* lagT[k];
        const su2double lsdlt =  lagS[j]*dLagT[k];

        for(int i=0; i<nDOFs1DGrid; ++i, ++ind)
        {
          const su2double x = mCoorNodalGridDOFs[0][ind];
          const su2double y = mCoorNodalGridDOFs[1][ind];
          const su2double z = mCoorNodalGridDOFs[2][ind];

          const su2double  lrlslt =  lagR[i]* lslt;
          const su2double dlrlslt = dLagR[i]* lslt;
          const su2double lrdlslt =  lagR[i]*dlslt;
          const su2double lrlsdlt =  lagR[i]*lsdlt;

          f0 -= x*lrlslt; f1 -= y*lrlslt; f2 -= z*lrlslt;

          a00 += x*dlrlslt; a01 += x*lrdlslt; a02 += x*lrlsdlt;
          a10 += y*dlrlslt; a11 += y*lrdlslt; a12 += y*lrlsdlt;
          a20 += z*dlrlslt; a21 += z*lrdlslt; a22 += z*lrlsdlt;
        }
      }
    }

    // Compute the updates of the parametric values. As minus the
    // Jacobian is computed, the updates should be added to parCoor.
    const su2double detInv = one/(a00*a11*a22 - a00*a12*a21 - a01*a10*a22
                           +      a01*a12*a20 + a02*a10*a21 - a02*a11*a20);
    const su2double dr =  detInv*(a01*a12*f2 - a01*a22*f1 - a02*a11*f2
                       +          a02*a21*f1 + a11*a22*f0 - a12*a21*f0);
    const su2double ds = -detInv*(a00*a12*f2 - a00*a22*f1 - a02*a10*f2
                       +          a02*a20*f1 + a10*a22*f0 - a12*a20*f0);
    const su2double dt =  detInv*(a00*a11*f2 - a00*a21*f1 - a01*a10*f2
                       +          a01*a20*f1 + a10*a21*f0 - a11*a20*f0);
    parCoor[0] += dr;
    parCoor[1] += ds;
    parCoor[2] += dt;

    // Check for convergence.
    if(FABS(dr) <= tolNewton && FABS(ds) <= tolNewton && FABS(dt) <= tolNewton) break;
  }

  // Terminate if the Newton algorithm did not converge.
  if(itCount == maxIt)
    Terminate("ElementClass::HighOrderContainmentSearch", __FILE__, __LINE__,
              "Newton did not converge");
}

//------------------------------------------------------------------------------

// Function, which initializes the solution in the DOFs of the element.
void ElementClass::InitSol(const InputParamClass      *inputParam,
                           const StandardElementClass *standardHex,
                           const su2double            *primVar)
{
  // Initialize the modal solution to zero.
  for(int m=0; m<nVar; ++m)
#pragma omp simd
    for(int l=0; l<standardHex->mNDOFsPad; ++l)
      mSol[m][l] = zero;

  // Easier storage of the primitive variables.
  const su2double rho = primVar[0];
  const su2double u   = primVar[1];
  const su2double v   = primVar[2];
  const su2double w   = primVar[3];
  const su2double p   = primVar[4];

  // Convert the primitive variables to the working variables. Make a
  // distinction between conservative and entropy variables.
  su2double workVar[nVar];
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      // The conservative variables are used as working variables.
      // Compute the conservative variables from the primitive ones.
      const su2double ovgm1 = one/(GamConstant - one);

      workVar[0] = rho;
      workVar[1] = rho*u;
      workVar[2] = rho*v;
      workVar[3] = rho*w;
      workVar[4] = p*ovgm1 + half*rho*(u*u + v*v + w*w);

      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      // The entropy variables are used as working variables.
      // Compute the entropy variables from the primitive ones.
      const su2double ovgm1 = one/(GamConstant - one);
      const su2double s     = LOG(p/POW(rho,GamConstant));
      const su2double pInv  = one/p;

      workVar[0] =  (GamConstant-s)*ovgm1 - half*pInv*rho*(u*u + v*v + w*w);
      workVar[1] =  rho*u*pInv;
      workVar[2] =  rho*v*pInv;
      workVar[3] =  rho*w*pInv;
      workVar[4] = -rho*pInv;

      break;
    }

    //--------------------------------------------------------------------------

    default:
    {
      // This is just to avoid a compiler warning.
      break;
    }
  }

  // The flow is initialized to a constant solution, i.e. in the modal form
  // only the first DOF is nonzero. However, it must be scaled, because the
  // basis functions are orthonormal over the standard element and therefore
  // the first basis funcion is not equal to 1. Compute this scale factor.
  const su2double V0        = NormJacobi(0, 0, 0, zero);
  const su2double scaleFact = one/(V0*V0*V0);

  // Set the first DOF of the modal solution to the working variables.
  for(int m=0; m<nVar; ++m)
    mSol[m][0] = workVar[m]*scaleFact;
}

//------------------------------------------------------------------------------

// Function, which determines the interpolation weights for the exchange location
// of the wall model of the wall boundary faces, if needed.
void ElementClass::InterpolationWeightsExchangeLocation(
                    const InputParamClass                                    *inputParam,
                    const StandardElementClass                               *standardHex,
                    std::vector<std::vector< std::vector<ElementClass *> > > &elements,
                    const int                                                nElemPerRankI,
                    const int                                                nElemPerRankJ,
                    const int                                                nElemPerRankK,
                    CADTElemClass                                            *localVolumeADT,
                    const std::vector<su2double>                             &rSample1D)
{
  // Check the 6 faces of the element for a wall model treatment and call
  // the function InterpolationWeightsExchangeLocationFace to carry out
  // the actual work. First the iMin boundary.
  if( mBCIMin )
  {
    if( mBCIMin->UseWallModel() )
      InterpolationWeightsExchangeLocationFace(mSurfMetricIntIMin.data(), one,
                                               inputParam, standardHex, elements,
                                               nElemPerRankI, nElemPerRankJ,
                                               nElemPerRankK, localVolumeADT,
                                               mExchangeDataIMin, rSample1D);
  }
 
  // The iMax boundary.
  if( mBCIMax )
  {
    if( mBCIMax->UseWallModel() )
      InterpolationWeightsExchangeLocationFace(mSurfMetricIntIMax.data(), -one,
                                               inputParam, standardHex, elements,
                                               nElemPerRankI, nElemPerRankJ,
                                               nElemPerRankK, localVolumeADT,
                                               mExchangeDataIMax, rSample1D);
  }

  // The jMin boundary.
  if( mBCJMin )
  {
    if( mBCJMin->UseWallModel() )
      InterpolationWeightsExchangeLocationFace(mSurfMetricIntJMin.data(), one,
                                               inputParam, standardHex, elements,
                                               nElemPerRankI, nElemPerRankJ,
                                               nElemPerRankK, localVolumeADT,
                                               mExchangeDataJMin, rSample1D);
  }
 
  // The jMax boundary.
  if( mBCJMax )
  {
    if( mBCJMax->UseWallModel() )
      InterpolationWeightsExchangeLocationFace(mSurfMetricIntJMax.data(), -one,
                                               inputParam, standardHex, elements,
                                               nElemPerRankI, nElemPerRankJ,
                                               nElemPerRankK, localVolumeADT,
                                               mExchangeDataJMax, rSample1D);
  }

  // The kMin boundary.
  if( mBCKMin )
  {
    if( mBCKMin->UseWallModel() )
      InterpolationWeightsExchangeLocationFace(mSurfMetricIntKMin.data(), one,
                                               inputParam, standardHex, elements,
                                               nElemPerRankI, nElemPerRankJ,
                                               nElemPerRankK, localVolumeADT,
                                               mExchangeDataKMin, rSample1D);
  }
 
  // The kMax boundary.
  if( mBCKMax )
  {
    if( mBCKMax->UseWallModel() )
      InterpolationWeightsExchangeLocationFace(mSurfMetricIntKMax.data(), -one,
                                               inputParam, standardHex, elements,
                                               nElemPerRankI, nElemPerRankJ,
                                               nElemPerRankK, localVolumeADT,
                                               mExchangeDataKMax, rSample1D);
  }
}

//------------------------------------------------------------------------------

// Function, which determines the interpolation weights for the exchange locations
// corresponding to the integration points of the given face.
void ElementClass::InterpolationWeightsExchangeLocationFace(
                    su2double                                                **faceMetric,
                    const su2double                                          factNorm,
                    const InputParamClass                                    *inputParam,
                    const StandardElementClass                               *standardHex,
                    std::vector<std::vector< std::vector<ElementClass *> > > &elements,
                    const int                                                nElemPerRankI,
                    const int                                                nElemPerRankJ,
                    const int                                                nElemPerRankK,
                    CADTElemClass                                            *localVolumeADT,
                    ExchangeDataWallModelClass                               &exchangeData,
                    const std::vector<su2double>                             &rSample1D)
{
  // Determine the inverse of the Vandermonde matrix in the 1D grid DOFs.
  const int nDOFs1DGrid = (int) standardHex->mRDOFsGrid1D.size();
  std::vector<su2double> VInv(nDOFs1DGrid*nDOFs1DGrid);

  Vandermonde1D(nDOFs1DGrid, standardHex->mRDOFsGrid1D, VInv);
  InverseMatrix(nDOFs1DGrid, VInv);

  // Store the number of integration points of the face a bit easier and determine
  // the number of elements per rank in the i-j plane.
  const int nInt2D         = standardHex->mNIntegration2D;
  const int nElemPerRankIJ = nElemPerRankI*nElemPerRankJ;

  // Define the vector to store the donor information.
  std::vector<CExchangeInfoWMClass> donorElements(nInt2D);

  // Loop over the number of integration points.
  for(int l=0; l<nInt2D; ++l)
  {
    // Determine the multiplication factor for the unit normal to obtain
    // the distance vector from the integration point.
    const su2double dist = factNorm*inputParam->mExchangeLocation;

    // Determine the coordinate of the exchange point.
    su2double coorExchange[] = {faceMetric[13][l] + dist*faceMetric[0][l],
                                faceMetric[14][l] + dist*faceMetric[1][l],
                                faceMetric[15][l] + dist*faceMetric[2][l]};

    // Search for the element, which contains the exchange location.
    unsigned short subElem;
    unsigned long  parElem;
    int            rank;
    su2double      parCoor[3], weightsInterpol[8];
    if( localVolumeADT->DetermineContainingElement(coorExchange, subElem,
                                                   parElem, rank, parCoor,
                                                   weightsInterpol) )
    {
      // Subelement found that contains the exchange location. However,
      // what is needed is the location in the high order parent element.
      // First determine the i-, j- and k-index of the parent element.
      int k   = parElem/nElemPerRankIJ;
      int rem = parElem - k*nElemPerRankIJ;
      int j   = rem/nElemPerRankI;
      int i   = rem - j*nElemPerRankI;

      // Increment i, j and k, because the owned elements start at 1.
      ++i; ++j; ++k;

      // Determine the location of the exchange coordinate in the
      // higher order element.
      elements[i][j][k]->HighOrderContainmentSearch(coorExchange, subElem, standardHex,
                                                    weightsInterpol, VInv, rSample1D,
                                                    parCoor);

      // Store the data for this integration point.
      donorElements[l].mI = i;
      donorElements[l].mJ = j;
      donorElements[l].mK = k;

      donorElements[l].mIntID = l;

      donorElements[l].mParCoor[0] = parCoor[0];
      donorElements[l].mParCoor[1] = parCoor[1];
      donorElements[l].mParCoor[2] = parCoor[2];
    }
    else
    {
      // No subelement in the local partition that contains the exchange point.
      // Create an error message and exit.
      std::ostringstream message;
      message << "Exchange location for the wall modeling not found in the current partition."
              << std::endl
              << "Reduce the number of MPI ranks or decrease the exchange location.";
      Terminate("ElementClass::InterpolationWeightsExchangeLocationFace",
                __FILE__, __LINE__, message.str());
    }
  }

  // Sort the donor info in increasing order.
  std::sort(donorElements.begin(), donorElements.end());

  // Determine the number of exchange points per donor element. The vector
  // mNExchangePointsPerDonorElement will be stored in cumulative storage
  // format, hence its first element is zero.
  exchangeData.mNExchangePointsPerDonorElement.resize(2,0);
  exchangeData.mNExchangePointsPerDonorElement[1] = 1;

  exchangeData.mIndicesDonorElements.push_back(donorElements[0].mI);
  exchangeData.mIndicesDonorElements.push_back(donorElements[0].mJ);
  exchangeData.mIndicesDonorElements.push_back(donorElements[0].mK);

  for(int l=1; l<nInt2D; ++l)
  {
    if((donorElements[l].mI != donorElements[l-1].mI) ||
       (donorElements[l].mJ != donorElements[l-1].mJ) ||
       (donorElements[l].mK != donorElements[l-1].mK))
    {
      exchangeData.mNExchangePointsPerDonorElement.push_back(l+1);

      exchangeData.mIndicesDonorElements.push_back(donorElements[l].mI);
      exchangeData.mIndicesDonorElements.push_back(donorElements[l].mJ);
      exchangeData.mIndicesDonorElements.push_back(donorElements[l].mK);
    }
    else
    {
      ++exchangeData.mNExchangePointsPerDonorElement.back();
    }
  }

  // Easier storage of the number of solution DOFs in 1D and 3D.
  const int nDOFs1DSol = standardHex->mNDOFs1D;
  const int nDOFs3DSol = standardHex->mNDOFs;

  // Allocate the memory to store the actual interpolation information.
  exchangeData.mExchangePointsPerDonorElement.resize(nInt2D);

  exchangeData.mWeightsExchangePoints.resize(nDOFs3DSol, NULL);
  for(int l=0; l<nDOFs3DSol; ++l)
  {
    exchangeData.mWeightsExchangePoints[l] = (su2double *) AllocateMemory(nInt2D*sizeof(su2double));
    if( !exchangeData.mWeightsExchangePoints[l] )
      Terminate("ElementClass::InterpolationWeightsExchangeLocationFace", __FILE__,
                __LINE__, "Memory allocation failure for mWeightsExchangePoints");
  }

  // Loop over the integration points to compute and store the interpolation weights.
  for(int l=0; l<nInt2D; ++l)
  {
    exchangeData.mExchangePointsPerDonorElement[l] = donorElements[l].mIntID;

    // Loop over the solution DOFs to create the interpolation data.
    int ind = 0;
    for(int k=0; k<nDOFs1DSol; ++k)
    {
      const su2double wK = NormJacobi(k, 0, 0, donorElements[l].mParCoor[2]);
      for(int j=0; j<nDOFs1DSol; ++j)
      {
        const su2double wJ = NormJacobi(j, 0, 0, donorElements[l].mParCoor[1]);
        for(int i=0; i<nDOFs1DSol; ++i, ++ind)
        {
          const su2double wI = NormJacobi(i, 0, 0, donorElements[l].mParCoor[0]);
          exchangeData.mWeightsExchangePoints[ind][l] = wK*wJ*wI;
        }
      }
    }
  }
}

//------------------------------------------------------------------------------

// Function, which computes the residuals of the i face(s) of this element.
void ElementClass::IFaceResiduals(
            const InputParamClass                                     *inputParam,
            const StandardElementClass                                *standardHex,
            std::vector<std::vector< std::vector<ElementClass *> > >  &elements,
            su2double                                                **workArray,
            const bool                                                 ComputeMonitoringData,
            su2double                                                 &EddyVisMax,
            su2double                                                 *forceCoef)
{
  // Easier storage of the indices of the current element.
  const int i = mLocalInd[0], j = mLocalInd[1], k = mLocalInd[2];

  // Initialize the length scale for the penalty terms to the spacing in i-direction
  // of the current element and initialize the LES length scale to the LES length
  // scale of the current element.
  su2double lenScale    = mLenScaleIDir;
  su2double lenScaleLES = mLenScaleLES;

  // Easier storage of the number of DOFs and integration points.
  const int nDOFs1D  = standardHex->mNDOFs1D;
  const int nDOFsPad = standardHex->mNDOFsPad;
  const int nInt1D   = standardHex->mNIntegration1D;
  const int nInt     = standardHex->mNIntegration2D;
  const int nIntPad  = standardHex->mNIntegration2DPad;

  // Set the pointers for the left and right solution and its gradients
  // in the integration points of the face.
  su2double **solL    = workArray;
  su2double **dSolDxL = solL    + nVar;
  su2double **dSolDyL = dSolDxL + nVar;
  su2double **dSolDzL = dSolDyL + nVar;

  su2double **solR    = dSolDzL + nVar;
  su2double **dSolDxR = solR    + nVar;
  su2double **dSolDyR = dSolDxR + nVar;
  su2double **dSolDzR = dSolDyR + nVar;

  // Set the pointer for the total normal flux at the interface. The total
  // flux contains the contributions from the inviscid flux, the viscous flux
  // and the penalty flux. It does not include the symmetrizing fluxes,
  // because these are distributed differently.
  su2double **fluxTot = dSolDzR + nVar;

	// Set the pointers for the left and right symmetrizing fluxes.
	su2double **dFluxSymDxL = fluxTot     + nVar;
	su2double **dFluxSymDyL = dFluxSymDxL + nVar;
	su2double **dFluxSymDzL = dFluxSymDyL + nVar;

	su2double **dFluxSymDxR = dFluxSymDzL + nVar;
	su2double **dFluxSymDyR = dFluxSymDxR + nVar;
	su2double **dFluxSymDzR = dFluxSymDyR + nVar;

  // Set the pointer for the eddy viscosities, which are needed when an
  // LES subgrid scale model is used.
  su2double *eddyVis = dFluxSymDzR[nVar];

  // Compute the right solution and its gradients in the integration points
  // of the face. The right solution corresponds to the iMin face of this
  // element, hence the Legendre basis functions and its derivatives must be
  // taken at the min boundary.
  TensorProductSolAndGradIFace(nInt1D, nVar, nDOFs1D,
                               standardHex->mLegendreInt1D,
                               standardHex->mDerLegendreInt1D,
                               standardHex->mLegendreMinFace1D,
                               standardHex->mDerLegendreMinFace1D,
                               mSol.data(), 
															 solR, dSolDxR, dSolDyR, dSolDzR);

  // Determine what type of face we are dealing with.
  if( elements[i-1][j][k] )
  {
    // This is an internal face. Modify the length scales.
    lenScale    = std::min(lenScale, elements[i-1][j][k]->mLenScaleIDir);
    lenScaleLES = half*(lenScaleLES + elements[i-1][j][k]->mLenScaleLES);

    // Compute the left solution and its gradients in the integration points
    // of the face. The left solution corresponds to the iMax face of the left
    // neighbor, hence the Legendre basis functions and its derivatives must
    // be taken at the max boundary.
    TensorProductSolAndGradIFace(nInt1D, nVar, nDOFs1D,
                                 standardHex->mLegendreInt1D,
                                 standardHex->mDerLegendreInt1D,
                                 standardHex->mLegendreMaxFace1D,
                                 standardHex->mDerLegendreMaxFace1D,
                                 elements[i-1][j][k]->mSol.data(), 
																 solL, dSolDxL, dSolDyL, dSolDzL);

    // Compute the fluxes for an internal face.
    FluxesInternalFace(inputParam, nInt, nIntPad, lenScale, lenScaleLES,
                       standardHex->mIntWeightsFace, 
											 solL, dSolDxL, dSolDyL, dSolDzL, 
											 elements[i-1][j][k]->mSurfMetricIntIMax.data(),
                       solR, dSolDxR, dSolDyR, dSolDzR, 
											 mSurfMetricIntIMin.data(),
                       eddyVis, fluxTot,
											 dFluxSymDxL, dFluxSymDyL, dFluxSymDzL,
											 dFluxSymDxR, dFluxSymDyR, dFluxSymDzR, 
											 ComputeMonitoringData, EddyVisMax);

    // Initialize the residual on the iMin boundary, which is the residual
    // for the left neighbor, to zero.
    for(int l=0; l<nVar; ++l)
#pragma omp simd
      for(int i=0; i<nDOFsPad; ++i)
        mResIMin[l][i] = zero;

    // Distribute the normal fluxes to the left element. The face corresponds
    // to the iMax face of the left element, hence the value of the Legendre
    // basis functions at the max face must be used in i-direction.
    TensorProductResIFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreMaxFace1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreInt1DTranspose,
                          fluxTot, mResIMin.data());

    // Distribute the symmetrizing fluxes in i-direction, stored in dFluxSymDxL
    // to the left element. These fluxes are multiplied by the derivative of
    // the basis functions in i-direction.
    TensorProductResIFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mDerLegendreMaxFace1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreInt1DTranspose,
                          dFluxSymDxL, mResIMin.data());

    // Distribute the symmetrizing fluxes in j-direction, stored in dFluxSymDyL,
    // to the left element. These fluxes are multiplied by the derivative of
    // the basis functions in j-direction.
    TensorProductResIFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreMaxFace1D,
                          standardHex->mDerLegendreInt1DTranspose,
                          standardHex->mLegendreInt1DTranspose,
                          dFluxSymDyL, mResIMin.data());

    // Distribute the symmetrizing fluxes in k-direction, stored in dFluxSymDzL,
    // to the left element. These fluxes are multiplied by the derivative of
    // the basis functions in k-direction.
    TensorProductResIFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreMaxFace1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mDerLegendreInt1DTranspose,
                          dFluxSymDzL, mResIMin.data());

    // Negate the normal fluxes, as these will have an opposite sign for
    // the current element.
    for(int l=0; l<nVar; ++l)
#pragma omp simd
      for(int i=0; i<nIntPad; ++i)
        fluxTot[l][i] = -fluxTot[l][i];

  }
  else
  {
    // This is a boundary face. Compute the fluxes for a boundary face. Note, it does not
		// matter if we switch the dFluxSym R and L state, since they are overwritten then the
		// values are returned the same for both states.
    FluxesBoundaryFace(inputParam, standardHex, elements, nInt, nIntPad, mBCIMin, lenScale,
                       lenScaleLES, standardHex->mIntWeightsFace, 
											 solR, dSolDxR, dSolDyR, dSolDzR, 
											-one, mSurfMetricIntIMin.data(), 
											 &mExchangeDataIMin, mPrescribedDataIMin.data(), 
											 solL, dSolDxL, dSolDyL, dSolDzL,
                       eddyVis, fluxTot, 
											 dFluxSymDxR, dFluxSymDyR, dFluxSymDzR,
											 dFluxSymDxL, dFluxSymDyL, dFluxSymDzL,
											 ComputeMonitoringData, EddyVisMax, forceCoef);
  }

  // Distribute the normal fluxes to the current element. The face corresponds
  // to the iMin face of the current element, hence the value of the Legendre
  // basis functions at the min face must be used in i-direction.
  TensorProductResIFace(nInt1D, nVar, nDOFs1D,
                        standardHex->mLegendreMinFace1D,
                        standardHex->mLegendreInt1DTranspose,
                        standardHex->mLegendreInt1DTranspose,
                        fluxTot, mRes.data());

  // Distribute the symmetrizing fluxes in i-direction, stored in dFluxSymDxR,
  // to the current element. These fluxes are multiplied by the derivative
  // of the basis functions in i-direction.
  TensorProductResIFace(nInt1D, nVar, nDOFs1D,
                        standardHex->mDerLegendreMinFace1D,
                        standardHex->mLegendreInt1DTranspose,
                        standardHex->mLegendreInt1DTranspose,
                        dFluxSymDxR, mRes.data());

  // Distribute the symmetrizing fluxes in j-direction, stored in dFluxSymDyR,
  // to the current element. These fluxes are multiplied by the derivative
  // of the basis functions in j-direction.
  TensorProductResIFace(nInt1D, nVar, nDOFs1D,
                        standardHex->mLegendreMinFace1D,
                        standardHex->mDerLegendreInt1DTranspose,
                        standardHex->mLegendreInt1DTranspose,
                        dFluxSymDyR, mRes.data());

  // Distribute the symmetrizing fluxes in k-direction, stored in dFluxSymDzR,
  // to the current element. These fluxes are multiplied by the derivative
  // of the basis functions in k-direction.
  TensorProductResIFace(nInt1D, nVar, nDOFs1D,
                        standardHex->mLegendreMinFace1D,
                        standardHex->mLegendreInt1DTranspose,
                        standardHex->mDerLegendreInt1DTranspose,
                        dFluxSymDzR, mRes.data());

  // Check if the iMax face is a boundary face.
  if( !elements[i+1][j][k] )
  {
    // Reset the length scales.
    lenScale    = mLenScaleIDir;
    lenScaleLES = mLenScaleLES;

    // Compute the left solution and its gradients in the integration points
    // of the face. The left solution corresponds to the iMax face of the current
    // element, hence the Legendre basis functions and its derivatives must be
    // taken at the max boundary.
    TensorProductSolAndGradIFace(nInt1D, nVar, nDOFs1D,
                                 standardHex->mLegendreInt1D,
                                 standardHex->mDerLegendreInt1D,
                                 standardHex->mLegendreMaxFace1D,
                                 standardHex->mDerLegendreMaxFace1D,
                                 mSol.data(), 
																 solL, dSolDxL, dSolDyL, dSolDzL);

    // Compute the fluxes for a boundary face.
    FluxesBoundaryFace(inputParam, standardHex, elements, nInt, nIntPad, mBCIMax, lenScale,
                       lenScaleLES, standardHex->mIntWeightsFace, 
											 solL, dSolDxL, dSolDyL, dSolDzL, 
											 one, mSurfMetricIntIMax.data(), 
											 &mExchangeDataIMax, mPrescribedDataIMax.data(), 
											 solR, dSolDxR, dSolDyR, dSolDzR,
                       eddyVis, fluxTot, 
											 dFluxSymDxL, dFluxSymDyL, dFluxSymDzL,
											 dFluxSymDxR, dFluxSymDyR, dFluxSymDzR,
											 ComputeMonitoringData, EddyVisMax, forceCoef);

    // Distribute the normal fluxes to the current element. The face corresponds
    // to the iMax face of the current element, hence the value of the Legendre
    // basis functions at the max face must be used in i-direction.
    TensorProductResIFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreMaxFace1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreInt1DTranspose,
                          fluxTot, mRes.data());

    // Distribute the symmetrizing fluxes in i-direction, stored in dFluxSymDxL,
    // to the current element. These fluxes are multiplied by the derivative of
    // the basis functions in i-direction.
    TensorProductResIFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mDerLegendreMaxFace1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreInt1DTranspose,
                          dFluxSymDxL, mRes.data());

    // Distribute the symmetrizing fluxes in j-direction, stored in dFluxSymDyL,
    // to the current element. These fluxes are multiplied by the derivative of
    // the basis functions in j-direction.
    TensorProductResIFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreMaxFace1D,
                          standardHex->mDerLegendreInt1DTranspose,
                          standardHex->mLegendreInt1DTranspose,
                          dFluxSymDyL, mRes.data());

    // Distribute the symmetrizing fluxes in k-direction, stored in dFluxSymDzL,
    // to the current element. These fluxes are multiplied by the derivative of
    // the basis functions in k-direction.
    TensorProductResIFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreMaxFace1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mDerLegendreInt1DTranspose,
                          dFluxSymDzL, mRes.data());
  }
}

//------------------------------------------------------------------------------

// Function, which computes the residuals of the j face(s) of this element.
void ElementClass::JFaceResiduals(
            const InputParamClass                                    *inputParam,
            const StandardElementClass                               *standardHex,
            std::vector<std::vector< std::vector<ElementClass *> > > &elements,
            su2double                                                **workArray,
            const bool                                               ComputeMonitoringData,
            su2double                                                &EddyVisMax,
            su2double                                                *forceCoef)
{
  // Easier storage of the indices of the current element.
  const int i = mLocalInd[0], j = mLocalInd[1], k = mLocalInd[2];

  // Initialize the length scale for the penalty terms to the spacing in j-direction
  // of the current element and initialize the LES length scale to the LES length
  // scale of the current element.
  su2double lenScale    = mLenScaleJDir;
  su2double lenScaleLES = mLenScaleLES;

  // Easier storage of the number of DOFs and integration points.
  const int nDOFs1D  = standardHex->mNDOFs1D;
  const int nDOFsPad = standardHex->mNDOFsPad;
  const int nInt1D   = standardHex->mNIntegration1D;
  const int nInt     = standardHex->mNIntegration2D;
  const int nIntPad  = standardHex->mNIntegration2DPad;

  // Set the pointers for the left and right solution and its gradients
  // in the integration points of the face.
  su2double **solL    = workArray;
  su2double **dSolDxL = solL    + nVar;
  su2double **dSolDyL = dSolDxL + nVar;
  su2double **dSolDzL = dSolDyL + nVar;

  su2double **solR    = dSolDzL + nVar;
  su2double **dSolDxR = solR    + nVar;
  su2double **dSolDyR = dSolDxR + nVar;
  su2double **dSolDzR = dSolDyR + nVar;

  // Set the pointer for the total normal flux at the interface. The total
  // flux contains the contributions from the inviscid flux, the viscous flux
  // and the penalty flux. It does not include the symmetrizing fluxes,
  // because these are distributed differently.
  su2double **fluxTot = dSolDzR + nVar;

	// Set the pointers for the left and right symmetrizing fluxes.
	su2double **dFluxSymDxL = fluxTot     + nVar;
	su2double **dFluxSymDyL = dFluxSymDxL + nVar;
	su2double **dFluxSymDzL = dFluxSymDyL + nVar;

	su2double **dFluxSymDxR = dFluxSymDzL + nVar;
	su2double **dFluxSymDyR = dFluxSymDxR + nVar;
	su2double **dFluxSymDzR = dFluxSymDyR + nVar;

  // Set the pointer for the eddy viscosities, which are needed when an
  // LES subgrid scale model is used.
  su2double *eddyVis = dFluxSymDzR[nVar];

  // Compute the right solution and its gradients in the integration points
  // of the face. The right solution corresponds to the jMin face of this
  // element, hence the Legendre basis functions and its derivatives must be
  // taken at the min boundary.
  TensorProductSolAndGradJFace(nInt1D, nVar, nDOFs1D,
                               standardHex->mLegendreInt1D,
                               standardHex->mDerLegendreInt1D,
                               standardHex->mLegendreMinFace1D,
                               standardHex->mDerLegendreMinFace1D,
                               mSol.data(), 
															 solR, dSolDxR, dSolDyR, dSolDzR);

  // Determine what type of face we are dealing with.
  if( elements[i][j-1][k] )
  {
    // This is an internal face. Modify the length scales.
    lenScale    = std::min(lenScale, elements[i][j-1][k]->mLenScaleJDir);
    lenScaleLES = half*(lenScaleLES + elements[i][j-1][k]->mLenScaleLES);

    // Compute the left solution and its gradients in the integration points
    // of the face. The left solution corresponds to the jMax face of the bottom
    // neighbor, hence the Legendre basis functions and its derivatives must
    // be taken at the max boundary.
    TensorProductSolAndGradJFace(nInt1D, nVar, nDOFs1D,
                                 standardHex->mLegendreInt1D,
                                 standardHex->mDerLegendreInt1D,
                                 standardHex->mLegendreMaxFace1D,
                                 standardHex->mDerLegendreMaxFace1D,
                                 elements[i][j-1][k]->mSol.data(), 
																 solL, dSolDxL, dSolDyL, dSolDzL);

    // Compute the fluxes for an internal face.
    FluxesInternalFace(inputParam, nInt, nIntPad, lenScale, lenScaleLES,
                       standardHex->mIntWeightsFace, 
											 solL, dSolDxL, dSolDyL, dSolDzL, 
											 elements[i][j-1][k]->mSurfMetricIntJMax.data(),
                       solR, dSolDxR, dSolDyR, dSolDzR, 
											 mSurfMetricIntJMin.data(),
                       eddyVis, fluxTot, 
											 dFluxSymDxL, dFluxSymDyL, dFluxSymDzL,
											 dFluxSymDxR, dFluxSymDyR, dFluxSymDzR,
											 ComputeMonitoringData, EddyVisMax);

    // Initialize the residual on the jMin boundary, which is the residual
    // for the bottom neighbor, to zero.
    for(int l=0; l<nVar; ++l)
#pragma omp simd
      for(int i=0; i<nDOFsPad; ++i)
        mResJMin[l][i] = zero;

    // Distribute the normal fluxes to the bottom element. The face corresponds
    // to the jMax face of the bottom element, hence the value of the Legendre
    // basis functions at the max face must be used in j-direction.
    TensorProductResJFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreMaxFace1D,
                          standardHex->mLegendreInt1DTranspose,
                          fluxTot, mResJMin.data());

    // Distribute the symmetrizing fluxes in i-direction, stored in dFluxSymDxL,
    // to the bottom element. These fluxes are multiplied by the derivative of
    // the basis functions in i-direction.
    TensorProductResJFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mDerLegendreInt1DTranspose,
                          standardHex->mLegendreMaxFace1D,
                          standardHex->mLegendreInt1DTranspose,
                          dFluxSymDxL, mResJMin.data());

    // Distribute the symmetrizing fluxes in j-direction, stored in dFluxSymDyL,
    // to the bottom element. These fluxes are multiplied by the derivative of
    // the basis functions in j-direction.
    TensorProductResJFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mDerLegendreMaxFace1D,
                          standardHex->mLegendreInt1DTranspose,
                          dFluxSymDyL, mResJMin.data());

    // Distribute the symmetrizing fluxes in k-direction, stored in dFluxSymDzL,
    // to the bottom element. These fluxes are multiplied by the derivative of
    // the basis functions in k-direction.
    TensorProductResJFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreMaxFace1D,
                          standardHex->mDerLegendreInt1DTranspose,
                          dFluxSymDzL, mResJMin.data());
    
		// Negate the normal fluxes, as these will have an opposite sign for
    // the current element.
    for(int l=0; l<nVar; ++l)
#pragma omp simd
      for(int i=0; i<nIntPad; ++i)
        fluxTot[l][i] = -fluxTot[l][i];
  }
  else
  {
    // This is a boundary face. Compute the fluxes for a boundary face. Note, it does not
		// matter if we switch the dFluxSym R and L state, since they are overwritten then the
		// values are returned the same for both states.
    FluxesBoundaryFace(inputParam, standardHex, elements, nInt, nIntPad, mBCJMin, lenScale,
                       lenScaleLES, standardHex->mIntWeightsFace, 
											 solR, dSolDxR, dSolDyR, dSolDzR, 
											-one, mSurfMetricIntJMin.data(), 
											 &mExchangeDataJMin, mPrescribedDataJMin.data(), 
											 solL, dSolDxL, dSolDyL, dSolDzL,
                       eddyVis, fluxTot,
											 dFluxSymDxR, dFluxSymDyR, dFluxSymDzR,
											 dFluxSymDxL, dFluxSymDyL, dFluxSymDzL,	 
											 ComputeMonitoringData, EddyVisMax, forceCoef);
  }

  // Distribute the normal fluxes to the current element. The face corresponds
  // to the jMin face of the current element, hence the value of the Legendre
  // basis functions at the min face must be used in j-direction.
  TensorProductResJFace(nInt1D, nVar, nDOFs1D,
                        standardHex->mLegendreInt1DTranspose,
                        standardHex->mLegendreMinFace1D,
                        standardHex->mLegendreInt1DTranspose,
                        fluxTot, mRes.data());

  // Distribute the symmetrizing fluxes in i-direction, stored in dFluxSymDxR,
  // to the current element. These fluxes are multiplied by the derivative
  // of the basis functions in i-direction.
  TensorProductResJFace(nInt1D, nVar, nDOFs1D,
                        standardHex->mDerLegendreInt1DTranspose,
                        standardHex->mLegendreMinFace1D,
                        standardHex->mLegendreInt1DTranspose,
                        dFluxSymDxR, mRes.data());

  // Distribute the symmetrizing fluxes in j-direction, stored in dFluxSymDyR,
  // to the current element. These fluxes are multiplied by the derivative
  // of the basis functions in j-direction.
  TensorProductResJFace(nInt1D, nVar, nDOFs1D,
                        standardHex->mLegendreInt1DTranspose,
                        standardHex->mDerLegendreMinFace1D,
                        standardHex->mLegendreInt1DTranspose,
                        dFluxSymDyR, mRes.data());

  // Distribute the symmetrizing fluxes in k-direction, stored in dFluxSymDzR,
  // to the current element. These fluxes are multiplied by the derivative
  // of the basis functions in k-direction.
  TensorProductResJFace(nInt1D, nVar, nDOFs1D,
                        standardHex->mLegendreInt1DTranspose,
                        standardHex->mLegendreMinFace1D,
                        standardHex->mDerLegendreInt1DTranspose,
                        dFluxSymDzR, mRes.data());

  // Check if the jMax face is a boundary face.
  if( !elements[i][j+1][k] )
  {
    // Reset the length scales.
    lenScale    = mLenScaleJDir;
    lenScaleLES = mLenScaleLES;

    // Compute the left solution and its gradients in the integration points
    // of the face. The left solution corresponds to the jMax face of the current
    // element, hence the Legendre basis functions and its derivatives must be
    // taken at the max boundary.
    TensorProductSolAndGradJFace(nInt1D, nVar, nDOFs1D,
                                 standardHex->mLegendreInt1D,
                                 standardHex->mDerLegendreInt1D,
                                 standardHex->mLegendreMaxFace1D,
                                 standardHex->mDerLegendreMaxFace1D,
                                 mSol.data(), 
																 solL, dSolDxL, dSolDyL, dSolDzL);

    // Compute the fluxes for a boundary face.
    FluxesBoundaryFace(inputParam, standardHex, elements, nInt, nIntPad, mBCJMax, lenScale,
                       lenScaleLES, standardHex->mIntWeightsFace, 
											 solL, dSolDxL, dSolDyL, dSolDzL, 
											 one, mSurfMetricIntJMax.data(), 
											 &mExchangeDataJMax, mPrescribedDataJMax.data(), 
											 solR, dSolDxR, dSolDyR, dSolDzR,
                       eddyVis, fluxTot, 
											 dFluxSymDxL, dFluxSymDyL, dFluxSymDzL,
											 dFluxSymDxR, dFluxSymDyR, dFluxSymDzR,
											 ComputeMonitoringData, EddyVisMax, forceCoef);

    // Distribute the normal fluxes to the current element. The face corresponds
    // to the jMax face of the current element, hence the value of the Legendre
    // basis functions at the max face must be used in j-direction.
    TensorProductResJFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreMaxFace1D,
                          standardHex->mLegendreInt1DTranspose,
                          fluxTot, mRes.data());

    // Distribute the symmetrizing fluxes in i-direction, stored in dFluxSymDxL,
    // to the current element. These fluxes are multiplied by the derivative of
    // the basis functions in i-direction.
    TensorProductResJFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mDerLegendreInt1DTranspose,
                          standardHex->mLegendreMaxFace1D,
                          standardHex->mLegendreInt1DTranspose,
                          dFluxSymDxL, mRes.data());

    // Distribute the symmetrizing fluxes in j-direction, stored in dFluxSymDyL,
    // to the current element. These fluxes are multiplied by the derivative of
    // the basis functions in j-direction.
    TensorProductResJFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mDerLegendreMaxFace1D,
                          standardHex->mLegendreInt1DTranspose,
                          dFluxSymDyL, mRes.data());

    // Distribute the symmetrizing fluxes in k-direction, stored in dFluxSymDzL,
    // to the current element. These fluxes are multiplied by the derivative of
    // the basis functions in k-direction.
    TensorProductResJFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreMaxFace1D,
                          standardHex->mDerLegendreInt1DTranspose,
                          dFluxSymDzL, mRes.data());
  }
}

//------------------------------------------------------------------------------

// Function, which computes the residuals of the k face(s) of this element.
void ElementClass::KFaceResiduals(
            const InputParamClass                                    *inputParam,
            const StandardElementClass                               *standardHex,
            std::vector<std::vector< std::vector<ElementClass *> > > &elements,
            su2double                                                **workArray,
            const bool                                               ComputeMonitoringData,
            su2double                                                &EddyVisMax,
            su2double                                                *forceCoef)
{
// Easier storage of the indices of the current element.
  const int i = mLocalInd[0], j = mLocalInd[1], k = mLocalInd[2];

  // Initialize the length scale for the penalty terms to the spacing in k-direction
  // of the current element and initialize the LES length scale to the LES length
  // scale of the current element.
  su2double lenScale    = mLenScaleKDir;
  su2double lenScaleLES = mLenScaleLES;

  // Easier storage of the number of DOFs and integration points.
  const int nDOFs1D  = standardHex->mNDOFs1D;
  const int nDOFsPad = standardHex->mNDOFsPad;
  const int nInt1D   = standardHex->mNIntegration1D;
  const int nInt     = standardHex->mNIntegration2D;
  const int nIntPad  = standardHex->mNIntegration2DPad;

  // Set the pointers for the left and right solution and its gradients
  // in the integration points of the face.
  su2double **solL    = workArray;
  su2double **dSolDxL = solL    + nVar;
  su2double **dSolDyL = dSolDxL + nVar;
  su2double **dSolDzL = dSolDyL + nVar;

  su2double **solR    = dSolDzL + nVar;
  su2double **dSolDxR = solR    + nVar;
  su2double **dSolDyR = dSolDxR + nVar;
  su2double **dSolDzR = dSolDyR + nVar;

  // Set the pointer for the total normal flux at the interface. The total
  // flux contains the contributions from the inviscid flux, the viscous flux
  // and the penalty flux. It does not include the symmetrizing fluxes,
  // because these are distributed differently.
  su2double **fluxTot = dSolDzR + nVar;

	// Set the pointers for the left and right symmetrizing fluxes.
	su2double **dFluxSymDxL = fluxTot     + nVar;
	su2double **dFluxSymDyL = dFluxSymDxL + nVar;
	su2double **dFluxSymDzL = dFluxSymDyL + nVar;

	su2double **dFluxSymDxR = dFluxSymDzL + nVar;
	su2double **dFluxSymDyR = dFluxSymDxR + nVar;
	su2double **dFluxSymDzR = dFluxSymDyR + nVar;

  // Set the pointer for the eddy viscosities, which are needed when an
  // LES subgrid scale model is used.
  su2double *eddyVis = dFluxSymDzR[nVar];

  // Compute the right solution and its gradients in the integration points
  // of the face. The right solution corresponds to the kMin face of this
  // element, hence the Legendre basis functions and its derivatives must be
  // taken at the min boundary.
  TensorProductSolAndGradKFace(nInt1D, nVar, nDOFs1D,
                               standardHex->mLegendreInt1D,
                               standardHex->mDerLegendreInt1D,
                               standardHex->mLegendreMinFace1D,
                               standardHex->mDerLegendreMinFace1D,
                               mSol.data(), 
															 solR, dSolDxR, dSolDyR, dSolDzR);

  // Determine what type of face we are dealing with.
  if( elements[i][j][k-1] )
  {
    // This is an internal face. Modify the length scales.
    lenScale    = std::min(lenScale, elements[i][j][k-1]->mLenScaleKDir);
    lenScaleLES = half*(lenScaleLES + elements[i][j][k-1]->mLenScaleLES);

    // Compute the left solution and its gradients in the integration points
    // of the face. The left solution corresponds to the kMax face of the back
    // neighbor, hence the Legendre basis functions and its derivatives must
    // be taken at the max boundary.
    TensorProductSolAndGradKFace(nInt1D, nVar, nDOFs1D,
                                 standardHex->mLegendreInt1D,
                                 standardHex->mDerLegendreInt1D,
                                 standardHex->mLegendreMaxFace1D,
                                 standardHex->mDerLegendreMaxFace1D,
                                 elements[i][j][k-1]->mSol.data(), 
																 solL, dSolDxL, dSolDyL, dSolDzL);

    // Compute the fluxes for an internal face.
    FluxesInternalFace(inputParam, nInt, nIntPad, lenScale, lenScaleLES,
                       standardHex->mIntWeightsFace, 
											 solL, dSolDxL, dSolDyL, dSolDzL, 
											 elements[i][j][k-1]->mSurfMetricIntKMax.data(),
                       solR, dSolDxR, dSolDyR, dSolDzR, 
											 mSurfMetricIntKMin.data(),
                       eddyVis, fluxTot, 
											 dFluxSymDxL, dFluxSymDyL, dFluxSymDzL,
											 dFluxSymDxR, dFluxSymDyR, dFluxSymDzR,
											 ComputeMonitoringData, EddyVisMax);

    // Initialize the residual on the kMin boundary, which is the residual
    // for the bottom neighbor, to zero.
    for(int l=0; l<nVar; ++l)
#pragma omp simd
      for(int i=0; i<nDOFsPad; ++i)
        mResKMin[l][i] = zero;

    // Distribute the normal fluxes to the back element. The face corresponds
    // to the kMax face of the back element, hence the value of the Legendre
    // basis functions at the max face must be used in k-direction.
    TensorProductResKFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreMaxFace1D,
                          fluxTot, mResKMin.data());

    // Distribute the symmetrizing fluxes in i-direction, stored in dFluxSymDxL,
    // to the back element. These fluxes are multiplied by the derivative of
    // the basis functions in i-direction.
    TensorProductResKFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mDerLegendreInt1DTranspose,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreMaxFace1D,
                          dFluxSymDxL, mResKMin.data());

    // Distribute the symmetrizing fluxes in j-direction, stored in dFluxSymDyL,
    // to the back element. These fluxes are multiplied by the derivative of
    // the basis functions in j-direction.
    TensorProductResKFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mDerLegendreInt1DTranspose,
                          standardHex->mLegendreMaxFace1D,
                          dFluxSymDyL, mResKMin.data());

    // Distribute the symmetrizing fluxes in k-direction, stored in dFluxSymDzL,
    // to the back element. These fluxes are multiplied by the derivative of
    // the basis functions in k-direction.
    TensorProductResKFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mDerLegendreMaxFace1D,
                          dFluxSymDzL, mResKMin.data());

    // Negate the normal fluxes, as these will have an opposite sign for
    // the current element.
    for(int l=0; l<nVar; ++l)
#pragma omp simd
      for(int i=0; i<nIntPad; ++i)
        fluxTot[l][i] = -fluxTot[l][i];
  }
  else
  {
    // This is a boundary face. Compute the fluxes for a boundary face. Note, it does not
		// matter if we switch the dFluxSym R and L state, since they are overwritten then the
		// values are returned the same for both states.
    FluxesBoundaryFace(inputParam, standardHex, elements, nInt, nIntPad, mBCKMin, lenScale,
                       lenScaleLES, standardHex->mIntWeightsFace, 
											 solR, dSolDxR, dSolDyR, dSolDzR, 
											-one, mSurfMetricIntKMin.data(), 
											 &mExchangeDataKMin, mPrescribedDataKMin.data(), 
											 solL, dSolDxL, dSolDyL, dSolDzL,
                       eddyVis, fluxTot, 
											 dFluxSymDxR, dFluxSymDyR, dFluxSymDzR,
											 dFluxSymDxL, dFluxSymDyL, dFluxSymDzL,
											 ComputeMonitoringData, EddyVisMax, forceCoef);
  }

  // Distribute the normal fluxes to the current element. The face corresponds
  // to the kMin face of the current element, hence the value of the Legendre
  // basis functions at the min face must be used in k-direction.
  TensorProductResKFace(nInt1D, nVar, nDOFs1D,
                        standardHex->mLegendreInt1DTranspose,
                        standardHex->mLegendreInt1DTranspose,
                        standardHex->mLegendreMinFace1D,
                        fluxTot, mRes.data());

  // Distribute the symmetrizing fluxes in i-direction, stored in dFluxSymDxR,
  // to the current element. These fluxes are multiplied by the derivative
  // of the basis functions in i-direction.
  TensorProductResKFace(nInt1D, nVar, nDOFs1D,
                        standardHex->mDerLegendreInt1DTranspose,
                        standardHex->mLegendreInt1DTranspose,
                        standardHex->mLegendreMinFace1D,
                        dFluxSymDxR, mRes.data());

  // Distribute the symmetrizing fluxes in j-direction, stored in dFluxSymDyR,
  // to the current element. These fluxes are multiplied by the derivative
  // of the basis functions in j-direction.
  TensorProductResKFace(nInt1D, nVar, nDOFs1D,
                        standardHex->mLegendreInt1DTranspose,
                        standardHex->mDerLegendreInt1DTranspose,
                        standardHex->mLegendreMinFace1D,
                        dFluxSymDyR, mRes.data());

  // Distribute the symmetrizing fluxes in k-direction, stored in dFluxSymDzR,
  // to the current element. These fluxes are multiplied by the derivative
  // of the basis functions in k-direction.
  TensorProductResKFace(nInt1D, nVar, nDOFs1D,
                        standardHex->mLegendreInt1DTranspose,
                        standardHex->mLegendreInt1DTranspose,
                        standardHex->mDerLegendreMinFace1D,
                        dFluxSymDzR, mRes.data());

  // Check if the kMax face is a boundary face.
  if( !elements[i][j][k+1] )
  {
    // Reset the length scales.
    lenScale    = mLenScaleKDir;
    lenScaleLES = mLenScaleLES;

    // Compute the left solution and its gradients in the integration points
    // of the face. The left solution corresponds to the kMax face of the current
    // element, hence the Legendre basis functions and its derivatives must be
    // taken at the max boundary.
    TensorProductSolAndGradKFace(nInt1D, nVar, nDOFs1D,
                                 standardHex->mLegendreInt1D,
                                 standardHex->mDerLegendreInt1D,
                                 standardHex->mLegendreMaxFace1D,
                                 standardHex->mDerLegendreMaxFace1D,
                                 mSol.data(), 
																 solL, dSolDxL, dSolDyL, dSolDzL);

    // Compute the fluxes for a boundary face.
    FluxesBoundaryFace(inputParam, standardHex, elements, nInt, nIntPad, mBCKMax, lenScale,
                       lenScaleLES, standardHex->mIntWeightsFace, 
											 solL, dSolDxL, dSolDyL, dSolDzL, 
											 one, mSurfMetricIntKMax.data(), 
											 &mExchangeDataKMax, mPrescribedDataKMax.data(), 
											 solR, dSolDxR, dSolDyR, dSolDzR,
                       eddyVis, fluxTot, 
											 dFluxSymDxL, dFluxSymDyL, dFluxSymDzL,
											 dFluxSymDxR, dFluxSymDyR, dFluxSymDzR,
											 ComputeMonitoringData, EddyVisMax, forceCoef);

    // Distribute the normal fluxes to the current element. The face corresponds
    // to the kMax face of the current element, hence the value of the Legendre
    // basis functions at the max face must be used in k-direction.
    TensorProductResKFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreMaxFace1D,
                          fluxTot, mRes.data());

    // Distribute the symmetrizing fluxes in i-direction, stored in dFluxSymDxL,
    // to the current element. These fluxes are multiplied by the derivative of
    // the basis functions in i-direction.
    TensorProductResKFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mDerLegendreInt1DTranspose,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreMaxFace1D,
                          dFluxSymDxL, mRes.data());

    // Distribute the symmetrizing fluxes in j-direction, stored in dFluxSymDyL,
    // to the current element. These fluxes are multiplied by the derivative of
    // the basis functions in j-direction.
    TensorProductResKFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mDerLegendreInt1DTranspose,
                          standardHex->mLegendreMaxFace1D,
                          dFluxSymDyL, mRes.data());

    // Distribute the symmetrizing fluxes in k-direction, stored in dFluxSymDzL,
    // to the current element. These fluxes are multiplied by the derivative of
    // the basis functions in k-direction.
    TensorProductResKFace(nInt1D, nVar, nDOFs1D,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mLegendreInt1DTranspose,
                          standardHex->mDerLegendreMaxFace1D,
                          dFluxSymDzL, mRes.data());
  }
}

//------------------------------------------------------------------------------

// Function, which determines whether or not an LGL distribution of the
// nodal grid DOFs is used.
bool ElementClass::LGLDistribution(const InputParamClass      *inputParam,
                                   const StandardElementClass *standardHex)
{
  // Define the tolerance in the distribution.
  const su2double tolDist = (su2double) 0.05;

  // Easier storage of the polynomial degree of the grid and allocate the
  // memory for the vector to store the length distribution of an edge.
  const int nPolyGrid = inputParam->mNPolyGridDOFs;
  std::vector<su2double> lenDistr(nPolyGrid+1);

  // Determine the point distribution for a standard edge with an equidistant
  // distribution and with an LGL distribution.
  std::vector<su2double> rDOFsEqui, rDOFsLGL;
  LocationDOFs1D(nPolyGrid+1, EQUIDISTANT, rDOFsEqui);
  LocationDOFs1D(nPolyGrid+1, LGL_POINTS,  rDOFsLGL);

  // Initialize the number of edges with an equidistant distribution and an LGL
  // distribution to zero. Als the number of bad edges are initialized to zero.
  int nEdgesEquiDistant = 0, nEdgesLGL = 0, nEdgesBad = 0;

  //-----------------------
  // Edges in i-direction.
  //-----------------------

  // Loop over the k- and j-direction.
  for(int k=0; k<=nPolyGrid; ++k)
  {
    for(int j=0; j<=nPolyGrid; ++j)
    {
      // Determine the 1D index of the starting node of the edge.
      const int indS = k*(nPolyGrid+1)*(nPolyGrid+1) + j*(nPolyGrid+1);

      // Set the length distribution of the first point to zero.
      lenDistr[0] = zero;

      // Loop over the remaining points to determine the distribution.
      for(int i=1; i<=nPolyGrid; ++i)
      {
        const int ind1 = indS + i;
        const int ind0 = ind1 - 1;
        const su2double dx = mCoorNodalGridDOFs[0][ind1] - mCoorNodalGridDOFs[0][ind0];
        const su2double dy = mCoorNodalGridDOFs[1][ind1] - mCoorNodalGridDOFs[1][ind0];
        const su2double dz = mCoorNodalGridDOFs[2][ind1] - mCoorNodalGridDOFs[2][ind0];

        lenDistr[i] = lenDistr[i-1] + SQRT(dx*dx + dy*dy + dz*dz);
      }

      // Create a scaled version, such that it corresponds to a standard
      // edge (-1..1).
      const su2double twoLenInv = two/lenDistr[nPolyGrid];
      for(int i=0; i<=nPolyGrid; ++i)
        lenDistr[i] = -one + twoLenInv*lenDistr[i];

      // Determine the maximum deviation from the equidistant and LGL
      // distribution. Only necessary to consider the internal points.
      su2double maxDevEqui = zero, maxDevLGL = zero;
      for(int i=1; i<nPolyGrid; ++i)
      {
        maxDevEqui = std::max(maxDevEqui, FABS(lenDistr[i]-rDOFsEqui[i]));
        maxDevLGL  = std::max(maxDevLGL,  FABS(lenDistr[i]-rDOFsLGL[i]));
      }

      // Check the situation and either update the number of equidistant edges
      // or the number of LGL edges. If the deviation from the particular
      // distribution is too high the number of bad edges is updated.
      if(maxDevEqui <= maxDevLGL)
      {
        ++nEdgesEquiDistant;
        if(maxDevEqui > tolDist) ++nEdgesBad;
      }
      else
      {
        ++nEdgesLGL;
        if(maxDevLGL > tolDist) ++nEdgesBad;
      }
    }
  }

  //-----------------------
  // Edges in j-direction.
  //-----------------------

  // Loop over the k- and i-direction.
  for(int k=0; k<=nPolyGrid; ++k)
  {
    for(int i=0; i<=nPolyGrid; ++i)
    {
      // Determine the 1D index of the starting node of the edge.
      const int indS = k*(nPolyGrid+1)*(nPolyGrid+1) + i;

      // Set the length distribution of the first point to zero.
      lenDistr[0] = zero;

      // Loop over the remaining points to determine the distribution.
      for(int j=1; j<=nPolyGrid; ++j)
      {
        const int ind1 = indS + j*(nPolyGrid+1);
        const int ind0 = ind1 - (nPolyGrid+1);
        const su2double dx = mCoorNodalGridDOFs[0][ind1] - mCoorNodalGridDOFs[0][ind0];
        const su2double dy = mCoorNodalGridDOFs[1][ind1] - mCoorNodalGridDOFs[1][ind0];
        const su2double dz = mCoorNodalGridDOFs[2][ind1] - mCoorNodalGridDOFs[2][ind0];

        lenDistr[j] = lenDistr[j-1] + SQRT(dx*dx + dy*dy + dz*dz);
      }

      // Create a scaled version, such that it corresponds to a standard
      // edge (-1..1).
      const su2double twoLenInv = two/lenDistr[nPolyGrid];
      for(int j=0; j<=nPolyGrid; ++j)
        lenDistr[j] = -one + twoLenInv*lenDistr[j];

      // Determine the maximum deviation from the equidistant and LGL
      // distribution. Only necessary to consider the internal points.
      su2double maxDevEqui = zero, maxDevLGL = zero;
      for(int j=1; j<nPolyGrid; ++j)
      {
        maxDevEqui = std::max(maxDevEqui, FABS(lenDistr[j]-rDOFsEqui[j]));
        maxDevLGL  = std::max(maxDevLGL,  FABS(lenDistr[j]-rDOFsLGL[j]));
      }

      // Check the situation and either update the number of equidistant edges
      // or the number of LGL edges. If the deviation from the particular
      // distribution is too high the number of bad edges is updated.
      if(maxDevEqui <= maxDevLGL)
      {
        ++nEdgesEquiDistant;
        if(maxDevEqui > tolDist) ++nEdgesBad;
      }
      else
      {
        ++nEdgesLGL;
        if(maxDevLGL > tolDist) ++nEdgesBad;
      }
    }
  }

  //-----------------------
  // Edges in k-direction.
  //-----------------------

  // Loop over the j- and i-direction.
  for(int j=0; j<=nPolyGrid; ++j)
  {
    for(int i=0; i<=nPolyGrid; ++i)
    {
      // Determine the 1D index of the starting node of the edge.
      const int indS = j*(nPolyGrid+1) + i;

      // Set the length distribution of the first point to zero.
      lenDistr[0] = zero;

      // Loop over the remaining points to determine the distribution.
      for(int k=1; k<=nPolyGrid; ++k)
      {
        const int ind1 = indS + k*(nPolyGrid+1)*(nPolyGrid+1);
        const int ind0 = ind1 - (nPolyGrid+1)*(nPolyGrid+1);
        const su2double dx = mCoorNodalGridDOFs[0][ind1] - mCoorNodalGridDOFs[0][ind0];
        const su2double dy = mCoorNodalGridDOFs[1][ind1] - mCoorNodalGridDOFs[1][ind0];
        const su2double dz = mCoorNodalGridDOFs[2][ind1] - mCoorNodalGridDOFs[2][ind0];

        lenDistr[k] = lenDistr[k-1] + SQRT(dx*dx + dy*dy + dz*dz);
      }

      // Create a scaled version, such that it corresponds to a standard
      // edge (-1..1).
      const su2double twoLenInv = two/lenDistr[nPolyGrid];
      for(int k=0; k<=nPolyGrid; ++k)
        lenDistr[k] = -one + twoLenInv*lenDistr[k];

      // Determine the maximum deviation from the equidistant and LGL
      // distribution. Only necessary to consider the internal points.
      su2double maxDevEqui = zero, maxDevLGL = zero;
      for(int k=1; k<nPolyGrid; ++k)
      {
        maxDevEqui = std::max(maxDevEqui, FABS(lenDistr[k]-rDOFsEqui[k]));
        maxDevLGL  = std::max(maxDevLGL,  FABS(lenDistr[k]-rDOFsLGL[k]));
      }

      // Check the situation and either update the number of equidistant edges
      // or the number of LGL edges. If the deviation from the particular
      // distribution is too high the number of bad edges is updated.
      if(maxDevEqui <= maxDevLGL)
      {
        ++nEdgesEquiDistant;
        if(maxDevEqui > tolDist) ++nEdgesBad;
      }
      else
      {
        ++nEdgesLGL;
        if(maxDevLGL > tolDist) ++nEdgesBad;
      }
    }
  }

  // Check if the distribution is consistent.
  if(nEdgesEquiDistant > nEdgesLGL)
  {
    if((nEdgesEquiDistant < 10*nEdgesLGL) ||
       (nEdgesEquiDistant < 10*nEdgesBad))
    {
      Terminate("ElementClass::LGLDistribution", __FILE__, __LINE__,
                "Not able to determine the point distribution");
    }
  }
  else
  {
    if((nEdgesLGL < 10*nEdgesEquiDistant) ||
       (nEdgesLGL < 10*nEdgesBad))
    {
      Terminate("ElementClass::LGLDistribution", __FILE__, __LINE__,
                "Not able to determine the point distribution");
    }
  }

  // Return the correct boolean.
  return nEdgesLGL > nEdgesEquiDistant;
}

//------------------------------------------------------------------------------

// Function, which multiplies the residual with the inverse of
// the mass matrix.
void ElementClass::MultiplyResInverseMassMatrix(const InputParamClass      *inputParam,
                                                const StandardElementClass *standardHex,
                                                su2double                  **workArray)
{
  // Make a distinction between the working variables and call the appropriate
  // function to carry out the actual work.
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      CG_ConservativeVarMassMatrix(standardHex, mSol.data(), mVolMetricInt[0],
                                   mVolElem, mRes.data(), workArray);
      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      CG_EntropyVarMassMatrix(standardHex, mSol.data(), mVolMetricInt[0],
                              mVolElem, mDUdVInt.data(), mRes.data(), workArray);
      break;
    }

    //--------------------------------------------------------------------------

    default:
    {
      // This is just to avoid a compiler warning.
      break;
    }
  }
}

//------------------------------------------------------------------------------

// Function, which determines the prescribed data in the integration
// points of the boundary faces.
void ElementClass::PrescribedDataIntegrationPoints(
                                 const InputParamClass      *inputParam,
                                 const StandardElementClass *standardHex,
                                 const su2double            *nonDimPrimVarFreeStream)
{
  // Call the function PrescribedDataIntegrationPointsFace, with the appropriate
  // arguments, to do the actual work.
  PrescribedDataIntegrationPointsFace(mSurfMetricIntIMin, mBCIMin, inputParam, standardHex,
                                      nonDimPrimVarFreeStream, mPrescribedDataIMin);
  PrescribedDataIntegrationPointsFace(mSurfMetricIntIMax, mBCIMax, inputParam, standardHex,
                                      nonDimPrimVarFreeStream, mPrescribedDataIMax);
  PrescribedDataIntegrationPointsFace(mSurfMetricIntJMin, mBCJMin, inputParam, standardHex,
                                      nonDimPrimVarFreeStream, mPrescribedDataJMin);
  PrescribedDataIntegrationPointsFace(mSurfMetricIntIMax, mBCJMax, inputParam, standardHex,
                                      nonDimPrimVarFreeStream, mPrescribedDataJMax);
  PrescribedDataIntegrationPointsFace(mSurfMetricIntKMin, mBCKMin, inputParam, standardHex,
                                      nonDimPrimVarFreeStream, mPrescribedDataKMin);
  PrescribedDataIntegrationPointsFace(mSurfMetricIntKMax, mBCKMax, inputParam, standardHex,
                                      nonDimPrimVarFreeStream, mPrescribedDataKMax);
}

//------------------------------------------------------------------------------

// Function, which determines the prescribed data in the integration
// points of the given boundary face.
void ElementClass::PrescribedDataIntegrationPointsFace(
                             const std::vector<su2double *> &faceMetric,
                             SubfaceBaseClass               *BC,
                             const InputParamClass          *inputParam,
                             const StandardElementClass     *standardHex,
                             const su2double                *nonDimPrimVarFreeStream,
                             std::vector<su2double *>       &prescribedData)
{
  // Return immediately if this is not a physical boundary.
  if( !BC ) return;

  // Determine the number of prescribed variables for this boundary.
  // Return if no variables are prescribed.
  const int nVarPrescribed = BC->GetNVarPrescribed();
  if(nVarPrescribed == 0) return;

  // Allocate the memory for prescribedData.
  const int nInt2D = standardHex->mNIntegration2DPad;
  prescribedData.resize(nVarPrescribed);
  for(int l=0; l<nVarPrescribed; ++l)
  {
    prescribedData[l] = (su2double *) AllocateMemory(nInt2D*sizeof(su2double));
    if( !prescribedData[l] )
      Terminate("ElementClass::PrescribedDataIntegrationPointsFace", __FILE__,
                __LINE__, "Memory allocation failure for prescribedData");
  }

  // Determine the, possibly interpolated, prescribed data
  // in the integration points.
  BC->PrescribedDataIntegrationPoints(nInt2D, faceMetric[13], faceMetric[14],
                                      faceMetric[15], nonDimPrimVarFreeStream,
                                      inputParam->mFEMVariables, prescribedData);
}

//------------------------------------------------------------------------------

// Function, which stores the local indices of the element.
void ElementClass::StoreLocalIndices(const int i,
                                     const int j,
                                     const int k)
{
  mLocalInd[0] = i;
  mLocalInd[1] = j;
  mLocalInd[2] = k;
}

//------------------------------------------------------------------------------

// Function, which updates the data for the averages.
void ElementClass::UpdateAverageData(const InputParamClass      *inputParam,
                                     const StandardElementClass *standardHex,
                                     const int                  nTimeSteps,
                                     su2double                  **primVar,
                                     su2double                  **gradVel)
{
  // Store the number of DOFs and padded DOFs per element a bit easier.
  const int nDOFs    = standardHex->mNDOFs;
  const int nDOFsPad = standardHex->mNDOFsPad;

  // Determine the factors for the averaging.
  const su2double factNew = one/nTimeSteps;
  const su2double factOld = one - factNew;

  // Compute the primitive variables in the nodal DOFs.
  ComputePrimitiveVariablesNodalDOFs(inputParam, standardHex, primVar);

  // Loop over the padded number of DOFs.
#pragma omp simd
  for(int l=0; l<nDOFsPad; ++l)
  {
    // Easier storage of the primitive variables.
    const su2double rho = primVar[0][l];
    const su2double u   = primVar[1][l];
    const su2double v   = primVar[2][l];
    const su2double w   = primVar[3][l];
    const su2double p   = primVar[4][l];

    // Update the average primitive variables.
    mAvePrim[0][l] = factOld*mAvePrim[0][l] + factNew*rho;
    mAvePrim[1][l] = factOld*mAvePrim[1][l] + factNew*u;
    mAvePrim[2][l] = factOld*mAvePrim[2][l] + factNew*v;
    mAvePrim[3][l] = factOld*mAvePrim[3][l] + factNew*w;
    mAvePrim[4][l] = factOld*mAvePrim[4][l] + factNew*p;

    // Update the average of the velocity products.
    mAveVelProd[0][l] = factOld*mAveVelProd[0][l] + factNew*u*u;
    mAveVelProd[1][l] = factOld*mAveVelProd[1][l] + factNew*v*v;
    mAveVelProd[2][l] = factOld*mAveVelProd[2][l] + factNew*w*w;
    mAveVelProd[3][l] = factOld*mAveVelProd[3][l] + factNew*u*v;
    mAveVelProd[4][l] = factOld*mAveVelProd[4][l] + factNew*u*w;
    mAveVelProd[5][l] = factOld*mAveVelProd[5][l] + factNew*v*w;
  }

// Check if a subgrid scale model is used.
  if(inputParam->mSGSModelType != NO_SGS_MODEL)
  {
    // Set the pointer for the eddy viscosity. Note that it is assumed that the first
    // index of primVar is allocated with dimension nVar+1.
    su2double *eddyVis = primVar[nVar];

    // Compute the velocity gradients in the nodal DOFs.
    ComputeVelocityGradientsNodalDOFs(inputParam, standardHex, primVar, gradVel);

    // Compute the eddy viscosities in the DOFs.
    inputParam->mSGSModel->EddyViscosity(nDOFs, mLenScaleLES, primVar,
                                         gradVel, eddyVis);

    // Update the average of the eddy viscosity.
#pragma omp simd
    for(int l=0; l<nDOFs; ++l)
      mAveEddyVis[l] = factOld*mAveEddyVis[l] + factNew*eddyVis[l];
  }
}  

//------------------------------------------------------------------------------

// Function, which computes the volume residuals of this element.
void ElementClass::VolumeResidual(const InputParamClass      *inputParam,
                                  const StandardElementClass *standardHex,
                                  su2double                  **workArray,
                                  const bool                 ComputeMonitoringData,
                                  su2double                  &Mach2Max,
                                  su2double                  &EddyVisMax)
{
  // Easier storage of the number of DOFs and integration points.
  const int nDOFs1D  = standardHex->mNDOFs1D;
  const int nDOFsPad = standardHex->mNDOFsPad;
  const int nInt1D   = standardHex->mNIntegration1D;
  const int nInt     = standardHex->mNIntegration;
  const int nIntPad  = standardHex->mNIntegrationPad;

  //--------------------------------------------------------------------------
  // Step 1: Compute the solution and the gradients in the integration points.
  //--------------------------------------------------------------------------

  // Set the pointers for the solution and its gradients.
  su2double **Var    = workArray;
  su2double **dVarDx = Var    + nVar;
  su2double **dVarDy = dVarDx + nVar;
  su2double **dVarDz = dVarDy + nVar;

  // Set the pointer for the eddy viscosity.
  su2double *eddyVis = dVarDz[nVar];

  // Call TensorProductSolAndGradVolume to carry out the actual tensor product
  // multiplication to compute both the solution and the gradients in the
  // integration points. Note that the gradients contain the derivatives w.r.t.
  // the parametric coordinates and not yet w.r.t. Cartesian coordinates.
  TensorProductSolAndGradVolume(nInt1D, nVar, nDOFs1D, standardHex->mLegendreInt1D,
                                standardHex->mDerLegendreInt1D, mSol.data(),
                                Var, dVarDx, dVarDy, dVarDz);

  // Initialize the padded values to avoid problems.
  for(int l=0; l<nVar; ++l)
  {
    for(int m=nInt; m<nIntPad; ++m)
    {
      Var[l][m]    = Var[l][0];
      dVarDx[l][m] = dVarDy[l][m] = dVarDz[l][m] = zero;
    }
  }

  // Convert the gradients to Cartesian values.
  for(int l=0; l<nVar; ++l)
  {
#pragma omp simd
    for(int m=0; m<nIntPad; ++m)
    {
      const su2double dvdr = dVarDx[l][m];
      const su2double dvds = dVarDy[l][m];
      const su2double dvdt = dVarDz[l][m];

      dVarDx[l][m] = dvdr*mVolMetricInt[1][l]
                   + dvds*mVolMetricInt[4][l]
                   + dvdt*mVolMetricInt[7][l];
      dVarDy[l][m] = dvdr*mVolMetricInt[2][l]
                   + dvds*mVolMetricInt[5][l]
                   + dvdt*mVolMetricInt[8][l];
      dVarDz[l][m] = dvdr*mVolMetricInt[3][l]
                   + dvds*mVolMetricInt[6][l]
                   + dvdt*mVolMetricInt[9][l];
    }
  }

  //--------------------------------------------------------------------------
  // Step 2: Convert the working variables and its gradients to primitive
  //         variables and gradients of velocity and internal energy.
  //--------------------------------------------------------------------------

  // Make a distinction between the working variables.
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      // The conservative variables are used as working variables.
      // Compute gam-1, as this appears in the expression for the pressure.
      const su2double gm1 = GamConstant - one;

      // Loop over the (padded) number of integration points.
#pragma omp simd
      for(int l=0; l<nIntPad; ++l)
      {
        // Compute the primitive variables from the conservative ones.
        const su2double rho    = Var[0][l];
        const su2double rhoInv = one/rho;
        const su2double u      = rhoInv*Var[1][l];
        const su2double v      = rhoInv*Var[2][l];
        const su2double w      = rhoInv*Var[3][l];
        const su2double p      = gm1*(Var[4][l]
                               -      half*(u*Var[1][l] + v*Var[2][l] + w*Var[3][l]));

        // Compute the velocity gradients from the gradients
        // of the conservative variables. Use the locations of the
        // momentum gradients to store these values.
        dVarDx[1][l] = rhoInv*(dVarDx[1][l] - u*dVarDx[0][l]);
        dVarDy[1][l] = rhoInv*(dVarDy[1][l] - u*dVarDy[0][l]);
        dVarDz[1][l] = rhoInv*(dVarDz[1][l] - u*dVarDz[0][l]);

        dVarDx[2][l] = rhoInv*(dVarDx[2][l] - v*dVarDx[0][l]);
        dVarDy[2][l] = rhoInv*(dVarDy[2][l] - v*dVarDy[0][l]);
        dVarDz[2][l] = rhoInv*(dVarDz[2][l] - v*dVarDz[0][l]);

        dVarDx[3][l] = rhoInv*(dVarDx[3][l] - w*dVarDx[0][l]);
        dVarDy[3][l] = rhoInv*(dVarDy[3][l] - w*dVarDy[0][l]);
        dVarDz[3][l] = rhoInv*(dVarDz[3][l] - w*dVarDz[0][l]);

        // Compute the gradients of the internal energy. Use the locations
        // of the total energy gradients to store them.
        const su2double ETot = rhoInv*Var[4][l];

        dVarDx[4][l] = rhoInv*(dVarDx[4][l] - ETot*dVarDx[0][l])
                     - u*dVarDx[1][l] - v*dVarDx[2][l] - w*dVarDx[3][l];
        dVarDy[4][l] = rhoInv*(dVarDy[4][l] - ETot*dVarDy[0][l])
                     - u*dVarDy[1][l] - v*dVarDy[2][l] - w*dVarDy[3][l];
        dVarDz[4][l] = rhoInv*(dVarDz[4][l] - ETot*dVarDz[0][l])
                     - u*dVarDz[1][l] - v*dVarDz[2][l] - w*dVarDz[3][l];

        // Store the velocities in Var1, Var2 and Var and use
        // dVar0Dx to store the pressure.
        Var[1][l] = u; Var[2][l] = v; Var[3][l] = w; dVarDx[0][l] = p;
      }

      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      // The entropy variables are used as working variables.
      // Easier storage of some expressions involving gamma.
      const su2double gm1   =  GamConstant - one;
      const su2double ov1mg = -one/gm1;

      // Loop over the (padded) number of integration points.
#pragma omp simd
      for(int l=0; l<nIntPad; ++l)
      {
        // Compute the primitive variables from the entropy ones.
        const su2double V4Inv =  one/Var[4][l];
        const su2double u     = -V4Inv*Var[1][l];
        const su2double v     = -V4Inv*Var[2][l];
        const su2double w     = -V4Inv*Var[3][l];
        const su2double eKin  =  half*(u*u + v*v + w*w);
        const su2double s     =  GamConstant - gm1*(Var[0][l] - Var[4][l]*eKin);
        const su2double tmp   = -Var[4][l]*EXP(s);
        const su2double rho   =  POW(tmp, ov1mg);
        const su2double p     = -rho*V4Inv;
        const su2double rE    = -p*ov1mg + rho*eKin;

        // Compute the velocity gradients, which are stored directly again
        // in dVar1Dx, dVar2Dx, etc.
        dVarDx[1][l] = -V4Inv*(dVarDx[1][l] + u*dVarDx[4][l]);
        dVarDy[1][l] = -V4Inv*(dVarDy[1][l] + u*dVarDy[4][l]);
        dVarDz[1][l] = -V4Inv*(dVarDz[1][l] + u*dVarDz[4][l]);

        dVarDx[2][l] = -V4Inv*(dVarDx[2][l] + v*dVarDx[4][l]);
        dVarDy[2][l] = -V4Inv*(dVarDy[2][l] + v*dVarDy[4][l]);
        dVarDz[2][l] = -V4Inv*(dVarDz[2][l] + v*dVarDz[4][l]);

        dVarDx[3][l] = -V4Inv*(dVarDx[3][l] + w*dVarDx[4][l]);
        dVarDy[3][l] = -V4Inv*(dVarDy[3][l] + w*dVarDy[4][l]);
        dVarDz[3][l] = -V4Inv*(dVarDz[3][l] + w*dVarDz[4][l]);

        // Compute the gradients of the internal energy.
        // Store them in dVar4Dx, etc.
        const su2double fact = -V4Inv*V4Inv*ov1mg;

        dVarDx[4][l] *= fact; dVarDy[4][l] *= fact; dVarDz[4][l] *= fact;

        // Store the density, velocities and total energy in Var0, etc.
        // and use dVar0Dx to store the pressure.
        Var[0][l] = rho; Var[1][l] = u; Var[2][l] = v; Var[3][l] = w; Var[4][l] = rE;
        dVarDx[0][l] = p;

        // Store the elements of the transformation matrix dUdV in this
        // integration point. Note that this matrix is symmetric and hence
        // only the upper-diagonal part (or lower diagonal part) is stored.
        const su2double rH = rE + p;
        const su2double ru = rho*u;
        const su2double rv = rho*v;
        const su2double rw = rho*w;

        mDUdVInt[0][l]  = rho;                // dUdV(0,0)
        mDUdVInt[1][l]  = ru;                 // dUdV(0,1) = dUdV(1,0)
        mDUdVInt[2][l]  = rv;                 // dUdV(0,2) = dUdV(2,0)
        mDUdVInt[3][l]  = rw;                 // dUdV(0,3) = dUdV(3,0)
        mDUdVInt[4][l]  = rE;                 // dUdV(0,4) = dUdV(4,0)
        mDUdVInt[5][l]  = ru*u + p;           // dUdV(1,1)
        mDUdVInt[6][l]  = ru*v;               // dUdV(1,2) = dUdV(2,1)
        mDUdVInt[7][l]  = ru*w;               // dUdV(1,3) = dUdV(3,1)
        mDUdVInt[8][l]  = rH*u;               // dUdV(1,4) = dUdV(4,1)
        mDUdVInt[9][l]  = rv*v + p;           // dUdV(2,2)
        mDUdVInt[10][l] = rv*w;               // dUdV(2,3) = dUdV(3,2)
        mDUdVInt[11][l] = rH*v;               // dUdV(2,4) = dUdV(4,2)
        mDUdVInt[12][l] = rw*w + p;           // dUdV(3,3)
        mDUdVInt[13][l] = rH*w;               // dUdV(3,4) = dUdV(4,3)
        mDUdVInt[14][l] = rE*rH/rho + p*eKin; // dUdV(4,4)
      }

      break;
    }

    //--------------------------------------------------------------------------

    default:
    {
      // This is just to avoid a compiler warning.
      break;
    }
  }

  //--------------------------------------------------------------------------
  // Step 3: Compute the eddy viscosity, if needed.
  //--------------------------------------------------------------------------

  if(inputParam->mSGSModelType != NO_SGS_MODEL)
  {
    // Set the pointers for the velocity gradients.
    su2double *velGrad[9];
    velGrad[0] = dVarDx[1]; velGrad[1] = dVarDx[2]; velGrad[2] = dVarDx[3];
    velGrad[3] = dVarDy[1]; velGrad[4] = dVarDy[2]; velGrad[5] = dVarDy[3];
    velGrad[6] = dVarDz[1]; velGrad[7] = dVarDz[2]; velGrad[8] = dVarDz[3];

    // Compute the eddy viscosities.
    inputParam->mSGSModel->EddyViscosity(nIntPad, mLenScaleLES, Var, velGrad, eddyVis);
  }
  else
  {
    // No SGS model is used. Initialize the eddy viscosity to zero.
#pragma omp simd
    for(int l=0; l<nIntPad; ++l)
      eddyVis[l] = zero;
  }

  //--------------------------------------------------------------------------
  // Step 4: Update the monitoring data, if needed.
  //--------------------------------------------------------------------------

  if( ComputeMonitoringData )
  {
    // Loop over the integration points. Note that the padded number can
    // be taken here, because the dummy values are initialized properly.
    for(int l=0; l<nIntPad; ++l)
    {
      // More readable storage of the variables.
      const su2double rho = Var[0][l], u = Var[1][l], v = Var[2][l], w = Var[3][l];
      const su2double p   = dVarDx[0][l];

      // Compute the square of the local Mach number and update Mach2Max.
      const su2double Mach2 = rho*(u*u + v*v + w*w)/(GamConstant*p);
      Mach2Max = std::max(Mach2Max, Mach2);

      // Update the maximum value of the eddy viscosity.
      EddyVisMax = std::max(EddyVisMax, eddyVis[l]);
    }
  }

  //--------------------------------------------------------------------------
  // Step 5: Compute the sum of the inviscid and viscous fluxes in the three
  //         index directions in the integration points.
  //--------------------------------------------------------------------------

  // Loop over the (padded) number of integration points.
#pragma omp simd
  for(int l=0; l<nIntPad; ++l)
  {
    // Determine the Jacobian times the integration weight for this integration
    // point. The minus sign is there, because in the weak formulation the
    // volume term appears with a minus sign (due to the integration by parts).
    const su2double weight = -mVolMetricInt[0][l]*standardHex->mIntWeights[l];

    // More readable storage of the variables.
    const su2double rho = Var[0][l], u = Var[1][l], v = Var[2][l], w = Var[3][l], rE = Var[4][l];
    const su2double p   = dVarDx[0][l];

    const su2double dudx = dVarDx[1][l], dudy = dVarDy[1][l], dudz = dVarDz[1][l];
    const su2double dvdx = dVarDx[2][l], dvdy = dVarDy[2][l], dvdz = dVarDz[2][l];
    const su2double dwdx = dVarDx[3][l], dwdy = dVarDy[3][l], dwdz = dVarDz[3][l];
    const su2double dedx = dVarDx[4][l], dedy = dVarDy[4][l], dedz = dVarDz[4][l];

    // Compute the momentum variables.
    const su2double ru = rho*u, rv = rho*v, rw = rho*w;

    // Compute the total viscosity and total factor for the heat flux.
    const su2double muTot           = mu + eddyVis[l];
    const su2double factHeatFluxTot = mu*factHeatFlux_Lam + eddyVis[l]*factHeatFlux_Turb;

    // Compute the components of the viscous stress tensor.
    const su2double divVelTerm = lambdaOverMu*(dudx + dvdy + dwdz);

    const su2double tauxx = muTot*(two*dudx + divVelTerm);
    const su2double tauyy = muTot*(two*dvdy + divVelTerm);
    const su2double tauzz = muTot*(two*dwdz + divVelTerm);

    const su2double tauxy = muTot*(dudy + dvdx);
    const su2double tauxz = muTot*(dudz + dwdx);
    const su2double tauyz = muTot*(dvdz + dwdy);

    // Compute the components of the heat flux vector.
    const su2double qx = -factHeatFluxTot*dedx;
    const su2double qy = -factHeatFluxTot*dedy;
    const su2double qz = -factHeatFluxTot*dedz;

    // Compute the sum of the inviscid and viscous fluxes in x-direction.
    const su2double fx0 = ru;
    const su2double fx1 = ru*u - tauxx + p;
    const su2double fx2 = rv*u - tauxy;
    const su2double fx3 = rw*u - tauxz;
    const su2double fx4 = (rE+p)*u - (u*tauxx + v*tauxy + w*tauxz - qx);

    // Compute the sum of the inviscid and viscous fluxes in y-direction.
    const su2double fy0 = rv;
    const su2double fy1 = ru*v - tauxy;
    const su2double fy2 = rv*v - tauyy + p;
    const su2double fy3 = rw*v - tauyz;
    const su2double fy4 = (rE+p)*v - (u*tauxy + v*tauyy + w*tauyz - qy);

    // Compute the sum of the inviscid and viscous fluxes in z-direction.
    const su2double fz0 = rw;
    const su2double fz1 = ru*w - tauxz;
    const su2double fz2 = rv*w - tauyz;
    const su2double fz3 = rw*w - tauzz + p;
    const su2double fz4 = (rE+p)*w - (u*tauxz + v*tauyz + w*tauzz - qz);

    // Compute the fluxes in i-direction, multiplied by the integration weight
    // and Jacobian. Fluxes in i-direction are stored in dVarDx.
    const su2double drdx = weight*mVolMetricInt[1][l];
    const su2double drdy = weight*mVolMetricInt[2][l];
    const su2double drdz = weight*mVolMetricInt[3][l];

    dVarDx[0][l] = fx0*drdx + fy0*drdy + fz0*drdz;
    dVarDx[1][l] = fx1*drdx + fy1*drdy + fz1*drdz;
    dVarDx[2][l] = fx2*drdx + fy2*drdy + fz2*drdz;
    dVarDx[3][l] = fx3*drdx + fy3*drdy + fz3*drdz;
    dVarDx[4][l] = fx4*drdx + fy4*drdy + fz4*drdz;

    // Compute the fluxes in j-direction, multiplied by the integration weight
    // and Jacobian. Fluxes in j-direction are stored in dVarDy.
    const su2double dsdx = weight*mVolMetricInt[4][l];
    const su2double dsdy = weight*mVolMetricInt[5][l];
    const su2double dsdz = weight*mVolMetricInt[6][l];

    dVarDy[0][l] = fx0*dsdx + fy0*dsdy + fz0*dsdz;
    dVarDy[1][l] = fx1*dsdx + fy1*dsdy + fz1*dsdz;
    dVarDy[2][l] = fx2*dsdx + fy2*dsdy + fz2*dsdz;
    dVarDy[3][l] = fx3*dsdx + fy3*dsdy + fz3*dsdz;
    dVarDy[4][l] = fx4*dsdx + fy4*dsdy + fz4*dsdz;

    // Compute the fluxes in k-direction, multiplied by the integration weight
    // and Jacobian. Fluxes in k-direction are stored in dVarDz.
    const su2double dtdx = weight*mVolMetricInt[7][l];
    const su2double dtdy = weight*mVolMetricInt[8][l];
    const su2double dtdz = weight*mVolMetricInt[9][l];

    dVarDz[0][l] = fx0*dtdx + fy0*dtdy + fz0*dtdz;
    dVarDz[1][l] = fx1*dtdx + fy1*dtdy + fz1*dtdz;
    dVarDz[2][l] = fx2*dtdx + fy2*dtdy + fz2*dtdz;
    dVarDz[3][l] = fx3*dtdx + fy3*dtdy + fz3*dtdz;
    dVarDz[4][l] = fx4*dtdx + fy4*dtdy + fz4*dtdz;
  }

  //--------------------------------------------------------------------------
  // Step 6: Compute the contribution of the volume integral to the residual.
  //         This also serves as an initialization of the residual.
  //--------------------------------------------------------------------------

  // Initialize the residual to zero.
  for(int l=0; l<nVar; ++l)
#pragma omp simd
    for(int m=0; m<nDOFsPad; ++m)
      mRes[l][m] = zero;

  // Call the function TensorProductVolumeResidual three times to accumulate
  // the residual contributions from the fluxes in the three directions.
  // First the i-direction.
  TensorProductVolumeResidual(nInt1D, nVar, nDOFs1D,
                              standardHex->mDerLegendreInt1DTranspose,
                              standardHex->mLegendreInt1DTranspose,
                              standardHex->mLegendreInt1DTranspose,
                              dVarDx, mRes.data());

  // The j-direction.
  TensorProductVolumeResidual(nInt1D, nVar, nDOFs1D,
                              standardHex->mLegendreInt1DTranspose,
                              standardHex->mDerLegendreInt1DTranspose,
                              standardHex->mLegendreInt1DTranspose,
                              dVarDy, mRes.data());

  // And the k-direction.
  TensorProductVolumeResidual(nInt1D, nVar, nDOFs1D,
                              standardHex->mLegendreInt1DTranspose,
                              standardHex->mLegendreInt1DTranspose,
                              standardHex->mDerLegendreInt1DTranspose,
                              dVarDz, mRes.data());

  //--------------------------------------------------------------------------
  // Step 7: Compute the contribution of the body force term to the residual.
  //--------------------------------------------------------------------------

  // Test if a body force was specified.
  if( inputParam->mFBodySpecified )
  {
    // Create the non-dimensional body force.
    const su2double pRefInv = one/pRef;
    const su2double fBodyX  = pRefInv*inputParam->mFBody[0];
    const su2double fBodyY  = pRefInv*inputParam->mFBody[1];
    const su2double fBodyZ  = pRefInv*inputParam->mFBody[2];

    // Loop over the integration points.
#pragma omp simd
    for(int l=0; l<nIntPad; ++l)
    {
      // Determine the multiplication weight for the source terms. This weight
      // is the weight of this integration point multiplied by the Jacobian.
      // The minus sign is present, because the source terms are defined on the
      // RHS of the Navier-Stokes equations, while the residual is defined on
      // the LHS.
      const su2double weight = -mVolMetricInt[0][l]*standardHex->mIntWeights[l];

      // Compute the source terms for the momentum and energy equations.
      // Use Var as storage.
      const su2double u = Var[1][l], v = Var[2][l], w = Var[3][l];

      Var[0][l] = weight*fBodyX;
      Var[1][l] = weight*fBodyY;
      Var[2][l] = weight*fBodyZ;
      Var[3][l] = weight*(u*fBodyX + v*fBodyY + w*fBodyZ);
    }

    // Add the source terms to the residuals of the momentum and energy
    // equations. The continuity equation does not have a source term,
    // hence no need to update it.
    TensorProductVolumeResidual(nInt1D, 4, nDOFs1D,
                                standardHex->mLegendreInt1DTranspose,
                                standardHex->mLegendreInt1DTranspose,
                                standardHex->mLegendreInt1DTranspose,
                                &Var[0], &mRes[1]);
  }
}
