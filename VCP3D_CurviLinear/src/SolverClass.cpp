//------------------------------------------------------------------------------
// File, which contains the implementation of the member functions of the
// class SolverClass.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"

//------------------------------------------------------------------------------

// Constructor. Initialize the member variables.
SolverClass::SolverClass(InputParamClass *inputParam)
{
  // Set the pointer to the input parameter.
  mInputParam = inputParam;

  // Initialize the current time step and the number of time steps for the
  // averaging to zero. The latter may get corrected in case of a restart.
  mCurrentTimeStep  = 0;
  mNTimeStepAverage = 0;

  // Initialize the averaged force coefficients. They may get
  // corrected in case of a restart.
  mCfxAve = zero;
  mCfyAve = zero;
  mCfzAve = zero;

  mCfxVisAve = zero;
  mCfyVisAve = zero;
  mCfzVisAve = zero;

	// Initialize the surface area values to zero.
	for(unsigned short iBoundary=0; iBoundary<6; iBoundary++)
		mSurfAreaBoundary[iBoundary] = zero;

  // Intialize the MPI variables, as some of them may not be needed.
#ifdef HAVE_MPI
  mCommBufIMin = NULL;
  mCommBufJMin = NULL;
  mCommBufKMin = NULL;

  mCommBufIMax = NULL;
  mCommBufJMax = NULL;
  mCommBufKMax = NULL;
#endif

  // Determine the data for the standard hexahedron.
  mStandardHex.DataStandardElement(mInputParam);

  // Distribute the grid over the MPI ranks.
  DistributeGrid();

  // Allocate the residual arrays and possibly the solution arrays
  // for the halo elements.
  AllocateResidualAndHaloSolutionArrays();

#ifdef HAVE_MPI
  // Create the MPI persistent communication pattern.
  CreatePersistentCommPattern();
#endif
}

//------------------------------------------------------------------------------

// Destructor. Release the memory.
SolverClass::~SolverClass()
{
  // Release the memory of allocated halo elements.
  for(unsigned int i=0; i<mElements.size(); ++i)
  {
    for(unsigned int j=0; j<mElements[i].size(); ++j)
    {
      for(unsigned int k=0; k<mElements[i][j].size(); ++k)
      {
        // Check if this is not a NULL pointer.
        if( mElements[i][j][k] )
        {
          // Check if this is an allocated halo element. If so, deallocate it.
          if(mElements[i][j][k]->mElemType != INTERNAL_ELEMENT)
            delete mElements[i][j][k];
        }
      }
    }
  }

  // Release the memory of the owned elements.
  for(int i=1; i<=mNElemPerRankI; ++i)
    for(int j=1; j<=mNElemPerRankJ; ++j)
      for(int k=1; k<=mNElemPerRankK; ++k)
        delete mElements[i][j][k];

  // Release the memory of the MPI communication buffers.
#ifdef HAVE_MPI
  FreeMemory((void **) &mCommBufIMin);
  FreeMemory((void **) &mCommBufJMin);
  FreeMemory((void **) &mCommBufKMin);

  FreeMemory((void **) &mCommBufIMax);
  FreeMemory((void **) &mCommBufJMax);
  FreeMemory((void **) &mCommBufKMax);
#endif
}

//-----------------------------------------------------------------------------

// Function, which adds fluctuations to the initial solution, if desired.
void SolverClass::AddFluctuations(void)
{
  // Return immediately if no fluctuations should be added.
  if(mInputParam->mFluctuationsInitSol == NO_FLUCTUATIONS) return;
  if(mCurrentTimeStep && (mInputParam->mFluctuationsInitSol != SYNTHETIC)) return;

  // Store the number of DOFs per element a bit easier.
  const int nDOFs = mStandardHex.mNDOFs;

  // Start of the parallel region.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    // Allocate the memory to store the dimensional primitive variables.
    std::vector<su2double *> primVar(nVar);
    for(int m=0; m<nVar; ++m)
    {
      primVar[m] = NULL;
      primVar[m] = (su2double *) AllocateMemory(nDOFs*sizeof(su2double));
      if( !primVar[m] )
        Terminate("SolverClass::AddFluctuations", __FILE__, __LINE__,
                  "Memory allocation failure for primVar[m].");
    }

    // Loop over the number of elements.
#ifdef HAVE_OPENMP
#pragma omp for collapse(3), schedule(static), nowait
#endif
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        for(int i=1; i<=mNElemPerRankI; ++i)
        {
          // It is easier to add the fluctuations for the dimensional primitive
          // variables. So convert the current variables first.
          mElements[i][j][k]->ComputePrimitiveVariablesNodalDOFs(mInputParam,
                                                                 &mStandardHex,
                                                                 primVar.data());

          // Add the fluctuations to the primitive variables.
          mElements[i][j][k]->AddFluctuationsPrimVar(mInputParam, &mStandardHex,
                                                     mCurrentTimeStep, primVar.data());

          // Convert the primitive variables back to the working variables.
          mElements[i][j][k]->ComputeWorkingVariables(mInputParam, &mStandardHex,
                                                      primVar.data());
        }
      }
    }

    // Release the memory of the primitive variables again.
    for(int m=0; m<nVar; ++m)
      FreeMemory((void **) &primVar[m]);

  } // End of the OpenMP parallel region.
}

//------------------------------------------------------------------------------

// Function, which controls the allocation of the residual and halo solution.
void SolverClass::AllocateResidualAndHaloSolutionArrays(void)
{
  // As also pointers are used for halo elements, the allocation of the
  // residual and halo solution arrays must happen in two loops. First allocate
  // these arrays for true halo elements.
  for(unsigned int i=0; i<mElements.size(); ++i)
  {
    for(unsigned int j=0; j<mElements[i].size(); ++j)
    {
      for(unsigned int k=0; k<mElements[i][j].size(); ++k)
      {
        // Check if this is not a NULL pointer.
        if( mElements[i][j][k] )
        {
          // Check if this is an allocated halo element. If so, allocate the arrays.
          if(mElements[i][j][k]->mElemType != INTERNAL_ELEMENT)
            mElements[i][j][k]->AllocateResidualAndHaloSolutionArrays(&mStandardHex);
        }
      }
    }
  }

  // The second loop is over the owned elements.
  for(int i=1; i<=mNElemPerRankI; ++i)
    for(int j=1; j<=mNElemPerRankJ; ++j)
      for(int k=1; k<=mNElemPerRankK; ++k)
        mElements[i][j][k]->AllocateResidualAndHaloSolutionArrays(&mStandardHex);
}

//------------------------------------------------------------------------------

// Function, which determines the bounding box for the fluctuations.
void SolverClass::BoundingBoxFluctuations(void)
{
  // Relative tolerance to be added to the bounding boxes.
  const su2double tolRel = (su2double) 1.e-5;

  // Determine the number of nodal grid DOFs per element.
  const int nDOFs1DGrid = mInputParam->mNPolyGridDOFs+1;
  const int nDOFs3DGrid = nDOFs1DGrid*nDOFs1DGrid*nDOFs1DGrid;

  // Initialize the bounding box of the local coordinates.
  su2double coorMin[] = { valLarge,  valLarge,  valLarge};
  su2double coorMax[] = {-valLarge, -valLarge, -valLarge};

  // Determine the bounding box of the local coordinates.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static), collapse(3), \
                         reduction(min: coorMin[:3]), \
                         reduction(max: coorMax[:3])
#endif
  for(int k=1; k<=mNElemPerRankK; ++k)
  {
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      for(int i=1; i<=mNElemPerRankI; ++i)
      {
        su2double **coor = mElements[i][j][k]->mCoorNodalGridDOFs.data();
        for(int l=0; l<nDOFs3DGrid; ++l)
        {
          coorMin[0] = std::min(coorMin[0], coor[0][l]);
          coorMin[1] = std::min(coorMin[1], coor[1][l]);
          coorMin[2] = std::min(coorMin[2], coor[2][l]);

          coorMax[0] = std::max(coorMax[0], coor[0][l]);
          coorMax[1] = std::max(coorMax[1], coor[1][l]);
          coorMax[2] = std::max(coorMax[2], coor[2][l]);
        }
      }
    }
  }

  // When MPI is used, the reduction over the mpi ranks must be carried out.
#ifdef HAVE_MPI
  su2double locBuf[] = {coorMin[0], coorMin[1], coorMin[2]};
  MPI_Allreduce(locBuf, coorMin, 3, MPI_SU2DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  locBuf[0] = coorMax[0]; locBuf[1] = coorMax[1]; locBuf[2] = coorMax[2];
  MPI_Allreduce(locBuf, coorMax, 3, MPI_SU2DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

  // Store the coordinates of the bounding box in mInputParam.
  mInputParam->mLowerCoorBoxFluctuations[0] = std::max(coorMin[0], mInputParam->mLowerCoorBoxFluctuations[0]);
  mInputParam->mLowerCoorBoxFluctuations[1] = std::max(coorMin[1], mInputParam->mLowerCoorBoxFluctuations[1]);
  mInputParam->mLowerCoorBoxFluctuations[2] = std::max(coorMin[2], mInputParam->mLowerCoorBoxFluctuations[2]);

  mInputParam->mUpperCoorBoxFluctuations[0] = std::min(coorMax[0], mInputParam->mUpperCoorBoxFluctuations[0]);
  mInputParam->mUpperCoorBoxFluctuations[1] = std::min(coorMax[1], mInputParam->mUpperCoorBoxFluctuations[1]);
  mInputParam->mUpperCoorBoxFluctuations[2] = std::min(coorMax[2], mInputParam->mUpperCoorBoxFluctuations[2]);

  // Add/subtract a tolerance to this bounding box.
  su2double tol = tolRel*(mInputParam->mUpperCoorBoxFluctuations[0]-mInputParam->mLowerCoorBoxFluctuations[0]);
  mInputParam->mLowerCoorBoxFluctuations[0] -= tol;
  mInputParam->mUpperCoorBoxFluctuations[0] += tol;

  tol = tolRel*(mInputParam->mUpperCoorBoxFluctuations[1]-mInputParam->mLowerCoorBoxFluctuations[1]);
  mInputParam->mLowerCoorBoxFluctuations[1] -= tol;
  mInputParam->mUpperCoorBoxFluctuations[1] += tol;

  tol = tolRel*(mInputParam->mUpperCoorBoxFluctuations[2]-mInputParam->mLowerCoorBoxFluctuations[2]);
  mInputParam->mLowerCoorBoxFluctuations[2] -= tol;
  mInputParam->mUpperCoorBoxFluctuations[2] += tol;
}

//------------------------------------------------------------------------------

// Function, which checks if the coordinates of the given faces match.
void SolverClass::CheckMatchingFaces(const int                    rankMin,
                                     const int                    rankMax,
                                     const int                    nElemMin,
                                     const int                    nElemMax,
                                     const std::vector<su2double> &coorMin,
                                     const std::vector<su2double> &coorMinComm,
                                     const std::vector<su2double> &coorMax,
                                     const std::vector<su2double> &coorMaxComm,
                                     const char                   *faceDir)
{
  std::ostringstream message;

  // Set the relative tolerance for matching.
  const su2double tolRel = (su2double) 1.e-4;

  // Determine the number of grid DOFs on a face.
  const int nDOFs1D = mInputParam->mNPolyGridDOFs + 1;
  const int nDOFs2D = nDOFs1D*nDOFs1D;

  // Make a distinction between parallel mode and sequential mode.
#ifdef HAVE_MPI
  // Parallel mode. Send the relevant buffers for the max boundary
  // to the neighboring rank. Note that the messages are sent even if the
  // ranks are the same. Performance is not an issue here.
  MPI_Request intRequest, doubleRequest;
  MPI_Isend(&nElemMax, 1, MPI_INT, rankMax, rank, MPI_COMM_WORLD, &intRequest);

  if(nElemMax > 0)
    MPI_Isend(coorMaxComm.data(), coorMaxComm.size(), MPI_SU2DOUBLE, rankMax,
              rank+1, MPI_COMM_WORLD, &doubleRequest);

  // Receive the integer message. Use a blocking receive.
  int nElemRecv;
  MPI_Recv(&nElemRecv, 1, MPI_INT, rankMin, rankMin, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  std::vector<su2double> coorRecv(nElemRecv*(1+3*nDOFs2D));
  if(nElemRecv > 0)
    MPI_Recv(coorRecv.data(), coorRecv.size(), MPI_SU2DOUBLE, rankMin, rankMin+1,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  // Complete the nonblocking sends.
  MPI_Wait(&intRequest, MPI_STATUS_IGNORE);
  if(nElemMax > 0) MPI_Wait(&doubleRequest, MPI_STATUS_IGNORE);

  // Check if the number of received elements is equal to the number of
  // elements on the min boundary. If not, something is wrong with the
  // provided subface information or there is a bug in the MPI. Just
  // write the first option.
  if(nElemMin != nElemRecv)
  {
    message << "Number of 1 to 1 matching elements in " << faceDir
            << " direction differs for the min and max boundaries.";
    TerminateAll("SolverClass::CheckMatchingFaces", __FILE__, __LINE__,
                  message.str());
  }

  // Check the coordinates on the min boundary. Loop over the number of elements.
  int ind = 0;
  for(int i=0; i<nElemMin; ++i)
  {
    // Determine the length scale, based on the length scale of the two elements.
    // Set the tolerance accordingly.
    const su2double tol = tolRel*std::min(coorMin[ind], coorRecv[ind]);
    ++ind;

    // Loop over the number of coordinates of this face.
    for(int j=0; j<nDOFs2D; ++j, ind+=3)
    {
      // Determine the difference in x-, y- and z-coordinates.
      const su2double dx = coorMin[ind]   - coorRecv[ind];
      const su2double dy = coorMin[ind+1] - coorRecv[ind+1];
      const su2double dz = coorMin[ind+2] - coorRecv[ind+2];

      // Check for matching coordinates. If not, terminate.
      if((FABS(dx) > tol) || (FABS(dy) > tol) || (FABS(dz) > tol))
      {
        message << "1 to 1 matching subfaces in " << faceDir
                << " direction do not match";
        TerminateAll("SolverClass::CheckMatchingFaces", __FILE__, __LINE__,
                  message.str());
      }
    }
  }

  // Send the relevant buffers for the min boundary to the neighboring rank.
  // Note that the messages are sent even if the ranks are the same.
  // Performance is not an issue here.
  MPI_Isend(&nElemMin, 1, MPI_INT, rankMin, rank, MPI_COMM_WORLD, &intRequest);

  if(nElemMin > 0)
    MPI_Isend(coorMinComm.data(), coorMinComm.size(), MPI_SU2DOUBLE, rankMin,
              rank+1, MPI_COMM_WORLD, &doubleRequest);

  // Receive the integer message. Use a blocking receive.
  MPI_Recv(&nElemRecv, 1, MPI_INT, rankMax, rankMax, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  coorRecv.resize(nElemRecv*(1+3*nDOFs2D));
  if(nElemRecv > 0)
    MPI_Recv(coorRecv.data(), coorRecv.size(), MPI_SU2DOUBLE, rankMax, rankMax+1,
             MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  // Complete the nonblocking sends.
  MPI_Wait(&intRequest, MPI_STATUS_IGNORE);
  if(nElemMin > 0) MPI_Wait(&doubleRequest, MPI_STATUS_IGNORE);

  // Check the coordinates on the max boundary. Loop over the number of elements.
  ind = 0;
  for(int i=0; i<nElemMax; ++i)
  {
    // Determine the length scale, based on the length scale of the two elements.
    // Set the tolerance accordingly.
    const su2double tol = tolRel*std::min(coorMax[ind], coorRecv[ind]);
    ++ind;

    // Loop over the number of coordinates of this face.
    for(int j=0; j<nDOFs2D; ++j, ind+=3)
    {
      // Determine the difference in x-, y- and z-coordinates.
      const su2double dx = coorMax[ind]   - coorRecv[ind];
      const su2double dy = coorMax[ind+1] - coorRecv[ind+1];
      const su2double dz = coorMax[ind+2] - coorRecv[ind+2];

      // Check for matching coordinates. If not, terminate.
      if((FABS(dx) > tol) || (FABS(dy) > tol) || (FABS(dz) > tol))
      {
        message << "1 to 1 matching subfaces in " << faceDir
                << " direction do not match";
        TerminateAll("SolverClass::CheckMatchingFaces", __FILE__, __LINE__,
                  message.str());
      }
    }
  }

#else
  // Sequential mode. Check if the number of elements are the same. If not,
  // something is wrong with the provided subface information.
  if(nElemMin != nElemMax)
  {
    message << "Number of 1 to 1 matching elements in " << faceDir
            << " direction differs for the min and max boundaries.";
    TerminateAll("SolverClass::CheckMatchingFaces", __FILE__, __LINE__,
                  message.str());
  }

  // Check the coordinates on the min boundary. Loop over the number of elements.
  int ind = 0;
  for(int i=0; i<nElemMin; ++i)
  {
    // Determine the length scale, based on the length scale of the two elements.
    // Set the tolerance accordingly.
    const su2double tol = tolRel*std::min(coorMin[ind], coorMaxComm[ind]);
    ++ind;

    // Loop over the number of coordinates of this face.
    for(int j=0; j<nDOFs2D; ++j, ind+=3)
    {
      // Determine the difference in x-, y- and z-coordinates.
      const su2double dx = coorMin[ind]   - coorMaxComm[ind];
      const su2double dy = coorMin[ind+1] - coorMaxComm[ind+1];
      const su2double dz = coorMin[ind+2] - coorMaxComm[ind+2];

      // Check for matching coordinates. If not, terminate.
      if((FABS(dx) > tol) || (FABS(dy) > tol) || (FABS(dz) > tol))
      {
        message << "1 to 1 matching subfaces in " << faceDir
                << " direction do not match";
        TerminateAll("SolverClass::CheckMatchingFaces", __FILE__, __LINE__,
                  message.str());
      }
    }
  }

  // Check the coordinates on the max boundary. Loop over the number of elements.
  ind = 0;
  for(int i=0; i<nElemMax; ++i)
  {
    // Determine the length scale, based on the length scale of the two elements.
    // Set the tolerance accordingly.
    const su2double tol = tolRel*std::min(coorMax[ind], coorMinComm[ind]);
    ++ind;

    // Loop over the number of coordinates of this face.
    for(int j=0; j<nDOFs2D; ++j, ind+=3)
    {
      // Determine the difference in x-, y- and z-coordinates.
      const su2double dx = coorMax[ind]   - coorMinComm[ind];
      const su2double dy = coorMax[ind+1] - coorMinComm[ind+1];
      const su2double dz = coorMax[ind+2] - coorMinComm[ind+2];

      // Check for matching coordinates. If not, terminate.
      if((FABS(dx) > tol) || (FABS(dy) > tol) || (FABS(dz) > tol))
      {
        message << "1 to 1 matching subfaces in " << faceDir
                << " direction do not match";
        TerminateAll("SolverClass::CheckMatchingFaces", __FILE__, __LINE__,
                  message.str());
      }
    }
  }

#endif
}

//------------------------------------------------------------------------------

// Function, which checks if matching subfaces actually match.
void SolverClass::CheckMatchingSubfaces(void)
{
  // Easier storage of the number of 1D and 2D grid DOFs.
  const int nDOFs1D = mInputParam->mNPolyGridDOFs + 1;
  const int nDOFs2D = nDOFs1D*nDOFs1D;

  //-----------------------------------------------------------------------
  // Store the coordinates on the 1-to-1 matching iMin faces of the local
  // iMin boundaries in the buffers coorMin and coorMinComm. The number
  // of elements on the iMin boundary for which information is copied
  // is also stored.
  //-----------------------------------------------------------------------

  int nElemMin = 0;
  std::vector<su2double> coorMin, coorMinComm;

  // Loop over the elements on the iMin boundary and test if halo's are present.
  for(int j=1; j<=mNElemPerRankJ; ++j)
  {
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      if( mElements[0][j][k] )
      {
        // Update nElemMin and set the pointers for the grid coordinates
        // and periodic translation to the donor face.
        ++nElemMin;
        su2double **coor = mElements[1][j][k]->mCoorNodalGridDOFs.data();
        su2double *trans = mElements[1][j][k]->mTransIMin;

        // Copy the relevant length scale of this element.
        coorMin.push_back(mElements[1][j][k]->mLenScaleIDir);
        coorMinComm.push_back(mElements[1][j][k]->mLenScaleIDir);

        // Loop over the grid DOFs of the iMin face of the element and
        // copy the coordinates to coorMin and coorMinComm.
        for(int kk=0; kk<nDOFs1D; ++kk)
        {
          for(int jj=0; jj<nDOFs1D; ++jj)
          {
            const int ind = kk*nDOFs2D + jj*nDOFs1D;
            coorMin.push_back(coor[0][ind]);
            coorMin.push_back(coor[1][ind]);
            coorMin.push_back(coor[2][ind]);

            coorMinComm.push_back(coor[0][ind]+trans[0]);
            coorMinComm.push_back(coor[1][ind]+trans[1]);
            coorMinComm.push_back(coor[2][ind]+trans[2]);
          }
        }
      }
    }
  }

  //-----------------------------------------------------------------------
  // Store the coordinates on the 1-to-1 matching iMax faces of the local
  // iMax boundaries in the buffers coorMax and coorMaxComm. The number
  // of elements on the iMax boundary for which information is copied
  // is also stored.
  //-----------------------------------------------------------------------

  int nElemMax = 0;
  std::vector<su2double> coorMax, coorMaxComm;

  // Loop over the elements on the iMax boundary and test if halo's are present.
  for(int j=1; j<=mNElemPerRankJ; ++j)
  {
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      if( mElements[mNElemPerRankI+1][j][k] )
      {
        // Update nElemMax and set the pointers for the grid coordinates
        // and periodic translation to the donor face.
        ++nElemMax;
        su2double **coor = mElements[mNElemPerRankI][j][k]->mCoorNodalGridDOFs.data();
        su2double *trans = mElements[mNElemPerRankI][j][k]->mTransIMax;

        // Copy the relevant length scale of this element.
        coorMax.push_back(mElements[mNElemPerRankI][j][k]->mLenScaleIDir);
        coorMaxComm.push_back(mElements[mNElemPerRankI][j][k]->mLenScaleIDir);

        // Loop over the grid DOFs of the iMax face of the element and
        // copy the coordinates to coorMax and coorMaxComm.
        for(int kk=0; kk<nDOFs1D; ++kk)
        {
          for(int jj=0; jj<nDOFs1D; ++jj)
          {
            const int ind = kk*nDOFs2D + jj*nDOFs1D + nDOFs1D-1;
            coorMax.push_back(coor[0][ind]);
            coorMax.push_back(coor[1][ind]);
            coorMax.push_back(coor[2][ind]);

            coorMaxComm.push_back(coor[0][ind]+trans[0]);
            coorMaxComm.push_back(coor[1][ind]+trans[1]);
            coorMaxComm.push_back(coor[2][ind]+trans[2]);
          }
        }
      }
    }
  }

  //---------------------------------------------------
  // Check for matching coordinates in i-direction.
  //---------------------------------------------------

  // Determine the ranks of the processors on the iMin and iMax boundaries.
  int rankNeig = mMyRankI-1;
  if(rankNeig == -1) rankNeig = mNRanksI-1;
  int rankMin = rankNeig + mMyRankJ*mNRanksI + mMyRankK*mNRanksI*mNRanksJ;

  rankNeig = mMyRankI+1;
  if(rankNeig == mNRanksI) rankNeig = 0;
  int rankMax = rankNeig + mMyRankJ*mNRanksI + mMyRankK*mNRanksI*mNRanksJ;

  // Call the function CheckMatchingFaces to do the actual checking.
  CheckMatchingFaces(rankMin, rankMax, nElemMin, nElemMax, coorMin,
                     coorMinComm, coorMax, coorMaxComm, "I");

  //-----------------------------------------------------------------------
  // Store the coordinates on the 1-to-1 matching jMin faces of the local
  // jMin boundaries in the buffers coorMin and coorMinComm. The number
  // of elements on the jMin boundary for which information is copied
  // is also stored.
  //-----------------------------------------------------------------------

  nElemMin = 0;
  coorMin.resize(0);
  coorMinComm.resize(0);

  // Loop over the elements on the jMin boundary and test if halo's are present.
  for(int i=1; i<=mNElemPerRankI; ++i)
  {
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      if( mElements[i][0][k] )
      {
        // Update nElemMin and set the pointers for the grid coordinates
        // and periodic translation to the donor face.
        ++nElemMin;
        su2double **coor = mElements[i][1][k]->mCoorNodalGridDOFs.data();
        su2double *trans = mElements[i][1][k]->mTransJMin;

        // Copy the relevant length scale of this element.
        coorMin.push_back(mElements[i][1][k]->mLenScaleJDir);
        coorMinComm.push_back(mElements[i][1][k]->mLenScaleJDir);

        // Loop over the grid DOFs of the jMin face of the element and
        // copy the coordinates to coorMin and coorMinComm.
        for(int kk=0; kk<nDOFs1D; ++kk)
        {
          for(int ii=0; ii<nDOFs1D; ++ii)
          {
            const int ind = kk*nDOFs2D + ii;
            coorMin.push_back(coor[0][ind]);
            coorMin.push_back(coor[1][ind]);
            coorMin.push_back(coor[2][ind]);

            coorMinComm.push_back(coor[0][ind]+trans[0]);
            coorMinComm.push_back(coor[1][ind]+trans[1]);
            coorMinComm.push_back(coor[2][ind]+trans[2]);
          }
        }
      }
    }
  }

  //-----------------------------------------------------------------------
  // Store the coordinates on the 1-to-1 matching jMax faces of the local
  // jMax boundaries in the buffers coorMax and coorMaxComm. The number
  // of elements on the jMax boundary for which information is copied
  // is also stored.
  //-----------------------------------------------------------------------

  nElemMax = 0;
  coorMax.resize(0);
  coorMaxComm.resize(0);

  // Loop over the elements on the jMax boundary and test if halo's are present.
  for(int i=1; i<=mNElemPerRankI; ++i)
  {
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      if( mElements[i][mNElemPerRankJ+1][k] )
      {
        // Update nElemMax and set the pointers for the grid coordinates
        // and periodic translation to the donor face.
        ++nElemMax;
        su2double **coor = mElements[i][mNElemPerRankJ][k]->mCoorNodalGridDOFs.data();
        su2double *trans = mElements[i][mNElemPerRankJ][k]->mTransJMax;

        // Copy the relevant length scale of this element.
        coorMax.push_back(mElements[i][mNElemPerRankJ][k]->mLenScaleJDir);
        coorMaxComm.push_back(mElements[i][mNElemPerRankJ][k]->mLenScaleJDir);

        // Loop over the grid DOFs of the jMax face of the element and
        // copy the coordinates to coorMax and coorMaxComm.
        for(int kk=0; kk<nDOFs1D; ++kk)
        {
          for(int ii=0; ii<nDOFs1D; ++ii)
          {
            const int ind = kk*nDOFs2D + ii + (nDOFs1D-1)*nDOFs1D;
            coorMax.push_back(coor[0][ind]);
            coorMax.push_back(coor[1][ind]);
            coorMax.push_back(coor[2][ind]);

            coorMaxComm.push_back(coor[0][ind]+trans[0]);
            coorMaxComm.push_back(coor[1][ind]+trans[1]);
            coorMaxComm.push_back(coor[2][ind]+trans[2]);
          }
        }
      }
    }
  }

  //---------------------------------------------------
  // Check for matching coordinates in j-direction.
  //---------------------------------------------------

  // Determine the ranks of the processors on the jMin and jMax boundaries.
  rankNeig = mMyRankJ-1;
  if(rankNeig == -1) rankNeig = mNRanksJ-1;
  rankMin = mMyRankI + rankNeig*mNRanksI + mMyRankK*mNRanksI*mNRanksJ;

  rankNeig = mMyRankJ+1;
  if(rankNeig == mNRanksJ) rankNeig = 0;
  rankMax = mMyRankI + rankNeig*mNRanksI + mMyRankK*mNRanksI*mNRanksJ;

  // Call the function CheckMatchingFaces to do the actual checking.
  CheckMatchingFaces(rankMin, rankMax, nElemMin, nElemMax, coorMin,
                     coorMinComm, coorMax, coorMaxComm, "J");

  //-----------------------------------------------------------------------
  // Store the coordinates on the 1-to-1 matching kMin faces of the local
  // kMin boundaries in the buffers coorMin and coorMinComm. The number
  // of elements on the kMin boundary for which information is copied
  // is also stored.
  //-----------------------------------------------------------------------

  nElemMin = 0;
  coorMin.resize(0);
  coorMinComm.resize(0);

  // Loop over the elements on the kMin boundary and test if halo's are present.
  for(int i=1; i<=mNElemPerRankI; ++i)
  {
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      if( mElements[i][j][0] )
      {
        // Update nElemMin and set the pointers for the grid coordinates
        // and periodic translation to the donor face.
        ++nElemMin;
        su2double **coor = mElements[i][j][1]->mCoorNodalGridDOFs.data();
        su2double *trans = mElements[i][j][1]->mTransKMin;

        // Copy the relevant length scale of this element.
        coorMin.push_back(mElements[i][j][1]->mLenScaleKDir);
        coorMinComm.push_back(mElements[i][j][1]->mLenScaleKDir);

        // Loop over the grid DOFs of the kMin face of the element and
        // copy the coordinates to coorMin and coorMinComm.
        for(int jj=0; jj<nDOFs1D; ++jj)
        {
          for(int ii=0; ii<nDOFs1D; ++ii)
          {
            const int ind = jj*nDOFs1D + ii;
            coorMin.push_back(coor[0][ind]);
            coorMin.push_back(coor[1][ind]);
            coorMin.push_back(coor[2][ind]);

            coorMinComm.push_back(coor[0][ind]+trans[0]);
            coorMinComm.push_back(coor[1][ind]+trans[1]);
            coorMinComm.push_back(coor[2][ind]+trans[2]);
          }
        }
      }
    }
  }

  //-----------------------------------------------------------------------
  // Store the coordinates on the 1-to-1 matching kMax faces of the local
  // kMax boundaries in the buffers coorMax and coorMaxComm. The number
  // of elements on the kMax boundary for which information is copied
  // is also stored.
  //-----------------------------------------------------------------------

  nElemMax = 0;
  coorMax.resize(0);
  coorMaxComm.resize(0);

  // Loop over the elements on the kMax boundary and test if halo's are present.
  for(int i=1; i<=mNElemPerRankI; ++i)
  {
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      if( mElements[i][j][mNElemPerRankK+1] )
      {
        // Update nElemMax and set the pointers for the grid coordinates
        // and periodic translation to the donor face.
        ++nElemMax;
        su2double **coor = mElements[i][j][mNElemPerRankK]->mCoorNodalGridDOFs.data();
        su2double *trans = mElements[i][j][mNElemPerRankK]->mTransKMax;

        // Copy the relevant length scale of this element.
        coorMax.push_back(mElements[i][j][mNElemPerRankK]->mLenScaleKDir);
        coorMaxComm.push_back(mElements[i][j][mNElemPerRankK]->mLenScaleKDir);

        // Loop over the grid DOFs of the kMax face of the element and
        // copy the coordinates to coorMax and coorMaxComm.
        for(int jj=0; jj<nDOFs1D; ++jj)
        {
          for(int ii=0; ii<nDOFs1D; ++ii)
          {
            const int ind = jj*nDOFs1D + ii + (nDOFs1D-1)*nDOFs2D;
            coorMax.push_back(coor[0][ind]);
            coorMax.push_back(coor[1][ind]);
            coorMax.push_back(coor[2][ind]);

            coorMaxComm.push_back(coor[0][ind]+trans[0]);
            coorMaxComm.push_back(coor[1][ind]+trans[1]);
            coorMaxComm.push_back(coor[2][ind]+trans[2]);
          }
        }
      }
    }
  }

  //---------------------------------------------------
  // Check for matching coordinates in k-direction.
  //---------------------------------------------------

  // Determine the ranks of the processors on the kMin and kMax boundaries.
  rankNeig = mMyRankK-1;
  if(rankNeig == -1) rankNeig = mNRanksK-1;
  rankMin = mMyRankI + mMyRankJ*mNRanksI + rankNeig*mNRanksI*mNRanksJ;

  rankNeig = mMyRankK+1;
  if(rankNeig == mNRanksK) rankNeig = 0;
  rankMax = mMyRankI + mMyRankJ*mNRanksI + rankNeig*mNRanksI*mNRanksJ;

  // Call the function CheckMatchingFaces to do the actual checking.
  CheckMatchingFaces(rankMin, rankMax, nElemMin, nElemMax, coorMin,
                     coorMinComm, coorMax, coorMaxComm, "K");
}

//------------------------------------------------------------------------------

// Function, which computes the metric terms.
void SolverClass::ComputeMetricTerms(void)
{
  // Define the variables to store the number of equidistant and LGL elements.
  int nElemEquiDistant = 0, nElemLGL = 0;

  // Start of the parallel region, if required.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    // Loop over the elements to determine the number of elements for which the
    // grid points are at the LGL location and the number of elements for which
    // an equidistant spacing is used.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3), \
                reduction(+: nElemEquiDistant), \
                reduction(+: nElemLGL)
#endif
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        for(int i=1; i<=mNElemPerRankI; ++i)
        {
          if( mElements[i][j][k]->LGLDistribution(mInputParam, &mStandardHex) )
            ++nElemLGL;
          else
            ++nElemEquiDistant;
        }
      }
    }

    // For a parallel executable, determine the global data for
    // nElemEquiDistant and nElemLGL.
#ifdef HAVE_MPI
#ifdef HAVE_OPENMP
#pragma omp single
#endif
    {
      int locBuf[] = {nElemEquiDistant, nElemLGL};
      int globBuf[2];
      MPI_Allreduce(locBuf, globBuf, 2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      nElemEquiDistant = globBuf[0];
      nElemLGL         = globBuf[1];
    }
#endif

    // Check if either all elements are equidistant or LGL. If not, terminate.
    if(nElemEquiDistant && nElemLGL)
      TerminateAll("SolverClass::ComputeMetricTerms", __FILE__, __LINE__,
                   "Elements with both equidistant and LGL spacing encountered");

    // Determine whether or not LGL spacing is used.
    const bool LGLSpacing = nElemLGL > 0;

    // Determine the location of the 1D grid DOFs in the standard element.
    // Note that these location are not necessarily equal to the member variable
    // mRDOFs1D of standardHex, because it is possible that a different polynomial
    // degree is used for the grid and solution, as well as a different positioning
    // of the DOFs.
    const int nPolyGrid = mInputParam->mNPolyGridDOFs;

    const int nDOFs1DGrid = nPolyGrid+1;

#ifdef HAVE_OPENMP
#pragma omp single
#endif
    {
      if( LGLSpacing ) LocationDOFs1D(nDOFs1DGrid, LGL_POINTS,  mStandardHex.mRDOFsGrid1D);
      else             LocationDOFs1D(nDOFs1DGrid, EQUIDISTANT, mStandardHex.mRDOFsGrid1D);
    }

    // Allocate the memory for the 1D Lagrangian basis functions and its derivatives
    // in both the integration points and on the min and max faces.
    const int sizeLag = nDOFs1DGrid*mStandardHex.mNIntegration1DPad;
    su2double *lagrangeInt1D    = (su2double *) AllocateMemory(sizeLag*sizeof(su2double));
    su2double *derLagrangeInt1D = (su2double *) AllocateMemory(sizeLag*sizeof(su2double));

    if(!lagrangeInt1D || !derLagrangeInt1D)
      Terminate("SolverClass::ComputeMetricTerms", __FILE__, __LINE__,
                "Memory allocation failure for interior Lagrangian basis functions");

    const int nDOFs1DGridPad = ((nDOFs1DGrid+vecLen1D-1)/vecLen1D)*vecLen1D;
    su2double *lagrangeMinFace1D    = (su2double *) AllocateMemory(nDOFs1DGridPad*sizeof(su2double));
    su2double *lagrangeMaxFace1D    = (su2double *) AllocateMemory(nDOFs1DGridPad*sizeof(su2double));
    su2double *derLagrangeMinFace1D = (su2double *) AllocateMemory(nDOFs1DGridPad*sizeof(su2double));
    su2double *derLagrangeMaxFace1D = (su2double *) AllocateMemory(nDOFs1DGridPad*sizeof(su2double));

    if(!lagrangeMinFace1D || !lagrangeMaxFace1D || !derLagrangeMinFace1D || !derLagrangeMaxFace1D)
      Terminate("SolverClass::ComputeMetricTerms", __FILE__, __LINE__,
                "Memory allocation failure for face Lagrangian basis functions");

    // Initialize these variables to zero, such that the padded values get a value.
    for(int i=0; i<sizeLag; ++i)
      lagrangeInt1D[i] = derLagrangeInt1D[i] = zero;

    for(int i=0; i<nDOFs1DGridPad; ++i)
      lagrangeMinFace1D[i]    = lagrangeMaxFace1D[i]    =
      derLagrangeMinFace1D[i] = derLagrangeMaxFace1D[i] = zero;

    // Set the corresponding entry of lagrangeMinFace1D and
    // lagrangeMaxFace1D to 1.
    lagrangeMinFace1D[0] = lagrangeMaxFace1D[nDOFs1DGrid-1] = one;

    // Compute the other Lagrangian basis functions.
    mStandardHex.LagrangianBasisFunctions(mStandardHex.mRDOFsGrid1D, mStandardHex.mRIntegration1D,
                                          lagrangeInt1D, derLagrangeInt1D,
                                          derLagrangeMinFace1D, derLagrangeMaxFace1D);

    // Allocate the memory for Lagrangian basis functions of the grid DOFs in the
    // solution DOFs and its derivatives.
    const int sizeLagSolDOFs     = nDOFs1DGrid*mStandardHex.mNDOFs1DPad;
    su2double *lagrangeDOFs1D    = (su2double *) AllocateMemory(sizeLagSolDOFs*sizeof(su2double));
    su2double *derLagrangeDOFs1D = (su2double *) AllocateMemory(sizeLagSolDOFs*sizeof(su2double));

    if(!lagrangeDOFs1D || !derLagrangeDOFs1D)
      Terminate("SolverClass::ComputeMetricTerms", __FILE__, __LINE__,
                "Memory allocation failure for lagrangeDOFs1D and derLagrangeDOFs1D");

    // Initialize these variables to zero, such that the padded values get a value.
    for(int i=0; i<sizeLagSolDOFs; ++i)
      lagrangeDOFs1D[i] = derLagrangeDOFs1D[i] = zero;

    // Compute the Lagrangian basis functions and its derivatives
    // for the nodal solution DOFs.
    mStandardHex.LagrangianBasisFunctions(mStandardHex.mRDOFsGrid1D, mStandardHex.mRDOFs1D,
                                          lagrangeDOFs1D, derLagrangeDOFs1D, NULL, NULL);

    // Loop over the elements to compute the metric terms, the length scales
    // and the coordinates of the nodal solution DOFs.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3)
#endif
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        for(int i=1; i<=mNElemPerRankI; ++i)
        {
          mElements[i][j][k]->ComputeMetricTerms(nDOFs1DGrid, mStandardHex.mNIntegration1D,
                                                 lagrangeInt1D, derLagrangeInt1D,
                                                 lagrangeMinFace1D, lagrangeMaxFace1D,
                                                 derLagrangeMinFace1D, derLagrangeMaxFace1D);

          mElements[i][j][k]->ComputeLengthScales(mInputParam, &mStandardHex);

          mElements[i][j][k]->ComputeCoorNodalSolDOFs(nDOFs1DGrid, mStandardHex.mNDOFs1D,
                                                      lagrangeDOFs1D, derLagrangeDOFs1D);
        }
      }
    }

    // Release the memory of the Lagrangian basis functions again.
    FreeMemory((void **) &lagrangeInt1D);
    FreeMemory((void **) &derLagrangeInt1D);
    FreeMemory((void **) &lagrangeMinFace1D);
    FreeMemory((void **) &lagrangeMaxFace1D);
    FreeMemory((void **) &derLagrangeMinFace1D);
    FreeMemory((void **) &derLagrangeMaxFace1D);
    FreeMemory((void **) &lagrangeDOFs1D);
    FreeMemory((void **) &derLagrangeDOFs1D);

  } // End of the OpenMP parallel region.

  //------------------------------------------------------
  // Communicate the metric terms of the max boundaries.
  //------------------------------------------------------

  // For a parallel computation the metric terms on the max
  // boundaries of the halo elements must be communicated.
#ifdef HAVE_MPI
  CommunicateMetricTermsMaxBoundaries();
#endif

  //------------------------------------------------------
  // Check for matching subfaces on the block boundaries.
  //------------------------------------------------------

  CheckMatchingSubfaces();

  //---------------------------------------------------------------
  // Determine the interpolation weights for the exchange location
  // for the wall model, if needed.
  //---------------------------------------------------------------

  InterpolationWeightsExchangeLocation();
}

//------------------------------------------------------------------------------

// Function, which determines the precision of the restart file.
short SolverClass::DeterminePrecisionRestartFile(const char *filename)
{
  // Open the restart file for binary reading.
  FILE *restartFile = fopen(filename, "rb");
  if( !restartFile )
    Terminate("SolverClass::DeterminePrecisionRestartFile", __FILE__, __LINE__,
              std::string(filename) + " could not be opened for binary reading");

  // Read the first 5 integer variables of the solution file, check if it
  // is in the correct format and determine the number of variables and
  // number of DOFs.
  int header[5];
  if(std::fread(header, sizeof(int), 5, restartFile) != 5)
    Terminate("SolverClass::DeterminePrecisionRestartFile", __FILE__, __LINE__,
              "File header could not be read");

  // Check if byte swapping must be applied.
  if(header[0] != SU2_MagicNumber)
    SwapBytes(header, sizeof(int), 5);

  // Check if this is indeed an SU2 restart file.
  if(header[0] != SU2_MagicNumber)
    Terminate("SolverClass::DeterminePrecisionRestartFile", __FILE__, __LINE__,
              "Restart file is not an SU2 restart file");

  // Easier storage of the number of variables and DOFS in the restart file.
  const int nVarRestart  = header[1];
  const int nDOFsRestart = header[2];

  // Store the current position in the file.
  long int currentPos = ftell(restartFile);

  // Jump to the end of the file and determine the number of bytes in between
  // the current position and the end.
  if(fseek(restartFile, 0L, SEEK_END))
    Terminate("SolverClass::DeterminePrecisionRestartFile", __FILE__, __LINE__,
              "Failed to jump to the end of the file");
  long int nBytesFloatingPoint = ftell(restartFile) - currentPos;
  fclose(restartFile);

  // Correct nBytesFloatingPoint for the variable names and iteration number.
  nBytesFloatingPoint -= nVarRestart*CGNS_STRING_SIZE + sizeof(int);

  // Determine the total number of floating point items stored in nBytesFloatingPoint,
  // which is the solution and 8 items for the meta data.
  const int nItemsFloatingPoint = nVarRestart*nDOFsRestart + 8;

  // Determine the number of bytes per floating point item, which is returned.
  short nBytesPerItem = (short) round(((double) nBytesFloatingPoint)
                      /               ((double) nItemsFloatingPoint));
  return nBytesPerItem;
}

//------------------------------------------------------------------------------

// Function, which determines the visualization data.
void SolverClass::DetermineVisualizationData(std::vector<std::string> &visNames,
                                             std::vector<float>       &visBuf)
{
  // Define the names of the variables to be visualized. Note that a vector
  // quantity must be defined contiguously, with first the x-component, followed
  // by the y-component and finally the z-component. Furthermore, the names of
  // the components of a vector field must end with _x, _y and _z, respectively.
  visNames.push_back("Density");
  visNames.push_back("Momentum_x");
  visNames.push_back("Momentum_y");
  visNames.push_back("Momentum_z");
  visNames.push_back("Energy");
  visNames.push_back("Pressure");
  visNames.push_back("Velocity_x");
  visNames.push_back("Velocity_y");
  visNames.push_back("Velocity_z");
  visNames.push_back("Temperature");
  visNames.push_back("Mach");
  visNames.push_back("Vorticity_x");
  visNames.push_back("Vorticity_y");
  visNames.push_back("Vorticity_z");
  visNames.push_back("QCriterion");

  // Add the ratio of eddy viscosity and laminar viscosity if a subgrid
  // scale model is used.
  if(mInputParam->mSGSModelType != NO_SGS_MODEL)
    visNames.push_back("RatioEddyViscosityLaminarViscosity");

  // Add the average variables if these have been computed.
  if( mInputParam->mComputeAverageSolution )
  {
    visNames.push_back("AverageDensity");
    visNames.push_back("AverageVelocity_x");
    visNames.push_back("AverageVelocity_y");
    visNames.push_back("AverageVelocity_z");
    visNames.push_back("AveragePressure");
    visNames.push_back("AverageTemperature");
    visNames.push_back("AverageMach");

    visNames.push_back("AverageR_uu");
    visNames.push_back("AverageR_vv");
    visNames.push_back("AverageR_ww");
    visNames.push_back("AverageR_uv");
    visNames.push_back("AverageR_uw");
    visNames.push_back("AverageR_vw");

    // Add the ratio of eddy viscosity and laminar viscosity if a subgrid
    // scale model is used.
    if(mInputParam->mSGSModelType != NO_SGS_MODEL)
      visNames.push_back("AverageRatioEddyViscosityLaminarViscosity");
  }

  // Determine the number of elements per MPI rank.
  const int nElemIJ  = mNElemPerRankI*mNElemPerRankJ;
  const int nElemIJK = nElemIJ*mNElemPerRankK;

  // Determine the number of DOFs per MPI rank.
  const int nDOFs    = mStandardHex.mNDOFs;
  const int nDOFsLoc = nElemIJK*nDOFs;

  // Allocate the memory for the visualization buffer and initialize
  // the data to zero.
  visBuf.resize(nDOFsLoc*visNames.size(), zero);

// Start of the parallel region.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    // Allocate the memory to store the dimensional primitive variables.
    // The additional entry is used for the eddy viscosity.
    std::vector<su2double *> primVar(nVar+1, NULL);
    for(int m=0; m<=nVar; ++m)
    {
      primVar[m] = (su2double *) AllocateMemory(nDOFs*sizeof(su2double));
      if( !primVar[m] )
        Terminate("SolverClass::DetermineVisualizationData", __FILE__, __LINE__,
                  "Memory allocation failure for primVar[m].");
    }

    // Allocate the memory to store the gradients of the velocities.
    // Note that it must be allocated with size 3*nVar, because it is also
    // used as temporary storage for the nodal gradients of the working
    // variables.
    std::vector<su2double *> gradVel(3*nVar, NULL);
    for(int m=0; m<(3*nVar); ++m)
    {
      gradVel[m] = (su2double *) AllocateMemory(nDOFs*sizeof(su2double));
      if( !gradVel[m] )
        Terminate("SolverClass::DetermineVisualizationData", __FILE__, __LINE__,
                  "Memory allocation failure for gradVel[m].");
    }

    // Loop over the owned elements.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), nowait
#endif
    for(int elem=0; elem<nElemIJK; ++elem)
    {
      // Retrieve the i-, j- and k-index of the element. Use a zero based
      // index here, because that is easier for the retrieval.
      int k   = elem/nElemIJ;
      int rem = elem - k*nElemIJ;
      int j   = rem/mNElemPerRankI;
      int i   = rem - j*mNElemPerRankI;

      // Increment i, j and k, because the owned elements start at 1.
      ++i; ++j; ++k;

      // Determine the dimensional primitive variables.
      mElements[i][j][k]->ComputePrimitiveVariablesNodalDOFs(mInputParam, &mStandardHex,
                                                             primVar.data());

      // Determine the gradients of the velocities in the DOFs.
      mElements[i][j][k]->ComputeVelocityGradientsNodalDOFs(mInputParam, &mStandardHex,
                                                            primVar.data(), gradVel.data());

      // Set the pointers for the average variables.
      su2double **avePrim    = mElements[i][j][k]->mAvePrim.data();
      su2double **aveVelProd = mElements[i][j][k]->mAveVelProd.data();
      su2double *aveEddyVis  = mElements[i][j][k]->mAveEddyVis;   

      // Loop over the number of variables to be determined.
      for(unsigned int visVar=0; visVar<visNames.size(); ++visVar)
      {
        // Check for the variable name and act accordingly.
        if(visNames[visVar] == "Density")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Store the density in the buffer.
          for(int l=0; l<nDOFs; ++l)
            buf[l] = (float) primVar[0][l];
        }
        else if(visNames[visVar] == "Momentum_x")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + 3*elem*nDOFs;

          // Momentum in x-direction. As momentum is a vector variable,
          // also the y- and z-momentum must be stored.
          for(int l=0; l<nDOFs; ++l)
          {
            const int l3 = 3*l;
            buf[l3]   = (float) (primVar[0][l]*primVar[1][l]);
            buf[l3+1] = (float) (primVar[0][l]*primVar[2][l]);
            buf[l3+2] = (float) (primVar[0][l]*primVar[3][l]);
          }

          // Update visVar by an additional 2, because the y- and z-momentum
          // have been stored as well.
          visVar += 2;
        }
        else if(visNames[visVar] == "Energy")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Abbreviate 1/(gam-1).
          const su2double ovgm1 = one/(GamConstant - one);

          // Store the total energy in the buffer.
          for(int l=0; l<nDOFs; ++l)
          {
            const su2double eKin = half*(primVar[1][l]*primVar[1][l]
                                 +       primVar[2][l]*primVar[2][l]
                                 +       primVar[3][l]*primVar[3][l]);
            const su2double ETot = ovgm1*primVar[4][l] + primVar[0][l]*eKin;
            buf[l] = (float) ETot;
          }
        }
        else if(visNames[visVar] == "Pressure")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Store the pressure in the buffer.
          for(int l=0; l<nDOFs; ++l)
            buf[l] = (float) primVar[4][l];
        }
        else if(visNames[visVar] == "Velocity_x")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + 3*elem*nDOFs;

          // Velocity in x-direction. As velocity is a vector variable,
          // also the y- and z-velocity must be stored.
          for(int l=0; l<nDOFs; ++l)
          {
            const int l3 = 3*l;
            buf[l3]   = (float) primVar[1][l];
            buf[l3+1] = (float) primVar[2][l];
            buf[l3+2] = (float) primVar[3][l];
          }

          // Update visVar by an additional 2, because the y- and z-velocity
          // have been stored as well.
          visVar += 2;
        }
        else if(visNames[visVar] == "Temperature")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Store the temperature in the buffer.
          for(int l=0; l<nDOFs; ++l)
          {
            const su2double T = primVar[4][l]/(RGas*primVar[0][l]);
            buf[l] = (float) T;
          }
        }
        else if(visNames[visVar] == "Mach")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Store the Mach number in the buffer.
          for(int l=0; l<nDOFs; ++l)
          {
            const su2double V2 = primVar[1][l]*primVar[1][l]
                               + primVar[2][l]*primVar[2][l]
                               + primVar[3][l]*primVar[3][l];
            const su2double a2 = GamConstant*primVar[4][l]/primVar[0][l];
            const su2double Ma = SQRT(V2/a2);
            buf[l] = (float) Ma;
          }
        }
        else if(visNames[visVar] == "Vorticity_x")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + 3*elem*nDOFs;

          // X-component of the vorticity. As vorticity is a vector variable,
          // also the y- and z-components must be stored.
          for(int l=0; l<nDOFs; ++l)
          {
            const int l3 = 3*l;
            buf[l3]   = (float) (gradVel[5][l] - gradVel[7][l]);  // dwdy - dvdz
            buf[l3+1] = (float) (gradVel[6][l] - gradVel[2][l]);  // dudz - dwdx
            buf[l3+2] = (float) (gradVel[1][l] - gradVel[3][l]);  // dvdx - dudy
          }

          // Update visVar by an additional 2, because the y- and z-vorticity
          // have been stored as well.
          visVar += 2;
        }
        else if(visNames[visVar] == "QCriterion")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Q-criterion must be stored in the buffer.
          for(int l=0; l<nDOFs; ++l)
          {
            // Compute the components of the strain and vorticity tensors.
            const su2double s11 = gradVel[0][l];                         // dudx
            const su2double s12 = half*(gradVel[3][l] + gradVel[1][l]);  // (dudy + dvdx)/2
            const su2double s13 = half*(gradVel[6][l] + gradVel[2][l]);  // (dudz + dwdx)/2
            const su2double s22 = gradVel[4][l];                         // dvdy
            const su2double s23 = half*(gradVel[7][l] + gradVel[5][l]);  // (dvdz + dwdy)/2
            const su2double s33 = gradVel[8][l];                         // dwdz

            const su2double o12 = half*(gradVel[3][l] - gradVel[1][l]);  // (dudy - dvdx)/2
            const su2double o13 = half*(gradVel[6][l] - gradVel[2][l]);  // (dudz - dwdx)/2
            const su2double o23 = half*(gradVel[7][l] - gradVel[5][l]);  // (dvdz - dwdy)/2

            // Compute the Q-criterion and store it in the buffer.
            const su2double Q = two*(o12*o12 + o13*o13 + o23*o23 - s12*s12 - s13*s13 - s23*s23)
                              - s11*s11 - s22*s22 - s33*s33;

            buf[l] = (float) Q;
          }
        }
        else if(visNames[visVar] == "RatioEddyViscosityLaminarViscosity")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // The ratio of the eddy viscosity and laminar viscosity must be stored in the buffer.
          // Set the pointer to the entry in primVar where the eddy viscosity should be stored
          // and compute it.
          su2double *eddyVis = primVar[nVar];

          mInputParam->mSGSModel->EddyViscosity(nDOFs, mElements[i][j][k]->mLenScaleLES,
                                                primVar.data(), gradVel.data(), eddyVis);
 
          // Loop over the DOFs and store the ratio of eddy viscosity and
          // laminar viscosity in buf.
          const su2double muDimInv = uRef/(mu*pRef);

          for(int l=0; l<nDOFs; ++l)
          {
            const su2double ratio = muDimInv*eddyVis[l];
            buf[l] = (float) ratio;
          }
        }
        else if(visNames[visVar] == "AverageDensity")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Store the average density in the buffer.
          for(int l=0; l<nDOFs; ++l)
            buf[l] = (float) avePrim[0][l];
        }
        else if(visNames[visVar] == "AverageVelocity_x")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + 3*elem*nDOFs;

          // Average velocity in x-direction. As velocity is a vector
          // variable, also the y- and z-velocity must be stored.
          for(int l=0; l<nDOFs; ++l)
          {
            const int l3 = 3*l;
            buf[l3]   = (float) avePrim[1][l];
            buf[l3+1] = (float) avePrim[2][l];
            buf[l3+2] = (float) avePrim[3][l];
          }

          // Update visVar by an additional 2, because the y- and z-velocity
          // have been stored as well.
          visVar += 2;
        }
        else if(visNames[visVar] == "AveragePressure")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Store the average pressure in the buffer.
          for(int l=0; l<nDOFs; ++l)
            buf[l] = (float) avePrim[4][l];
        }
        else if(visNames[visVar] == "AverageTemperature")
        {
          // Check if the average temperature can be computed.
          if(mNTimeStepAverage > 0)
          {
            // Set the pointer to the correction location in visBuf.
            float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

            // Store the average temperature in the buffer.
            for(int l=0; l<nDOFs; ++l)
            {
              const su2double T = avePrim[4][l]/(RGas*avePrim[0][l]);
              buf[l] = (float) T;
            }
          }
        }
        else if(visNames[visVar] == "AverageMach")
        {
          // Check if the average Mach number can be computed.
          if(mNTimeStepAverage > 0)
          {
            // Set the pointer to the correction location in visBuf.
            float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

            // Store the average Mach number in the buffer.
            for(int l=0; l<nDOFs; ++l)
            {
              const su2double V2 = avePrim[1][l]*avePrim[1][l]
                                 + avePrim[2][l]*avePrim[2][l]
                                 + avePrim[3][l]*avePrim[3][l];
              const su2double a2 = GamConstant*avePrim[4][l]/avePrim[0][l];
              const su2double Ma = SQRT(V2/a2);
              buf[l] = (float) Ma;
            }
          }
        }
        else if(visNames[visVar] == "AverageR_uu")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Store the average u'u' Reynolds stress.
          for(int l=0; l<nDOFs; ++l)
            buf[l] = (float) (aveVelProd[0][l] - avePrim[1][l]*avePrim[1][l]);
        }
        else if(visNames[visVar] == "AverageR_vv")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Store the average v'v' Reynolds stress.
          for(int l=0; l<nDOFs; ++l)
            buf[l] = (float) (aveVelProd[1][l] - avePrim[2][l]*avePrim[2][l]);
        }
        else if(visNames[visVar] == "AverageR_ww")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Store the average w'w' Reynolds stress.
          for(int l=0; l<nDOFs; ++l)
            buf[l] = (float) (aveVelProd[2][l] - avePrim[3][l]*avePrim[3][l]);
        }
        else if(visNames[visVar] == "AverageR_uv")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Store the average u'v' Reynolds stress.
          for(int l=0; l<nDOFs; ++l)
            buf[l] = (float) (aveVelProd[3][l] - avePrim[1][l]*avePrim[2][l]);
        }
        else if(visNames[visVar] == "AverageR_uw")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Store the average u'w' Reynolds stress.
          for(int l=0; l<nDOFs; ++l)
            buf[l] = (float) (aveVelProd[4][l] - avePrim[1][l]*avePrim[3][l]);
        }
        else if(visNames[visVar] == "AverageR_vw")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Store the average v'w' Reynolds stress.
          for(int l=0; l<nDOFs; ++l)
            buf[l] = (float) (aveVelProd[5][l] - avePrim[2][l]*avePrim[3][l]);
        }
        else if(visNames[visVar] == "AverageRatioEddyViscosityLaminarViscosity")
        {
          // Set the pointer to the correction location in visBuf.
          float *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs;

          // Store the average ratio of eddy and laminar viscosity.
          const su2double muDimInv = uRef/(mu*pRef);
          for(int l=0; l<nDOFs; ++l)
            buf[l] = (float) (muDimInv*aveEddyVis[l]);
        }
        else
        {
          // Unknown visualization variable encountered.
#ifdef HAVE_OPENMP
          if(omp_get_thread_num() == 0)
#endif
            TerminateAll("SolverClass::DetermineVisualizationData",
                         __FILE__, __LINE__, "Unknown visualization name");
        }
      }
    }

    // Release the memory of the primitive variables
    // and velocity gradients again.
    for(int m=0; m<=nVar; ++m)
      FreeMemory((void **) &primVar[m]);

    for(int m=0; m<(3*nVar); ++m)
      FreeMemory((void **) &gradVel[m]);

  } // End of the OpenMP parallel region.
}

//------------------------------------------------------------------------------

// Function, which distributes the grid over the MPI ranks.
void SolverClass::DistributeGrid(void)
{
  // Determine the number of elements in the three index direction. It
  // suffices to look at the i-index of the iMax subface, etc.
  const int nElemI = mInputParam->mSubfaces[1][0]->mElemIBeg;
  const int nElemJ = mInputParam->mSubfaces[3][0]->mElemJBeg;
  const int nElemK = mInputParam->mSubfaces[5][0]->mElemKBeg;

  //----------------------------------------------------------------------------
  // Determine the number of ranks in the three coordinate directions.
  //----------------------------------------------------------------------------

  // Initialize the number of ranks in the three directions to 1.
  mNRanksI = mNRanksJ = mNRanksK = 1;

  // Initialize the number of elements in the three directions per rank to
  // the total number of elements in these directions.
  mNElemPerRankI = nElemI;
  mNElemPerRankJ = nElemJ;
  mNElemPerRankK = nElemK;

  // Infinite loop to determine the number of ranks in the three directions.
  for(;;)
  {
    // Criterion to exit the loop
    if(mNRanksI*mNRanksJ*mNRanksK == nRanks) break;

    // Determine the new number of elements in each direction if the number of
    // ranks would be doubled in this direction. If this leads to an uneven split
    // in the number of elements, nI, nJ and nK are set to -1 to indicate this.
    int nI = mNElemPerRankI%2 ? -1 : mNElemPerRankI/2;
    int nJ = mNElemPerRankJ%2 ? -1 : mNElemPerRankJ/2;
    int nK = mNElemPerRankK%2 ? -1 : mNElemPerRankK/2;

    // Overrule nI, nJ and nK if it is not allowed to split the block in
    // the corresponding direction.
    if( !mInputParam->mAllowSplitIDir ) nI = -1;
    if( !mInputParam->mAllowSplitJDir ) nJ = -1;
    if( !mInputParam->mAllowSplitKDir ) nK = -1;

    // Check if a division is possible at all. If not, terminate.
    if((nI == -1) && (nJ == -1) && (nK == -1))
      TerminateAll("SolverClass::DistributeGrid", __FILE__, __LINE__,
                   "Not possible to evenly distribute the grid over the ranks. Change the number of ranks.");

    // Determine in which direction the number of ranks must be increased.
    if((nI > nJ) && (nI > nK)) {mNElemPerRankI /= 2; mNRanksI *= 2;}
    else if(nJ > nK)           {mNElemPerRankJ /= 2; mNRanksJ *= 2;}
    else                       {mNElemPerRankK /= 2; mNRanksK *= 2;}
  }

  // Determine the i-, j- and k-indices of this rank.
  mMyRankK = rank/(mNRanksI*mNRanksJ);
  const int remainder = rank - mMyRankK*mNRanksI*mNRanksJ;
  mMyRankJ = remainder/mNRanksI;
  mMyRankI = remainder - mMyRankJ*mNRanksI;

  // Determine the global indices of 1st owned element of this rank.
  const int iLowRank = mMyRankI*mNElemPerRankI;
  const int jLowRank = mMyRankJ*mNElemPerRankJ;
  const int kLowRank = mMyRankK*mNElemPerRankK;

  // Allocate the memory for mElements. In mElements a halo layer is stored,
  // which possibly contains the information of the neighboring element. Hence
  // the size in each coordinate direction is the number of elements in this
  // direction plus 2 (one halo on each side).
  mElements.resize(mNElemPerRankI+2);
  for(int i=0; i<(mNElemPerRankI+2); ++i)
  {
    mElements[i].resize(mNElemPerRankJ+2);
    for(int j=0; j<(mNElemPerRankJ+2); ++j)
    {
      mElements[i][j].resize(mNElemPerRankK+2);
      for(int k=0; k<(mNElemPerRankK+2); ++k)
        mElements[i][j][k] = NULL;
    }
  }

  // Allocate the memory for the owned elements. Note that the starting
  // index of owned elements is 1.
  for(int i=1; i<=mNElemPerRankI; ++i)
  {
    int globIndElem[3];
    globIndElem[0] = iLowRank + i-1;
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      globIndElem[1] = jLowRank + j-1;
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        globIndElem[2] = kLowRank + k-1;
        mElements[i][j][k] = new ElementClass(mInputParam, &mStandardHex,
                                              globIndElem, INTERNAL_ELEMENT);
        mElements[i][j][k]->StoreLocalIndices(i, j, k);
      }
    }
  }

  //----------------------------------------
  // Elements on the iMin boundary.
  //----------------------------------------

  // Create the double int vector to store the subface index for the cells.
  std::vector<std::vector<int> > subfaceIDBoundElem;
  subfaceIDBoundElem.resize(nElemJ);
  for(int j=0; j<nElemJ; ++j)
    subfaceIDBoundElem[j].assign(nElemK, -1);

  // Loop over the subfaces of the iMin boundary and fill subfaceIDBoundElem
  // with the correct data.
  for(unsigned int l=0; l<mInputParam->mSubfaces[0].size(); ++l)
  {
    // Determine the minimum and maximum of the j- and k-ranges.
    const int jBeg = std::min(mInputParam->mSubfaces[0][l]->mElemJBeg,
                              mInputParam->mSubfaces[0][l]->mElemJEnd);
    const int jEnd = std::max(mInputParam->mSubfaces[0][l]->mElemJBeg,
                              mInputParam->mSubfaces[0][l]->mElemJEnd);
    const int kBeg = std::min(mInputParam->mSubfaces[0][l]->mElemKBeg,
                              mInputParam->mSubfaces[0][l]->mElemKEnd);
    const int kEnd = std::max(mInputParam->mSubfaces[0][l]->mElemKBeg,
                              mInputParam->mSubfaces[0][l]->mElemKEnd);

    // Loop over the element range of this subface and set subface ID.
    for(int j=jBeg; j<jEnd; ++j)
      for(int k=kBeg; k<kEnd; ++k)
        subfaceIDBoundElem[j][k] = l;
  }

  // Determine the situation.
  if(mMyRankI == 0)
  {
    // The iMin boundary is a true block boundary. Loop over the owned
    // elements adjacent to this boundary.
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        // Determine the global j and k indices of this element and
        // the corresponding subface ID..
        const int jG  = mElements[1][j][k]->mGlobalInd[1];
        const int kG  = mElements[1][j][k]->mGlobalInd[2];
        const int sID = subfaceIDBoundElem[jG][kG];

        // Check if this is a physical boundary condition.
        if( mInputParam->mSubfaces[0][sID]->BCIsPhysical() )
        {
          // Set the iMin boundary to a physical boundary.
          mElements[1][j][k]->mBCIMin = mInputParam->mSubfaces[0][sID];
        }
        else
        {
          // The iMin boundary is an internal boundary. Check if there
          // is only one rank in i-direction.
          if(mNRanksI == 1)
          {
            // It suffices to set the pointer to the halo element.
            mElements[0][j][k] = mElements[mNElemPerRankI][j][k];
          }
          else
          {
            // Multiple ranks in i-direction, so the halo must be allocated.
            int globIndElem[] = {nElemI-1, jLowRank+j-1, kLowRank+k-1};
            mElements[0][j][k] = new ElementClass(mInputParam, &mStandardHex,
                                                  globIndElem, HALO_ELEMENT_IMIN);
          }

          // Store the possible periodic transformation.
          mInputParam->mSubfaces[0][sID]->GetPeriodicTransformation(mElements[1][j][k]->mTransIMin);
        }
      }
    }
  }
  else
  {
    // The iMin boundary is an internal boundary, which is created by the
    // distribution over the MPI ranks. The elements must be allocated.
    int globIndElem[3];
    globIndElem[0] = iLowRank-1;

    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      globIndElem[1] = jLowRank + j-1;
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        globIndElem[2] = kLowRank + k-1;
        mElements[0][j][k] = new ElementClass(mInputParam, &mStandardHex,
                                              globIndElem, HALO_ELEMENT_IMIN);
      }
    }
  }

  //----------------------------------------
  // Elements on the iMax boundary.
  //----------------------------------------

  // Loop over the subfaces of the iMax boundary and fill subfaceIDBoundElem
  // with the correct data.
  for(unsigned int l=0; l<mInputParam->mSubfaces[1].size(); ++l)
  {
    // Determine the minimum and maximum of the j- and k-ranges.
    const int jBeg = std::min(mInputParam->mSubfaces[1][l]->mElemJBeg,
                              mInputParam->mSubfaces[1][l]->mElemJEnd);
    const int jEnd = std::max(mInputParam->mSubfaces[1][l]->mElemJBeg,
                              mInputParam->mSubfaces[1][l]->mElemJEnd);
    const int kBeg = std::min(mInputParam->mSubfaces[1][l]->mElemKBeg,
                              mInputParam->mSubfaces[1][l]->mElemKEnd);
    const int kEnd = std::max(mInputParam->mSubfaces[1][l]->mElemKBeg,
                              mInputParam->mSubfaces[1][l]->mElemKEnd);

    // Loop over the element range of this subface and set subface ID.
    for(int j=jBeg; j<jEnd; ++j)
      for(int k=kBeg; k<kEnd; ++k)
        subfaceIDBoundElem[j][k] = l;
  }

  // Determine the situation.
  if(mMyRankI == (mNRanksI-1))
  {
    // The iMax boundary is a true block boundary. Loop over the owned
    // elements adjacent to this boundary.
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        // Determine the global j and k indices of this element and
        // the corresponding subface ID..
        const int jG  = mElements[mNElemPerRankI][j][k]->mGlobalInd[1];
        const int kG  = mElements[mNElemPerRankI][j][k]->mGlobalInd[2];
        const int sID = subfaceIDBoundElem[jG][kG];

        // Check if this is a physical boundary condition.
        if( mInputParam->mSubfaces[1][sID]->BCIsPhysical() )
        {
          // Set the iMax boundary to a physical boundary.
          mElements[mNElemPerRankI][j][k]->mBCIMax = mInputParam->mSubfaces[1][sID];
        }
        else
        {
          // The iMax boundary is an internal boundary. Check if there
          // is only one rank in i-direction.
          if(mNRanksI == 1)
          {
            // It suffices to set the pointer to the halo element.
            mElements[mNElemPerRankI+1][j][k] = mElements[1][j][k];
          }
          else
          {
            // Multiple ranks in i-direction, so the halo must be allocated.
            int globIndElem[] = {0, jLowRank+j-1, kLowRank+k-1};
            mElements[mNElemPerRankI+1][j][k] = new ElementClass(mInputParam, &mStandardHex,
                                                                 globIndElem, HALO_ELEMENT_IMAX);
          }

          // Store the possible periodic transformation.
          mInputParam->mSubfaces[1][sID]->GetPeriodicTransformation(mElements[mNElemPerRankI][j][k]->mTransIMax);
        }
      }
    }
  }
  else
  {
    // The iMax boundary is an internal boundary, which is created by the
    // distribution over the MPI ranks. The elements must be allocated.
    int globIndElem[3];
    globIndElem[0] = iLowRank + mNElemPerRankI;

    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      globIndElem[1] = jLowRank + j-1;
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        globIndElem[2] = kLowRank + k-1;
        mElements[mNElemPerRankI+1][j][k] = new ElementClass(mInputParam, &mStandardHex,
                                                             globIndElem, HALO_ELEMENT_IMAX);
      }
    }
  }

  //----------------------------------------
  // Elements on the jMin boundary.
  //----------------------------------------

  // Resize the double vector subfaceIDBoundElem, such that it can be used
  // on the jMin and jMax boundaries.
  subfaceIDBoundElem.resize(nElemI);
  for(int i=0; i<nElemI; ++i)
    subfaceIDBoundElem[i].assign(nElemK, -1);

  // Loop over the subfaces of the jMin boundary and fill subfaceIDBoundElem
  // with the correct data.
  for(unsigned int l=0; l<mInputParam->mSubfaces[2].size(); ++l)
  {
    // Determine the minimum and maximum of the i- and k-ranges.
    const int iBeg = std::min(mInputParam->mSubfaces[2][l]->mElemIBeg,
                              mInputParam->mSubfaces[2][l]->mElemIEnd);
    const int iEnd = std::max(mInputParam->mSubfaces[2][l]->mElemIBeg,
                              mInputParam->mSubfaces[2][l]->mElemIEnd);
    const int kBeg = std::min(mInputParam->mSubfaces[2][l]->mElemKBeg,
                              mInputParam->mSubfaces[2][l]->mElemKEnd);
    const int kEnd = std::max(mInputParam->mSubfaces[2][l]->mElemKBeg,
                              mInputParam->mSubfaces[2][l]->mElemKEnd);

    // Loop over the element range of this subface and set subface ID.
    for(int i=iBeg; i<iEnd; ++i)
      for(int k=kBeg; k<kEnd; ++k)
        subfaceIDBoundElem[i][k] = l;
  }

  // Determine the situation.
  if(mMyRankJ == 0)
  {
    // The jMin boundary is a true block boundary. Loop over the owned
    // elements adjacent to this boundary.
    for(int i=1; i<=mNElemPerRankI; ++i)
    {
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        // Determine the global i and k indices of this element and
        // the corresponding subface ID..
        const int iG  = mElements[i][1][k]->mGlobalInd[0];
        const int kG  = mElements[i][1][k]->mGlobalInd[2];
        const int sID = subfaceIDBoundElem[iG][kG];

        // Check if this is a physical boundary condition.
        if( mInputParam->mSubfaces[2][sID]->BCIsPhysical() )
        {
          // Set the jMin boundary to a physical boundary.
          mElements[i][1][k]->mBCJMin = mInputParam->mSubfaces[2][sID];
        }
        else
        {
          // The jMin boundary is an internal boundary. Check if there
          // is only one rank in j-direction.
          if(mNRanksJ == 1)
          {
            // It suffices to set the pointer to the halo element.
            mElements[i][0][k] = mElements[i][mNElemPerRankJ][k];
          }
          else
          {
            // Multiple ranks in j-direction, so the halo must be allocated.
            int globIndElem[] = {iLowRank+i-1, nElemJ-1, kLowRank+k-1};
            mElements[i][0][k] = new ElementClass(mInputParam, &mStandardHex,
                                                  globIndElem, HALO_ELEMENT_JMIN);
          }

          // Store the possible periodic transformation.
          mInputParam->mSubfaces[2][sID]->GetPeriodicTransformation(mElements[i][1][k]->mTransJMin);
        }
      }
    }
  }
  else
  {
    // The jMin boundary is an internal boundary, which is created by the
    // distribution over the MPI ranks. The elements must be allocated.
    int globIndElem[3];
    globIndElem[1] = jLowRank-1;

    for(int i=1; i<=mNElemPerRankI; ++i)
    {
      globIndElem[0] = iLowRank + i-1;
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        globIndElem[2] = kLowRank + k-1;
        mElements[i][0][k] = new ElementClass(mInputParam, &mStandardHex,
                                                  globIndElem, HALO_ELEMENT_JMIN);
      }
    }
  }

  //----------------------------------------
  // Elements on the jMax boundary.
  //----------------------------------------

  // Loop over the subfaces of the jMax boundary and fill subfaceIDBoundElem
  // with the correct data.
  for(unsigned int l=0; l<mInputParam->mSubfaces[3].size(); ++l)
  {
    // Determine the minimum and maximum of the i- and k-ranges.
    const int iBeg = std::min(mInputParam->mSubfaces[3][l]->mElemIBeg,
                              mInputParam->mSubfaces[3][l]->mElemIEnd);
    const int iEnd = std::max(mInputParam->mSubfaces[3][l]->mElemIBeg,
                              mInputParam->mSubfaces[3][l]->mElemIEnd);
    const int kBeg = std::min(mInputParam->mSubfaces[3][l]->mElemKBeg,
                              mInputParam->mSubfaces[3][l]->mElemKEnd);
    const int kEnd = std::max(mInputParam->mSubfaces[3][l]->mElemKBeg,
                              mInputParam->mSubfaces[3][l]->mElemKEnd);

    // Loop over the element range of this subface and set subface ID.
    for(int i=iBeg; i<iEnd; ++i)
      for(int k=kBeg; k<kEnd; ++k)
        subfaceIDBoundElem[i][k] = l;
  }

  // Determine the situation.
  if(mMyRankJ == (mNRanksJ-1))
  {
    // The jMax boundary is a true block boundary. Loop over the owned
    // elements adjacent to this boundary.
    for(int i=1; i<=mNElemPerRankI; ++i)
    {
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        // Determine the global i and k indices of this element and
        // the corresponding subface ID..
        const int iG  = mElements[i][mNElemPerRankJ][k]->mGlobalInd[0];
        const int kG  = mElements[i][mNElemPerRankJ][k]->mGlobalInd[2];
        const int sID = subfaceIDBoundElem[iG][kG];

        // Check if this is a physical boundary condition.
        if( mInputParam->mSubfaces[3][sID]->BCIsPhysical() )
        {
          // Set the jMax boundary to a physical boundary.
          mElements[i][mNElemPerRankJ][k]->mBCJMax = mInputParam->mSubfaces[3][sID];
        }
        else
        {
          // The jMax boundary is an internal boundary. Check if there
          // is only one rank in j-direction.
          if(mNRanksJ == 1)
          {
            // It suffices to set the pointer to the halo element.
            mElements[i][mNElemPerRankJ+1][k] = mElements[i][1][k];
          }
          else
          {
            // Multiple ranks in j-direction, so the halo must be allocated.
            int globIndElem[] = {iLowRank+i-1, 0, kLowRank+k-1};
            mElements[i][mNElemPerRankJ+1][k] = new ElementClass(mInputParam, &mStandardHex,
                                                                 globIndElem, HALO_ELEMENT_JMAX);
          }

          // Store the possible periodic transformation.
          mInputParam->mSubfaces[3][sID]->GetPeriodicTransformation(mElements[i][mNElemPerRankJ][k]->mTransJMax);
        }
      }
    }
  }
  else
  {
    // The jMax boundary is an internal boundary, which is created by the
    // distribution over the MPI ranks. The elements must be allocated.
    int globIndElem[3];
    globIndElem[1] = jLowRank + mNElemPerRankJ;

    for(int i=1; i<=mNElemPerRankI; ++i)
    {
      globIndElem[0] = iLowRank + i-1;
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        globIndElem[2] = kLowRank + k-1;
        mElements[i][mNElemPerRankJ+1][k] = new ElementClass(mInputParam, &mStandardHex,
                                                             globIndElem, HALO_ELEMENT_JMAX);
      }
    }
  }

  //----------------------------------------
  // Elements on the kMin boundary.
  //----------------------------------------

  // Resize the double vector subfaceIDBoundElem, such that it can be used
  // on the kMin and kMax boundaries.
  for(int i=0; i<nElemI; ++i)
    subfaceIDBoundElem[i].assign(nElemJ, -1);

  // Loop over the subfaces of the kMin boundary and fill subfaceIDBoundElem
  // with the correct data.
  for(unsigned int l=0; l<mInputParam->mSubfaces[4].size(); ++l)
  {
    // Determine the minimum and maximum of the i- and k-ranges.
    const int iBeg = std::min(mInputParam->mSubfaces[4][l]->mElemIBeg,
                              mInputParam->mSubfaces[4][l]->mElemIEnd);
    const int iEnd = std::max(mInputParam->mSubfaces[4][l]->mElemIBeg,
                              mInputParam->mSubfaces[4][l]->mElemIEnd);
    const int jBeg = std::min(mInputParam->mSubfaces[4][l]->mElemJBeg,
                              mInputParam->mSubfaces[4][l]->mElemJEnd);
    const int jEnd = std::max(mInputParam->mSubfaces[4][l]->mElemJBeg,
                              mInputParam->mSubfaces[4][l]->mElemJEnd);

    // Loop over the element range of this subface and set subface ID.
    for(int i=iBeg; i<iEnd; ++i)
      for(int j=jBeg; j<jEnd; ++j)
        subfaceIDBoundElem[i][j] = l;
  }

  // Determine the situation.
  if(mMyRankK == 0)
  {
    // The kMin boundary is a true block boundary. Loop over the owned
    // elements adjacent to this boundary.
    for(int i=1; i<=mNElemPerRankI; ++i)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        // Determine the global i and j indices of this element and
        // the corresponding subface ID..
        const int iG  = mElements[i][j][1]->mGlobalInd[0];
        const int jG  = mElements[i][j][1]->mGlobalInd[1];
        const int sID = subfaceIDBoundElem[iG][jG];

        // Check if this is a physical boundary condition.
        if( mInputParam->mSubfaces[4][sID]->BCIsPhysical() )
        {
          // Set the kMin boundary to a physical boundary.
          mElements[i][j][1]->mBCKMin = mInputParam->mSubfaces[4][sID];
        }
        else
        {
          // The kMin boundary is an internal boundary. Check if there
          // is only one rank in k-direction.
          if(mNRanksK == 1)
          {
            // It suffices to set the pointer to the halo element.
            mElements[i][j][0] = mElements[i][j][mNElemPerRankK];
          }
          else
          {
            // Multiple ranks in k-direction, so the halo must be allocated.
            int globIndElem[] = {iLowRank+i-1, jLowRank+j-1, nElemK-1};
            mElements[i][j][0] = new ElementClass(mInputParam, &mStandardHex,
                                                  globIndElem, HALO_ELEMENT_KMIN);
          }

          // Store the possible periodic transformation.
          mInputParam->mSubfaces[4][sID]->GetPeriodicTransformation(mElements[i][j][1]->mTransKMin);
        }
      }
    }
  }
  else
  {
    // The kMin boundary is an internal boundary, which is created by the
    // distribution over the MPI ranks. The elements must be allocated.
    int globIndElem[3];
    globIndElem[2] = kLowRank-1;

    for(int i=1; i<=mNElemPerRankI; ++i)
    {
      globIndElem[0] = iLowRank + i-1;
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        globIndElem[1] = jLowRank + j-1;
        mElements[i][j][0] = new ElementClass(mInputParam, &mStandardHex,
                                              globIndElem, HALO_ELEMENT_KMIN);
      }
    }
  }

  //----------------------------------------
  // Elements on the kMax boundary.
  //----------------------------------------

  // Loop over the subfaces of the kMax boundary and fill subfaceIDBoundElem
  // with the correct data.
  for(unsigned int l=0; l<mInputParam->mSubfaces[5].size(); ++l)
  {
    // Determine the minimum and maximum of the i- and k-ranges.
    const int iBeg = std::min(mInputParam->mSubfaces[5][l]->mElemIBeg,
                              mInputParam->mSubfaces[5][l]->mElemIEnd);
    const int iEnd = std::max(mInputParam->mSubfaces[5][l]->mElemIBeg,
                              mInputParam->mSubfaces[5][l]->mElemIEnd);
    const int jBeg = std::min(mInputParam->mSubfaces[5][l]->mElemJBeg,
                              mInputParam->mSubfaces[5][l]->mElemJEnd);
    const int jEnd = std::max(mInputParam->mSubfaces[5][l]->mElemJBeg,
                              mInputParam->mSubfaces[5][l]->mElemJEnd);

    // Loop over the element range of this subface and set subface ID.
    for(int i=iBeg; i<iEnd; ++i)
      for(int j=jBeg; j<jEnd; ++j)
        subfaceIDBoundElem[i][j] = l;
  }

  // Determine the situation.
  if(mMyRankK == (mNRanksK-1))
  {
    // The kMax boundary is a true block boundary. Loop over the owned
    // elements adjacent to this boundary.
    for(int i=1; i<=mNElemPerRankI; ++i)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        // Determine the global i and j indices of this element and
        // the corresponding subface ID..
        const int iG  = mElements[i][j][mNElemPerRankK]->mGlobalInd[0];
        const int jG  = mElements[i][j][mNElemPerRankK]->mGlobalInd[1];
        const int sID = subfaceIDBoundElem[iG][jG];

        // Check if this is a physical boundary condition.
        if( mInputParam->mSubfaces[5][sID]->BCIsPhysical() )
        {
          // Set the kMax boundary to a physical boundary.
          mElements[i][j][mNElemPerRankK]->mBCKMax = mInputParam->mSubfaces[5][sID];
        }
        else
        {
          // The kMax boundary is an internal boundary. Check if there
          // is only one rank in k-direction.
          if(mNRanksK == 1)
          {
            // It suffices to set the pointer to the halo element.
            mElements[i][j][mNElemPerRankK+1] = mElements[i][j][1];
          }
          else
          {
            // Multiple ranks in k-direction, so the halo must be allocated.
            int globIndElem[] = {iLowRank+i-1, jLowRank+j-1, 0};
            mElements[i][j][mNElemPerRankK+1] = new ElementClass(mInputParam, &mStandardHex,
                                                                 globIndElem, HALO_ELEMENT_KMAX);
          }

          // Store the possible periodic transformation.
          mInputParam->mSubfaces[5][sID]->GetPeriodicTransformation(mElements[i][j][mNElemPerRankK]->mTransKMax);
        }
      }
    }
  }
  else
  {
    // The kMax boundary is an internal boundary, which is created by the
    // distribution over the MPI ranks. The elements must be allocated.
    int globIndElem[3];
    globIndElem[2] = kLowRank + mNElemPerRankK;

    for(int i=1; i<=mNElemPerRankI; ++i)
    {
      globIndElem[0] = iLowRank + i-1;
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        globIndElem[1] = jLowRank + j-1;
        mElements[i][j][mNElemPerRankK+1] = new ElementClass(mInputParam, &mStandardHex,
                                                             globIndElem, HALO_ELEMENT_KMAX);
      }
    }
  }
}

//------------------------------------------------------------------------------

// Function, which carries out one time synchronization step.
void SolverClass::ExecuteTimeSynchronizationStep(su2double *monitoringData)
{
  // Definition of the functions f1 and f2, which give the multiplication factor
  // for the inviscid and viscous terms in the time step definition as a function
  // of the polynomial degree.
  const su2double f1[] = {(su2double)  2.0, (su2double)  6.0, (su2double) 12.0,
                          (su2double) 20.0, (su2double) 28.0, (su2double) 38.0,
                          (su2double) 50.0, (su2double) 64.0, (su2double) 80.0,
                          (su2double) 98.0};
  const su2double f2[] = {(su2double)     6.0, (su2double)   36.0, (su2double)   150.0,
                          (su2double)   420.0, (su2double)  980.0, (su2double)  1975.0,
                          (su2double)  3575.0, (su2double) 7000.0, (su2double) 14000.0,
                          (su2double) 28000.0};

  // Set the values of the Runge Kutta stages.
  const int nRKStages = 3;
  const su2double alphaRK[] = {one, three/four, one/three};
  const su2double betaRK[]  = {one, one/four,   two/three};

  // Determine the number of DOFs per element and its padded value.
  const int nDOFs    = mStandardHex.mNDOFs;
  const int nDOFsPad = mStandardHex.mNDOFsPad;

  // While loop until the synchronization time has been reached.
  bool synchronizationTimeReached = false, firstTimeStep = true;
  su2double timeElapsed = zero;
  while( !synchronizationTimeReached )
  {
    // Update the current time step.
    ++mCurrentTimeStep;

    // Definition of the shared variables needed in the parallel region.
    su2double deltaT;

#ifdef HAVE_MPI
    MPI_Request reduceRequest;
#endif

    // Start of the parallel region.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
    {
      //----------------------------------------------------------------
      // Determine the time step dictated by the stability requirements.
      //----------------------------------------------------------------

      // Allocate the memory to store the solution and eddy viscosity
      // in the nodal DOFs.
      std::vector<su2double *> solNodal(nVar+1);
      for(int m=0; m<=nVar; ++m)
      {
        solNodal[m] = NULL;
        solNodal[m] = (su2double *) AllocateMemory(nDOFs*sizeof(su2double));
        if( !solNodal[m] )
          Terminate("SolverClass::ExecuteTimeSynchronizationStep", __FILE__,
                    __LINE__, "Memory allocation failure for solNodal[m].");
      }

      // Allocate the memory to store the gradients of the solution, if needed.
      std::vector<su2double *> gradSolNodal;
      if(mInputParam->mSGSModelType != NO_SGS_MODEL)
      {
        gradSolNodal.resize(3*nVar);
        for(int m=0; m<(3*nVar); ++m)
        {
          gradSolNodal[m] = NULL;
          gradSolNodal[m] = (su2double *) AllocateMemory(nDOFs*sizeof(su2double));
          if( !gradSolNodal[m] )
            Terminate("SolverClass::ExecuteTimeSynchronizationStep", __FILE__,
                      __LINE__, "Memory allocation failure for gradSolNodal[m].");
        }
      }

      // Initialize the value of deltaT to a large value.
#ifdef HAVE_OPENMP
#pragma omp single
#endif
      deltaT = (su2double) 1.e+10;

      // Loop over the local elements to determine the local value of dtMin.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3), reduction(min:deltaT)
#endif
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        for(int j=1; j<=mNElemPerRankJ; ++j)
        {
          for(int i=1; i<=mNElemPerRankI; ++i)
          {
            // Compute the timestep of this element. If it is less than the
            // currently stored value, update the latter.
            const su2double dtElem = mElements[i][j][k]->ComputeTimeStep(mInputParam, &mStandardHex,
                                                                         f1, f2, solNodal.data(),
                                                                         gradSolNodal.data());
            deltaT = std::min(deltaT, dtElem);
          }
        }
      }

      // Determine the global value of deltaT if MPI is used. If also OpenMP
      // is used, only one thread needs to post this communication request.
      // No need for the other threads to wait until the MPI call is finished.
#ifdef HAVE_MPI
#ifdef HAVE_OPENMP
#pragma omp single nowait
#endif
      {
        su2double locBuf = deltaT;
        MPI_Iallreduce(&locBuf, &deltaT, 1, MPI_SU2DOUBLE, MPI_MIN,
                       MPI_COMM_WORLD, &reduceRequest);
      }
#endif

      // Release the memory of the vectors used to compute the time step.
      for(unsigned long i=0; i<solNodal.size(); ++i)
        FreeMemory((void **) &solNodal[i]);

      for(unsigned long i=0; i<gradSolNodal.size(); ++i)
        FreeMemory((void **) &gradSolNodal[i]);

      // Store the values of the zeroth Runge Kutta stage. The nowait clause
      // is added because the parallel region ends after this loop.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3), nowait
#endif
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        for(int j=1; j<=mNElemPerRankJ; ++j)
        {
          for(int i=1; i<=mNElemPerRankI; ++i)
          {
            for(int l=0; l<nVar; ++l)
            {
              const su2double *sol = mElements[i][j][k]->mSol[l];
              su2double *solOld    = mElements[i][j][k]->mSolOld[l];
#pragma omp simd
              for(int m=0; m<nDOFsPad; ++m)
                solOld[m] = sol[m];
            }
          }
        }
      }
    } // End of the OpenMP parallel region.

    // In case MPI is used, the non-blocking collective communication
    // must be completed first.
#ifdef HAVE_MPI
    MPI_Wait(&reduceRequest, MPI_STATUS_IGNORE);
#endif
    // If this is the first time step, check if the synchronization step
    // is significantly lower than the actual time step. If this happens
    // print a warning. The division by uRef is there, because deltaT is
    // non-dimensional. Keep in mind that the reference length is 1.0
    // and is therefore ignored.
    if( firstTimeStep )
    {
      if(mInputParam->mDeltaTSynchr < half*deltaT/uRef)
      {
        if(rank == 0)
        {
          std::cout << "#" << std::endl;
          std::cout << "#                 WARNING" << std::endl;
          std::cout << "#-----------------------------------------------------" << std::endl;
          std::cout << "# The time synchronization step, " << mInputParam->mDeltaTSynchr
                    << "," << std::endl;
          std::cout << "# is significantly lower than the time step needed for" << std::endl;
          std::cout << "# stability, " << deltaT/uRef << ". This is not efficient." << std::endl;
          std::cout << "#-----------------------------------------------------" << std::endl;
          std::cout << "#" << std::endl << std::flush;
        }
      }

      // Set firstTimeStep to false.
      firstTimeStep = false;
    }

    // Compute the non-dimensional synchronization step and check if it has been reached.
    // If so, adapt the time step and set synchronizationTimeReached to true.
    const su2double deltaTSynchr = mInputParam->mDeltaTSynchr*uRef;
    const su2double tNew         = timeElapsed + deltaT;
    if(tNew >= ((su2double) 0.99)*deltaTSynchr)
    {
      deltaT = deltaTSynchr - timeElapsed;
      synchronizationTimeReached = true;
    }

    // Update the value of timeElapsed.
    timeElapsed += deltaT;
    
    // Loop over the number of RK stages.
    for(int stage=0; stage<nRKStages; ++stage)
    {
      // Easier storage of some of the coefficients of this stage.
      const su2double alpha       = alphaRK[stage];
      const su2double oneMinAlpha = one - alpha;
      const su2double betaDt      = betaRK[stage]*deltaT;

      // Loop over the local elements to determine the new solution.
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(3), schedule(static)
#endif
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        for(int j=1; j<=mNElemPerRankJ; ++j)
        {
          for(int i=1; i<=mNElemPerRankI; ++i)
          {
            for(int l=0; l<nVar; ++l)
            {
              su2double *sol = mElements[i][j][k]->mSol[l];
              const su2double *solOld = mElements[i][j][k]->mSolOld[l];
              const su2double *res    = mElements[i][j][k]->mRes[l];
#pragma omp simd
              for(int m=0; m<nDOFsPad; ++m)
                sol[m] = oneMinAlpha*sol[m] + alpha*solOld[m] - betaDt*res[m];
            }
          }
        }
      }

      // Determine whether or not the monitoringData must be
      // computed inside the Residual function.
      bool ComputeMonitoringData = false;
      if(synchronizationTimeReached && (stage == (nRKStages-1)))
        ComputeMonitoringData = true;

			// Compute the average data required in the NSCBC, in case imposed.
			if( mInputParam->mNSCBC_Specified ) AverageDataNSCBC();

      // Compute the residual of the new solution.
      Residual(monitoringData, ComputeMonitoringData);
    }

    // Possible add the fluctuations for the turbulence.
    AddFluctuations();

  } // End of the while loop for the synchronization time step.

  // If needed, update the average data.
  if( mInputParam->mComputeAverageSolution )
  {
    ++mNTimeStepAverage;
    UpdateAverageData(monitoringData);
  }
}

//------------------------------------------------------------------------------

// Function, which initializes the solution.
void SolverClass::InitSol(void)
{
  // Determine the non-dimensional free-stream primitive variables and
  // store them a bit easier.
  mInputParam->DetermineNonDimPrimVarFreeStream(mNonDimPrimVarFreeStream);

  // Loop over the owned elements.
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(3), schedule(static)
#endif
  for(int k=1; k<=mNElemPerRankK; ++k)
    for(int j=1; j<=mNElemPerRankJ; ++j)
      for(int i=1; i<=mNElemPerRankI; ++i)
        mElements[i][j][k]->InitSol(mInputParam, &mStandardHex,
                                    mNonDimPrimVarFreeStream);
}

//-----------------------------------------------------------------------------

// Function, which initializes the NSCBC parameters.
void SolverClass::InitNSCBC(void)
{
	// Initialize surface area of each boundary face.
	su2double LocalSurfAreaBoundary[] = {zero, zero, zero, zero, zero, zero};
	
  // Start of the parallel region.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {

#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3), reduction(+: LocalSurfAreaBoundary[:6])
#endif
		for(int k=1; k<=mNElemPerRankK; ++k)
			for(int j=1; j<=mNElemPerRankJ; ++j)
				for(int i=1; i<=mNElemPerRankI; ++i)
					mElements[i][j][k]->ComputeSurfaceAreaBoundary(mInputParam, &mStandardHex,
																												 LocalSurfAreaBoundary);

	} // End of OpenMP parallel loop.

#ifdef HAVE_MPI
#ifdef HAVE_OPENMP
#pragma omp single nowait
#endif
	{
		su2double locBuf[] = { LocalSurfAreaBoundary[0], LocalSurfAreaBoundary[1],
													 LocalSurfAreaBoundary[2], LocalSurfAreaBoundary[3],
													 LocalSurfAreaBoundary[4], LocalSurfAreaBoundary[5] };

		// Accumulate all the surface areas across all processors per each face. 
		// Thus, the surface boundary area per face is computed and accessible on 
		// all processes.
		MPI_Allreduce(locBuf, mSurfAreaBoundary, 6, MPI_SU2DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
	}
#else 
	// This is a single rank simulation, thus use the local area as the entire surface area.
	for(unsigned short iBoundary=0; iBoundary<6; iBoundary++)
		mSurfAreaBoundary[iBoundary] = LocalSurfAreaBoundary[iBoundary];
#endif

	// Assign respective tuning parameters.
#ifdef HAVE_OPENMP	
#pragma omp for schedule(static), collapse(3)
#endif
	for(int k=1; k<=mNElemPerRankK; ++k)
		for(int j=1; j<=mNElemPerRankJ; ++j)
			for(int i=1; i<=mNElemPerRankI; ++i)
				mElements[i][j][k]->TuningParamNSCBC(mInputParam, &mStandardHex);
}

//-----------------------------------------------------------------------------

// Function, which computes the average data required for the NSCBC.
void SolverClass::AverageDataNSCBC(void)
{
	// Initialize the weighted local Mach number of each boundary face.
	su2double LocalMachWeighted[]   = {zero, zero, zero, zero, zero, zero};
	// Initialize the average Mach number on each boundary face.
	su2double AverageBoundaryMach[] = {zero, zero, zero, zero, zero, zero};

  // Easier storage of the padded number of integration points.
  const int nInt2DPad = mStandardHex.mNIntegration2DPad;

  // Determine the size of the indices of the work array.
  const int sizeWorkArrayInd1 = nVar;
  const int sizeWorkArrayInd2 = nInt2DPad;

	// Start of the parallel region.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    // Allocate the memory for the work arrays, which is needed to compute
    // the residuals.
    std::vector<su2double *> workArray(sizeWorkArrayInd1, NULL);
    for(int i=0; i<sizeWorkArrayInd1; ++i)
    {
      workArray[i] = (su2double *) AllocateMemory(sizeWorkArrayInd2*sizeof(su2double));
      if( !workArray[i] )
        Terminate("SolverClass::AverageDataNSCBC", __FILE__, __LINE__,
                  "Memory allocation failure for workArray");
    }


#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3), reduction(+: LocalMachWeighted[:6])
#endif
		for(int k=1; k<=mNElemPerRankK; ++k)
			for(int j=1; j<=mNElemPerRankJ; ++j)
				for(int i=1; i<=mNElemPerRankI; ++i)
					mElements[i][j][k]->AverageMachNSCBC(mInputParam, &mStandardHex, 
																							 workArray.data(), LocalMachWeighted);

    // Free the local memory.
    for(int i=0; i<sizeWorkArrayInd1; ++i)
      FreeMemory((void **) &workArray[i]);

  } // End of the OpenMP parallel region.

#ifdef HAVE_MPI
#ifdef HAVE_OPENMP
#pragma omp single nowait
#endif
	{
		su2double locBuf[] = { LocalMachWeighted[0], LocalMachWeighted[1],
													 LocalMachWeighted[2], LocalMachWeighted[3],
													 LocalMachWeighted[4], LocalMachWeighted[5] };

		// Accumulate all the weighted Mach number across all processors per each face. 
		// Thus, the weighted Mach number per face is computed and accessible on 
		// all processes.
		MPI_Allreduce(locBuf, AverageBoundaryMach, 6, MPI_SU2DOUBLE, MPI_SUM, MPI_COMM_WORLD); 
	}
#else 
	// This is a single rank simulation, thus use the local Mach as the total Mach.
	for(unsigned short iBoundary=0; iBoundary<6; iBoundary++)
		AverageBoundaryMach[iBoundary] = LocalMachWeighted[iBoundary];
#endif

	// Normalize by the true surface area per boundary face to obtain the actual average.
	for(unsigned short iBoundary=0; iBoundary<6; iBoundary++)
		if( FABS(mSurfAreaBoundary[iBoundary]) > epsSmall ) 
			AverageBoundaryMach[iBoundary] /= mSurfAreaBoundary[iBoundary];

	// Assign averaged Mach number on respective characteristic BCs.
#ifdef HAVE_OPENMP	
#pragma omp for schedule(static), collapse(3)
#endif
	for(int k=1; k<=mNElemPerRankK; ++k)
		for(int j=1; j<=mNElemPerRankJ; ++j)
			for(int i=1; i<=mNElemPerRankI; ++i)
				mElements[i][j][k]->SetAverageBoundaryMachNSCBC(AverageBoundaryMach);
}

//------------------------------------------------------------------------------

// Function, which determines the interpolation weights for the
// exchange location of the wall model, if needed.
void SolverClass::InterpolationWeightsExchangeLocation(void)
{
  // Return immediately if no wall model is used.
  if(mInputParam->mWallModelType == NO_WALL_MODEL) return;

  //----------------------------------------------------------------------------
  // Determine the 1D parametric coordinates used for the subdivision into
  // trilinear sub-elements. Usually an oversampling is needed compared to the
  // grid DOFs to account for curved elements. Several choices are possible,
  // but a practical choice are the integration points augmented with the
  // boundaries.
  //----------------------------------------------------------------------------

  // Determine the sampling coordinates.
  const int nSample1D = mStandardHex.mNIntegration1D + 2;
  std::vector<su2double> rSample1D(nSample1D);

  rSample1D[0] = -one;
  for(int i=0; i<mStandardHex.mNIntegration1D; ++i)
    rSample1D[i+1] = mStandardHex.mRIntegration1D[i];
  rSample1D[nSample1D-1] = one;

  // Allocate the memory for the 1D Lagrangian basis functions in the
  // sampling points. Initialize its values to zero, such that the
  // padded values are initialized.
  const int nDOFs1DGrid  = mInputParam->mNPolyGridDOFs +1;
  const int nSample1DPad = ((nSample1D+vecLen1D-1)/vecLen1D)*vecLen1D;
  const int sizeLag      = nDOFs1DGrid*nSample1DPad;

  su2double *lagrangeSample1D = (su2double *) AllocateMemory(sizeLag*sizeof(su2double));

  for(int i=0; i<sizeLag; ++i)
    lagrangeSample1D[i] = zero;

  // Compute the other Lagrangian basis functions.
  mStandardHex.LagrangianBasisFunctions(mStandardHex.mRDOFsGrid1D, rSample1D,
                                        lagrangeSample1D, NULL, NULL, NULL);

  // Determine the number of sample points in two and three space dimensions.
  const int nSample2D = nSample1D*nSample1D;
  const int nSample3D = nSample2D*nSample1D;

  // Allocate the memory to store the coordinates of the sampling points
  // of an element.
  std::vector<su2double *> coorSample(3, NULL);
  for(unsigned int m=0; m<coorSample.size(); ++m)
  {
    coorSample[m] = (su2double *) AllocateMemory(nSample3D*sizeof(su2double));
    if( !coorSample[m] )
      Terminate("SolverClass::InterpolationWeightsExchangeLocation", __FILE__, __LINE__,
                "Memory allocation failure for coorSample");
  }

  //----------------------------------------------------------------------------
  // Determine the local ADT of the volume sub-elements. It is assumed that the
  // donor information is stored in the local partitioning.
  //----------------------------------------------------------------------------

  // Determine the number of elements per rank and the number
  // of sub-elements per element.
  const int nElemPerRank = mNElemPerRankI*mNElemPerRankJ*mNElemPerRankK;
  const int nSubElem = (nSample1D-1)*(nSample1D-1)*(nSample1D-1);

  // Define the vectors needed to build the local ADT.
  std::vector<unsigned long>  parentElement;
  std::vector<unsigned short> subElementIDInParent;
  std::vector<unsigned short> VTK_TypeElem;
  std::vector<unsigned long>  elemConn;
  std::vector<su2double>      volCoor; 

  // Reserve the memory for these vectors.
  parentElement.reserve(nElemPerRank*nSubElem);
  subElementIDInParent.reserve(nElemPerRank*nSubElem);
  VTK_TypeElem.reserve(nElemPerRank*nSubElem);
  elemConn.reserve(8*nElemPerRank*nSubElem);
  volCoor.reserve(3*nElemPerRank*nSample3D);

  // Loop over the local owned elements.
  int elemID = 0;
  int offsetConn = 0;
  for(int k=1; k<=mNElemPerRankK; ++k)
  {
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      for(int i=1; i<=mNElemPerRankI; ++i, ++elemID)
      {
        // Determine the coordinates of the sample points.
        TensorProductSolAndGradVolume(nSample1D, 3, nDOFs1DGrid, lagrangeSample1D, NULL,
                                      mElements[i][j][k]->mCoorNodalGridDOFs.data(),
                                      coorSample.data(), NULL, NULL, NULL);

        // Store the coordinates of the sample DOFs in volCoor.
        for(int l=0; l<nSample3D; ++l)
        {
          volCoor.push_back(coorSample[0][l]);
          volCoor.push_back(coorSample[1][l]);
          volCoor.push_back(coorSample[2][l]);
        }

        // Loop over the sub-elements to create the connectivity of
        // these linear sub-elements, as the ADT is not able to handle
        // high order elements.
        int subElemID = 0;
        for(int kk=1; kk<nSample1D; ++kk)
        {
         const int kOff = offsetConn + (kk-1)*nSample2D;
          for(int jj=1; jj<nSample1D; ++jj)
          {
            const int jOff = kOff + (jj-1)*nSample1D;
            for(int ii=1; ii<nSample1D; ++ii, ++subElemID)
            {
              // Store the connectivity.
              const int n0 = jOff + ii-1;

              elemConn.push_back(n0);
              elemConn.push_back(n0 + 1);
              elemConn.push_back(n0 + 1 + nSample1D);
              elemConn.push_back(n0 + nSample1D);
              elemConn.push_back(n0 + nSample2D);
              elemConn.push_back(n0 + 1 + nSample2D);
              elemConn.push_back(n0 + 1 + nSample1D + nSample2D);
              elemConn.push_back(n0 + nSample1D + nSample2D);

              // Store the parent element ID and the sub-element ID
              // inside the parent element.
              parentElement.push_back(elemID);
              subElementIDInParent.push_back(subElemID);

              // Store the VTK element type, which is a hexahedron.
              VTK_TypeElem.push_back(HEXAHEDRON);
            }
          }
        }

        // Increment the offset for the connectivity.
        offsetConn += nSample3D;
      }
    }
  }

  // Build the local volume ADT.
  CADTElemClass localVolumeADT(3, volCoor, elemConn, VTK_TypeElem,
                               subElementIDInParent, parentElement);

  // Release the memory of the vectors used to build the ADT. To make sure
  // that all the memory is deleted, the swap function is used.
  std::vector<unsigned short>().swap(subElementIDInParent);
  std::vector<unsigned short>().swap(VTK_TypeElem);
  std::vector<unsigned long>().swap(parentElement);
  std::vector<unsigned long>().swap(elemConn);
  std::vector<su2double>().swap(volCoor);

  //----------------------------------------------------------------------------
  // Search for the donor elements for the corresponding exchange locations of
  // the integration points of faces that belong to viscous walls for which
  // a wall model must be used.
  //----------------------------------------------------------------------------

  // Parallel loop of the local elements. This is not the most efficient way
  // of doing this, but efficiency is not an issue here.
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(3), schedule(static)
#endif
  for(int k=1; k<=mNElemPerRankK; ++k)
    for(int j=1; j<=mNElemPerRankJ; ++j)
      for(int i=1; i<=mNElemPerRankI; ++i)
        mElements[i][j][k]->InterpolationWeightsExchangeLocation(mInputParam, &mStandardHex,
                                                                 mElements, mNElemPerRankI,
                                                                 mNElemPerRankJ, mNElemPerRankK,
                                                                 &localVolumeADT, rSample1D);

  // Release the memory of locally allocated arrays.
  FreeMemory((void **) &lagrangeSample1D);

  for(unsigned int m=0; m<coorSample.size(); ++m)
    FreeMemory((void **) &coorSample[m]);
}

//------------------------------------------------------------------------------

// Function, which determines the prescribed data in the integration
// points of the boundary faces.
void SolverClass::PrescribedDataIntegrationPoints(void)
{
  // Parallel loop of the local elements. This is not the most efficient way
  // of doing this, but efficiency is not an issue here.
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(3), schedule(static)
#endif
  for(int k=1; k<=mNElemPerRankK; ++k)
    for(int j=1; j<=mNElemPerRankJ; ++j)
      for(int i=1; i<=mNElemPerRankI; ++i)
        mElements[i][j][k]->PrescribedDataIntegrationPoints(mInputParam,
                                                            &mStandardHex,
                                                            mNonDimPrimVarFreeStream);
}

//------------------------------------------------------------------------------

// Function to read the grid.
void SolverClass::ReadGrid(void)
{
  // Determine the global number of elements in the three index direction.
  // It suffices to look at the i-index of the iMax subface, etc.
  const int nElemI = mInputParam->mSubfaces[1][0]->mElemIBeg;
  const int nElemJ = mInputParam->mSubfaces[3][0]->mElemJBeg;
  const int nElemK = mInputParam->mSubfaces[5][0]->mElemKBeg;

  // Easier storage of the polynomial degree of the grid DOFs.
  const int nPolyGrid = mInputParam->mNPolyGridDOFs;

  // Determine the number of grid points that must be present in the grid file.
  const int ilG = nElemI*nPolyGrid + 1;
  const int jlG = nElemJ*nPolyGrid + 1;
  const int klG = nElemK*nPolyGrid + 1;

  // Define the vector, which stores the coordinates as one block
  // for the current rank.
  std::vector<su2double> coorBuf;

#ifdef HAVE_MPI
  // Parallel mode. Open the grid file for reading and check if it went OK.
  MPI_File fh;
  if(MPI_File_open(MPI_COMM_WORLD, mInputParam->mGridFile.c_str(),
                   MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) != MPI_SUCCESS)
    TerminateAll("SolverClass::ReadGrid", __FILE__, __LINE__,
                 "Grid file could not be opened for reading");

  // Define the variable, which keeps track of the size of the header.
  MPI_Offset sizeHeader = 0;

  // Read the first integer of the grid file. All processors
  // read it, because a collective read is nonsense here.
  int recordSize;
  MPI_File_read(fh, &recordSize, 1, MPI_INT, MPI_STATUS_IGNORE);
  sizeHeader += sizeof(int);

  // Check if byte swapping must be applied.
  bool byteSwap = false;
  if((recordSize != sizeof(int)) && (recordSize != 3*sizeof(int)))
  {
    byteSwap = true;
    SwapBytes(&recordSize, sizeof(int), 1);
  }

  // If a multiblock format is specified, read the number of blocks.
  if(recordSize == sizeof(int))
  {
    // Read the next three integers and apply byte swapping, if needed.
    int recordNBlocks[3];
    MPI_File_read(fh, recordNBlocks, 3, MPI_INT, MPI_STATUS_IGNORE);
    if( byteSwap ) SwapBytes(recordNBlocks, sizeof(int), 3);
    sizeHeader += 3*sizeof(int);

    // Check if only one block is present in the grid.
    if(recordNBlocks[0] != 1)
      TerminateAll("SolverClass::ReadGrid", __FILE__, __LINE__,
                   "More than one block present in the grid.");

    // Set the record size for the next record.
    recordSize = recordNBlocks[2];
  }

  // Check that the next record consists of three integers.
  if(recordSize != 3*sizeof(int))
    TerminateAll("SolverClass::ReadGrid", __FILE__, __LINE__,
                 "Unexpected record size for the block size");

  // Read the block sizes, the closing integer of that record
  // and the leading integer of the record of the coordinates.
  int blockSizes[5];
  MPI_File_read(fh, blockSizes, 5, MPI_INT, MPI_STATUS_IGNORE);
  if( byteSwap ) SwapBytes(blockSizes, sizeof(int), 5);
  sizeHeader += 5*sizeof(int);

  // Check if the block size is correct.
  if((ilG != blockSizes[0]) || (jlG != blockSizes[1]) || (klG != blockSizes[2]))
    TerminateAll("SolverClass::ReadGrid", __FILE__, __LINE__,
                 "Block size does not correspond to the number of elements and polynomial degree");

  // Determine the precision of the coordinates.
  const int precisionCoor = blockSizes[4]/(3*ilG*jlG*klG);

  // Make a distinction between the precision and call ReadGridCoordinates
  // with the appropriate template variable to carry out the reading.
  switch( precisionCoor )
  {
    case 4:
    {
      ReadGridCoordinates<float>(fh, coorBuf, byteSwap, sizeHeader, ilG, jlG, klG);
      break;
    }

    case 8:
    {
      ReadGridCoordinates<double>(fh, coorBuf, byteSwap, sizeHeader, ilG, jlG, klG);
      break;
    }

    default:
      Terminate("SolverClass::ReadGrid", __FILE__, __LINE__,
                "Unknown precision encountered in the grid file");
  }

  // Close the file.
  MPI_File_close(&fh);  

#else

  // No MPI is used, so the reading of the file can be done sequentially.
  // Open the grid file for binary reading.
  FILE *gridFile = std::fopen(mInputParam->mGridFile.c_str(), "rb");
  if( !gridFile )
    Terminate("SolverClass::ReadGrid", __FILE__, __LINE__,
              "Grid file could not be opened for binary reading");

  // Read the first integer of the grid file.
  int recordSize;
  if(std::fread(&recordSize, sizeof(int), 1, gridFile) != 1)
    Terminate("SolverClass::ReadGrid", __FILE__, __LINE__,
              "Something wrong when reading the first record size");

  // Check if byte swapping must be applied.
  bool byteSwap = false;
  if((recordSize != sizeof(int)) && (recordSize != 3*sizeof(int)))
  {
    byteSwap = true;
    SwapBytes(&recordSize, sizeof(int), 1);
  }

  // If a multiblock format is specified, read the number of blocks.
  if(recordSize == sizeof(int))
  {
    // Read the next three integers and apply byte swapping, if needed.
    int recordNBlocks[3];
    if(std::fread(recordNBlocks, sizeof(int), 3, gridFile) != 3)
      Terminate("SolverClass::ReadGrid", __FILE__, __LINE__,
                "Something wrong when reading the number of blocks.");
    if( byteSwap ) SwapBytes(recordNBlocks, sizeof(int), 3);

    // Check if only one block is present in the grid.
    if(recordNBlocks[0] != 1)
      TerminateAll("SolverClass::ReadGrid", __FILE__, __LINE__,
                   "More than one block present in the grid.");

    // Set the record size for the next record.
    recordSize = recordNBlocks[2];
  }

  // Check that the next record consists of three integers.
  if(recordSize != 3*sizeof(int))
    TerminateAll("SolverClass::ReadGrid", __FILE__, __LINE__,
                 "Unexpected record size for the block size");

  // Read the block sizes, the closing integer of that record
  // and the leading integer of the record of the coordinates.
  int blockSizes[5];
  if(std::fread(blockSizes, sizeof(int), 5, gridFile) != 5)
    Terminate("SolverClass::ReadGrid", __FILE__, __LINE__,
              "Something wrong when reading the block sizes.");
  if( byteSwap ) SwapBytes(blockSizes, sizeof(int), 5);

  // Check if the block size is correct.
  if((ilG != blockSizes[0]) || (jlG != blockSizes[1]) || (klG != blockSizes[2]))
    TerminateAll("SolverClass::ReadGrid", __FILE__, __LINE__,
                 "Block size does not correspond to the number of elements and polynomial degree");

  // Determine the precision of the coordinates.
  const int precisionCoor = blockSizes[4]/(3*ilG*jlG*klG);

  // Make a distinction between the precision and call ReadGridCoordinates
  // with the appropriate template variable to carry out the reading.
  switch( precisionCoor )
  {
    case 4:
    {
      ReadGridCoordinates<float>(gridFile, coorBuf, byteSwap, ilG, jlG, klG);
      break;
    }

    case 8:
    {
      ReadGridCoordinates<double>(gridFile, coorBuf, byteSwap, ilG, jlG, klG);
      break;
    }

    default:
      Terminate("SolverClass::ReadGrid", __FILE__, __LINE__,
                "Unknown precision encountered in the grid file");
  }

  // Close the file again.
  std::fclose(gridFile);
#endif

  // Determine the number of grid DOFs per element.
  const int nGridDOFs = (nPolyGrid+1)*(nPolyGrid+1)*(nPolyGrid+1);

  // Determine the local number of grid points in each direction
  // present in coorBuf.
  const int il = mNElemPerRankI*nPolyGrid + 1;
  const int jl = mNElemPerRankJ*nPolyGrid + 1;
  const int kl = mNElemPerRankK*nPolyGrid + 1;  

  // Loop over the elements. This loop can be parallelized with OpenMP.
#ifdef HAVE_OPENMP
#pragma omp parallel for collapse(3), schedule(static)
#endif
  for(int k=1; k<=mNElemPerRankK; ++k)
  {
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      for(int i=1; i<=mNElemPerRankI; ++i)
      {
        // Set the pointers for the x-, y- and z-coordinates of the
        // lower left corner of the element.
        const su2double *x = coorBuf.data() + ((k-1)*il*jl + (j-1)*il + (i-1))*nPolyGrid;
        const su2double *y = x + il*jl*kl;
        const su2double *z = y + il*jl*kl;

        // Allocate the memory for the coordinates of the grid DOFs.
        mElements[i][j][k]->mCoorNodalGridDOFs.resize(3, NULL);
        for(int ii=0; ii<3; ++ii)
        {
          mElements[i][j][k]->mCoorNodalGridDOFs[ii] = (su2double *) AllocateMemory(nGridDOFs*sizeof(su2double));
          if( !mElements[i][j][k]->mCoorNodalGridDOFs[ii] )
            Terminate("SolverClass::ReadGrid", __FILE__, __LINE__,
                      "Memory allocation failure for mCoorNodalGridDOFs");
        }

        // Loop over the three directions of the element and copy the data.
        int ind = 0;
        for(int kk=0; kk<=nPolyGrid; ++kk)
        {
          for(int jj=0; jj<=nPolyGrid; ++jj)
          {
            for(int ii=0; ii<=nPolyGrid; ++ii, ++ind)
            {
              const int jnd = kk*il*jl + jj*il + ii;
              mElements[i][j][k]->mCoorNodalGridDOFs[0][ind] = x[jnd];
              mElements[i][j][k]->mCoorNodalGridDOFs[1][ind] = y[jnd];
              mElements[i][j][k]->mCoorNodalGridDOFs[2][ind] = z[jnd];
            }
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------

// Function to read the restart solution.
void SolverClass::ReadRestartSolution(void)
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
    precisionRestartFile = DeterminePrecisionRestartFile(mInputParam->mRestartFile.c_str());

#ifdef HAVE_MPI
  // Broadcast the precision to all ranks.
  MPI_Bcast(&precisionRestartFile, 1, MPI_UNSIGNED_SHORT, 0, MPI_COMM_WORLD);

  // Open the restart file for reading and check if it went OK.
  MPI_File fh;
  if(MPI_File_open(MPI_COMM_WORLD, mInputParam->mRestartFile.c_str(),
                   MPI_MODE_RDONLY, MPI_INFO_NULL, &fh) != MPI_SUCCESS)
    TerminateAll("SolverClass::ReadRestartSolution", __FILE__, __LINE__,
                 "Restart file could not be opened for reading");

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
    TerminateAll("SolverClass::ReadRestartSolution", __FILE__, __LINE__,
                 "Restart file is not an SU2 restart file");

  // Easier storage of the number of variables and DOFS in the restart file.
  const int nVarRestart  = header[1];
  const int nDOFsRestart = header[2];

  // Check if the number of DOFs is correct.
  if(nDOFsRestart != nDOFsTot)
    TerminateAll("SolverClass::ReadRestartSolution", __FILE__, __LINE__,
                 "Wrong number of DOFs in the restart file");

  // Determine the byte offset where the data for the solution
  // variables starts.
  MPI_Offset offset = 5*sizeof(int) + nVarRestart*CGNS_STRING_SIZE;

  // Make a distinction between the precision and call ReadSolutionFromRestartFile
  // with the appropriate template parameter to read the restart file accordingly.
  switch( precisionRestartFile )
  {
    case 4:
    {
      ReadSolutionFromRestartFile<float>(fh, nVarRestart, byteSwap, offset);
      break;
    }

    case 8:
    {
      ReadSolutionFromRestartFile<double>(fh, nVarRestart, byteSwap, offset);
      break;
    }

    default:
      TerminateAll("SolverClass::ReadRestartSolution", __FILE__, __LINE__,
                   "Unknown precision encountered in the restart file");
  }

  // Check if the average solution is also read. In that case the number of
  // time steps for which the averaging has been applied must also be read
  // as well as the time averaged force coefficients.
  if( mInputParam->mContinueAveragingFromRestart )
  {
    // Determine the new offset and reset the file view.
    offset += nDOFsTot*nVarRestart*precisionRestartFile;

    char datarep[] = "native";
    MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, datarep, MPI_INFO_NULL);

    // Read the number of time steps for which the averaging has been applied.
    // All processors read it, because a collective read is nonsense here.
    MPI_File_read_at(fh, offset, &mNTimeStepAverage, 1, MPI_INT, MPI_STATUS_IGNORE);

    // Apply byte swapping, if needed.
    if( byteSwap ) SwapBytes(&mNTimeStepAverage, sizeof(int), 1);

    // Read the restart meta data, which contains the time averaged force coefficients.
    // All processors read it, because a collective read is nonsense here.
    offset += sizeof(int);
    su2double restartMetaData[8];
    MPI_File_read_at(fh, offset, restartMetaData, 8, MPI_SU2DOUBLE,  MPI_STATUS_IGNORE);

    if( byteSwap ) SwapBytes(restartMetaData, sizeof(su2double), 8);

    // Copy the data into the member variables to store the averaged force coefficients.
    mCfxAve    = restartMetaData[0];
    mCfyAve    = restartMetaData[1];
    mCfzAve    = restartMetaData[2];
    mCfxVisAve = restartMetaData[3];
    mCfyVisAve = restartMetaData[4];
    mCfzVisAve = restartMetaData[5];
  }

  // Close the file.
  MPI_File_close(&fh);

#else
  // No MPI is used, so the reading of the file can be done sequentially.
  // Open the restart file for binary reading.
  FILE *restartFile = std::fopen(mInputParam->mRestartFile.c_str(), "rb");
  if( !restartFile )
    Terminate("SolverClass::ReadRestartSolution", __FILE__, __LINE__,
              "Restart file could not be opened for binary reading");

  // Read the first 5 integer variables of the solution file.
  int header[5];
  if(std::fread(header, sizeof(int), 5, restartFile) != 5)
    Terminate("SolverClass::ReadRestartSolution", __FILE__, __LINE__,
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
    Terminate("SolverClass::ReadRestartSolution", __FILE__, __LINE__,
              "Restart file is not an SU2 restart file");

  // Easier storage of the number of variables and DOFS in the restart file.
  const int nVarRestart  = header[1];
  const int nDOFsRestart = header[2];

  // Check if the number of DOFs is correct.
  if(nDOFsRestart != nDOFsTot)
    Terminate("SolverClass::ReadRestartSolution", __FILE__, __LINE__,
              "Wrong number of DOFs in the restart file");

  // Read the names of the solution variables.
  for(int i=0; i<nVarRestart; ++i)
  {
    char varName[CGNS_STRING_SIZE];
    if(std::fread(varName, sizeof(char), CGNS_STRING_SIZE, restartFile) != CGNS_STRING_SIZE)
      Terminate("SolverClass::ReadRestartSolution", __FILE__, __LINE__,
                "Variable name could not be read");
  }

  // Make a distinction between the precision and call ReadSolutionFromRestartFile
  // with the appropriate template parameter to read the restart file accordingly.
  switch( precisionRestartFile )
  {
    case 4:
    {
      ReadSolutionFromRestartFile<float>(restartFile, nVarRestart, byteSwap);
      break;
    }

    case 8:
    {
      ReadSolutionFromRestartFile<double>(restartFile, nVarRestart, byteSwap);
      break;
    }

    default:
      Terminate("SolverClass::ReadRestartSolution", __FILE__, __LINE__,
                "Unknown precision encountered in the restart file");
  }

  // Check if the average solution is also read. In that case the number of
  // time steps for which the averaging has been applied must also be read
  // as well as the time averaged force coefficients.
  // Apply byte swapping, if needed.
  if( mInputParam->mContinueAveragingFromRestart )
  {
    if(std::fread(&mNTimeStepAverage, sizeof(int), 1, restartFile) != 1)
      Terminate("SolverClass::ReadRestartSolution", __FILE__, __LINE__,
                "Number of time steps for averaging could not be read");

    if( byteSwap ) SwapBytes(&mNTimeStepAverage, sizeof(int), 1);

    su2double restartMetaData[8];
    if(std::fread(restartMetaData, sizeof(su2double), 8, restartFile) != 8)
      Terminate("SolverClass::ReadRestartSolution", __FILE__, __LINE__,
                "Restart meta data could not be read");

    if( byteSwap ) SwapBytes(restartMetaData, sizeof(su2double), 8);

    // Copy the data into the member variables to store the averaged force coefficients.
    mCfxAve    = restartMetaData[0];
    mCfyAve    = restartMetaData[1];
    mCfzAve    = restartMetaData[2];
    mCfxVisAve = restartMetaData[3];
    mCfyVisAve = restartMetaData[4];
    mCfzVisAve = restartMetaData[5];
  }

  // Close the file again.
  std::fclose(restartFile);
#endif

  // Start of the parallel region.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    // Allocate the memory to store the nodal solution.
    std::vector<su2double *> solNode(nVar, NULL);
    for(int i=0; i<nVar; ++i)
    {
      solNode[i] = (su2double *) AllocateMemory(nDOFs*sizeof(su2double));
      if( !solNode[i] )
        Terminate("SolverClass::ReadRestartSolution", __FILE__, __LINE__,
                  "Memory allocation failure for solNode");
    }

    // Loop over the owned elements.
#ifdef HAVE_OPENMP
#pragma omp for collapse(3), schedule(static)
#endif
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        for(int i=1; i<=mNElemPerRankI; ++i)
        {
          // Set the pointer to the solution vector.
          su2double **sol = mElements[i][j][k]->mSol.data();

          // Make a distinction between the working variables.
          switch( mInputParam->mFEMVariables )
          {
            case CONSERVATIVE_VARIABLES:
            {
              // Conservative variables are used as working variables.
              // Non-dimensionalize them.
              const su2double rhoRefInv = one/rhoRef;
              const su2double momRefInv = one/(rhoRef*uRef);
              const su2double pRefInv   = one/pRef;

              // Loop over the DOFs.
#pragma omp simd
              for(int DOF=0; DOF<nDOFs; ++DOF)
              {
                solNode[0][DOF] = sol[0][DOF]*rhoRefInv;
                solNode[1][DOF] = sol[1][DOF]*momRefInv;
                solNode[2][DOF] = sol[2][DOF]*momRefInv;
                solNode[3][DOF] = sol[3][DOF]*momRefInv;
                solNode[4][DOF] = sol[4][DOF]*pRefInv;
              }

              break;
            }

            //------------------------------------------------------------------------

            case ENTROPY_VARIABLES:
            {
              // The entropy variables are used as working variables.
              // Compute the non-dimensional primitive variables from the conservative
              // ones and compute the entropy variables.
              const su2double rhoRefInv = one/rhoRef;
              const su2double uRefInv   = one/uRef;
              const su2double pRefInv   = one/pRef;

              const su2double gm1   = GamConstant - one;
              const su2double ovgm1 = one/gm1;

              // Loop over the DOFs.
#pragma omp simd
              for(int DOF=0; DOF<nDOFs; ++DOF)
              {
                su2double rho    = sol[0][DOF];
                su2double rhoInv = one/rho;
                su2double u      = rhoInv*sol[1][DOF];
                su2double v      = rhoInv*sol[2][DOF];
                su2double w      = rhoInv*sol[3][DOF];
                su2double p      = gm1*(sol[4][DOF] - half*(u*sol[1][DOF] + v*sol[2][DOF] + w*sol[3][DOF]));

                rho *= rhoRefInv;
                u   *= uRefInv;
                v   *= uRefInv;
                w   *= uRefInv;
                p   *= pRefInv;

                const su2double s    = LOG(p/POW(rho,GamConstant));
                const su2double pInv = one/p;

                solNode[0][DOF] =  (GamConstant-s)*ovgm1 - half*rho*(u*u + v*v + w*w)*pInv;
                solNode[1][DOF] =  rho*u*pInv;
                solNode[2][DOF] =  rho*v*pInv;
                solNode[3][DOF] =  rho*w*pInv;
                solNode[4][DOF] = -rho*pInv;
              }

              break;
            }

            //------------------------------------------------------------------------

            default:
            {
               // This is just to avoid a compiler warning.
               break;
            }
          }

          // Convert the nodal form to the modal form.
          const int nDOFs1D = mStandardHex.mNDOFs1D;
          TensorProductSolAndGradVolume(nDOFs1D, nVar, nDOFs1D,
                                        mStandardHex.mVandermonde1DInverse,
                                        mStandardHex.mVandermonde1DInverse,
                                        solNode.data(), sol, NULL, NULL, NULL);
        }
      }
    }

    // Release the memory of the nodal solution again.
    for(int i=0; i<nVar; ++i)
      FreeMemory((void **) &solNode[i]);

  } // End of the OpenMP parallel region.
}

//------------------------------------------------------------------------------

// Function, which computes the spatial residual in the DOFs of the owned elements.
void SolverClass::Residual(su2double  *monitoringData,
                           const bool ComputeMonitoringData)
{
  // Define the variables for the monitoring.
  su2double Mach2Max = zero, EddyVisMax = -one;
  su2double forceCoef[] = {zero, zero, zero, zero, zero, zero};

  // Easier storage of the padded number of DOFs and integration points.
  const int nDOFsPad = mStandardHex.mNDOFsPad;
  const int nIntPad  = mStandardHex.mNIntegrationPad;

  // Determine the size needed for the wall model, if needed.
  int nSizeWallModel = 0;
  if(mInputParam->mWallModelType == EQUILIBRIUM_WALL_MODEL)
    nSizeWallModel = mInputParam->mWallModel->GetNWallModelPoints();

  // Determine the size needed for the preconditioned CG algorithm,
  // which is for vectors for the actual preconditioned CG algorithm
  // and one vector to store interpolated data in the integration points,
  // needed for the matrix free evaluation.
  const int sizeCG = 4*nVar + nVar;

  // Determine the size of the indices of the work array. 
	// The size is 15*nVar+1, because in the numerical surfaces, the indices are:
	// [00-03]: left -state solution and gradient,
	// [04-07]: right-state solution and gradient,
	// [08   ]: flux total,
	// [09-11]: left -state symmetrizing fluxes,
	// [12-14]: right-state symmetrizing fluxes,
	// [15   ]: eddy viscosity.
  const int sizeWorkArrayInd1 = std::max(15*nVar+1, sizeCG);
  const int sizeWorkArrayInd2 = std::max(std::max(nIntPad, nDOFsPad), nSizeWallModel);

  // Start of the parallel region.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    // Allocate the memory for the work arrays, which is needed to compute
    // the residuals.
    std::vector<su2double *> workArray(sizeWorkArrayInd1, NULL);
    for(int i=0; i<sizeWorkArrayInd1; ++i)
    {
      workArray[i] = (su2double *) AllocateMemory(sizeWorkArrayInd2*sizeof(su2double));
      if( !workArray[i] )
        Terminate("SolverClass::Residual", __FILE__, __LINE__,
                  "Memory allocation failure for workArray");
    }

    // Start the communication of the solution if MPI is used. One OpenMP thread
    // handles this per direction if there is something to be communicated.
    // No need for the other threads to wait.
#ifdef HAVE_MPI
    if(mNRanksI > 1)
    {
#ifdef HAVE_OPENMP
#pragma omp single nowait
#endif
      StartCommunicationSolutionIDir();
    }

    if(mNRanksJ > 1)
    {
#ifdef HAVE_OPENMP
#pragma omp single nowait
#endif
      StartCommunicationSolutionJDir();
    }

    if(mNRanksK > 1)
    {
#ifdef HAVE_OPENMP
#pragma omp single nowait
#endif
      StartCommunicationSolutionKDir();
    }
#endif  // HAVE_MPI

    // Compute the contributions from the volume integral to the residual.
    // This also serves as an initialization of the residual.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3), nowait, \
                reduction(max: Mach2Max), \
                reduction(max: EddyVisMax)
#endif
    for(int k=1; k<=mNElemPerRankK; ++k)
      for(int j=1; j<=mNElemPerRankJ; ++j)
        for(int i=1; i<=mNElemPerRankI; ++i)
          mElements[i][j][k]->VolumeResidual(mInputParam, &mStandardHex,
                                             workArray.data(), ComputeMonitoringData,
                                             Mach2Max, EddyVisMax);

    // Complete the communication of the solution in i-direction if MPI is used.
    // One thread handles this, while the others wait.
#ifdef HAVE_MPI
    if(mNRanksI > 1)
    {
#ifdef HAVE_OPENMP
#pragma omp single
#endif
      CompleteCommunicationSolutionIDir();
    }
#endif  // HAVE_MPI

    // Compute the contributions from the faces in i-direction to the residual.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3), nowait, \
                reduction(max: EddyVisMax), \
                reduction(+: forceCoef[:6])
#endif
    for(int k=1; k<=mNElemPerRankK; ++k)
      for(int j=1; j<=mNElemPerRankJ; ++j)
        for(int i=1; i<=mNElemPerRankI; ++i)
          mElements[i][j][k]->IFaceResiduals(mInputParam, &mStandardHex,
                                             mElements, workArray.data(),
                                             ComputeMonitoringData, EddyVisMax,
                                             forceCoef);

    // Start the communication of the face residuals in i-direction if
    // MPI is used. One OpenMP thread handles this. However, the other threads
    // should have reached this point before the message can be sent. No need
    // for the other threads to wait when the actual MPI sending takes place.
#ifdef HAVE_MPI
    if(mNRanksI > 1)
    {
#ifdef HAVE_OPENMP
#pragma omp barrier
#pragma omp single nowait
#endif
      StartCommunicationIFaceResiduals();
    }
#endif  // HAVE_MPI

    // Complete the communication of the solution in j-direction if MPI is used.
    // One thread handles this, while the others wait.
#ifdef HAVE_MPI
    if(mNRanksJ > 1)
    {
#ifdef HAVE_OPENMP
#pragma omp single
#endif
      CompleteCommunicationSolutionJDir();
    }
#endif  // HAVE_MPI

    // Compute the contributions from the faces in j-direction to the residual.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3), nowait, \
                reduction(max: EddyVisMax), \
                reduction(+: forceCoef[:6])
#endif
    for(int k=1; k<=mNElemPerRankK; ++k)
      for(int j=1; j<=mNElemPerRankJ; ++j)
        for(int i=1; i<=mNElemPerRankI; ++i)
          mElements[i][j][k]->JFaceResiduals(mInputParam, &mStandardHex,
                                             mElements, workArray.data(),
                                             ComputeMonitoringData, EddyVisMax,
                                             forceCoef);

    // Start the communication of the face residuals in j-direction if
    // MPI is used. One OpenMP thread handles this. However, the other threads
    // should have reached this point before the message can be sent. No need
    // for the other threads to wait when the actual MPI sending takes place.
#ifdef HAVE_MPI
    if(mNRanksJ > 1)
    {
#ifdef HAVE_OPENMP
#pragma omp barrier
#pragma omp single nowait
#endif
      StartCommunicationJFaceResiduals();
    }
#endif  // HAVE_MPI

    // Complete the communication of the solution in k-direction if MPI is used.
    // One thread handles this, while the others wait.
#ifdef HAVE_MPI
    if(mNRanksK > 1)
    {
#ifdef HAVE_OPENMP
#pragma omp single
#endif
      CompleteCommunicationSolutionKDir();
    }
#endif  // HAVE_MPI

    // Compute the contributions from the faces in k-direction to the residual.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3), nowait, \
                reduction(max: EddyVisMax), \
                reduction(+: forceCoef[:6])
#endif
    for(int k=1; k<=mNElemPerRankK; ++k)
      for(int j=1; j<=mNElemPerRankJ; ++j)
        for(int i=1; i<=mNElemPerRankI; ++i)
          mElements[i][j][k]->KFaceResiduals(mInputParam, &mStandardHex,
                                             mElements, workArray.data(),
                                             ComputeMonitoringData, EddyVisMax,
                                             forceCoef);

    // Start the communication of the face residuals in k-direction if
    // MPI is used. One OpenMP thread handles this. However, the other threads
    // should have reached this point before the message can be sent. No need
    // for the other threads to wait when the actual MPI sending takes place.
#ifdef HAVE_MPI
    if(mNRanksK > 1)
    {
#ifdef HAVE_OPENMP
#pragma omp barrier
#pragma omp single nowait
#endif
      StartCommunicationKFaceResiduals();
    }
#endif  // HAVE_MPI

    // Complete the communications of the face residuals in i-direction,
    // if MPI is used.
#ifdef HAVE_MPI
    if(mNRanksI > 1)
    {
#ifdef HAVE_OPENMP
#pragma omp single
#endif
      CompleteCommunicationIFaceResiduals();
    }
    else if(nRanks == 1)
    {
#ifdef HAVE_OPENMP
#pragma omp barrier
#endif
    }
#elif HAVE_OPENMP
#pragma omp barrier
#endif
    
    // Loop over the owned elements to add the face residuals in
    // i-direction, which were computed in the neighboring element.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3), nowait
#endif
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        for(int i=1; i<=mNElemPerRankI; ++i)
        {
          if( mElements[i+1][j][k] )
          {
            for(int l=0; l<nVar; ++l)
            {
              const su2double *resIMin = mElements[i+1][j][k]->mResIMin[l];
              su2double *res           = mElements[i][j][k]->mRes[l];
#pragma omp simd
              for(int m=0; m<nDOFsPad; ++m)
                res[m] += resIMin[m];
            }
          }
        }
      }
    }

    // Complete the communications of the face residuals in j-direction,
    // if MPI is used.
#ifdef HAVE_MPI
    if(mNRanksJ > 1)
    {
#ifdef HAVE_OPENMP
#pragma omp single
#endif
      CompleteCommunicationJFaceResiduals();
    }
    else
    {
#ifdef HAVE_OPENMP
#pragma omp barrier
#endif
    }
#elif HAVE_OPENMP
#pragma omp barrier
#endif

    // Loop over the owned elements to add the face residuals in
    // j-direction, which were computed in the neighboring element.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3), nowait
#endif
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        for(int i=1; i<=mNElemPerRankI; ++i)
        {
          if( mElements[i][j+1][k] )
          {
            for(int l=0; l<nVar; ++l)
            {
              const su2double *resJMin = mElements[i][j+1][k]->mResJMin[l];
              su2double *res           = mElements[i][j][k]->mRes[l];
#pragma omp simd
              for(int m=0; m<nDOFsPad; ++m)
                res[m] += resJMin[m];
            }
          }
        }
      }
    }

    // Complete the communications of the face residuals in k-direction,
    // if MPI is used.
#ifdef HAVE_MPI
    if(mNRanksK > 1)
    {
#ifdef HAVE_OPENMP
#pragma omp single
#endif
      CompleteCommunicationKFaceResiduals();
    }
    else
    {
#ifdef HAVE_OPENMP
#pragma omp barrier
#endif
    }
#elif HAVE_OPENMP
#pragma omp barrier
#endif

    // Loop over the owned elements to add the face residuals in
    // k-direction, which were computed in the neighboring element.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3), nowait
#endif
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        for(int i=1; i<=mNElemPerRankI; ++i)
        {
          if( mElements[i][j][k+1] )
          {
            for(int l=0; l<nVar; ++l)
            {
              const su2double *resKMin = mElements[i][j][k+1]->mResKMin[l];
              su2double *res           = mElements[i][j][k]->mRes[l];
#pragma omp simd
              for(int m=0; m<nDOFsPad; ++m)
                res[m] += resKMin[m];
            }
          }
        }
      }
    }

    // Compute the final residual by multiplying the current residual with the
    // inverse of the mass matrix (or the generalized mass matrix in case of
    // entropy variables). The nowait is there, because the parallel region
    // ends after this loop.
#ifdef HAVE_OPENMP
#pragma omp for schedule(static), collapse(3) nowait
#endif
    for(int k=1; k<=mNElemPerRankK; ++k)
      for(int j=1; j<=mNElemPerRankJ; ++j)
        for(int i=1; i<=mNElemPerRankI; ++i)
          mElements[i][j][k]->MultiplyResInverseMassMatrix(mInputParam, &mStandardHex,
                                                           workArray.data());

    // Free the local memory.
    for(int i=0; i<sizeWorkArrayInd1; ++i)
      FreeMemory((void **) &workArray[i]);

  } // End of the OpenMP parallel region.

  // Store the data to be monitored in monitoringData.
  monitoringData[0] = forceCoef[0];
  monitoringData[1] = forceCoef[1];
  monitoringData[2] = forceCoef[2];
  monitoringData[3] = forceCoef[3];
  monitoringData[4] = forceCoef[4];
  monitoringData[5] = forceCoef[5];
  monitoringData[6] = Mach2Max;
  monitoringData[7] = EddyVisMax;

  // For an MPI computation the global values must be obtained. The result
  // only needs to be known on the root rank, hence no need for Allreduce.
#ifdef HAVE_MPI
  if( ComputeMonitoringData )
  {
    MPI_Reduce(forceCoef, monitoringData, 6, MPI_SU2DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    su2double locBuf[] = {Mach2Max, EddyVisMax};
    MPI_Reduce(locBuf, &monitoringData[6], 2, MPI_SU2DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  }
#endif

  // Non-dimensionalize the force coefficients appropriately.
  const su2double factCoef = two/(GamConstant*mInputParam->mMach
                           *      mInputParam->mMach*mInputParam->mReferenceArea);
  monitoringData[0] *= factCoef;
  monitoringData[1] *= factCoef;
  monitoringData[2] *= factCoef;
  monitoringData[3] *= factCoef;
  monitoringData[4] *= factCoef;
  monitoringData[5] *= factCoef;

  //----------------------------------------------------------------------------
  // Part meant for debugging. The residuals of all DOFs are written to file,
  // such that they can be compared when modifications have been made.
  //----------------------------------------------------------------------------

  /*
  WriteRestartSolution(0, true);
  TerminateAll("SolverClass::Residual", __FILE__, __LINE__,
               "Residuals written to the restart file for testing purposes");
  */
}

//------------------------------------------------------------------------------

// Function, which updates the data for the averages.
void SolverClass::UpdateAverageData(const su2double *monitoringData)
{
  // Store the number of DOFs and padded DOFs per element a bit easier.
  const int nDOFs    = mStandardHex.mNDOFs;
  const int nDOFsPad = mStandardHex.mNDOFsPad;

  // Start of the parallel region.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    // Allocate the memory to store the dimensional primitive variables
    // and the eddy viscosity, if needed.
    std::vector<su2double *> primVar(nVar+1);
    for(int m=0; m<=nVar; ++m)
    {
      primVar[m] = NULL;
      primVar[m] = (su2double *) AllocateMemory(nDOFsPad*sizeof(su2double));
      if( !primVar[m] )
        Terminate("SolverClass::UpdateAverageData", __FILE__, __LINE__,
                  "Memory allocation failure for primVar[m].");
    }

    // Allocate the memory to store the gradients of the velocities, if needed.
    std::vector<su2double *> gradVel(3*nVar, NULL);
    if(mInputParam->mSGSModelType != NO_SGS_MODEL)
    {
      for(int m=0; m<(3*nVar); ++m)
      {
        gradVel[m] = (su2double *) AllocateMemory(nDOFs*sizeof(su2double));
        if( !gradVel[m] )
          Terminate("SolverClass::UpdateAverageData", __FILE__,
                    __LINE__, "Memory allocation failure for gradVel[m].");
      }
    }

    // Loop over the number of elements and call the member function
    // of the elements to do the actual work.
#ifdef HAVE_OPENMP
#pragma omp for collapse(3), schedule(static), nowait
#endif
    for(int k=1; k<=mNElemPerRankK; ++k)
      for(int j=1; j<=mNElemPerRankJ; ++j)
        for(int i=1; i<=mNElemPerRankI; ++i)
          mElements[i][j][k]->UpdateAverageData(mInputParam, &mStandardHex,
                                                mNTimeStepAverage, primVar.data(),
                                                gradVel.data());

    // Release the memory of the primitive variables and velocity gradients again.
    for(unsigned long i=0; i<primVar.size(); ++i)
      FreeMemory((void **) &primVar[i]);

    for(unsigned long i=0; i<gradVel.size(); ++i)
      FreeMemory((void **) &gradVel[i]);


  } // End of the OpenMP parallel region.

  // Update the average values of the force coefficients.
  const su2double factNew = one/mNTimeStepAverage;
  const su2double factOld = one - factNew;

  mCfxAve = factOld*mCfxAve + factNew*monitoringData[0];
  mCfyAve = factOld*mCfyAve + factNew*monitoringData[1];
  mCfzAve = factOld*mCfzAve + factNew*monitoringData[2];

  mCfxVisAve = factOld*mCfxVisAve + factNew*monitoringData[3];
  mCfyVisAve = factOld*mCfyVisAve + factNew*monitoringData[4];
  mCfzVisAve = factOld*mCfzVisAve + factNew*monitoringData[5];
}

//------------------------------------------------------------------------------

// Function to write the binary paraview solution file.
void SolverClass::WriteParaviewSolution(const int timeStep)
{
  // Define the maximum string length for the writing.
  const int MAX_STRING_LENGTH = 255;

  // Check for big endian. If not, the bytes must be swapped when
  // writing a binary VTK file.
  bool BigEndian;
  union {int i; char c[4];} val;
  val.i = 0x76543210;
  if (val.c[0] == 0x10) BigEndian = false;
  else                  BigEndian = true;

  // Determine the name of the file to be written.
  std::ostringstream fileName;
  fileName << mInputParam->mParaviewFile << "_TimeStep_" << timeStep << ".vtk";

  // Easier storage of the total number of elements in the three index directions.
  const int nElemI = mInputParam->mSubfaces[1][0]->mElemIBeg;
  const int nElemJ = mInputParam->mSubfaces[3][0]->mElemJBeg;
  const int nElemK = mInputParam->mSubfaces[5][0]->mElemKBeg;

  // Easier storage of the number of DOFs in one direction, the number of DOFs
  // in 2D (i.e. a face) and the total number of DOFs.
  const int nDOFs1D = mStandardHex.mNDOFs1D;
  const int nDOFs2D = nDOFs1D*nDOFs1D;
  const int nDOFs   = mStandardHex.mNDOFs;

  // Determine the total number of DOFs in the simulation and
  // the number of DOFs stored on one MPI rank.
  const int nElemPerRankIJ  = mNElemPerRankI*mNElemPerRankJ;
  const int nElemPerRankIJK = mNElemPerRankK*nElemPerRankIJ;

  const int nDOFsTot = nElemI*nElemJ*nElemK*nDOFs;
  const int nDOFsLoc = nElemPerRankIJK*nDOFs;

  // Allocate the memory for the buffers to store the coordinates and
  // the connectivities. Note that single precision is used for the
  // visualization. The element types to be written are all the same,
  // namely a standard hexahedron, VTK number 12. Hence the data for
  // the vector elemTypeBuf can be set directly.
  const int nSubElem = (nDOFs1D-1)*(nDOFs1D-1)*(nDOFs1D-1);
  std::vector<float> coorBuf(3*nDOFsLoc);
  std::vector<int>   connBuf(9*nElemPerRankIJK*nSubElem);
  std::vector<int>   elemTypeBuf(nElemPerRankIJK*nSubElem, 12);

  // Loop over the owned elements. A 1D loop is used instead of a 3D loop,
  // such that it is easier to retrieve the position in the arrays
  // used for the writing of the actual data.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
  for(int elem=0; elem<nElemPerRankIJK; ++elem)
  {
    // Retrieve the i-, j- and k-index of the element. Use a zero based
    // index here, because that is easier for the retrieval.
    int k   = elem/nElemPerRankIJ;
    int rem = elem - k*nElemPerRankIJ;
    int j   = rem/mNElemPerRankI;
    int i   = rem - j*mNElemPerRankI;

    // Determine the global i-, j- and k-index of this element.
    const int elemIndK = mMyRankK*mNElemPerRankK + k;
    const int elemIndJ = mMyRankJ*mNElemPerRankJ + j;
    const int elemIndI = mMyRankI*mNElemPerRankI + i;

    // Increment i, j and k, because the owned elements start at 1.
    ++i; ++j; ++k;

    // Determine the location in the buffers where the data of
    // this element should be stored.
    float *coorElem = coorBuf.data() + 3*elem*nDOFs;
    int   *connElem = connBuf.data() + 9*elem*nSubElem;

    // Loop over the DOFs of the element to copy the coordinates to coorElem.
    int ind = 0;
    for(int ii=0; ii<nDOFs; ++ii)
    {
      coorElem[ind++] = (float) mElements[i][j][k]->mCoorNodalSolDOFs[0][ii];
      coorElem[ind++] = (float) mElements[i][j][k]->mCoorNodalSolDOFs[1][ii];
      coorElem[ind++] = (float) mElements[i][j][k]->mCoorNodalSolDOFs[2][ii];
    }

    // The current element has global indices (elemIndI, elemIndJ, elemIndK).
    // Convert this to a 1D index and determine the global offset where the
    // connectivities of this element start.
    const long offsetConn = (elemIndK*nElemI*nElemJ + elemIndJ*nElemI + elemIndI)*nDOFs;

    // Loop over the subelements of the current element and determine the
    // connectivities. Note that the first entry in the connectivity array
    // is the number of points per element, which is 8.
    ind = 0;
    for(int kk=1; kk<nDOFs1D; ++kk)
    {
      const int kOff = offsetConn + (kk-1)*nDOFs2D;
      for(int jj=1; jj<nDOFs1D; ++jj)
      {
        const int jOff = kOff + (jj-1)*nDOFs1D;
        for(int ii=1; ii<nDOFs1D; ++ii)
        {
          const int n0 = jOff + ii-1;
          connElem[ind++] = 8;
          connElem[ind++] = n0;
          connElem[ind++] = n0 + 1;
          connElem[ind++] = n0 + 1 + nDOFs1D;
          connElem[ind++] = n0 + nDOFs1D;
          connElem[ind++] = n0 + nDOFs2D;
          connElem[ind++] = n0 + 1 + nDOFs2D;
          connElem[ind++] = n0 + 1 + nDOFs1D + nDOFs2D;
          connElem[ind++] = n0 + nDOFs1D + nDOFs2D;
        }
      }
    }
  }

  // Determine the visualization data to be written to the paraview file
  std::vector<std::string> visNames;
  std::vector<float> visBuf;
  DetermineVisualizationData(visNames, visBuf);

  // Convert the buffers to big endian, if needed.
  if( !BigEndian )
  {
    SwapBytes(coorBuf.data(),     sizeof(float), coorBuf.size());
    SwapBytes(connBuf.data(),     sizeof(int),   connBuf.size());
    SwapBytes(elemTypeBuf.data(), sizeof(int),   elemTypeBuf.size());
    SwapBytes(visBuf.data(),      sizeof(float), visBuf.size());
  }

#ifdef HAVE_MPI
  // Wipe out any previous output file. The reason for doing this is that there
  // is no guarantee that an old file will be clobbered by the MPI write routines
  if(rank == 0) MPI_File_delete(fileName.str().c_str(), MPI_INFO_NULL);
  MPI_Barrier(MPI_COMM_WORLD);

  // Open the visualization file for writing and check if it went OK.
  MPI_File fh;
  if(MPI_File_open(MPI_COMM_WORLD, fileName.str().c_str(),
                   MPI_MODE_WRONLY | MPI_MODE_CREATE,
                   MPI_INFO_NULL, &fh) != MPI_SUCCESS)
    TerminateAll("SolverClass::WriteParaviewSolution", __FILE__, __LINE__,
                 "Paraview visualization file could not be opened for writing");

  // Initialize the file offset to zero.
  MPI_Offset offset = 0;

  // Write the header of the visualization file.
  char str_buf[MAX_STRING_LENGTH];
  std::strcpy(str_buf, "# vtk DataFile Version 3.0\n");
  if(rank == 0)
    MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  offset += strlen(str_buf)*sizeof(char);

  std::strcpy(str_buf, "vtk output\n");
  if(rank == 0)
    MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  offset += strlen(str_buf)*sizeof(char);

  std::strcpy(str_buf, "BINARY\n");
  if(rank == 0)
    MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  offset += strlen(str_buf)*sizeof(char);

  std::strcpy(str_buf, "DATASET UNSTRUCTURED_GRID\n");
  if(rank == 0)
    MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  offset += strlen(str_buf)*sizeof(char);

  std::sprintf(str_buf, "POINTS %i float\n", nDOFsTot);
  if(rank == 0)
    MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  offset += strlen(str_buf)*sizeof(char);

  // Define the variables to create the subarray data type for the actual
  // writing of the coordinates.
  int gsizes[]   = {(int) (3*nDOFs), nElemI, nElemJ, nElemK};
  int lsizes[]   = {(int) (3*nDOFs), mNElemPerRankI, mNElemPerRankJ, mNElemPerRankK};
  int startInd[] = {0, mMyRankI*mNElemPerRankI, mMyRankJ*mNElemPerRankJ,
                    mMyRankK*mNElemPerRankK};

  // Create the derived datatype for the writing of the coordinates.
  MPI_Datatype writeType;
  MPI_Type_create_subarray(4, gsizes, lsizes, startInd, MPI_ORDER_FORTRAN,
                           MPI_FLOAT, &writeType);
  MPI_Type_commit(&writeType);

  // Create the file view needed for the collective write. The offset is set
  // to the current position.
  char datarep[] = "native";
  MPI_File_set_view(fh, offset, MPI_FLOAT, writeType, datarep, MPI_INFO_NULL);

  // Write the coordinate data in a single call.
  MPI_File_write_all(fh, coorBuf.data(), coorBuf.size(), MPI_FLOAT,
                     MPI_STATUS_IGNORE);

  // Release the memory of the derived data type again.
  MPI_Type_free(&writeType);

  // Determine the new offset and reset the file view.
  offset += 3*nDOFsTot*sizeof(float);

  MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, datarep, MPI_INFO_NULL);

  // Write the ASCII line for the connectivity info.
  std::sprintf(str_buf, "\nCELLS %i %i\n", nElemPerRankIJK*nRanks*nSubElem,
                                         9*nElemPerRankIJK*nRanks*nSubElem);
  if(rank == 0)
    MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  offset += strlen(str_buf)*sizeof(char);

  // Define the variables to create the subarray data type for the actual
  // writing of the connectivity. Note that only the first element of
  // gsizes and lsizes must be changed.
  lsizes[0] = 9*nSubElem;
  gsizes[0] = lsizes[0];

  // Create the derived datatype for the writing of the connectivities.
  MPI_Type_create_subarray(4, gsizes, lsizes, startInd, MPI_ORDER_FORTRAN,
                           MPI_INT, &writeType);
  MPI_Type_commit(&writeType);

  // Create the file view needed for the collective write. The offset is set
  // to the current position.
  MPI_File_set_view(fh, offset, MPI_INT, writeType, datarep, MPI_INFO_NULL);

  // Write the connectivity data in a single call.
  MPI_File_write_all(fh, connBuf.data(), connBuf.size(), MPI_INT,
                     MPI_STATUS_IGNORE);

  // Release the memory of the derived data type again.
  MPI_Type_free(&writeType);

  // Determine the new offset and reset the file view.
  offset += 9*nElemPerRankIJK*nRanks*nSubElem*sizeof(int);

  MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, datarep, MPI_INFO_NULL);

  // Write the ASCII line for the cell types.
  std::sprintf(str_buf, "\nCELL_TYPES %i\n", nElemPerRankIJK*nRanks*nSubElem);

  if(rank == 0)
    MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  offset += strlen(str_buf)*sizeof(char);

  // Define the variables to create the subarray data type for the actual
  // writing of the element type. Note that only the first element of
  // gsizes and lsizes must be changed.
  lsizes[0] = nSubElem;
  gsizes[0] = lsizes[0];

  // Create the derived datatype for the writing of the element types.
  MPI_Type_create_subarray(4, gsizes, lsizes, startInd, MPI_ORDER_FORTRAN,
                           MPI_INT, &writeType);
  MPI_Type_commit(&writeType);

  // Create the file view needed for the collective write. The offset is set
  // to the current position.
  MPI_File_set_view(fh, offset, MPI_INT, writeType, datarep, MPI_INFO_NULL);

  // Write the element type data in a single call.
  MPI_File_write_all(fh, elemTypeBuf.data(), elemTypeBuf.size(), MPI_INT,
                     MPI_STATUS_IGNORE);

  // Release the memory of the derived data type again.
  MPI_Type_free(&writeType);

  // Determine the new offset and reset the file view.
  offset += nElemPerRankIJK*nRanks*nSubElem*sizeof(int);

  MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, datarep, MPI_INFO_NULL);

  // Write the ASCII line for the point data.
  std::sprintf(str_buf, "\nPOINT_DATA %i\n", nDOFsTot);

  if(rank == 0)
    MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                      MPI_CHAR, MPI_STATUS_IGNORE);
  offset += strlen(str_buf)*sizeof(char);

  // Loop over the variables to be stored.
  for(unsigned int var=0; var<visNames.size(); ++var)
  {
    // Check whether this is a scalar or vector variable.
    // Note that the y- and z-components of a vector are skipped.
    bool writeVar = true, isVector = false;
    size_t found = visNames[var].find("_x");
    if(found != std::string::npos) isVector = true;

    found = visNames[var].find("_y");
    if(found != std::string::npos) writeVar = false;

    found = visNames[var].find("_z");
    if(found != std::string::npos) writeVar = false;
 
    // Check for a vector field.
    if( isVector )
    {
      // Vector variable. Remove the trailing "_x".
      visNames[var].erase(visNames[var].end()-2, visNames[var].end());

      // Set the file view for writing the ASCII line with information.
      MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, datarep, MPI_INFO_NULL);

      // Write the ASCII line with the information.
      std::sprintf(str_buf, "\nVECTORS %s float\n", visNames[var].c_str());

      if(rank == 0)
        MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                          MPI_CHAR, MPI_STATUS_IGNORE);
      offset += strlen(str_buf)*sizeof(char);

      // Define the variables to create the subarray data type for the actual
      // writing of the vector data. Note that only the first element of
      // gsizes and lsizes must be changed.
      lsizes[0] = 3*nDOFs;
      gsizes[0] = lsizes[0];

      // Create the derived datatype for the writing of the vector data.
      MPI_Type_create_subarray(4, gsizes, lsizes, startInd, MPI_ORDER_FORTRAN,
                               MPI_FLOAT, &writeType);
      MPI_Type_commit(&writeType);

      // Create the file view needed for the collective write. The offset is set
      // to the current position.
      MPI_File_set_view(fh, offset, MPI_FLOAT, writeType, datarep, MPI_INFO_NULL);

      // Write the vector data in a single call.
      MPI_File_write_all(fh, visBuf.data() + var*nDOFsLoc, 3*nDOFsLoc,
                         MPI_FLOAT, MPI_STATUS_IGNORE);

      // Release the memory of the derived data type again.
      MPI_Type_free(&writeType);

      // Determine the new offset.
      offset += 3*nDOFsTot*sizeof(float);
    }
    else if( writeVar )
    {
      // Scalar variable. Set the file view for writing the ASCII line with information.
      MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, datarep, MPI_INFO_NULL);

      // Write the first ASCII line with the information.
      std::sprintf(str_buf, "\nSCALARS %s float 1\n", visNames[var].c_str());

      if(rank == 0)
        MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                          MPI_CHAR, MPI_STATUS_IGNORE);
      offset += strlen(str_buf)*sizeof(char);

      // Reset the file view again.
      MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, datarep, MPI_INFO_NULL);

      // Write the second ASCII line with the information.
      std::sprintf(str_buf, "LOOKUP_TABLE default\n");

      if(rank == 0)
        MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                          MPI_CHAR, MPI_STATUS_IGNORE);
      offset += strlen(str_buf)*sizeof(char);

      // Define the variables to create the subarray data type for the actual
      // writing of the scalar data. Note that only the first element of
      // gsizes and lsizes must be changed.
      lsizes[0] = nDOFs;
      gsizes[0] = lsizes[0];

      // Create the derived datatype for the writing of the scalar data.
      MPI_Type_create_subarray(4, gsizes, lsizes, startInd, MPI_ORDER_FORTRAN,
                               MPI_FLOAT, &writeType);
      MPI_Type_commit(&writeType);

      // Create the file view needed for the collective write. The offset is set
      // to the current position.
      MPI_File_set_view(fh, offset, MPI_FLOAT, writeType, datarep, MPI_INFO_NULL);

      // Write the scalar data in a single call.
      MPI_File_write_all(fh, visBuf.data() + var*nDOFsLoc, nDOFsLoc,
                         MPI_FLOAT, MPI_STATUS_IGNORE);

      // Release the memory of the derived data type again.
      MPI_Type_free(&writeType);

      // Determine the new offset.
      offset += nDOFsTot*sizeof(float);
    }
  }

  // Close the file.
  MPI_File_close(&fh);

#else

  // No MPI is used, so the writing of the file can be done sequentially.
  // Open the visualization file for binary writing.
  FILE *fh = std::fopen(fileName.str().c_str(), "wb");
  if( !fh )
    TerminateAll("SolverClass::WriteParaviewSolution", __FILE__, __LINE__,
                 "Paraview visualization file could not be opened for writing");

  // Write the header of the visualization file.
  char str_buf[MAX_STRING_LENGTH];
  std::strcpy(str_buf, "# vtk DataFile Version 3.0\n");
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

  std::strcpy(str_buf, "vtk output\n");
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

  std::strcpy(str_buf, "BINARY\n");
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

  std::strcpy(str_buf, "DATASET UNSTRUCTURED_GRID\n");
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

  // Write the coordinates.
  std::sprintf(str_buf, "POINTS %i float\n", nDOFsTot);
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);
  std::fwrite(coorBuf.data(), sizeof(float), coorBuf.size(), fh);

  // Write the connectivity data.
  std::sprintf(str_buf, "\nCELLS %i %i\n", nElemPerRankIJK*nRanks*nSubElem,
                                         9*nElemPerRankIJK*nRanks*nSubElem);
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);
  std::fwrite(connBuf.data(), sizeof(int), connBuf.size(), fh);

  // Write the element type data.
  std::sprintf(str_buf, "\nCELL_TYPES %i\n", nElemPerRankIJK*nRanks*nSubElem);
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);
  std::fwrite(elemTypeBuf.data(), sizeof(int), elemTypeBuf.size(), fh);

  // Write the ASCII line for the point data.
  std::sprintf(str_buf, "\nPOINT_DATA %i\n", nDOFsTot);
  std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

  // Loop over the variables to be stored.
  for(unsigned int var=0; var<visNames.size(); ++var)
  {
    // Check whether this is a scalar or vector variable.
    // Note that the y- and z-components of a vector are skipped.
    bool writeVar = true, isVector = false;
    size_t found = visNames[var].find("_x");
    if(found != std::string::npos) isVector = true;

    found = visNames[var].find("_y");
    if(found != std::string::npos) writeVar = false;

    found = visNames[var].find("_z");
    if(found != std::string::npos) writeVar = false;
 
    // Check for a vector field.
    if( isVector )
    {
      // Vector variable. Remove the trailing "_x".
      visNames[var].erase(visNames[var].end()-2, visNames[var].end());

      // Write the ASCII line with the information.
      std::sprintf(str_buf, "\nVECTORS %s float\n", visNames[var].c_str());
      std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

      // Write the vector data.
      std::fwrite(visBuf.data() + var*nDOFsLoc, sizeof(float),
                  3*nDOFsLoc, fh);
    }
    else if( writeVar )
    {
      // Scalar variable. Write the ASCII lines with the information.
      std::sprintf(str_buf, "\nSCALARS %s float 1\n", visNames[var].c_str());
      std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

      std::sprintf(str_buf, "LOOKUP_TABLE default\n");
      std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

      // Write the scalar data.
      std::fwrite(visBuf.data() + var*nDOFsLoc, sizeof(float),
                  nDOFsLoc, fh);
    }
  }

  // Close the file again.
  std::fclose(fh);

#endif

  // Write a message that a paraview visualization file has been written.
  if(rank == 0)
  {
    std::cout << "#" << std::endl;
    std::cout << "# Paraview visualization file " << fileName.str() << " written." << std::endl;
    std::cout << "#" << std::endl << std::flush;
  }
}

//-----------------------------------------------------------------------------

// Function to write the restart solution.
void SolverClass::WriteRestartSolution(const int  timeStep,
                                       const bool writeRes)
{
  // Easier storage of the total number of elements in the three index directions.
  const int nElemI = mInputParam->mSubfaces[1][0]->mElemIBeg;
  const int nElemJ = mInputParam->mSubfaces[3][0]->mElemJBeg;
  const int nElemK = mInputParam->mSubfaces[5][0]->mElemKBeg;

  // Determine the name of the restart file.
  std::ostringstream fileName;
  fileName << mInputParam->mSolFile << "_TimeStep_" << timeStep;

  // Easier storage of the number of DOFs in one direction
  // and the total number of DOFs.
  const int nDOFs1D = mStandardHex.mNDOFs1D;
  const int nDOFs   = mStandardHex.mNDOFs;

  // Determine the total number of DOFs in the simulation and
  // the number of DOFs stored on one MPI rank.
  const int nElemPerRankIJ  = mNElemPerRankI*mNElemPerRankJ;
  const int nElemPerRankIJK = mNElemPerRankK*nElemPerRankIJ;

  const int nDOFsTot = nElemI*nElemJ*nElemK*nDOFs;
  const int nDOFsLoc = nElemPerRankIJK*nDOFs;

  // Determine the total number of variables to write.
  int nVarWrite = nVar + 3;
  if(mInputParam->mComputeAverageSolution && !writeRes) nVarWrite += nVar + 7;

  // Determine the 5 entries of the integer information.
  int header[] = {SU2_MagicNumber, nVarWrite, nDOFsTot, 1, 8};

  // Set the restart meta data. The first six entries contain the
  // time averaged force coefficients.
  su2double restartMetaData[] = {mCfxAve, mCfyAve, mCfzAve, mCfxVisAve, mCfyVisAve, mCfzVisAve, zero, zero};

  // Determine the time step to be written to the restart file. If averaging is applied,
  // the number of time steps for which averaging is applied is written. Otherwise
  // it is the current time step.
  int timeStepWrite = mInputParam->mComputeAverageSolution ? mNTimeStepAverage : timeStep;

  // Allocate the memory for the write buffer.
  std::vector<su2double> writeBuf(nDOFsLoc*nVarWrite);

  // Start of the parallel region.
#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    // Allocate the memory to store the solution in the nodal DOFs.
    std::vector<su2double *> solNodal(nVar, NULL);
    for(int i=0; i<nVar; ++i)
    {
      solNodal[i] = (su2double *) AllocateMemory(nDOFs*sizeof(su2double));
      if( !solNodal[i] )
        Terminate("SolverClass::WriteRestartSolution", __FILE__, __LINE__,
                  "Memory allocation failure for solNodal");
    }

    // Loop over the owned elements. A 1D loop is used, such that
    // the location in the write buffer can be easily determined.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int elem=0; elem<nElemPerRankIJK; ++elem)
    {
      // Retrieve the i-, j- and k-index of the element. Use a zero based
      // index here, because that is easier for the retrieval.
      int k   = elem/nElemPerRankIJ;
      int rem = elem - k*nElemPerRankIJ;
      int j   = rem/mNElemPerRankI;
      int i   = rem - j*mNElemPerRankI;

      // Increment i, j and k, because the owned elements start at 1.
      ++i; ++j; ++k;

      // Determine the location in writeBuf where the data of
      // this element should be stored.
      su2double *writeBufElem = writeBuf.data() + elem*nDOFs*nVarWrite;

      // Check whether residuals are written.
      if( writeRes )
      {
        // Loop over the DOFs of the element.
        int ind = 0;
        for(int DOF=0; DOF<nDOFs; ++DOF)
        {
          // Write the coordinates of the nodal solution DOFs.
          writeBufElem[ind++] = mElements[i][j][k]->mCoorNodalSolDOFs[0][DOF];
          writeBufElem[ind++] = mElements[i][j][k]->mCoorNodalSolDOFs[1][DOF];
          writeBufElem[ind++] = mElements[i][j][k]->mCoorNodalSolDOFs[2][DOF];

          // Copy the residuals to the write buffer.
          writeBufElem[ind++] = mElements[i][j][k]->mRes[0][DOF];
          writeBufElem[ind++] = mElements[i][j][k]->mRes[1][DOF];
          writeBufElem[ind++] = mElements[i][j][k]->mRes[2][DOF];
          writeBufElem[ind++] = mElements[i][j][k]->mRes[3][DOF];
          writeBufElem[ind++] = mElements[i][j][k]->mRes[4][DOF];
        }
      }
      else
      {
        // The solution is written. First determine the solution in
        // the nodal DOFs.
        TensorProductSolAndGradVolume(nDOFs1D, nVar, nDOFs1D, mStandardHex.mLegendreDOFs1D,
                                      mStandardHex.mDerLegendreDOFs1D,
                                      mElements[i][j][k]->mSol.data(), solNodal.data(),
                                      NULL, NULL, NULL);

        // Set the double pointers for the solution.
        su2double **sol        = solNodal.data();
        su2double **primAvg    = mElements[i][j][k]->mAvePrim.data();
        su2double **velProdAvg = mElements[i][j][k]->mAveVelProd.data();
        su2double *eddyVis     = mElements[i][j][k]->mAveEddyVis;

        // Convert the nodal solution to conservative variables, if needed.
        if(mInputParam->mFEMVariables == ENTROPY_VARIABLES)
        {
          // Abbreviate some expressions involving gamma.
          const su2double gm1    =  GamConstant - one;
          const su2double ov1mg  = -one/gm1;

          // Loop over the DOFs of the element.
#pragma omp simd
          for(int DOF=0; DOF<nDOFs; ++DOF)
          {
            // Compute the non-dimensional primitive variables from
            // the entropy variables.
            const su2double V4Inv =  one/sol[4][DOF];
            const su2double u     = -V4Inv*sol[1][DOF];
            const su2double v     = -V4Inv*sol[2][DOF];
            const su2double w     = -V4Inv*sol[3][DOF];
            const su2double eKin  =  half*(u*u + v*v + w*w);
            const su2double s     =  GamConstant - gm1*(sol[0][DOF] - sol[4][DOF]*eKin);
            const su2double tmp   = -sol[4][DOF]*EXP(s);
            const su2double rho   =  POW(tmp, ov1mg);
            const su2double p     = -rho*V4Inv;

            // Store the conservative variables again in sol.
            sol[0][DOF] = rho;
            sol[1][DOF] = rho*u;
            sol[2][DOF] = rho*v;
            sol[3][DOF] = rho*w;
            sol[4][DOF] = rho*eKin - ov1mg*p;
          }
        }

        // Reference value for the momentum variables.
        const su2double momRef = rhoRef*uRef;

        // Loop over the DOFs and write the coordinates and the
        // dimensional conservative variables to the write buffer.
        int ind = 0;
        for(int DOF=0; DOF<nDOFs; ++DOF)
        {
          // Write the coordinates of the nodal solution DOFs.
          writeBufElem[ind++] = mElements[i][j][k]->mCoorNodalSolDOFs[0][DOF];
          writeBufElem[ind++] = mElements[i][j][k]->mCoorNodalSolDOFs[1][DOF];
          writeBufElem[ind++] = mElements[i][j][k]->mCoorNodalSolDOFs[2][DOF];

          // Write the dimensional conservative variables.
          writeBufElem[ind++] = rhoRef*sol[0][DOF];
          writeBufElem[ind++] = momRef*sol[1][DOF];
          writeBufElem[ind++] = momRef*sol[2][DOF];
          writeBufElem[ind++] = momRef*sol[3][DOF];
          writeBufElem[ind++] = pRef  *sol[4][DOF];

          // Add the average solution, if needed.
          if( mInputParam->mComputeAverageSolution )
          {
            writeBufElem[ind++] = primAvg[0][DOF];
            writeBufElem[ind++] = primAvg[1][DOF];
            writeBufElem[ind++] = primAvg[2][DOF];
            writeBufElem[ind++] = primAvg[3][DOF];
            writeBufElem[ind++] = primAvg[4][DOF];

            writeBufElem[ind++] = velProdAvg[0][DOF];
            writeBufElem[ind++] = velProdAvg[1][DOF];
            writeBufElem[ind++] = velProdAvg[2][DOF];
            writeBufElem[ind++] = velProdAvg[3][DOF];
            writeBufElem[ind++] = velProdAvg[4][DOF];
            writeBufElem[ind++] = velProdAvg[5][DOF];

            writeBufElem[ind++] = eddyVis[DOF];
          }
        }
      }
    }

    // Release the memory of solNodal again.
    for(int i=0; i<nVar; ++i)
      FreeMemory((void **) &solNodal[i]);

  } // End of the OpenMP parallel region.

#ifdef HAVE_MPI
  // Wipe out any previous output file. The reason for doing this is that there
  // is no guarantee that an old file will be clobbered by the MPI write routines
  if(rank == 0) MPI_File_delete(fileName.str().c_str(), MPI_INFO_NULL);
  MPI_Barrier(MPI_COMM_WORLD);

  // Open the restart file for writing and check if it went OK.
  MPI_File fh;
  if(MPI_File_open(MPI_COMM_WORLD, fileName.str().c_str(),
                   MPI_MODE_WRONLY | MPI_MODE_CREATE,
                   MPI_INFO_NULL, &fh) != MPI_SUCCESS)
    TerminateAll("SolverClass::WriteRestartSolution", __FILE__, __LINE__,
                 "Restart file could not be opened for writing");

  // Define the variables to create the subarray data type for the actual IO.
  const int gsizes[]   = {(int) (nVarWrite*nDOFs), nElemI, nElemJ, nElemK};
  const int lsizes[]   = {(int) (nVarWrite*nDOFs), mNElemPerRankI, mNElemPerRankJ,
                          mNElemPerRankK};
  const int startInd[] = {0, mMyRankI*mNElemPerRankI, mMyRankJ*mNElemPerRankJ,
                          mMyRankK*mNElemPerRankK};

  // Create the derived datatype for the writing.
  MPI_Datatype writeType;
  MPI_Type_create_subarray(4, gsizes, lsizes, startInd, MPI_ORDER_FORTRAN,
                           MPI_SU2DOUBLE, &writeType);
  MPI_Type_commit(&writeType);

  // Write the header of the file. As this is a limited amount of data
  // the writing of the header is carried out by rank 0.
  if(rank == 0)
  {
    // Write the integer header.
    MPI_File_write(fh, header, 5, MPI_INT, MPI_STATUS_IGNORE);

    // Write the names of 3 coordinates and 5 conservative solution variables.
    char varName[CGNS_STRING_SIZE];
    strncpy(varName, "x", CGNS_STRING_SIZE);
    MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

    strncpy(varName, "y", CGNS_STRING_SIZE);
    MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

    strncpy(varName, "z", CGNS_STRING_SIZE);
    MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

    strncpy(varName, "Density", CGNS_STRING_SIZE);
    MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

    strncpy(varName, "Momentum_x", CGNS_STRING_SIZE);
    MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

    strncpy(varName, "Momentum_y", CGNS_STRING_SIZE);
    MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

    strncpy(varName, "Momentum_z", CGNS_STRING_SIZE);
    MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

    strncpy(varName, "Energy", CGNS_STRING_SIZE);
    MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

        // Write the names of the average variables, if needed.
    if( mInputParam->mComputeAverageSolution )
    {
      strncpy(varName, "Density_Average", CGNS_STRING_SIZE);
      MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

      strncpy(varName, "VelocityX_Average", CGNS_STRING_SIZE);
      MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

      strncpy(varName, "VelocityY_Average", CGNS_STRING_SIZE);
      MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

      strncpy(varName, "VelocityZ_Average", CGNS_STRING_SIZE);
      MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

      strncpy(varName, "Pressure_Average", CGNS_STRING_SIZE);
      MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

      strncpy(varName, "VelocityX_VelocityX_Average", CGNS_STRING_SIZE);
      MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

      strncpy(varName, "VelocityY_VelocityY_Average", CGNS_STRING_SIZE);
      MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

      strncpy(varName, "VelocityZ_VelocityZ_Average", CGNS_STRING_SIZE);
      MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

      strncpy(varName, "VelocityX_VelocityY_Average", CGNS_STRING_SIZE);
      MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

      strncpy(varName, "VelocityX_VelocityZ_Average", CGNS_STRING_SIZE);
      MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

      strncpy(varName, "VelocityY_VelocityZ_Average", CGNS_STRING_SIZE);
      MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);

      strncpy(varName, "EddyVis_Average", CGNS_STRING_SIZE);
      MPI_File_write(fh, varName, CGNS_STRING_SIZE, MPI_CHAR, MPI_STATUS_IGNORE);
    }
  }

  // Determine the byte offset where the data for the solution
  // variables starts.
  MPI_Offset offset = 5*sizeof(int) + nVarWrite*CGNS_STRING_SIZE;

  // Create the file view needed for the collective write. The offset is set
  // to the position where the data of the element starts;
  char datarep[] = "native";
  MPI_File_set_view(fh, offset, MPI_SU2DOUBLE, writeType, datarep, MPI_INFO_NULL);

  // Write the solution data in a single call.
  MPI_File_write_all(fh, writeBuf.data(), writeBuf.size(), MPI_SU2DOUBLE,
                     MPI_STATUS_IGNORE);

  // Determine the new offset and reset the file view.
  offset += nDOFsTot*nVarWrite*sizeof(su2double);

  MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, datarep, MPI_INFO_NULL);

  // Rank 0 writes the meta data.
  if(rank == 0)
  {
    // Write the number of time steps.
    MPI_File_write_at(fh, offset, &timeStepWrite, 1, MPI_INT, MPI_STATUS_IGNORE);
    offset += sizeof(int);

    // Write the additional floating point data.
    MPI_File_write_at(fh, offset, restartMetaData, 8, MPI_SU2DOUBLE,  MPI_STATUS_IGNORE);
  }

  // Close the file.
  MPI_File_close(&fh);

  // Release the memory of the derived data type again.
  MPI_Type_free(&writeType);

#else
  // No MPI is used, so the writing of the file can be done sequentially.
  // Open the restart file for binary writing.
  FILE *fh = std::fopen(fileName.str().c_str(), "wb");
  if( !fh )
    Terminate("SolverClass::WriteRestartSolution", __FILE__, __LINE__,
              "Solution file could not be opened for binary writing.");

  // Write the integer header.
  std::fwrite(header, 5, sizeof(int), fh);

  // Write the names of 3 coordinates and 5 conservative solution variables.
  char varName[CGNS_STRING_SIZE];
  strncpy(varName, "x", CGNS_STRING_SIZE);
  std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

  strncpy(varName, "y", CGNS_STRING_SIZE);
  std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

  strncpy(varName, "z", CGNS_STRING_SIZE);
  std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

  strncpy(varName, "Density", CGNS_STRING_SIZE);
  std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

  strncpy(varName, "Momentum_x", CGNS_STRING_SIZE);
  std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

  strncpy(varName, "Momentum_y", CGNS_STRING_SIZE);
  std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

  strncpy(varName, "Momentum_z", CGNS_STRING_SIZE);
  std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

  strncpy(varName, "Energy", CGNS_STRING_SIZE);
  std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

  // Write the names of the average variables, if needed.
  if( mInputParam->mComputeAverageSolution )
  {
    strncpy(varName, "Density_Average", CGNS_STRING_SIZE);
    std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

    strncpy(varName, "VelocityX_Average", CGNS_STRING_SIZE);
    std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

    strncpy(varName, "VelocityY_Average", CGNS_STRING_SIZE);
    std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

    strncpy(varName, "VelocityZ_Average", CGNS_STRING_SIZE);
    std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

    strncpy(varName, "Pressure_Average", CGNS_STRING_SIZE);
    std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

    strncpy(varName, "VelocityX_VelocityX_Average", CGNS_STRING_SIZE);
    std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

    strncpy(varName, "VelocityY_VelocityY_Average", CGNS_STRING_SIZE);
    std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

    strncpy(varName, "VelocityZ_VelocityZ_Average", CGNS_STRING_SIZE);
    std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

    strncpy(varName, "VelocityX_VelocityY_Average", CGNS_STRING_SIZE);
    std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

    strncpy(varName, "VelocityX_VelocityZ_Average", CGNS_STRING_SIZE);
    std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

    strncpy(varName, "VelocityY_VelocityZ_Average", CGNS_STRING_SIZE);
    std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);

    strncpy(varName, "EddyVis_Average", CGNS_STRING_SIZE);
    std::fwrite(varName, CGNS_STRING_SIZE, sizeof(char), fh);
  }

  // Write the solution buffer.
  std::fwrite(writeBuf.data(), writeBuf.size(), sizeof(su2double), fh);

  // Write the time step.
  std::fwrite(&timeStepWrite, 1, sizeof(int), fh);

  // Write the restart meta data.
  std::fwrite(restartMetaData, 8, sizeof(double), fh);

  // Close the file again.
  std::fclose(fh);

#endif

  // Write a message that a restart file has been written.
  if(rank == 0)
  {
    std::cout << "#" << std::endl;
    std::cout << "# Restart file " << fileName.str() << " written." << std::endl;
    std::cout << "#" << std::endl << std::flush;
  }
}

//-----------------------------------------------------------------------------

// Function, which writes the space time average data if averaging is applied.
void SolverClass::WriteSpaceTimeAverageData(void)
{
  // Return immediately if no average is computed.
  if( !mInputParam->mComputeAverageSolution ) return;

  // Define the maximum string length for the writing.
  const int MAX_STRING_LENGTH = 255;

  // Check for big endian. If not, the bytes must be swapped when
  // writing a binary VTK file.
  bool BigEndian;
  union {int i; char c[4];} val;
  val.i = 0x76543210;
  if (val.c[0] == 0x10) BigEndian = false;
  else                  BigEndian = true;

  // Define the names of the averaged variables to be visualized. Note that a vector
  // quantity must be defined contiguously, with first the x-component, followed
  // by the y-component and finally the z-component. Furthermore, the names of
  // the components of a vector field must end with _x, _y and _z, respectively.
  std::vector<std::string> visNames;

  visNames.push_back("AverageDensity");
  visNames.push_back("AverageVelocity_x");
  visNames.push_back("AverageVelocity_y");
  visNames.push_back("AverageVelocity_z");
  visNames.push_back("AveragePressure");
  visNames.push_back("AverageTemperature");
  visNames.push_back("AverageMach");

  visNames.push_back("AverageR_uu");
  visNames.push_back("AverageR_vv");
  visNames.push_back("AverageR_ww");
  visNames.push_back("AverageR_uv");
  visNames.push_back("AverageR_uw");
  visNames.push_back("AverageR_vw");

  if(mInputParam->mSGSModelType != NO_SGS_MODEL)
    visNames.push_back("AverageRatioEddyViscosityLaminarViscosity");

  // Easier storage of the number of DOFs and the integration weights
  // of the 1D integration rule based on the DOFs.
  const int       nDOFs1D  = mStandardHex.mNDOFs1D;
  const int       nDOFs2D  = nDOFs1D*nDOFs1D;
  const su2double *weights = mStandardHex.mIntegrationWeights1DDOFs.data();

  // Easier storage of the total number of elements in the i-, j- and k-direction.
  const int nElemI = mInputParam->mSubfaces[1][0]->mElemIBeg;
  const int nElemJ = mInputParam->mSubfaces[3][0]->mElemJBeg;
  const int nElemK = mInputParam->mSubfaces[5][0]->mElemKBeg;

  // Determine the scaling factor for the weights, such that the weights of
  // all integration points of all elements in k-direction sum up to 1.
  const su2double scaleWeights = one/(two*nElemK);

  // Determine the number of local elements on a k-plane and the
  // number of local DOFs per MPI rank.
  const int nElemPerRankIJ = mNElemPerRankI*mNElemPerRankJ;
  const int nDOFsLoc       = nElemPerRankIJ*nDOFs2D;

  // Determine the total number of DOFs in a k-plane
  const int nDOFsTot = nElemI*nElemJ*nDOFs2D;

  // Allocate the memory for the visualization buffer of
  // the averaged solution in k-direction.
  const int nVarAvg = (int) visNames.size();
  std::vector<su2double> visBuf(mNElemPerRankI*mNElemPerRankJ*nDOFs2D*nVarAvg, zero);

  // Loop over the local elements in k-direction.
  for(int k=1; k<=mNElemPerRankK; ++k)
  {
    // Loop over the 1D DOFs in k-direction.
    for(int kk=0; kk<nDOFs1D; ++kk)
    {
      // Determine the weight of the current integration point, corrected for
      // the number of elements in k-direction.
      const su2double weight = weights[kk]*scaleWeights;

      // Determine the k-offset for the DOFs in the elements that
      // are considered for the averaging.
      const int kkOff = kk*nDOFs2D;

      // Loop over the owned elements in j- and i-direction.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
      for(int elem=0; elem<nElemPerRankIJ; ++elem)
      {
        // Retrieve the i- and j-index of the element. Use a zero based
        // index here, because that is easier for the retrieval.
        int j = elem/mNElemPerRankI;
        int i = elem - j*mNElemPerRankI;

        // Increment i and j, because the owned elements start at 1.
        ++i; ++j;

        // Set the pointers for the average variables.
        su2double **avePrim    = mElements[i][j][k]->mAvePrim.data();
        su2double **aveVelProd = mElements[i][j][k]->mAveVelProd.data();
        su2double *aveEddyVis  = mElements[i][j][k]->mAveEddyVis;

        // Loop over the number of variables to be determined.
        for(unsigned int visVar=0; visVar<visNames.size(); ++visVar)
        {
          // Check for the variable name and act accordingly.
          if(visNames[visVar] == "AverageDensity")
          {
            // Set the pointer to the correction location in visBuf.
            su2double *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs2D;

            // Update the average density in the buffer.
            for(int l=0; l<nDOFs2D; ++l)
              buf[l] += weight*avePrim[0][l+kkOff];
          }
          else if(visNames[visVar] == "AverageVelocity_x")
          {
            // Set the pointer to the correction location in visBuf.
            su2double *buf = visBuf.data() + visVar*nDOFsLoc + 3*elem*nDOFs2D;

            // Average velocity in x-direction. As velocity is a vector
            // variable, also the y- and z-velocity must be stored.
            for(int l=0; l<nDOFs2D; ++l)
            {
              const int l3 = 3*l;
              const int ll = l+kkOff;
              buf[l3]   += weight*avePrim[1][ll];
              buf[l3+1] += weight*avePrim[2][ll];
              buf[l3+2] += weight*avePrim[3][ll];
            }

            // Update visVar by an additional 2, because the y- and z-velocity
            // have been stored as well.
            visVar += 2;
          }
          else if(visNames[visVar] == "AveragePressure")
          {
            // Set the pointer to the correction location in visBuf.
            su2double *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs2D;

            // Update the average pressure in the buffer.
            for(int l=0; l<nDOFs2D; ++l)
              buf[l] += weight*avePrim[4][l+kkOff];
          }
          else if(visNames[visVar] == "AverageTemperature")
          {
            // Check if the average temperature can be computed.
            if(mNTimeStepAverage > 0)
            {
              // Set the pointer to the correction location in visBuf.
              su2double *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs2D;

              // Update the average temperature in the buffer.
              for(int l=0; l<nDOFs2D; ++l)
              {
                const int ll = l+kkOff;
                const su2double T = avePrim[4][ll]/(RGas*avePrim[0][ll]);
                buf[l] += weight*T;
              }
            }
          }
          else if(visNames[visVar] == "AverageMach")
          {
            // Check if the average Mach number can be computed.
            if(mNTimeStepAverage > 0)
            {
              // Set the pointer to the correction location in visBuf.
              su2double *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs2D;

              // Update the average Mach number in the buffer.
              for(int l=0; l<nDOFs2D; ++l)
              {
                const int ll = l+kkOff;
                const su2double V2 = avePrim[1][ll]*avePrim[1][ll]
                                   + avePrim[2][ll]*avePrim[2][ll]
                                   + avePrim[3][ll]*avePrim[3][ll];
                const su2double a2 = GamConstant*avePrim[4][ll]/avePrim[0][ll];
                const su2double Ma = SQRT(V2/a2);
                buf[l] += weight*Ma;
              }
            }
          }
          else if(visNames[visVar] == "AverageR_uu")
          {
            // Set the pointer to the correction location in visBuf.
            su2double *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs2D;

            // Update the average u'u' Reynolds stress in the buffer.
            for(int l=0; l<nDOFs2D; ++l)
            {
              const int ll = l+kkOff;
              buf[l] += weight*(aveVelProd[0][ll] - avePrim[1][ll]*avePrim[1][ll]);
            }
          }
          else if(visNames[visVar] == "AverageR_vv")
          {
            // Set the pointer to the correction location in visBuf.
            su2double *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs2D;

            // Update the average v'v' Reynolds stress in the buffer.
            for(int l=0; l<nDOFs2D; ++l)
            {
              const int ll = l+kkOff;
              buf[l] += weight*(aveVelProd[1][ll] - avePrim[2][ll]*avePrim[2][ll]);
            }
          }
          else if(visNames[visVar] == "AverageR_ww")
          {
            // Set the pointer to the correction location in visBuf.
            su2double *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs2D;

            // Update the average w'w' Reynolds stress in the buffer.
            for(int l=0; l<nDOFs2D; ++l)
            {
              const int ll = l+kkOff;
              buf[l] += weight*(aveVelProd[2][ll] - avePrim[3][ll]*avePrim[3][ll]);
            }
          }
          else if(visNames[visVar] == "AverageR_uv")
          {
            // Set the pointer to the correction location in visBuf.
            su2double *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs2D;

            // Update the average u'v' Reynolds stress in the buffer.
            for(int l=0; l<nDOFs2D; ++l)
            {
              const int ll = l+kkOff;
              buf[l] += weight*(aveVelProd[3][ll] - avePrim[1][ll]*avePrim[2][ll]);
            }
          }
          else if(visNames[visVar] == "AverageR_uw")
          {
            // Set the pointer to the correction location in visBuf.
            su2double *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs2D;

            // Update the average u'w' Reynolds stress in the buffer.
            for(int l=0; l<nDOFs2D; ++l)
            {
              const int ll = l+kkOff;
              buf[l] += weight*(aveVelProd[4][ll] - avePrim[1][ll]*avePrim[3][ll]);
            }
          }
          else if(visNames[visVar] == "AverageR_vw")
          {
            // Set the pointer to the correction location in visBuf.
            su2double *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs2D;

            // Update the average v'w' Reynolds stress in the buffer.
            for(int l=0; l<nDOFs2D; ++l)
            {
              const int ll = l+kkOff;
              buf[l] += weight*(aveVelProd[5][ll] - avePrim[2][ll]*avePrim[3][ll]);
            }
          }
          else if(visNames[visVar] == "AverageRatioEddyViscosityLaminarViscosity")
          {
            // Set the pointer to the correction location in visBuf.
            su2double *buf = visBuf.data() + visVar*nDOFsLoc + elem*nDOFs2D;

            // Update the average ratio of eddy and laminar viscosity.
            const su2double muDimInv = uRef/(mu*pRef);
            for(int l=0; l<nDOFs2D; ++l)
              buf[l] += weight*(muDimInv*aveEddyVis[l+kkOff]);
          }
          else
          {
            // Unknown visualization variable encountered.
#ifdef HAVE_OPENMP
            if(omp_get_thread_num() == 0)
#endif
              TerminateAll("SolverClass::WriteSpaceTimeAverageData",
                           __FILE__, __LINE__, "Unknown visualization name");
          }
        }
      }
    } 
  }

  // When MPI is used, it may be necessary to carry out a reduce in
  // the k-direction to obtain the correct data.
#ifdef HAVE_MPI

  // The reduction in k-direction is only needed when there is more than
  // 1 rank in the k-direction. Test this.
  if(mNRanksK > 1)
  {
    // Create the subcommunicator needed to carry out the reduction
    // in k-direction.
    int color = mMyRankJ*mNRanksI + mMyRankI;

    MPI_Comm commK;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &commK);

    // Make a copy of visBuf and sum the data in k-direction.
    // The result only needs to be known on the root rank of
    // commK, because only this rank is involved in the actual
    // writing of the file with the averaged data.
    std::vector<su2double> tmpVisBuf = visBuf;
    MPI_Reduce(tmpVisBuf.data(), visBuf.data(), visBuf.size(),
               MPI_SU2DOUBLE, MPI_SUM, 0, commK);

    // Release the memory of the subcommunicator.
    MPI_Comm_free(&commK);
  }

  // Create the subcommunicator for the ranks in the plane k = constant.
  // This communicator is needed for the actual writing of the file.
  MPI_Comm commIJ;
  MPI_Comm_split(MPI_COMM_WORLD, mMyRankK, rank, &commIJ);

#endif
 
  // Check if this rank participates in the writing of the spatial
  // averaged file.
  if(mMyRankK == 0)
  {
    // Allocate the memory for the buffers to store the coordinates and
    // the connectivities. The element types to be written are all the same,
    // namely a standard quadrilateral, VTK number 9. Hence the data for
    // the vector elemTypeBuf can be set directly.
    const int nSubElem = (nDOFs1D-1)*(nDOFs1D-1);
    std::vector<su2double> coorBuf(3*nDOFsLoc);
    std::vector<int>       connBuf(5*nElemPerRankIJ*nSubElem);
    std::vector<int>       elemTypeBuf(nElemPerRankIJ*nSubElem, 9);

    // Loop over the owned elements. A 1D loop is used instead of a 2D loop,
    // such that it is easier to retrieve the position in the arrays
    // used for the writing of the actual data.
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
    for(int elem=0; elem<nElemPerRankIJ; ++elem)
    {
      // Retrieve the i- and j-index of the element. Use a zero based
      // index here, because that is easier for the retrieval.
      int j = elem/mNElemPerRankI;
      int i = elem - j*mNElemPerRankI;

      // Determine the global i- and j-index of this element.
      const int elemIndJ = mMyRankJ*mNElemPerRankJ + j;
      const int elemIndI = mMyRankI*mNElemPerRankI + i;

      // Increment i and j, because the owned elements start at 1.
      ++i; ++j;

      // Determine the location in the buffers where the data of
      // this element should be stored.
      su2double *coorElem = coorBuf.data() + 3*elem*nDOFs2D;
      int       *connElem = connBuf.data() + 5*elem*nSubElem;

      // Loop over the 2D DOFs of the element to copy the coordinates to coorElem.
      int ind = 0;
      for(int ii=0; ii<nDOFs2D; ++ii)
      {
        coorElem[ind++] = mElements[i][j][1]->mCoorNodalSolDOFs[0][ii];
        coorElem[ind++] = mElements[i][j][1]->mCoorNodalSolDOFs[1][ii];
        coorElem[ind++] = mElements[i][j][1]->mCoorNodalSolDOFs[2][ii];
      }

      // The current element has global indices (elemIndI, elemIndJ).
      // Convert this to a 1D index and determine the global offset where the
      // connectivities of this element start.
      const long offsetConn = (elemIndJ*nElemI + elemIndI)*nDOFs2D;

      // Loop over the subelements of the current element and determine the
      // connectivities. Note that the first entry in the connectivity array
      // is the number of points per element, which is 4.
      ind = 0;
      for(int jj=1; jj<nDOFs1D; ++jj)
      {
        const int jOff = offsetConn + (jj-1)*nDOFs1D;
        for(int ii=1; ii<nDOFs1D; ++ii)
        {
          const int n0 = jOff + ii-1;
          connElem[ind++] = 4;
          connElem[ind++] = n0;
          connElem[ind++] = n0 + 1;
          connElem[ind++] = n0 + 1 + nDOFs1D;
          connElem[ind++] = n0 + nDOFs1D;
        }
      }
    }

    // Convert the buffers to big endian, if needed.
    if( !BigEndian )
    {
      SwapBytes(coorBuf.data(),     sizeof(su2double), coorBuf.size());
      SwapBytes(connBuf.data(),     sizeof(int),       connBuf.size());
      SwapBytes(elemTypeBuf.data(), sizeof(int),       elemTypeBuf.size());
      SwapBytes(visBuf.data(),      sizeof(su2double), visBuf.size());
    }

#ifdef HAVE_MPI
    // Wipe out any previous output file. The reason for doing this is that there
    // is no guarantee that an old file will be clobbered by the MPI write routines
    if(rank == 0) MPI_File_delete(mInputParam->mSpaceAverageFile.c_str(),
                                  MPI_INFO_NULL);
    MPI_Barrier(commIJ);

    // Open the averaging file for writing and check if it went OK.
    MPI_File fh;
    if(MPI_File_open(commIJ, mInputParam->mSpaceAverageFile.c_str(),
                     MPI_MODE_WRONLY | MPI_MODE_CREATE,
                     MPI_INFO_NULL, &fh) != MPI_SUCCESS)
      TerminateAll("SolverClass::WriteSpaceTimeAverageData", __FILE__, __LINE__,
                   "Space average file could not be opened for writing");

    // Initialize the file offset to zero.
    MPI_Offset offset = 0;

    // Write the header of the visualization file.
    char str_buf[MAX_STRING_LENGTH];
    std::strcpy(str_buf, "# vtk DataFile Version 3.0\n");
    if(rank == 0)
      MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                        MPI_CHAR, MPI_STATUS_IGNORE);
    offset += strlen(str_buf)*sizeof(char);

    std::strcpy(str_buf, "vtk output\n");
    if(rank == 0)
      MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                        MPI_CHAR, MPI_STATUS_IGNORE);
    offset += strlen(str_buf)*sizeof(char);

    std::strcpy(str_buf, "BINARY\n");
    if(rank == 0)
      MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                        MPI_CHAR, MPI_STATUS_IGNORE);
    offset += strlen(str_buf)*sizeof(char);

    std::strcpy(str_buf, "DATASET UNSTRUCTURED_GRID\n");
    if(rank == 0)
      MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                        MPI_CHAR, MPI_STATUS_IGNORE);
    offset += strlen(str_buf)*sizeof(char);

    if(sizeof(su2double) == sizeof(double))
      std::sprintf(str_buf, "POINTS %i double\n", nDOFsTot);
    else
      std::sprintf(str_buf, "POINTS %i float\n", nDOFsTot);

    if(rank == 0)
      MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                        MPI_CHAR, MPI_STATUS_IGNORE);
    offset += strlen(str_buf)*sizeof(char);

    // Define the variables to create the subarray data type for the actual
    // writing of the coordinates.
    int gsizes[]   = {(int) (3*nDOFs2D), nElemI, nElemJ};
    int lsizes[]   = {(int) (3*nDOFs2D), mNElemPerRankI, mNElemPerRankJ};
    int startInd[] = {0, mMyRankI*mNElemPerRankI, mMyRankJ*mNElemPerRankJ};

    // Create the derived datatype for the writing of the coordinates.
    MPI_Datatype writeType;
    MPI_Type_create_subarray(3, gsizes, lsizes, startInd, MPI_ORDER_FORTRAN,
                             MPI_SU2DOUBLE, &writeType);
    MPI_Type_commit(&writeType);

    // Create the file view needed for the collective write. The offset is set
    // to the current position.
    char datarep[] = "native";
    MPI_File_set_view(fh, offset, MPI_SU2DOUBLE, writeType, datarep, MPI_INFO_NULL);

    // Write the coordinate data in a single call.
    MPI_File_write_all(fh, coorBuf.data(), coorBuf.size(), MPI_SU2DOUBLE,
                       MPI_STATUS_IGNORE);

    // Release the memory of the derived data type again.
    MPI_Type_free(&writeType);

    // Determine the new offset and reset the file view.
    offset += 3*nDOFsTot*sizeof(su2double);

    MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, datarep, MPI_INFO_NULL);

    // Write the ASCII line for the connectivity info.
    std::sprintf(str_buf, "\nCELLS %i %i\n", nElemI*nElemJ*nSubElem,
                                           5*nElemI*nElemJ*nSubElem);
    if(rank == 0)
      MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                        MPI_CHAR, MPI_STATUS_IGNORE);
    offset += strlen(str_buf)*sizeof(char);

    // Define the variables to create the subarray data type for the actual
    // writing of the connectivity. Note that only the first element of
    // gsizes and lsizes must be changed.
    lsizes[0] = 5*nSubElem;
    gsizes[0] = lsizes[0];

    // Create the derived datatype for the writing of the connectivities.
    MPI_Type_create_subarray(3, gsizes, lsizes, startInd, MPI_ORDER_FORTRAN,
                             MPI_INT, &writeType);
    MPI_Type_commit(&writeType);

    // Create the file view needed for the collective write. The offset is set
    // to the current position.
    MPI_File_set_view(fh, offset, MPI_INT, writeType, datarep, MPI_INFO_NULL);

    // Write the connectivity data in a single call.
    MPI_File_write_all(fh, connBuf.data(), connBuf.size(), MPI_INT,
                       MPI_STATUS_IGNORE);

    // Release the memory of the derived data type again.
    MPI_Type_free(&writeType);

    // Determine the new offset and reset the file view.
    offset += 5*nElemI*nElemJ*nSubElem*sizeof(int);

    MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, datarep, MPI_INFO_NULL);

    // Write the ASCII line for the cell types.
    std::sprintf(str_buf, "\nCELL_TYPES %i\n", nElemI*nElemJ*nSubElem);

    if(rank == 0)
      MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                        MPI_CHAR, MPI_STATUS_IGNORE);
    offset += strlen(str_buf)*sizeof(char);

    // Define the variables to create the subarray data type for the actual
    // writing of the element type. Note that only the first element of
    // gsizes and lsizes must be changed.
    lsizes[0] = nSubElem;
    gsizes[0] = lsizes[0];
 
    // Create the derived datatype for the writing of the element types.
    MPI_Type_create_subarray(3, gsizes, lsizes, startInd, MPI_ORDER_FORTRAN,
                             MPI_INT, &writeType);
    MPI_Type_commit(&writeType);

    // Create the file view needed for the collective write. The offset is set
    // to the current position.
    MPI_File_set_view(fh, offset, MPI_INT, writeType, datarep, MPI_INFO_NULL);

    // Write the element type data in a single call.
    MPI_File_write_all(fh, elemTypeBuf.data(), elemTypeBuf.size(), MPI_INT,
                       MPI_STATUS_IGNORE);

    // Release the memory of the derived data type again.
    MPI_Type_free(&writeType);

    // Determine the new offset and reset the file view.
    offset += nElemI*nElemJ*nSubElem*sizeof(int);

    MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, datarep, MPI_INFO_NULL);

    // Write the ASCII line for the point data.
    std::sprintf(str_buf, "\nPOINT_DATA %i\n", nDOFsTot);

    if(rank == 0)
      MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                        MPI_CHAR, MPI_STATUS_IGNORE);
    offset += strlen(str_buf)*sizeof(char);

    // Loop over the variables to be stored.
    for(unsigned int var=0; var<visNames.size(); ++var)
    {
      // Check whether this is a scalar or vector variable.
      // Note that the y- and z-components of a vector are skipped.
      bool writeVar = true, isVector = false;
      size_t found = visNames[var].find("_x");
      if(found != std::string::npos) isVector = true;

      found = visNames[var].find("_y");
      if(found != std::string::npos) writeVar = false;

      found = visNames[var].find("_z");
      if(found != std::string::npos) writeVar = false;

      // Check for a vector field.
      if( isVector )
      {
        // Vector variable. Remove the trailing "_x".
        visNames[var].erase(visNames[var].end()-2, visNames[var].end());

        // Set the file view for writing the ASCII line with information.
        MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, datarep, MPI_INFO_NULL);

        // Write the ASCII line with the information.
        if(sizeof(su2double) == sizeof(double))
          std::sprintf(str_buf, "\nVECTORS %s double\n", visNames[var].c_str());
        else
          std::sprintf(str_buf, "\nVECTORS %s float\n", visNames[var].c_str());

        if(rank == 0)
          MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                            MPI_CHAR, MPI_STATUS_IGNORE);
        offset += strlen(str_buf)*sizeof(char);

        // Define the variables to create the subarray data type for the actual
        // writing of the vector data. Note that only the first element of
        // gsizes and lsizes must be changed.
        lsizes[0] = 3*nDOFs2D;
        gsizes[0] = lsizes[0];

        // Create the derived datatype for the writing of the vector data.
        MPI_Type_create_subarray(3, gsizes, lsizes, startInd, MPI_ORDER_FORTRAN,
                                 MPI_SU2DOUBLE, &writeType);
        MPI_Type_commit(&writeType);

        // Create the file view needed for the collective write. The offset is set
        // to the current position.
        MPI_File_set_view(fh, offset, MPI_SU2DOUBLE, writeType, datarep, MPI_INFO_NULL);

        // Write the vector data in a single call.
        MPI_File_write_all(fh, visBuf.data() + var*nDOFsLoc, 3*nDOFsLoc,
                           MPI_SU2DOUBLE, MPI_STATUS_IGNORE);

        // Release the memory of the derived data type again.
        MPI_Type_free(&writeType);

        // Determine the new offset.
        offset += 3*nDOFsTot*sizeof(su2double);
      }
      else if( writeVar )
      {
        // Scalar variable. Set the file view for writing the ASCII line with information.
        MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, datarep, MPI_INFO_NULL);

        // Write the first ASCII line with the information.
        if(sizeof(su2double) == sizeof(double))
          std::sprintf(str_buf, "\nSCALARS %s double 1\n", visNames[var].c_str());
        else
          std::sprintf(str_buf, "\nSCALARS %s float 1\n", visNames[var].c_str());

        if(rank == 0)
          MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                            MPI_CHAR, MPI_STATUS_IGNORE);
        offset += strlen(str_buf)*sizeof(char);

        // Reset the file view again.
        MPI_File_set_view(fh, 0, MPI_BYTE, MPI_BYTE, datarep, MPI_INFO_NULL);

        // Write the second ASCII line with the information.
        std::sprintf(str_buf, "LOOKUP_TABLE default\n");

        if(rank == 0)
          MPI_File_write_at(fh, offset, str_buf, strlen(str_buf),
                            MPI_CHAR, MPI_STATUS_IGNORE);
        offset += strlen(str_buf)*sizeof(char);

        // Define the variables to create the subarray data type for the actual
        // writing of the scalar data. Note that only the first element of
        // gsizes and lsizes must be changed.
        lsizes[0] = nDOFs2D;
        gsizes[0] = lsizes[0];

        // Create the derived datatype for the writing of the vector data.
        MPI_Type_create_subarray(3, gsizes, lsizes, startInd, MPI_ORDER_FORTRAN,
                                 MPI_SU2DOUBLE, &writeType);
        MPI_Type_commit(&writeType);

        // Create the file view needed for the collective write. The offset is set
        // to the current position.
        MPI_File_set_view(fh, offset, MPI_SU2DOUBLE, writeType, datarep, MPI_INFO_NULL);

        // Write the scalar data in a single call.
        MPI_File_write_all(fh, visBuf.data() + var*nDOFsLoc, nDOFsLoc,
                           MPI_SU2DOUBLE, MPI_STATUS_IGNORE);

        // Release the memory of the derived data type again.
        MPI_Type_free(&writeType);

        // Determine the new offset.
        offset += nDOFsTot*sizeof(su2double);
      }
    }
    
    // Close the file.
    MPI_File_close(&fh);

#else
    // No MPI is used, so the writing of the file can be done sequentially.
    // Open the averaging file for binary writing
    FILE *fh = std::fopen(mInputParam->mSpaceAverageFile.c_str(), "wb");
    if( !fh )
      TerminateAll("SolverClass::WriteSpaceTimeAverageData", __FILE__, __LINE__,
                   "Space average file could not be opened for writing");

    // Write the header of the averaging file.
    char str_buf[MAX_STRING_LENGTH];
    std::strcpy(str_buf, "# vtk DataFile Version 3.0\n");
    std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

    std::strcpy(str_buf, "vtk output\n");
    std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

    std::strcpy(str_buf, "BINARY\n");
    std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

    std::strcpy(str_buf, "DATASET UNSTRUCTURED_GRID\n");
    std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

    // Write the coordinates.
    if(sizeof(su2double) == sizeof(double))
      std::sprintf(str_buf, "POINTS %i double\n", nDOFsTot);
    else
      std::sprintf(str_buf, "POINTS %i float\n", nDOFsTot);
    std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);
    std::fwrite(coorBuf.data(), sizeof(su2double), coorBuf.size(), fh);
    
    // Write the connectivity data.
    std::sprintf(str_buf, "\nCELLS %i %i\n", nElemI*nElemJ*nSubElem,
                                           5*nElemI*nElemJ*nSubElem);
    std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);
    std::fwrite(connBuf.data(), sizeof(int), connBuf.size(), fh);

    // Write the element type data.
    std::sprintf(str_buf, "\nCELL_TYPES %i\n", nElemI*nElemJ*nSubElem);
    std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);
    std::fwrite(elemTypeBuf.data(), sizeof(int), elemTypeBuf.size(), fh);

    // Write the ASCII line for the point data.
    std::sprintf(str_buf, "\nPOINT_DATA %i\n", nDOFsTot);
    std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);
 
    // Loop over the variables to be stored.
    for(unsigned int var=0; var<visNames.size(); ++var)
    {
      // Check whether this is a scalar or vector variable.
      // Note that the y- and z-components of a vector are skipped.
      bool writeVar = true, isVector = false;
      size_t found = visNames[var].find("_x");
      if(found != std::string::npos) isVector = true;

      found = visNames[var].find("_y");
      if(found != std::string::npos) writeVar = false;

      found = visNames[var].find("_z");
      if(found != std::string::npos) writeVar = false;

      // Check for a vector field.
      if( isVector )
      {
        // Vector variable. Remove the trailing "_x".
        visNames[var].erase(visNames[var].end()-2, visNames[var].end());

        // Write the ASCII line with the information.
        if(sizeof(su2double) == sizeof(double))
          std::sprintf(str_buf, "\nVECTORS %s double\n", visNames[var].c_str());
        else
          std::sprintf(str_buf, "\nVECTORS %s float\n", visNames[var].c_str());
        std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

        // Write the vector data.
        std::fwrite(visBuf.data() + var*nDOFsLoc, sizeof(su2double),
                    3*nDOFsLoc, fh);
      }
      else if( writeVar )
      {
        // Scalar variable. Write the ASCII lines with the information.
        if(sizeof(su2double) == sizeof(double))
          std::sprintf(str_buf, "\nSCALARS %s double 1\n", visNames[var].c_str());
        else
          std::sprintf(str_buf, "\nSCALARS %s float 1\n", visNames[var].c_str());
        std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

        std::sprintf(str_buf, "LOOKUP_TABLE default\n");
        std::fwrite(str_buf, sizeof(char), strlen(str_buf), fh);

        // Write the scalar data.
        std::fwrite(visBuf.data() + var*nDOFsLoc, sizeof(su2double),
                    nDOFsLoc, fh);
      }
    }

    // Close the file again.
    std::fclose(fh);

#endif
  }

  // Write a message that a file containing the space time averaged
  // data has been written. Also write the averaged data for the
  // lift and drag coefficients.
  if(rank == 0)
  {
    std::cout << "#" << std::endl;
    std::cout << "# Space time averaged data file "
              << mInputParam->mSpaceAverageFile << " written." << std::endl;
    std::cout << "#" << std::endl << std::flush;

    if( mInputParam->mMonitorForces )
    {
      const su2double ClInv = mCfxAve   *mInputParam->mLiftDirection[0]
                            + mCfyAve   *mInputParam->mLiftDirection[1]
                            + mCfzAve   *mInputParam->mLiftDirection[2];
      const su2double ClVis = mCfxVisAve*mInputParam->mLiftDirection[0]
                            + mCfyVisAve*mInputParam->mLiftDirection[1]
                            + mCfzVisAve*mInputParam->mLiftDirection[2];
      const su2double CdInv = mCfxAve   *mInputParam->mVelDirFreeStream[0]
                            + mCfyAve   *mInputParam->mVelDirFreeStream[1]
                            + mCfzAve   *mInputParam->mVelDirFreeStream[2];
      const su2double CdVis = mCfxVisAve*mInputParam->mVelDirFreeStream[0]
                            + mCfyVisAve*mInputParam->mVelDirFreeStream[1]
                            + mCfzVisAve*mInputParam->mVelDirFreeStream[2];

      std::cout << "# Time averaged Cl:    " << std::setw(13) << std::scientific << ClInv + ClVis << std::endl;
      std::cout << "# Time averaged ClVis: " << std::setw(13) << std::scientific << ClVis         << std::endl;
      std::cout << "# Time averaged Cd:    " << std::setw(13) << std::scientific << CdInv + CdVis << std::endl;
      std::cout << "# Time averaged CdVis: " << std::setw(13) << std::scientific << CdVis         << std::endl;
      std::cout << "#" << std::endl << std::flush;
    }
  }

  // For MPI the memory of the subcommunicator must be released.
  // Also synchronize the MPI ranks, because not all MPI ranks
  // may have contributed to the writing.
#ifdef HAVE_MPI
  MPI_Comm_free(&commIJ);
  MPI_Barrier(MPI_COMM_WORLD);
#endif 

}

//-----------------------------------------------------------------------------

// Function, which writes the grid of the solution DOFs in SU2 format.
void SolverClass::WriteSU2GridFile(void)
{
  // The writing of the SU2 grid file only works in sequential mode.
  if(nRanks > 1)
    TerminateAll("SolverClass::WriteSU2GridFile", __FILE__, __LINE__,
                 "This function only works for one MPI rank");

  // Easier storage of the polynomial degree of the grid and solution.
  const int nPolyGrid = mInputParam->mNPolyGridDOFs;
  const int nPolySol  = mInputParam->mNPolySolDOFs;

  // Determine the element type in SU2 format.
  int elemType  = 12 + 100*(nPolySol-1);
  int boundType =  9 + 100*(nPolySol-1);

  if(nPolyGrid != nPolySol)
  {
    elemType  += 10000*(nPolySol+1);
    boundType += 10000*(nPolySol+1);
  }

  // Determine the number of grid points in the three directions.
  const int il = nPolyGrid*mNElemPerRankI + 1;
  const int jl = nPolyGrid*mNElemPerRankJ + 1;
  const int kl = nPolyGrid*mNElemPerRankK + 1;

  // Allocate the memory for the mapping of the point numbers and
  // initialize it.
  std::vector<int> mapPoints(il*jl*kl);
  for(int i=0; i<(il*jl*kl); ++i)
    mapPoints[i] = i;

  // Determine the actual mapping by taking the 1 to 1 connectivities into account.
  // Use a while loop in case multiple 1 to 1 subfaces are present.
  bool situationChanged = true;
  while( situationChanged )
  {
    // Initialize situationChanged to false.
    situationChanged = false;

    // Loop over the subfaces.
    for(unsigned int i=0; i<mInputParam->mSubfaces.size(); ++i)
    {
      for(unsigned int j=0; j<mInputParam->mSubfaces[i].size(); ++j)
      {
        // Check for a 1 to 1 matching subface.
        if( mInputParam->mSubfaces[i][j]->Is1To1Matching() )
        {
          // Determine the actual subface range where the polynomial degree
          // is taken into account.
          const int iBeg = nPolyGrid*mInputParam->mSubfaces[i][j]->mElemIBeg;
          const int iEnd = nPolyGrid*mInputParam->mSubfaces[i][j]->mElemIEnd;
          const int jBeg = nPolyGrid*mInputParam->mSubfaces[i][j]->mElemJBeg;
          const int jEnd = nPolyGrid*mInputParam->mSubfaces[i][j]->mElemJEnd;
          const int kBeg = nPolyGrid*mInputParam->mSubfaces[i][j]->mElemKBeg;
          const int kEnd = nPolyGrid*mInputParam->mSubfaces[i][j]->mElemKEnd;

          // Same for the donor range, but only for the begin indices.
          const int dIBeg = nPolyGrid*mInputParam->mSubfaces[i][j]->GetElemDonorIBeg();
          const int dJBeg = nPolyGrid*mInputParam->mSubfaces[i][j]->GetElemDonorJBeg();
          const int dKBeg = nPolyGrid*mInputParam->mSubfaces[i][j]->GetElemDonorKBeg();

          // Loop over the nodal range of the subface.
          for(int kk=kBeg; kk<=kEnd; ++kk)
          {
            for(int jj=jBeg; jj<=jEnd; ++jj)
            {
              for(int ii=iBeg; ii<=iEnd; ++ii)
              {
                // Determine the donor indices. Note that the transformation matrix
                // between subfaces is the identity matrix.
                const int dii = ii - iBeg + dIBeg;
                const int djj = jj - jBeg + dJBeg;
                const int dkk = kk - kBeg + dKBeg;

                // Convert the i,j,k indices to one dimensional indices.
                const int thisInd  =  kk*jl*il +  jj*il +  ii;
                const int donorInd = dkk*jl*il + djj*il + dii;

                // Check if the donor index is larger than the index
                // of the current node.
                if(donorInd > thisInd)
                {
                  // If the mapping of the donor index is not equal to
                  // the mapping of the current index, make it so. In that
                  // case also set situationChanged to true.
                  if(mapPoints[donorInd] != mapPoints[thisInd])
                  {
                    mapPoints[donorInd] = mapPoints[thisInd];
                    situationChanged = true;
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  // Determine the number of different points in the grid
  // and the mapping from old to new points.
  std::vector<int> old2NewPoints(il*jl*kl);
  int nPoints = 0;

  for(int i=0; i<(il*jl*kl); ++i)
  {
    // Check if this point maps onto itself.
    if(mapPoints[i] == i)
    {
      // This point remains in the grid.
      old2NewPoints[i] = nPoints;
      ++nPoints;
    }
    else
    {
      // The point is duplicated. As the mapping point has already been
      // treated old2NewPoints can be set.
      old2NewPoints[i] = old2NewPoints[mapPoints[i]];
    }
  }

  // Open the file for writing. If not successful, print an error message and exit.
  std::ostringstream message;
  std::ofstream gridFile(mInputParam->mSU2GridFile.c_str());
  if( !gridFile )
  {
    message << "SU2 grid file " << mInputParam->mSU2GridFile
            << " could not be opened for writing";
    Terminate("SolverClass::WriteSU2GridFile", __FILE__, __LINE__, message.str());
  }

  // Write the dimension of the problem, which is 3.
  gridFile << "%" << std::endl
           << "% Problem dimension" << std::endl
           << "%" << std::endl
           << "NDIME= 3" << std::endl;

  // Write the number of elements.
  gridFile << "%" << std::endl
           << "% Inner element connectivity" << std::endl
           << "%" << std::endl
           << "NELEM= " << mNElemPerRankI*mNElemPerRankJ*mNElemPerRankK << std::endl;

  // Loop over the elements.
  for(int k=0; k<mNElemPerRankK; ++k)
  {
    for(int j=0; j<mNElemPerRankJ; ++j)
    {
      for(int i=0; i<mNElemPerRankI; ++i)
      {
        gridFile << " " << elemType;

        // Determine the index of the back lower left index of this hexahedron.
        const int indBLL = nPolyGrid*(k*jl*il + j*il + i);

        // Make a distinction between linear and high order elements.
        if(nPolyGrid == 1)
        {
          const int indFLL = indBLL + jl*il;
          const int indBUL = indBLL + il;
          const int indFUL = indFLL + il;

          gridFile << " " << old2NewPoints[indBLL]   << " " << old2NewPoints[indBLL+1]
                   << " " << old2NewPoints[indBUL+1] << " " << old2NewPoints[indBUL]
                   << " " << old2NewPoints[indFLL]   << " " << old2NewPoints[indFLL+1]
                   << " " << old2NewPoints[indFUL+1] << " " << old2NewPoints[indFUL];
        }
        else
        {
          for(int kk=0; kk<=nPolyGrid; ++kk)
          {
            const int indLL = indBLL + kk*jl*il;
            for(int jj=0; jj<=nPolyGrid; ++jj)
            {
              const int indL = indLL + jj*il;
              for(int ii=0; ii<=nPolyGrid; ++ii)
                gridFile << " " << old2NewPoints[indL + ii];
            }
          }
        }

        gridFile << " " << k*mNElemPerRankI*mNElemPerRankJ + j*mNElemPerRankI + i << std::endl;
      }
    }
  }

  // Write the number of coordinates.
  gridFile << "%" << std::endl
           << "% Node coordinates" << std::endl
           << "%" << std::endl
           << "NPOIN= " << nPoints << std::endl;

  // Loop over the elements and polynomial degree.
  int indOld = 0, indNew = 0;
  for(int k=1; k<=mNElemPerRankK; ++k)
  {
    const int kkStart = (k==1) ? 0 : 1;
    for(int kk=kkStart; kk<=nPolyGrid; ++kk)
    {
      const int offK = kk*(nPolyGrid+1)*(nPolyGrid+1);
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        const int jjStart = (j==1) ? 0 : 1;
        for(int jj=jjStart; jj<=nPolyGrid; ++jj)
        {
          const int offJ = offK + jj*(nPolyGrid+1);
          for(int i=1; i<=mNElemPerRankI; ++i)
          {
            const int iiStart = (i==1) ? 0 : 1;
            for(int ii=iiStart; ii<=nPolyGrid; ++ii, ++indOld)
            {
              // Check if this coordinate must be written.
              if(mapPoints[indOld] == indOld)
              {
                const int DOFNode = offJ + ii;
                gridFile.precision(15);
                gridFile << std::setw(25) << std::scientific << mElements[i][j][k]->mCoorNodalGridDOFs[0][DOFNode]
                         << std::setw(25) << std::scientific << mElements[i][j][k]->mCoorNodalGridDOFs[1][DOFNode]
                         << std::setw(25) << std::scientific << mElements[i][j][k]->mCoorNodalGridDOFs[2][DOFNode]
                         << " " << indNew << std::endl;
                ++indNew;
              }
            }
          }
        }
      }
    }
  }

  // Determine the total number of boundary subfaces, i.e. markers.
  int nMarkers = 0;
  for(unsigned int i=0; i<mInputParam->mSubfaces.size(); ++i)
    for(unsigned int j=0; j<mInputParam->mSubfaces[i].size(); ++j)
      if( !(mInputParam->mSubfaces[i][j]->Is1To1Matching()) )
        ++nMarkers;

  // Write the number of markers.
  gridFile << "%" << std::endl
           << "% Boundary elements" << std::endl
           << "%" << std::endl
           << "NMARK= " << nMarkers << std::endl;

  // Loop over the subfaces of the iMin boundary.
  for(unsigned int l=0; l<mInputParam->mSubfaces[0].size(); ++l)
  {
    if( !(mInputParam->mSubfaces[0][l]->Is1To1Matching()) )
    {
      // Determine the name for this subface.
      std::ostringstream nameSubface;
      nameSubface << "iMin_Subface_" << l;

      // Write the name of this subface.
      gridFile << "MARKER_TAG= " << nameSubface.str() << std::endl;

      // Determine the number of boundary elements and write this number.
      const int nBoundElem = (mInputParam->mSubfaces[0][l]->mElemJEnd - mInputParam->mSubfaces[0][l]->mElemJBeg)
                           * (mInputParam->mSubfaces[0][l]->mElemKEnd - mInputParam->mSubfaces[0][l]->mElemKBeg);
      gridFile << "MARKER_ELEMS= " << nBoundElem << std::endl;

      // Loop over the boundary elements of this iMin subface.
      for(int k=mInputParam->mSubfaces[0][l]->mElemKBeg; k<mInputParam->mSubfaces[0][l]->mElemKEnd; ++k)
      {
        for(int j=mInputParam->mSubfaces[0][l]->mElemJBeg; j<mInputParam->mSubfaces[0][l]->mElemJEnd; ++j)
        {
          // Determine the indices in the i,j,k grid and write the boundary element type.
          const int indLL = nPolyGrid*(k*jl*il + j*il);
          gridFile << " " << boundType;

          // Write the true connectivity of this boundary element.
          if(nPolyGrid == 1)
          {
            const int indUL = indLL + jl*il;
            gridFile << " " << old2NewPoints[indLL]    << " " << old2NewPoints[indLL+il]
                     << " " << old2NewPoints[indUL+il] << " " << old2NewPoints[indUL];
          }
          else
          {
            for(int kk=0; kk<=nPolyGrid; ++kk)
            {
              const int indL = indLL + kk*jl*il;
              for(int jj=0; jj<=nPolyGrid; ++jj)
                 gridFile << " " << old2NewPoints[indL + jj*il];
            }
          }

          gridFile << std::endl;
        }
      }
    }
  }

  // Loop over the subfaces of the iMax boundary.
  for(unsigned int l=0; l<mInputParam->mSubfaces[1].size(); ++l)
  {
    if( !(mInputParam->mSubfaces[1][l]->Is1To1Matching()) )
    {
      // Determine the name for this subface.
      std::ostringstream nameSubface;
      nameSubface << "iMax_Subface_" << l;

      // Write the name of this subface.
      gridFile << "MARKER_TAG= " << nameSubface.str() << std::endl;

      // Determine the number of boundary elements and write this number.
      const int nBoundElem = (mInputParam->mSubfaces[1][l]->mElemJEnd - mInputParam->mSubfaces[1][l]->mElemJBeg)
                           * (mInputParam->mSubfaces[1][l]->mElemKEnd - mInputParam->mSubfaces[1][l]->mElemKBeg);
      gridFile << "MARKER_ELEMS= " << nBoundElem << std::endl;

      // Loop over the boundary elements of this iMax subface.
      for(int k=mInputParam->mSubfaces[1][l]->mElemKBeg; k<mInputParam->mSubfaces[1][l]->mElemKEnd; ++k)
      {
        for(int j=mInputParam->mSubfaces[1][l]->mElemJBeg; j<mInputParam->mSubfaces[1][l]->mElemJEnd; ++j)
        {
          // Determine the indices in the i,j,k grid and write the boundary element type.
          const int indLL = il-1 + nPolyGrid*(k*jl*il + j*il);
          gridFile << " " << boundType;

          // Write the true connectivity of this boundary element.
          if(nPolyGrid == 1)
          {
            const int indUL = indLL + jl*il;
            gridFile << " " << old2NewPoints[indLL]    << " " << old2NewPoints[indLL+il]
                     << " " << old2NewPoints[indUL+il] << " " << old2NewPoints[indUL];
          }
          else
          {
            for(int kk=0; kk<=nPolyGrid; ++kk)
            {
              const int indL = indLL + kk*jl*il;
              for(int jj=0; jj<=nPolyGrid; ++jj)
                 gridFile << " " << old2NewPoints[indL + jj*il];
            }
          }

          gridFile << std::endl;
        }
      }
    }
  }

  // Loop over the subfaces of the jMin boundary.
  for(unsigned int l=0; l<mInputParam->mSubfaces[2].size(); ++l)
  {
    if( !(mInputParam->mSubfaces[2][l]->Is1To1Matching()) )
    {
      // Determine the name for this subface.
      std::ostringstream nameSubface;
      nameSubface << "jMin_Subface_" << l;

      // Write the name of this subface.
      gridFile << "MARKER_TAG= " << nameSubface.str() << std::endl;

      // Determine the number of boundary elements and write this number.
      const int nBoundElem = (mInputParam->mSubfaces[2][l]->mElemIEnd - mInputParam->mSubfaces[2][l]->mElemIBeg)
                           * (mInputParam->mSubfaces[2][l]->mElemKEnd - mInputParam->mSubfaces[2][l]->mElemKBeg);
      gridFile << "MARKER_ELEMS= " << nBoundElem << std::endl;

      // Loop over the boundary elements of this jMin subface.
      for(int k=mInputParam->mSubfaces[2][l]->mElemKBeg; k<mInputParam->mSubfaces[2][l]->mElemKEnd; ++k)
      {
        for(int i=mInputParam->mSubfaces[2][l]->mElemIBeg; i<mInputParam->mSubfaces[2][l]->mElemIEnd; ++i)
        {
          // Determine the indices in the i,j,k grid and write the boundary element type.
          const int indLL = nPolyGrid*(k*jl*il + i);
          gridFile << " " << boundType;

          // Write the true connectivity of this boundary element.
          if(nPolyGrid == 1)
          {
            const int indUL = indLL + jl*il;
            gridFile << " " << old2NewPoints[indLL]   << " " << old2NewPoints[indLL+1]
                     << " " << old2NewPoints[indUL+1] << " " << old2NewPoints[indUL];
          }
          else
          {
            for(int kk=0; kk<=nPolyGrid; ++kk)
            {
              const int indL = indLL + kk*jl*il;
              for(int ii=0; ii<=nPolyGrid; ++ii)
                 gridFile << " " << old2NewPoints[indL + ii];
            }
          }

          gridFile << std::endl;
        }
      }
    }
  }

  // Loop over the subfaces of the jMax boundary.
  for(unsigned int l=0; l<mInputParam->mSubfaces[3].size(); ++l)
  {
    if( !(mInputParam->mSubfaces[3][l]->Is1To1Matching()) )
    {
      // Determine the name for this subface.
      std::ostringstream nameSubface;
      nameSubface << "jMax_Subface_" << l;

      // Write the name of this subface.
      gridFile << "MARKER_TAG= " << nameSubface.str() << std::endl;

      // Determine the number of boundary elements and write this number.
      const int nBoundElem = (mInputParam->mSubfaces[3][l]->mElemIEnd - mInputParam->mSubfaces[3][l]->mElemIBeg)
                           * (mInputParam->mSubfaces[3][l]->mElemKEnd - mInputParam->mSubfaces[3][l]->mElemKBeg);
      gridFile << "MARKER_ELEMS= " << nBoundElem << std::endl;

      // Loop over the boundary elements of this jMax subface.
      for(int k=mInputParam->mSubfaces[3][l]->mElemKBeg; k<mInputParam->mSubfaces[3][l]->mElemKEnd; ++k)
      {
        for(int i=mInputParam->mSubfaces[3][l]->mElemIBeg; i<mInputParam->mSubfaces[3][l]->mElemIEnd; ++i)
        {
          // Determine the indices in the i,j,k grid and write the boundary element type.
          const int indLL = (jl-1)*il + nPolyGrid*(k*jl*il + i);
          gridFile << " " << boundType;

          // Write the true connectivity of this boundary element.
          if(nPolyGrid == 1)
          {
            const int indUL = indLL + jl*il;
            gridFile << " " << old2NewPoints[indLL]   << " " << old2NewPoints[indLL+1]
                     << " " << old2NewPoints[indUL+1] << " " << old2NewPoints[indUL];
          }
          else
          {
            for(int kk=0; kk<=nPolyGrid; ++kk)
            {
              const int indL = indLL + kk*jl*il;
              for(int ii=0; ii<=nPolyGrid; ++ii)
                 gridFile << " " << old2NewPoints[indL + ii];
            }
          }

          gridFile << std::endl;
        }
      }
    }
  }

  // Loop over the subfaces of the kMin boundary.
  for(unsigned int l=0; l<mInputParam->mSubfaces[4].size(); ++l)
  {
    if( !(mInputParam->mSubfaces[4][l]->Is1To1Matching()) )
    {
      // Determine the name for this subface.
      std::ostringstream nameSubface;
      nameSubface << "kMin_Subface_" << l;

      // Write the name of this subface.
      gridFile << "MARKER_TAG= " << nameSubface.str() << std::endl;

      // Determine the number of boundary elements and write this number.
      const int nBoundElem = (mInputParam->mSubfaces[4][l]->mElemIEnd - mInputParam->mSubfaces[4][l]->mElemIBeg)
                           * (mInputParam->mSubfaces[4][l]->mElemJEnd - mInputParam->mSubfaces[4][l]->mElemJBeg);
      gridFile << "MARKER_ELEMS= " << nBoundElem << std::endl;

      // Loop over the boundary elements of this kMin subface.
      for(int j=mInputParam->mSubfaces[4][l]->mElemJBeg; j<mInputParam->mSubfaces[4][l]->mElemJEnd; ++j)
      {
        for(int i=mInputParam->mSubfaces[4][l]->mElemIBeg; i<mInputParam->mSubfaces[4][l]->mElemIEnd; ++i)
        {
          // Determine the indices in the i,j,k grid and write the boundary element type.
          const int indLL = nPolyGrid*(j*il + i);
          gridFile << " " << boundType;

          // Write the true connectivity of this boundary element.
          if(nPolyGrid == 1)
          {
            const int indUL = indLL + il;
            gridFile << " " << old2NewPoints[indLL]   << " " << old2NewPoints[indLL+1]
                     << " " << old2NewPoints[indUL+1] << " " << old2NewPoints[indUL];
          }
          else
          {
            for(int jj=0; jj<=nPolyGrid; ++jj)
            {
              const int indL = indLL + jj*il;
              for(int ii=0; ii<=nPolyGrid; ++ii)
                 gridFile << " " << old2NewPoints[indL + ii];
            }
          }

          gridFile << std::endl;
        }
      }
    }
  }

  // Loop over the subfaces of the kMax boundary.
  for(unsigned int l=0; l<mInputParam->mSubfaces[5].size(); ++l)
  {
    if( !(mInputParam->mSubfaces[5][l]->Is1To1Matching()) )
    {
      // Determine the name for this subface.
      std::ostringstream nameSubface;
      nameSubface << "kMax_Subface_" << l;

      // Write the name of this subface.
      gridFile << "MARKER_TAG= " << nameSubface.str() << std::endl;

      // Determine the number of boundary elements and write this number.
      const int nBoundElem = (mInputParam->mSubfaces[5][l]->mElemIEnd - mInputParam->mSubfaces[5][l]->mElemIBeg)
                           * (mInputParam->mSubfaces[5][l]->mElemJEnd - mInputParam->mSubfaces[5][l]->mElemJBeg);
      gridFile << "MARKER_ELEMS= " << nBoundElem << std::endl;

      // Loop over the boundary elements of this kMax subface.
      for(int j=mInputParam->mSubfaces[5][l]->mElemJBeg; j<mInputParam->mSubfaces[5][l]->mElemJEnd; ++j)
      {
        for(int i=mInputParam->mSubfaces[5][l]->mElemIBeg; i<mInputParam->mSubfaces[5][l]->mElemIEnd; ++i)
        {
          // Determine the indices in the i,j,k grid and write the boundary element type.
          const int indLL = (kl-1)*jl*il + nPolyGrid*(j*il + i);
          gridFile << " " << boundType;

          // Write the true connectivity of this boundary element.
          if(nPolyGrid == 1)
          {
            const int indUL = indLL + il;
            gridFile << " " << old2NewPoints[indLL]   << " " << old2NewPoints[indLL+1]
                     << " " << old2NewPoints[indUL+1] << " " << old2NewPoints[indUL];
          }
          else
          {
            for(int jj=0; jj<=nPolyGrid; ++jj)
            {
              const int indL = indLL + jj*il;
              for(int ii=0; ii<=nPolyGrid; ++ii)
                 gridFile << " " << old2NewPoints[indL + ii];
            }
          }

          gridFile << std::endl;
        }
      }
    }
  }

  // Close the grid file.
  gridFile.close();

  // Write a message that the SU2 grid file has been written.
  std::cout << "#" << std::endl;
  std::cout << "# SU2 grid file " << mInputParam->mSU2GridFile << " written." << std::endl;
  std::cout << "#" << std::endl << std::flush;
}

//-----------------------------------------------------------------------------

// Function, which writes the time step data to std::cout.
void SolverClass::WriteTimeStepData(const int       timeStep,
                                    const double    time0,
                                    const su2double *monitoringData)
{
  // Only rank 0 performs this task.
  if(rank == 0)
  {
    // Compute the maximum Mach number, as currently the square of
    // this value is stored.
    const su2double MachMax = SQRT(monitoringData[6]);

    // Write the time step and CPU time to stdout.
    std::cout.precision(5);
    std::cout << std::setw(10) << timeStep
              << std::setw(13) << std::scientific << GetWallTime() - time0;

    // Write the dimensional time and maximum Mach number.
    std::cout << std::setw(13) << std::scientific << timeStep*mInputParam->mDeltaTSynchr;
    std::cout << std::setw(13) << std::scientific << MachMax;

    // Write the ratio of eddy viscosity and laminar viscosity if a
    // subgrid scale model is used.
    if( mInputParam->mSGSModel )
      std::cout << std::setw(13) << std::scientific << monitoringData[7]/mu;

    // Write the lift and drag coefficients if the force coefficients
    // must be monitored.
    if( mInputParam->mMonitorForces )
    {
      // Compute the inviscid and viscous lift and drag coefficients.
      const su2double ClInv = monitoringData[0]*mInputParam->mLiftDirection[0]
                            + monitoringData[1]*mInputParam->mLiftDirection[1]
                            + monitoringData[2]*mInputParam->mLiftDirection[2];
      const su2double ClVis = monitoringData[3]*mInputParam->mLiftDirection[0]
                            + monitoringData[4]*mInputParam->mLiftDirection[1]
                            + monitoringData[5]*mInputParam->mLiftDirection[2];
      const su2double CdInv = monitoringData[0]*mInputParam->mVelDirFreeStream[0]
                            + monitoringData[1]*mInputParam->mVelDirFreeStream[1]
                            + monitoringData[2]*mInputParam->mVelDirFreeStream[2];
      const su2double CdVis = monitoringData[3]*mInputParam->mVelDirFreeStream[0]
                            + monitoringData[4]*mInputParam->mVelDirFreeStream[1]
                            + monitoringData[5]*mInputParam->mVelDirFreeStream[2];

      // Write the data.
      std::cout << std::setw(13) << std::scientific << ClInv + ClVis
                << std::setw(13) << std::scientific << CdInv + CdVis
                << std::setw(13) << std::scientific << ClVis
                << std::setw(13) << std::scientific << CdVis;
    }

    // Write the new line.
    std::cout << std::endl << std::flush;
  }
}

//-----------------------------------------------------------------------------
// Below this line there are only functions that are needed when MPI is used.
//-----------------------------------------------------------------------------

#ifdef HAVE_MPI

// Function, which communicates the metric terms of the max boundaries.
void SolverClass::CommunicateMetricTermsMaxBoundaries(void)
{
  // Determine the number of metric items that must be communicated per element.
  const int nItemsPerElem = 5 + mElements[1][1][1]->mSurfMetricIntIMin.size()*mStandardHex.mNIntegration2DPad;

  // Determine for how many items per element the communication buffers
  // have been allocated.
  const int nItemsAlloc = nVar*mStandardHex.mNDOFs;

  // Determine the number of communication cycles, which are needed to
  // communicate all the metric data.
  const int nCommCycles = (nItemsPerElem%nItemsAlloc) ? (1 + nItemsPerElem/nItemsAlloc)
                                                      : (nItemsPerElem/nItemsAlloc);

  //------------------------------------------------------------
  // Communication of the metric data in i-direction, if needed.
  //------------------------------------------------------------

  // Check if data must be communicated in i-direction.
  if(mNRanksI > 1)
  {
    // Allocate the memory for the buffers to store the metric data.
    std::vector<su2double> sendBuf(nItemsPerElem*mNElemCommIMax);
    std::vector<su2double> recvBuf(nItemsPerElem*mNElemCommIMin);

    // Copy the metric data of the iMax boundary for the elements on
    // the iMax boundary of the local domain into sendBuf.
    unsigned int ii = 0;
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        if( mElements[mNElemPerRankI+1][j][k] )
        {
          if(mElements[mNElemPerRankI+1][j][k]->mElemType == HALO_ELEMENT_IMAX)
          {
            sendBuf[ii++] = mElements[mNElemPerRankI][j][k]->mLenScaleLES;
            sendBuf[ii++] = mElements[mNElemPerRankI][j][k]->mLenScaleVolume;
            sendBuf[ii++] = mElements[mNElemPerRankI][j][k]->mLenScaleIDir;
            sendBuf[ii++] = mElements[mNElemPerRankI][j][k]->mLenScaleJDir;
            sendBuf[ii++] = mElements[mNElemPerRankI][j][k]->mLenScaleKDir;

            for(unsigned l=0; l<mElements[mNElemPerRankI][j][k]->mSurfMetricIntIMax.size(); ++l)
              for(int m=0; m<mStandardHex.mNIntegration2DPad; ++m)
                sendBuf[ii++] = mElements[mNElemPerRankI][j][k]->mSurfMetricIntIMax[l][m];
          }
        }
      }
    }

    // Loop over the number of communication cycles.
    for(int comm=0; comm<nCommCycles; ++comm)
    {
      // Copy the required range from sendBuf into mCommBufIMax.
      unsigned int indBeg = comm*nItemsAlloc*mNElemCommIMax;
      unsigned int indEnd = (comm+1)*nItemsAlloc*mNElemCommIMax;
      if(indEnd > sendBuf.size()) indEnd = sendBuf.size();

      ii = 0;
      for(unsigned int i=indBeg; i<indEnd; ++i, ++ii)
        mCommBufIMax[ii] = sendBuf[i];

      // Start and complete the persistent communication.
      MPI_Startall(mCommRequestSolIDir.size(), mCommRequestSolIDir.data());
      MPI_Waitall(mCommRequestSolIDir.size(), mCommRequestSolIDir.data(),
                  MPI_STATUSES_IGNORE);

      // Copy the required range from mCommBufIMin into recvBuf.
      indBeg = comm*nItemsAlloc*mNElemCommIMin;
      indEnd = (comm+1)*nItemsAlloc*mNElemCommIMin;
      if(indEnd > recvBuf.size()) indEnd = recvBuf.size();

      ii = 0;
      for(unsigned int i=indBeg; i<indEnd; ++i, ++ii)
        recvBuf[i] = mCommBufIMin[ii];
    }

    // Copy the metric data from recvBuf into into the iMax boundary of the
    // elements on the iMin boundary of the local domain. The memory of
    // mSurfMetricIntIMax must be allocated for the halo elements.
    ii = 0;
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        if( mElements[0][j][k] )
        {
          if(mElements[0][j][k]->mElemType == HALO_ELEMENT_IMIN)
          {
            mElements[0][j][k]->mLenScaleLES    = recvBuf[ii++];
            mElements[0][j][k]->mLenScaleVolume = recvBuf[ii++];
            mElements[0][j][k]->mLenScaleIDir   = recvBuf[ii++];
            mElements[0][j][k]->mLenScaleJDir   = recvBuf[ii++];
            mElements[0][j][k]->mLenScaleKDir   = recvBuf[ii++];

            mElements[0][j][k]->mSurfMetricIntIMax.resize(mElements[1][j][k]->mSurfMetricIntIMax.size());
            for(unsigned l=0; l<mElements[0][j][k]->mSurfMetricIntIMax.size(); ++l)
            {
              mElements[0][j][k]->mSurfMetricIntIMax[l] = (su2double *) AllocateMemory(mStandardHex.mNIntegration2DPad*sizeof(su2double));
              if( !(mElements[0][j][k]->mSurfMetricIntIMax[l]) )
                Terminate("SolverClass::CommunicateMetricTermsMaxBoundaries", __FILE__, __LINE__,
                          "Memory allocation failure for mSurfMetricIntIMax");

              for(int m=0; m<mStandardHex.mNIntegration2DPad; ++m)
                mElements[0][j][k]->mSurfMetricIntIMax[l][m] = recvBuf[ii++];
            }
          }
        }
      }
    }
  }

  //------------------------------------------------------------
  // Communication of the metric data in j-direction, if needed.
  //------------------------------------------------------------

  // Check if data must be communicated in j-direction.
  if(mNRanksJ > 1)
  {
    // Allocate the memory for the buffers to store the metric data.
    std::vector<su2double> sendBuf(nItemsPerElem*mNElemCommJMax);
    std::vector<su2double> recvBuf(nItemsPerElem*mNElemCommJMin);

    // Copy the metric data of the jMax boundary for the elements on
    // the jMax boundary of the local domain into sendBuf.
    unsigned int ii = 0;
    for(int i=1; i<=mNElemPerRankI; ++i)
    {
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        if( mElements[i][mNElemPerRankJ+1][k] )
        {
          if(mElements[i][mNElemPerRankJ+1][k]->mElemType == HALO_ELEMENT_JMAX)
          {
            sendBuf[ii++] = mElements[i][mNElemPerRankJ][k]->mLenScaleLES;
            sendBuf[ii++] = mElements[i][mNElemPerRankJ][k]->mLenScaleVolume;
            sendBuf[ii++] = mElements[i][mNElemPerRankJ][k]->mLenScaleIDir;
            sendBuf[ii++] = mElements[i][mNElemPerRankJ][k]->mLenScaleJDir;
            sendBuf[ii++] = mElements[i][mNElemPerRankJ][k]->mLenScaleKDir;

            for(unsigned l=0; l<mElements[i][mNElemPerRankJ][k]->mSurfMetricIntJMax.size(); ++l)
              for(int m=0; m<mStandardHex.mNIntegration2DPad; ++m)
                sendBuf[ii++] = mElements[i][mNElemPerRankJ][k]->mSurfMetricIntJMax[l][m];
          }
        }
      }
    }

    // Loop over the number of communication cycles.
    for(int comm=0; comm<nCommCycles; ++comm)
    {
      // Copy the required range from sendBuf into mCommBufJMax.
      unsigned int indBeg = comm*nItemsAlloc*mNElemCommJMax;
      unsigned int indEnd = (comm+1)*nItemsAlloc*mNElemCommJMax;
      if(indEnd > sendBuf.size()) indEnd = sendBuf.size();

      ii = 0;
      for(unsigned int i=indBeg; i<indEnd; ++i, ++ii)
        mCommBufJMax[ii] = sendBuf[i];

      // Start and complete the persistent communication.
      MPI_Startall(mCommRequestSolJDir.size(), mCommRequestSolJDir.data());
      MPI_Waitall(mCommRequestSolJDir.size(), mCommRequestSolJDir.data(),
                  MPI_STATUSES_IGNORE);

      // Copy the required range from mCommBufJMin into recvBuf.
      indBeg = comm*nItemsAlloc*mNElemCommJMin;
      indEnd = (comm+1)*nItemsAlloc*mNElemCommJMin;
      if(indEnd > recvBuf.size()) indEnd = recvBuf.size();

      ii = 0;
      for(unsigned int i=indBeg; i<indEnd; ++i, ++ii)
        recvBuf[i] = mCommBufJMin[ii];
    }

    // Copy the metric data from recvBuf into into the jMax boundary of the
    // elements on the jMin boundary of the local domain. The memory of
    // mSurfMetricIntJMax must be allocated for the halo elements.
    ii = 0;
    for(int i=1; i<=mNElemPerRankI; ++i)
    {
      for(int k=1; k<=mNElemPerRankK; ++k)
      {
        if( mElements[i][0][k] )
        {
          if(mElements[i][0][k]->mElemType == HALO_ELEMENT_JMIN)
          {
            mElements[i][0][k]->mLenScaleLES    = recvBuf[ii++];
            mElements[i][0][k]->mLenScaleVolume = recvBuf[ii++];
            mElements[i][0][k]->mLenScaleIDir   = recvBuf[ii++];
            mElements[i][0][k]->mLenScaleJDir   = recvBuf[ii++];
            mElements[i][0][k]->mLenScaleKDir   = recvBuf[ii++];

            mElements[i][0][k]->mSurfMetricIntJMax.resize(mElements[i][1][k]->mSurfMetricIntJMax.size());
            for(unsigned l=0; l<mElements[i][0][k]->mSurfMetricIntJMax.size(); ++l)
            {
              mElements[i][0][k]->mSurfMetricIntJMax[l] = (su2double *) AllocateMemory(mStandardHex.mNIntegration2DPad*sizeof(su2double));
              if( !(mElements[i][0][k]->mSurfMetricIntJMax[l]) )
                Terminate("SolverClass::CommunicateMetricTermsMaxBoundaries", __FILE__, __LINE__,
                          "Memory allocation failure for mSurfMetricIntJMax");

              for(int m=0; m<mStandardHex.mNIntegration2DPad; ++m)
                mElements[i][0][k]->mSurfMetricIntJMax[l][m] = recvBuf[ii++];
            }
          }
        }
      }
    }
  }

  //------------------------------------------------------------
  // Communication of the metric data in k-direction, if needed.
  //------------------------------------------------------------

  // Check if data must be communicated in k-direction.
  if(mNRanksK > 1)
  {
    // Allocate the memory for the buffers to store the metric data.
    std::vector<su2double> sendBuf(nItemsPerElem*mNElemCommKMax);
    std::vector<su2double> recvBuf(nItemsPerElem*mNElemCommKMin);

    // Copy the metric data of the kMax boundary for the elements on
    // the kMax boundary of the local domain into sendBuf.
    unsigned int ii = 0;
    for(int i=1; i<=mNElemPerRankI; ++i)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        if( mElements[i][j][mNElemPerRankK+1] )
        {
          if(mElements[i][j][mNElemPerRankK+1]->mElemType == HALO_ELEMENT_KMAX)
          {
            sendBuf[ii++] = mElements[i][j][mNElemPerRankK]->mLenScaleLES;
            sendBuf[ii++] = mElements[i][j][mNElemPerRankK]->mLenScaleVolume;
            sendBuf[ii++] = mElements[i][j][mNElemPerRankK]->mLenScaleIDir;
            sendBuf[ii++] = mElements[i][j][mNElemPerRankK]->mLenScaleJDir;
            sendBuf[ii++] = mElements[i][j][mNElemPerRankK]->mLenScaleKDir;

            for(unsigned l=0; l<mElements[i][j][mNElemPerRankK]->mSurfMetricIntKMax.size(); ++l)
              for(int m=0; m<mStandardHex.mNIntegration2DPad; ++m)
                sendBuf[ii++] = mElements[i][j][mNElemPerRankK]->mSurfMetricIntKMax[l][m];
          }
        }
      }
    }

    // Loop over the number of communication cycles.
    for(int comm=0; comm<nCommCycles; ++comm)
    {
      // Copy the required range from sendBuf into mCommBufKMax.
      unsigned int indBeg = comm*nItemsAlloc*mNElemCommKMax;
      unsigned int indEnd = (comm+1)*nItemsAlloc*mNElemCommKMax;
      if(indEnd > sendBuf.size()) indEnd = sendBuf.size();

      ii = 0;
      for(unsigned int i=indBeg; i<indEnd; ++i, ++ii)
        mCommBufKMax[ii] = sendBuf[i];

      // Start and complete the persistent communication.
      MPI_Startall(mCommRequestSolKDir.size(), mCommRequestSolKDir.data());
      MPI_Waitall(mCommRequestSolKDir.size(), mCommRequestSolKDir.data(),
                  MPI_STATUSES_IGNORE);

      // Copy the required range from mCommBufKMin into recvBuf.
      indBeg = comm*nItemsAlloc*mNElemCommKMin;
      indEnd = (comm+1)*nItemsAlloc*mNElemCommKMin;
      if(indEnd > recvBuf.size()) indEnd = recvBuf.size();

      ii = 0;
      for(unsigned int i=indBeg; i<indEnd; ++i, ++ii)
        recvBuf[i] = mCommBufKMin[ii];
    }

    // Copy the metric data from recvBuf into into the kMax boundary of the
    // elements on the kMin boundary of the local domain. The memory of
    // mSurfMetricIntKMax must be allocated for the halo elements.
    ii = 0;
    for(int i=1; i<=mNElemPerRankI; ++i)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        if( mElements[i][j][0] )
        {
          if(mElements[i][j][0]->mElemType == HALO_ELEMENT_KMIN)
          {
            mElements[i][j][0]->mLenScaleLES    = recvBuf[ii++];
            mElements[i][j][0]->mLenScaleVolume = recvBuf[ii++];
            mElements[i][j][0]->mLenScaleIDir   = recvBuf[ii++];
            mElements[i][j][0]->mLenScaleJDir   = recvBuf[ii++];
            mElements[i][j][0]->mLenScaleKDir   = recvBuf[ii++];

            mElements[i][j][0]->mSurfMetricIntKMax.resize(mElements[i][j][1]->mSurfMetricIntKMax.size());
            for(unsigned l=0; l<mElements[i][j][0]->mSurfMetricIntKMax.size(); ++l)
            {
              mElements[i][j][0]->mSurfMetricIntKMax[l] = (su2double *) AllocateMemory(mStandardHex.mNIntegration2DPad*sizeof(su2double));
              if( !(mElements[i][j][0]->mSurfMetricIntKMax[l]) )
                Terminate("SolverClass::CommunicateMetricTermsMaxBoundaries", __FILE__, __LINE__,
                          "Memory allocation failure for mSurfMetricIntKMax");

              for(int m=0; m<mStandardHex.mNIntegration2DPad; ++m)
                mElements[i][j][0]->mSurfMetricIntKMax[l][m] = recvBuf[ii++];
            }
          }
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------

// Function, which completes the communication of the solution in i-direction.
void SolverClass::CompleteCommunicationSolutionIDir(void)
{
  // Easier storage of the number of DOFs.
  const int nDOFs = mStandardHex.mNDOFs;

  // Complete the communication.
  MPI_Waitall(mCommRequestSolIDir.size(), mCommRequestSolIDir.data(),
              MPI_STATUSES_IGNORE);

  // Check if there actually are halo elements on the iMin boundary.  
  if( mNElemCommIMin )
  {
    // Loop over the elements on the iMin boundary and copy the data.
    int ind = 0; 
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        if( mElements[0][j][k] )
        {
          for(int l=0; l<nVar; ++l, ind+=nDOFs)
            std::memcpy(mElements[0][j][k]->mSol[l], mCommBufIMin+ind,
                        nDOFs*sizeof(su2double));
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------

// Function, which completes the communication of the solution in j-direction.
void SolverClass::CompleteCommunicationSolutionJDir(void)
{
  // Easier storage of the number of DOFs.
  const int nDOFs = mStandardHex.mNDOFs;

  // Complete the communication.
  MPI_Waitall(mCommRequestSolJDir.size(), mCommRequestSolJDir.data(),
              MPI_STATUSES_IGNORE);

  // Check if there actually are halo elements on the jMin boundary.  
  if( mNElemCommJMin )
  {
    int ind = 0; 
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int i=1; i<=mNElemPerRankI; ++i)
      {
        if( mElements[i][0][k] )
        {
          for(int l=0; l<nVar; ++l, ind+=nDOFs)
            std::memcpy(mElements[i][0][k]->mSol[l], mCommBufJMin+ind,
                        nDOFs*sizeof(su2double));
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------

// Function, which completes the communication of the solution in k-direction.
void SolverClass::CompleteCommunicationSolutionKDir(void)
{
  // Easier storage of the number of DOFs.
  const int nDOFs = mStandardHex.mNDOFs;

  // Complete the communication.
  MPI_Waitall(mCommRequestSolKDir.size(), mCommRequestSolKDir.data(),
              MPI_STATUSES_IGNORE);

  // Check if there actually are halo elements on the kMin boundary.  
  if( mNElemCommKMin )
  {
    int ind = 0; 
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      for(int i=1; i<=mNElemPerRankI; ++i)
      {
        if( mElements[i][j][0] )
        {
          for(int l=0; l<nVar; ++l, ind+=nDOFs)
            std::memcpy(mElements[i][j][0]->mSol[l], mCommBufKMin+ind,
                        nDOFs*sizeof(su2double));
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------

// Function, which completes the communication of the i-face residuals.
void SolverClass::CompleteCommunicationIFaceResiduals(void)
{
  // Easier storage of the number of DOFs.
  const int nDOFs = mStandardHex.mNDOFs;

  // Complete the communication.
  MPI_Waitall(mCommRequestResIDir.size(), mCommRequestResIDir.data(),
              MPI_STATUSES_IGNORE);

  // Check if there actually are halo elements on the iMax boundary.  
  if( mNElemCommIMax )
  {
    // Loop over the elements on the iMax boundary and copy the data.
    int ind = 0; 
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        if( mElements[mNElemPerRankI+1][j][k] )
        {
          for(int l=0; l<nVar; ++l, ind+=nDOFs)
            std::memcpy(mElements[mNElemPerRankI+1][j][k]->mResIMin[l],
                        mCommBufIMax+ind, nDOFs*sizeof(su2double));
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------

// Function, which completes the communication of the j-face residuals.
void SolverClass::CompleteCommunicationJFaceResiduals(void)
{
  // Easier storage of the number of DOFs.
  const int nDOFs = mStandardHex.mNDOFs;

  // Complete the communication.
  MPI_Waitall(mCommRequestResJDir.size(), mCommRequestResJDir.data(),
              MPI_STATUSES_IGNORE);

  // Check if there actually are halo elements on the jMax boundary.  
  if( mNElemCommJMax )
  {
    // Loop over the elements on the jMax boundary and copy the data.
    int ind = 0; 
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int i=1; i<=mNElemPerRankI; ++i)
      {
        if( mElements[i][mNElemPerRankJ+1][k] )
        {
          for(int l=0; l<nVar; ++l, ind+=nDOFs)
            std::memcpy(mElements[i][mNElemPerRankJ+1][k]->mResJMin[l],
                        mCommBufJMax+ind, nDOFs*sizeof(su2double));
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------

// Function, which completes the communication of the k-face residuals.
void SolverClass::CompleteCommunicationKFaceResiduals(void)
{
  // Easier storage of the number of DOFs.
  const int nDOFs = mStandardHex.mNDOFs;

  // Complete the communication.
  MPI_Waitall(mCommRequestResKDir.size(), mCommRequestResKDir.data(),
              MPI_STATUSES_IGNORE);

  // Check if there actually are halo elements on the kMax boundary.  
  if( mNElemCommKMax )
  {
    // Loop over the elements on the kMax boundary and copy the data.
    int ind = 0; 
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      for(int i=1; i<=mNElemPerRankI; ++i)
      {
        if( mElements[i][j][mNElemPerRankK+1] )
        {
          for(int l=0; l<nVar; ++l, ind+=nDOFs)
            std::memcpy(mElements[i][j][mNElemPerRankK+1]->mResKMin[l],
                        mCommBufKMax+ind, nDOFs*sizeof(su2double));
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------

// Function, which creates the MPI persistent communication pattern.
void SolverClass::CreatePersistentCommPattern(void)
{
  //---------------------------------------------------------------------------
  //   Communication pattern of the solution and residual in i-direction.
  //---------------------------------------------------------------------------

  // Determine the index of the neighboring rank on the iMin and iMax side.
  // Take the possible periodicity in i-direction into account for these neighbors.
  int rankMinNeigh = mMyRankI - 1;
  if(rankMinNeigh == -1) rankMinNeigh = mNRanksI - 1;

  int rankMaxNeigh = mMyRankI + 1;
  if(rankMaxNeigh == mNRanksI) rankMaxNeigh = 0;

  // Determine the actual MPI rank of these neighbors.
  rankMinNeigh += mMyRankJ*mNRanksI + mMyRankK*mNRanksI*mNRanksJ;
  rankMaxNeigh += mMyRankJ*mNRanksI + mMyRankK*mNRanksI*mNRanksJ;

  // Determine the number of elements on the iMin and iMax boundaries for which
  // data must be communicated.
  mNElemCommIMin = mNElemCommIMax = 0;
  for(int j=1; j<=mNElemPerRankJ; ++j)
  {
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      if( mElements[0][j][k] )
        if(mElements[0][j][k]->mElemType == HALO_ELEMENT_IMIN) ++mNElemCommIMin;

      if( mElements[mNElemPerRankI+1][j][k] )
        if(mElements[mNElemPerRankI+1][j][k]->mElemType == HALO_ELEMENT_IMAX) ++mNElemCommIMax;
    }
  }

  // Determine the number of communication requests for the solution and residuals
  // and allocate the memory for the communication requests.
  int nCommRequests = 0;
  if( mNElemCommIMin ) ++nCommRequests;
  if( mNElemCommIMax ) ++nCommRequests;

  mCommRequestSolIDir.resize(nCommRequests);
  mCommRequestResIDir.resize(nCommRequests);

  // Check if the solution must be sent on the iMax boundary. If so,
  // allocate the communication buffer and create the persistent request.
  nCommRequests = 0;
  if( mNElemCommIMax )
  {
    const int sizeComm = nVar*mStandardHex.mNDOFs*mNElemCommIMax;
    mCommBufIMax = (su2double *) AllocateMemory(sizeComm*sizeof(su2double));
    if( !mCommBufIMax )
      Terminate("SolverClass::CreatePersistentCommPattern", __FILE__, __LINE__,
                "Memory allocation failure for mCommBufIMax");

    MPI_Send_init(mCommBufIMax, sizeComm, MPI_SU2DOUBLE, rankMaxNeigh,
                  rank, MPI_COMM_WORLD, &mCommRequestSolIDir[nCommRequests]);
    ++nCommRequests;
  }

  // Check if the solution must be received on the iMin boundary. If so,
  // allocate the communication buffer and create the persistent request.
  if( mNElemCommIMin )
  {
    const int sizeComm = nVar*mStandardHex.mNDOFs*mNElemCommIMin;
    mCommBufIMin = (su2double *) AllocateMemory(sizeComm*sizeof(su2double));
    if( !mCommBufIMin )
      Terminate("SolverClass::CreatePersistentCommPattern", __FILE__, __LINE__,
                "Memory allocation failure for mCommBufIMin");

    MPI_Recv_init(mCommBufIMin, sizeComm, MPI_SU2DOUBLE, rankMinNeigh,
                  rankMinNeigh, MPI_COMM_WORLD, &mCommRequestSolIDir[nCommRequests]);
  }

  // Check if the residuals must be sent on the iMin boundary. If so,
  // create the persistent request.
  nCommRequests = 0;
  if( mNElemCommIMin )
  {
    const int sizeComm = nVar*mStandardHex.mNDOFs*mNElemCommIMin;
    MPI_Send_init(mCommBufIMin, sizeComm, MPI_SU2DOUBLE, rankMinNeigh,
                  rank, MPI_COMM_WORLD, &mCommRequestResIDir[nCommRequests]);
    ++nCommRequests;
  }

  // Check if the residuals must be received on the iMax boundary. If so,
  // create the persistent request.
  if( mNElemCommIMax )
  {
    const int sizeComm = nVar*mStandardHex.mNDOFs*mNElemCommIMax;
    MPI_Recv_init(mCommBufIMax, sizeComm, MPI_SU2DOUBLE, rankMaxNeigh,
                  rankMaxNeigh, MPI_COMM_WORLD, &mCommRequestResIDir[nCommRequests]);
  }

  //---------------------------------------------------------------------------
  //   Communication pattern of the solution and residual in j-direction.
  //---------------------------------------------------------------------------

  // Determine the index of the neighboring rank on the jMin and jMax side.
  // Take the possible periodicity in j-direction into account for these neighbors.
  rankMinNeigh = mMyRankJ - 1;
  if(rankMinNeigh == -1) rankMinNeigh = mNRanksJ - 1;

  rankMaxNeigh = mMyRankJ + 1;
  if(rankMaxNeigh == mNRanksJ) rankMaxNeigh = 0;

  // Determine the actual MPI rank of these neighbors.
  rankMinNeigh = mMyRankI + rankMinNeigh*mNRanksI + mMyRankK*mNRanksI*mNRanksJ;
  rankMaxNeigh = mMyRankI + rankMaxNeigh*mNRanksI + mMyRankK*mNRanksI*mNRanksJ;

  // Determine the number of elements on the jMin and jMax boundaries for which
  // data must be communicated.
  mNElemCommJMin = mNElemCommJMax = 0;
  for(int i=1; i<=mNElemPerRankI; ++i)
  {
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      if( mElements[i][0][k] )
        if(mElements[i][0][k]->mElemType == HALO_ELEMENT_JMIN) ++mNElemCommJMin;

      if( mElements[i][mNElemPerRankJ+1][k] )
        if(mElements[i][mNElemPerRankJ+1][k]->mElemType == HALO_ELEMENT_JMAX) ++mNElemCommJMax;
    }
  }

  // Determine the number of communication requests for the solution and residuals
  // and allocate the memory for the communication requests.
  nCommRequests = 0;
  if( mNElemCommJMin ) ++nCommRequests;
  if( mNElemCommJMax ) ++nCommRequests;

  mCommRequestSolJDir.resize(nCommRequests);
  mCommRequestResJDir.resize(nCommRequests);

  // Check if the solution must be sent on the jMax boundary. If so,
  // allocate the communication buffer and create the persistent request.
  nCommRequests = 0;
  if( mNElemCommJMax )
  {
    const int sizeComm = nVar*mStandardHex.mNDOFs*mNElemCommJMax;
    mCommBufJMax = (su2double *) AllocateMemory(sizeComm*sizeof(su2double));
    if( !mCommBufJMax )
      Terminate("SolverClass::CreatePersistentCommPattern", __FILE__, __LINE__,
                "Memory allocation failure for mCommBufJMax");

    MPI_Send_init(mCommBufJMax, sizeComm, MPI_SU2DOUBLE, rankMaxNeigh,
                  rank, MPI_COMM_WORLD, &mCommRequestSolJDir[nCommRequests]);
    ++nCommRequests;
  }

  // Check if the solution must be received on the jMin boundary. If so,
  // allocate the communication buffer and create the persistent request.
  if( mNElemCommJMin )
  {
    const int sizeComm = nVar*mStandardHex.mNDOFs*mNElemCommJMin;
    mCommBufJMin = (su2double *) AllocateMemory(sizeComm*sizeof(su2double));
    if( !mCommBufJMin )
      Terminate("SolverClass::CreatePersistentCommPattern", __FILE__, __LINE__,
                "Memory allocation failure for mCommBufJMin");

    MPI_Recv_init(mCommBufJMin, sizeComm, MPI_SU2DOUBLE, rankMinNeigh,
                  rankMinNeigh, MPI_COMM_WORLD, &mCommRequestSolJDir[nCommRequests]);
  }

  // Check if the residuals must be sent on the jMin boundary. If so,
  // create the persistent request.
  nCommRequests = 0;
  if( mNElemCommJMin )
  {
    const int sizeComm = nVar*mStandardHex.mNDOFs*mNElemCommJMin;
    MPI_Send_init(mCommBufJMin, sizeComm, MPI_SU2DOUBLE, rankMinNeigh,
                  rank, MPI_COMM_WORLD, &mCommRequestResJDir[nCommRequests]);
    ++nCommRequests;
  }

  // Check if the residuals must be received on the jMax boundary. If so,
  // create the persistent request.
  if( mNElemCommJMax )
  {
    const int sizeComm = nVar*mStandardHex.mNDOFs*mNElemCommJMax;
    MPI_Recv_init(mCommBufJMax, sizeComm, MPI_SU2DOUBLE, rankMaxNeigh,
                  rankMaxNeigh, MPI_COMM_WORLD, &mCommRequestResJDir[nCommRequests]);
  }

  //---------------------------------------------------------------------------
  //   Communication pattern of the solution and residual in k-direction.
  //---------------------------------------------------------------------------

  // Determine the index of the neighboring rank on the kMin and kMax side.
  // Take the possible periodicity in k-direction into account for these neighbors.
  rankMinNeigh = mMyRankK - 1;
  if(rankMinNeigh == -1) rankMinNeigh = mNRanksK - 1;

  rankMaxNeigh = mMyRankK + 1;
  if(rankMaxNeigh == mNRanksK) rankMaxNeigh = 0;

  // Determine the actual MPI rank of these neighbors.
  rankMinNeigh = mMyRankI + mMyRankJ*mNRanksI + rankMinNeigh*mNRanksI*mNRanksJ;
  rankMaxNeigh = mMyRankI + mMyRankJ*mNRanksI + rankMaxNeigh*mNRanksI*mNRanksJ;

  // Determine the number of elements on the kMin and kMax boundaries for which
  // data must be communicated.
  mNElemCommKMin = mNElemCommKMax = 0;
  for(int i=1; i<=mNElemPerRankI; ++i)
  {
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      if( mElements[i][j][0] )
        if(mElements[i][j][0]->mElemType == HALO_ELEMENT_KMIN) ++mNElemCommKMin;

      if( mElements[i][j][mNElemPerRankK+1] )
        if(mElements[i][j][mNElemPerRankK+1]->mElemType == HALO_ELEMENT_KMAX) ++mNElemCommKMax;
    }
  }

  // Determine the number of communication requests for the solution and residuals
  // and allocate the memory for the communication requests.
  nCommRequests = 0;
  if( mNElemCommKMin ) ++nCommRequests;
  if( mNElemCommKMax ) ++nCommRequests;

  mCommRequestSolKDir.resize(nCommRequests);
  mCommRequestResKDir.resize(nCommRequests);

  // Check if the solution must be sent on the kMax boundary. If so,
  // allocate the communication buffer and create the persistent request.
  nCommRequests = 0;
  if( mNElemCommKMax )
  {
    const int sizeComm = nVar*mStandardHex.mNDOFs*mNElemCommKMax;
    mCommBufKMax = (su2double *) AllocateMemory(sizeComm*sizeof(su2double));
    if( !mCommBufKMax )
      Terminate("SolverClass::CreatePersistentCommPattern", __FILE__, __LINE__,
                "Memory allocation failure for mCommBufKMax");

    MPI_Send_init(mCommBufKMax, sizeComm, MPI_SU2DOUBLE, rankMaxNeigh,
                  rank, MPI_COMM_WORLD, &mCommRequestSolKDir[nCommRequests]);
    ++nCommRequests;
  }

  // Check if the solution must be received on the kMin boundary. If so,
  // allocate the communication buffer and create the persistent request.
  if( mNElemCommKMin )
  {
    const int sizeComm = nVar*mStandardHex.mNDOFs*mNElemCommKMin;
    mCommBufKMin = (su2double *) AllocateMemory(sizeComm*sizeof(su2double));
    if( !mCommBufKMin )
      Terminate("SolverClass::CreatePersistentCommPattern", __FILE__, __LINE__,
                "Memory allocation failure for mCommBufKMin");

    MPI_Recv_init(mCommBufKMin, sizeComm, MPI_SU2DOUBLE, rankMinNeigh,
                  rankMinNeigh, MPI_COMM_WORLD, &mCommRequestSolKDir[nCommRequests]);
  }

  // Check if the residuals must be sent on the kMin boundary. If so,
  // create the persistent request.
  nCommRequests = 0;
  if( mNElemCommKMin )
  {
    const int sizeComm = nVar*mStandardHex.mNDOFs*mNElemCommKMin;
    MPI_Send_init(mCommBufKMin, sizeComm, MPI_SU2DOUBLE, rankMinNeigh,
                  rank, MPI_COMM_WORLD, &mCommRequestResKDir[nCommRequests]);
    ++nCommRequests;
  }

  // Check if the residuals must be received on the kMax boundary. If so,
  // create the persistent request.
  if( mNElemCommKMax )
  {
    const int sizeComm = nVar*mStandardHex.mNDOFs*mNElemCommKMax;
    MPI_Recv_init(mCommBufKMax, sizeComm, MPI_SU2DOUBLE, rankMaxNeigh,
                  rankMaxNeigh, MPI_COMM_WORLD, &mCommRequestResKDir[nCommRequests]);
  }
}

//------------------------------------------------------------------------------

// Function, which releases the memory of the persistent communication requests.
void SolverClass::FreeCommRequests(void)
{
  for(unsigned long i=0; i<mCommRequestSolIDir.size(); ++i)
    MPI_Request_free(&mCommRequestSolIDir[i]);

  for(unsigned long i=0; i<mCommRequestSolJDir.size(); ++i)
    MPI_Request_free(&mCommRequestSolJDir[i]);

  for(unsigned long i=0; i<mCommRequestSolKDir.size(); ++i)
    MPI_Request_free(&mCommRequestSolKDir[i]);

  for(unsigned long i=0; i<mCommRequestResIDir.size(); ++i)
    MPI_Request_free(&mCommRequestResIDir[i]);

  for(unsigned long i=0; i<mCommRequestResJDir.size(); ++i)
    MPI_Request_free(&mCommRequestResJDir[i]);

  for(unsigned long i=0; i<mCommRequestResKDir.size(); ++i)
    MPI_Request_free(&mCommRequestResKDir[i]);
}

//------------------------------------------------------------------------------

// Function, which starts the communication of the solution in i-direction.
void SolverClass::StartCommunicationSolutionIDir(void)
{
  // Easier storage of the number of DOFs.
  const int nDOFs = mStandardHex.mNDOFs;

  // Copy the solution of the elements on the iMax boundary,
  // which are adjacent to a halo element, to mCommBufIMax.
  if( mNElemCommIMax )
  {
    int ind = 0; 
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        if( mElements[mNElemPerRankI+1][j][k] )
        {
          for(int l=0; l<nVar; ++l, ind+=nDOFs)
            std::memcpy(mCommBufIMax+ind, mElements[mNElemPerRankI][j][k]->mSol[l],
                        nDOFs*sizeof(su2double));
        }
      }
    }
  }

  // Start the communication.
  MPI_Startall(mCommRequestSolIDir.size(), mCommRequestSolIDir.data());
}

//------------------------------------------------------------------------------

// Function, which starts the communication of the solution in j-direction.
void SolverClass::StartCommunicationSolutionJDir(void)
{
  // Easier storage of the number of DOFs.
  const int nDOFs = mStandardHex.mNDOFs;

  // Copy the solution of the elements on the jMax boundary,
  // which are adjacent to a halo element, to mCommBufJMax.
  if( mNElemCommJMax )
  {
    int ind = 0; 
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int i=1; i<=mNElemPerRankI; ++i)
      {
        if( mElements[i][mNElemPerRankJ+1][k] )
        {
          for(int l=0; l<nVar; ++l, ind+=nDOFs)
            std::memcpy(mCommBufJMax+ind, mElements[i][mNElemPerRankJ][k]->mSol[l],
                        nDOFs*sizeof(su2double));
        }
      }
    }
  }

  // Start the communication.
  MPI_Startall(mCommRequestSolJDir.size(), mCommRequestSolJDir.data());
}

//------------------------------------------------------------------------------

// Function, which starts the communication of the solution in k-direction.
void SolverClass::StartCommunicationSolutionKDir(void)
{
  // Easier storage of the number of DOFs.
  const int nDOFs = mStandardHex.mNDOFs;

  // Copy the solution of the elements on the kMax boundary,
  // which are adjacent to a halo element, to mCommBufKMax.
  if( mNElemCommKMax )
  {
    int ind = 0; 
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      for(int i=1; i<=mNElemPerRankI; ++i)
      {
        if( mElements[i][j][mNElemPerRankK+1] )
        {
          for(int l=0; l<nVar; ++l, ind+=nDOFs)
            std::memcpy(mCommBufKMax+ind, mElements[i][j][mNElemPerRankK]->mSol[l],
                        nDOFs*sizeof(su2double));
        }
      }
    }
  }

  // Start the communication.
  MPI_Startall(mCommRequestSolKDir.size(), mCommRequestSolKDir.data());
}

//------------------------------------------------------------------------------

// Function, which starts the communication of the i-face residuals.
void SolverClass::StartCommunicationIFaceResiduals(void)
{
  // Easier storage of the number of DOFs.
  const int nDOFs = mStandardHex.mNDOFs;

  // Copy the residuals of the elements on the iMin boundary,
  // which are adjacent to a halo element, to mCommBufIMin.
  if( mNElemCommIMin )
  {
    int ind = 0; 
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int j=1; j<=mNElemPerRankJ; ++j)
      {
        if( mElements[0][j][k] )
        {
          for(int l=0; l<nVar; ++l, ind+=nDOFs)
            std::memcpy(mCommBufIMin+ind, mElements[1][j][k]->mResIMin[l],
                        nDOFs*sizeof(su2double));
        }
      }
    }
  }

  // Start the communication.
  MPI_Startall(mCommRequestResIDir.size(), mCommRequestResIDir.data());
}

//------------------------------------------------------------------------------

// Function, which starts the communication of the j-face residuals.
void SolverClass::StartCommunicationJFaceResiduals(void)
{
  // Easier storage of the number of DOFs.
  const int nDOFs = mStandardHex.mNDOFs;

  // Copy the residuals of the elements on the jMin boundary,
  // which are adjacent to a halo element, to mCommBufJMin.
  if( mNElemCommJMin )
  {
    int ind = 0; 
    for(int k=1; k<=mNElemPerRankK; ++k)
    {
      for(int i=1; i<=mNElemPerRankI; ++i)
      {
        if( mElements[i][0][k] )
        {
          for(int l=0; l<nVar; ++l, ind+=nDOFs)
            std::memcpy(mCommBufJMin+ind, mElements[i][1][k]->mResJMin[l],
                        nDOFs*sizeof(su2double));
        }
      }
    }
  }

  // Start the communication.
  MPI_Startall(mCommRequestResJDir.size(), mCommRequestResJDir.data());
}

//------------------------------------------------------------------------------

// Function, which starts the communication of the k-face residuals.
void SolverClass::StartCommunicationKFaceResiduals(void)
{
  // Easier storage of the number of DOFs.
  const int nDOFs = mStandardHex.mNDOFs;

  // Copy the residuals of the elements on the kMin boundary,
  // which are adjacent to a halo element, to mCommBufKMin.
  if( mNElemCommKMin )
  {
    int ind = 0; 
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      for(int i=1; i<=mNElemPerRankI; ++i)
      {
        if( mElements[i][j][0] )
        {
          for(int l=0; l<nVar; ++l, ind+=nDOFs)
            std::memcpy(mCommBufKMin+ind, mElements[i][j][1]->mResKMin[l],
                        nDOFs*sizeof(su2double));
        }
      }
    }
  }

  // Start the communication.
  MPI_Startall(mCommRequestResKDir.size(), mCommRequestResKDir.data());
}

#endif  // HAVE_MPI
