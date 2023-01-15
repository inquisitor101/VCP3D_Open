//------------------------------------------------------------------------------
// SolverClass.hpp: Definition of the class, which contains the variables and
//                  functions to solve the 3D Navier-Stokes equations for a
//                  single block grid. This include file should not be included
//                  by any source code file directly. These files should just
//                  include VCP3D_CurviLinear.hpp.
//------------------------------------------------------------------------------

#pragma once

// Include files needed.
#include "StandardElementClass.hpp"
#include "ElementClass.hpp"

// Definition of the class SolverClass.
class SolverClass
{
public:
  //------------------------------------------
  // Constructor and destructor.
  //------------------------------------------

  // Constructor. Initialize the member variables.
  SolverClass(InputParamClass *inputParam);

  // Destructor. Release the memory.
  ~SolverClass();

  //--------------------------------------
  // Public member variables.
  //--------------------------------------

  // Pointer to the object, that contains the input parameters.
  InputParamClass *mInputParam;

  //--------------------------------------
  // Public member functions.
  //--------------------------------------

  // Function, which adds fluctuations to the initial solution, if desired.
  void AddFluctuations(void);

  // Function, which determines the bounding box for the fluctuations.
  void BoundingBoxFluctuations(void);

  // Function, which computes the metric terms.
  void ComputeMetricTerms(void);

  // Function, which carries out one time synchronization step.
  void ExecuteTimeSynchronizationStep(su2double *monitoringData);

  // Function, which initializes the solution.
  void InitSol(void);

	// Function, which initializes the NSCBC parameters.
	void InitNSCBC(void);

	// Function, which computes the average data required for the NSCBC.
	void AverageDataNSCBC(void);

  // Function, which determines the prescribed data in the integration
  // points of the boundary faces.
  void PrescribedDataIntegrationPoints(void);

  // Function to read the grid.
  void ReadGrid(void);

  // Function to read the restart solution.
  void ReadRestartSolution(void);

  // Function, which computes the spatial residual in the DOFs of the owned elements.
  void Residual(su2double  *monitoringData,
                const bool ComputeMonitoringData);

  // Function, which updates the data for the averages.
  void UpdateAverageData(const su2double *monitoringData);

  // Function to write the binary paraview solution file.
  void WriteParaviewSolution(const int timeStep);

  // Function to write the restart solution.
  void WriteRestartSolution(const int  timeStep,
                            const bool writeRes = false);

  // Function, which writes the space time average data if averaging is applied.
  void WriteSpaceTimeAverageData(void);

  // Function, which writes the grid of the solution DOFs in SU2 format.
  void WriteSU2GridFile(void);

  // Function, which writes the time step data to std::cout.
  void WriteTimeStepData(const int       timeStep,
                         const double    time0,
                         const su2double *monitoringData);

#ifdef HAVE_MPI
  // Function, which releases the memory of the persistent communication requests.
  void FreeCommRequests(void);
#endif

private:
  //--------------------------------------
  // Private member variables.
  //--------------------------------------

  // Triple vector of the pointers to the elements. Note that a halo layer
  // is stored as well, so the owned elements start at index 1.
  std::vector<std::vector< std::vector<ElementClass *> > > mElements;

  // Rank index in i-direction of the this rank.
  int mMyRankI;

  // Rank index in j-direction of the this rank.
  int mMyRankJ;

  // Rank index in k-direction of the this rank.
  int mMyRankK;

  // Number of elements per rank in i-direction.
  int mNElemPerRankI;

  // Number of elements per rank in j-direction.
  int mNElemPerRankJ;

  // Number of elements per rank in k-direction.
  int mNElemPerRankK;

  // Number of ranks in i-direction.
  int mNRanksI;

  // Number of ranks in j-direction.
  int mNRanksJ;

  // Number of ranks in k-direction.
  int mNRanksK;

  // Counter, which keeps track of the number of time synchronization
  // steps used in the averaging. It is an integer to be consistent
  // with the format of the SU2 restart file.
  int mNTimeStepAverage;

  // Counter, which keeps track of the time step.
  unsigned long mCurrentTimeStep;

  // Time averaged force coefficient in x-direction.
  su2double mCfxAve;

  // Time averaged force coefficient in y-direction.
  su2double mCfyAve;

  // Time averaged force coefficient in z-direction.
  su2double mCfzAve;

  // Time averaged viscous force coefficient in x-direction.
  su2double mCfxVisAve;

  // Time averaged viscous force coefficient in y-direction.
  su2double mCfyVisAve;

  // Time averaged viscous force coefficient in z-direction.
  su2double mCfzVisAve;

  // Non-dimensional primitive variables of the free stream.
  su2double mNonDimPrimVarFreeStream[nVar];

	// Surface area of (physical) boundary faces.
	su2double mSurfAreaBoundary[6];

  // Standard hexahedron. 
  StandardElementClass mStandardHex;

#ifdef HAVE_MPI
  // The number of elements that must be communicated on the iMin boundary.
  int mNElemCommIMin;

  // The number of elements that must be communicated on the iMax boundary.
  int mNElemCommIMax;

  // The number of elements that must be communicated on the jMin boundary.
  int mNElemCommJMin;

  // The number of elements that must be communicated on the jMax boundary.
  int mNElemCommJMax;

  // The number of elements that must be communicated on the kMin boundary.
  int mNElemCommKMin;

  // The number of elements that must be communicated on the kMax boundary.
  int mNElemCommKMax;

  // The communication requests for the solution in i-direction.
  std::vector<MPI_Request> mCommRequestSolIDir;

  // The communication requests for the solution in j-direction.
  std::vector<MPI_Request> mCommRequestSolJDir;

  // The communication requests for the solution in k-direction.
  std::vector<MPI_Request> mCommRequestSolKDir;

  // The communication requests for the residuals in i-direction.
  std::vector<MPI_Request> mCommRequestResIDir;

  // The communication requests for the residuals in j-direction.
  std::vector<MPI_Request> mCommRequestResJDir;

  // The communication requests for the residuals in k-direction.
  std::vector<MPI_Request> mCommRequestResKDir;

  // MPI communication buffer for the iMin boundary.
  su2double *mCommBufIMin;

  // MPI communication buffer for the iMax boundary.
  su2double *mCommBufIMax;

  // MPI communication buffer for the jMin boundary.
  su2double *mCommBufJMin;

  // MPI communication buffer for the jMax boundary.
  su2double *mCommBufJMax;

  // MPI communication buffer for the kMin boundary.
  su2double *mCommBufKMin;

  // MPI communication buffer for the kMax boundary.
  su2double *mCommBufKMax;
#endif

  //------------------------------------------------
  // Private member functions.
  //------------------------------------------------

#ifdef HAVE_MPI
  // Function, which communicates the metric terms of the max boundaries.
  void CommunicateMetricTermsMaxBoundaries(void);

  // Function, which completes the communication of the solution in i-direction.
  void CompleteCommunicationSolutionIDir(void);

  // Function, which completes the communication of the solution in j-direction.
  void CompleteCommunicationSolutionJDir(void);

  // Function, which completes the communication of the solution in k-direction.
  void CompleteCommunicationSolutionKDir(void);

  // Function, which completes the communication of the i-face residuals.
  void CompleteCommunicationIFaceResiduals(void);

  // Function, which completes the communication of the j-face residuals.
  void CompleteCommunicationJFaceResiduals(void);

  // Function, which completes the communication of the k-face residuals.
  void CompleteCommunicationKFaceResiduals(void);

  // Function, which creates the MPI persistent communication pattern.
  void CreatePersistentCommPattern(void);

  // Function, which starts the communication of the solution in i-direction.
  void StartCommunicationSolutionIDir(void);

  // Function, which starts the communication of the solution in j-direction.
  void StartCommunicationSolutionJDir(void);

  // Function, which starts the communication of the solution in k-direction.
  void StartCommunicationSolutionKDir(void);

  // Function, which starts the communication of the i-face residuals.
  void StartCommunicationIFaceResiduals(void);

  // Function, which starts the communication of the j-face residuals.
  void StartCommunicationJFaceResiduals(void);

  // Function, which starts the communication of the k-face residuals.
  void StartCommunicationKFaceResiduals(void);
#endif

  // Function, which controls the allocation of the residual and halo solution.
  void AllocateResidualAndHaloSolutionArrays(void);

  // Function, which checks if the coordinates of the given faces match.
  void CheckMatchingFaces(const int                    rankMin,
                          const int                    rankMax,
                          const int                    nElemMin,
                          const int                    nElemMax,
                          const std::vector<su2double> &coorMin,
                          const std::vector<su2double> &coorMinComm,
                          const std::vector<su2double> &coorMax,
                          const std::vector<su2double> &coorMaxComm,
                          const char                   *faceDir);

  // Function, which checks if matching subfaces actually match.
  void CheckMatchingSubfaces(void);

  // Function, which determines the precision of the restart file.
  short DeterminePrecisionRestartFile(void);

  // Function, which determines the visualization data.
  void DetermineVisualizationData(std::vector<std::string> &visNames,
                                  std::vector<float>       &visBuf);

  // Function, which distributes the grid over the MPI ranks.
  void DistributeGrid(void);

  // Function, which determines the interpolation weights for the
  // exchange location of the wall model, if needed.
  void InterpolationWeightsExchangeLocation(void);

#ifdef HAVE_MPI
  // Function to read the grid coordinates from the grid file. It is a
  // template function, such that different precisions can be handled.
  template <class T>
  void ReadGridCoordinates(MPI_File               &fh,
                           std::vector<su2double> &coorBuf,
                           const bool             byteSwap,
                           const MPI_Offset       sizeHeader,
                           const int              ilG,
                           const int              jlG,
                           const int              klG);

  // Function to read the variables from the restart file. It is a template
  // function, such that different precisions can be handled.
  template <class T>
  void ReadSolutionFromRestartFile(MPI_File         &fh,
                                   const int        nVarRestart,
                                   const bool       byteSwap,
                                   const MPI_Offset sizeHeader);
#else
  // Function to read the grid coordinates from the grid file. It is a
  // template function, such that different precisions can be handled.
  template <class T>
  void ReadGridCoordinates(FILE                   *&gridFile,
                           std::vector<su2double> &coorBuf,
                           const bool             byteSwap,
                           const int              ilG,
                           const int              jlG,
                           const int              klG);

  // Function to read the variables from the restart file. It is a template
  // function, such that different precisions can be handled.
  template <class T>
  void ReadSolutionFromRestartFile(FILE       *&restartFile,
                                   const int  nVarRestart,
                                   const bool byteSwap);
#endif

  //------------------------------------------------
  // Disabled constructors and assignment operators.
  //------------------------------------------------

  // Default constructor, disabled.
  SolverClass();

  // Copy constructor, disabled.
  SolverClass(const SolverClass&);

  // Assignment operator, disabled.
  SolverClass& operator=(const SolverClass&);
};

//------------------------------------------------------------------------------

#ifdef HAVE_MPI
// Implementation of the template function ReadGridCoordinates.
// This function reads the grid using the parallel MPI IO functionality.
template <class T>
void SolverClass::ReadGridCoordinates(MPI_File               &fh,
                                      std::vector<su2double> &coorBuf,
                                      const bool             byteSwap,
                                      const MPI_Offset       sizeHeader,
                                      const int              ilG,
                                      const int              jlG,
                                      const int              klG)
{
  // Determine the corresponding MPI datatype.
  T dummyT = 0;  // It is initialized to avoid a compiler warning.
  MPI_Datatype dataType = DetermineMPI_Type(dummyT);

  // Determine the local number of grid points in the three directions.
  const int il = mNElemPerRankI*mInputParam->mNPolyGridDOFs + 1;
  const int jl = mNElemPerRankJ*mInputParam->mNPolyGridDOFs + 1;
  const int kl = mNElemPerRankK*mInputParam->mNPolyGridDOFs + 1;

  // Define the variables to create the subarray data type for the actual IO.
  const int gsizes[]   = {ilG, jlG, klG, 3};
  const int lsizes[]   = {il,  jl,  kl,  3};
  const int startInd[] = {mMyRankI*mNElemPerRankI*mInputParam->mNPolyGridDOFs,
                          mMyRankJ*mNElemPerRankJ*mInputParam->mNPolyGridDOFs,
                          mMyRankK*mNElemPerRankK*mInputParam->mNPolyGridDOFs, 0};

  // Create the derived datatype for the reading.
  MPI_Datatype readType;
  MPI_Type_create_subarray(4, gsizes, lsizes, startInd, MPI_ORDER_FORTRAN,
                           dataType, &readType);
  MPI_Type_commit(&readType);

  // Create the file view needed for the collective read. The offset is set
  // to the position where the data of the element starts;
  char datarep[] = "native";
  MPI_File_set_view(fh, sizeHeader, dataType, readType, datarep, MPI_INFO_NULL);

  // Allocate the memory for the read buffer.
  std::vector<T> readBuf(3*il*jl*kl);

  // Read the data using one collective read.
  if(MPI_File_read_all(fh, readBuf.data(), readBuf.size(),
                       dataType, MPI_STATUS_IGNORE) != MPI_SUCCESS)
    TerminateAll("SolverClass::ReadGridCoordinates", __FILE__, __LINE__,
                 "MPI_File_read_all was not able to read data from file.");

  // Carry out the byte swapping, if needed.
  if( byteSwap ) SwapBytes(readBuf.data(), sizeof(T), readBuf.size());

   // Release the memory of the derived data type used for the reading.
  MPI_Type_free(&readType);

  // Copy the data from the readBuf into coorBuf.
  coorBuf.resize(readBuf.size());
  for(unsigned int i=0; i<readBuf.size(); ++i)
    coorBuf[i] = (su2double) readBuf[i];
}

//------------------------------------------------------------------------------

// Implementation of the template function ReadSolutionFromRestartFile.
// This function reads the solution using the parallel MPI IO functionality.
template <class T>
void SolverClass::ReadSolutionFromRestartFile(MPI_File         &fh,
                                              const int        nVarRestart,
                                              const bool       byteSwap,
                                              const MPI_Offset sizeHeader)
{
  // Easier storage of the total number of elements in the three index directions.
  const int nElemI = mInputParam->mSubfaces[1][0]->mElemIBeg;
  const int nElemJ = mInputParam->mSubfaces[3][0]->mElemJBeg;
  const int nElemK = mInputParam->mSubfaces[5][0]->mElemKBeg;

  // Determine the number of variables to read from the restart file.
  int nVarRead = nVar;
  if( mInputParam->mContinueAveragingFromRestart ) nVarRead += nVar + 7;

  // Check if enough variables are present in the restart file. The additional 3
  // comes from the coordinates that are also stored in the restart file.
  if(nVarRestart < (nVarRead+3))
    TerminateAll("SolverClass::ReadSolutionFromRestartFile", __FILE__, __LINE__,
                 "Not enough variables present in the restart file.");

  // Determine the corresponding MPI datatype.
  T dummyT = 0;  // It is initialized to avoid a compiler warning.
  MPI_Datatype dataType = DetermineMPI_Type(dummyT);

  // Define the variables to create the subarray data type for the actual IO.
  const int gsizes[]   = {nVarRestart, (int) mStandardHex.mNDOFs, nElemI, nElemJ, nElemK};
  const int lsizes[]   = {nVarRead, (int) mStandardHex.mNDOFs, mNElemPerRankI, mNElemPerRankJ,
                          mNElemPerRankK};
  const int startInd[] = {3, 0, mMyRankI*mNElemPerRankI, mMyRankJ*mNElemPerRankJ,
                          mMyRankK*mNElemPerRankK};

  // Create the derived datatype for the reading.
  MPI_Datatype readType;
  MPI_Type_create_subarray(5, gsizes, lsizes, startInd, MPI_ORDER_FORTRAN,
                           dataType, &readType);
  MPI_Type_commit(&readType);

  // Create the file view needed for the collective read. The offset is set
  // to the position where the data of the element starts;
  char datarep[] = "native";
  MPI_File_set_view(fh, sizeHeader, dataType, readType, datarep, MPI_INFO_NULL);

  // Allocate the memory for the read buffer.
  std::vector<T> readBuf(lsizes[0]*lsizes[1]*lsizes[2]*lsizes[3]*lsizes[4]);

  // Read the data using one collective read.
  if(MPI_File_read_all(fh, readBuf.data(), readBuf.size(),
                       dataType, MPI_STATUS_IGNORE) != MPI_SUCCESS)
    TerminateAll("SolverClass::ReadSolutionFromRestartFile", __FILE__, __LINE__,
                 "MPI_File_read_all was not able to read data from file.");

  // Carry out the byte swapping, if needed.
  if( byteSwap ) SwapBytes(readBuf.data(), sizeof(T), readBuf.size());

  // Loop over the elements and copy the data from the read buffer into
  // the data structures of the element.
  long ii = 0;
  for(int k=1; k<=mNElemPerRankK; ++k)
  {
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      for(int i=1; i<=mNElemPerRankI; ++i)
      {
        for(int m=0; m<mStandardHex.mNDOFs; ++m)
        {
          for(int l=0; l<nVar; ++l, ++ii)
            mElements[i][j][k]->mSol[l][m] = (su2double) readBuf[ii];

          if( mInputParam->mContinueAveragingFromRestart )
          {
            for(int l=0; l<nVar; ++l, ++ii)
              mElements[i][j][k]->mAvePrim[l][m] = (su2double) readBuf[ii];

            for(int l=0; l<6; ++l, ++ii)
              mElements[i][j][k]->mAveVelProd[l][m] = (su2double) readBuf[ii];

            mElements[i][j][k]->mAveEddyVis[m] = (su2double) readBuf[ii]; ++ii;
          }
        }
      }
    }
  }

  // Release the memory of the derived data type used for the reading.
  MPI_Type_free(&readType);
}

#else
// Implementation of the template function ReadGridCoordinates.
// This function reads the coordinates from the current position in the grid file.
template <class T>
void SolverClass::ReadGridCoordinates(FILE                   *&gridFile,
                                      std::vector<su2double> &coorBuf,
                                      const bool             byteSwap,
                                      const int              ilG,
                                      const int              jlG,
                                      const int              klG)
{
  // Allocate the memory for the read buffer.
  std::vector<T> readBuf(3*ilG*jlG*klG);

  // Read the data from file.
  if(std::fread(readBuf.data(), sizeof(T), readBuf.size(), gridFile) != readBuf.size())
    Terminate("SolverClass::ReadGridCoordinates", __FILE__, __LINE__,
              "Coordinates could not be read");

  // Carry out the byte swapping, if needed.
  if( byteSwap ) SwapBytes(readBuf.data(), sizeof(T), readBuf.size());

  // Copy the data from the readBuf into coorBuf.
  coorBuf.resize(readBuf.size());
  for(unsigned int i=0; i<readBuf.size(); ++i)
    coorBuf[i] = (su2double) readBuf[i];
}

//------------------------------------------------------------------------------

// Implementation of the template function ReadSolutionFromRestartFile.
// This function reads the solution from the current position in the restart file.
template <class T>
void SolverClass::ReadSolutionFromRestartFile(FILE       *&restartFile,
                                              const int  nVarRestart,
                                              const bool byteSwap)
{
  // Determine the number of variables to read from the restart file.
  int nVarRead = nVar;
  if( mInputParam->mContinueAveragingFromRestart ) nVarRead += nVar + 7;

  // Check if enough variables are present in the restart file. The additional 3
  // comes from the coordinates that are also stored in the restart file.
  if(nVarRestart < (nVarRead+3))
    TerminateAll("SolverClass::ReadSolutionFromRestartFile", __FILE__, __LINE__,
                 "Not enough variables present in the restart file.");

  // Allocate the memory for the read buffer.
  const int sizeRead = mStandardHex.mNDOFs*nVarRestart;
  std::vector<T> readBuf(sizeRead);

  // Loop over the owned element range.
  for(int k=1; k<=mNElemPerRankK; ++k)
  {
    for(int j=1; j<=mNElemPerRankJ; ++j)
    {
      for(int i=1; i<=mNElemPerRankI; ++i)
      {
        // Read the buffer from file and check if it went OK.
        if(std::fread(readBuf.data(), sizeof(T), sizeRead, restartFile) != (unsigned int) sizeRead)
          Terminate("SolverClass::ReadSolutionFromRestartFile", __FILE__, __LINE__,
                    "Solution could not be read");

        // Carry out the byte swapping, if needed.
        if( byteSwap ) SwapBytes(readBuf.data(), sizeof(T), sizeRead);

        // Loop over the DOFs.
        for(int m=0; m<mStandardHex.mNDOFs; ++m)
        {
          // Set the pointer in readBuf where the solution data of this DOF starts.
          // The additional 3 is there to skip the three coordinates.
          T *buf = readBuf.data() + m*nVarRestart + 3;

          // Copy the data from the read buffer into the data structures of the element.
          int ii = 0; 
          for(int l=0; l<nVar; ++l, ++ii)
            mElements[i][j][k]->mSol[l][m] = (su2double) buf[ii];

          if( mInputParam->mContinueAveragingFromRestart )
          {
            for(int l=0; l<nVar; ++l, ++ii)
              mElements[i][j][k]->mAvePrim[l][m] = (su2double) buf[ii];

            for(int l=0; l<6; ++l, ++ii)
              mElements[i][j][k]->mAveVelProd[l][m] = (su2double) buf[ii];

            mElements[i][j][k]->mAveEddyVis[m] = (su2double) buf[ii];
          }
        }
      }
    }
  }
}

#endif
