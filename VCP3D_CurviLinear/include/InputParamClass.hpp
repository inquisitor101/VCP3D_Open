//------------------------------------------------------------------------------
// InputParamClass.hpp: Definition of the class, which contains the input
//                      parameters. This include file should not be included by
//                      any source code file directly. These files should just
//                      include VCP3D_CurviLinear.hpp.
//------------------------------------------------------------------------------

#pragma once

// Include the classes, which contain the information of the boundary subfaces.
#include "SubfaceClasses.hpp"

// Forward definition of CSGSModel and CWallModel.
class CSGSModel;
class CWallModel;

// Definition of the class InputParamClass.
class InputParamClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

  // Constructor. Set the default values.
  InputParamClass();

  // Destructor. Release the memory, if needed.
  ~InputParamClass();

  //------------------------------------------
  // Public member functions.
  //------------------------------------------

  // Function, which computes the reference conditions.
  void ComputeReferenceConditions(void);

  // Function, which determines the non-dimensional primitive variables
  // of the free-stream.
  void DetermineNonDimPrimVarFreeStream(su2double *primVar) const;

  // Function to read the prescribed boundary data.
  void ReadBoundaryData(const char *fileName);

  // Function to read the input parameters from the given file.
  void ReadInputParameters(const char *fileName);

  // Function to read the subface information.
  void ReadSubfaceInformation(const char *fileName);

  // Function, which writes the header for the time step history.
  void WriteHeaderTimeStepHistory(void);

  // Function to read the sponge layer data.
  void ReadSpongeLayerData(const char *fileName);


  //------------------------------------------
  // Public member variables.
  //------------------------------------------

  // Double vector of pointers containing the subface information.
  std::vector<std::vector<SubfaceBaseClass *> > mSubfaces;

  // Location of the DOFs inside the hexahedral elements.
  // Used for visualization and averaging, not for the actual discretization.
  ENUM_DOF_LOCATION mDOFLocation;

  // Set of variables that vary polynomially inside an element.
  ENUM_FEM_VARIABLES mFEMVariables;

  // Riemann solver for the inviscid fluxes.
  ENUM_RIEMANN mRiemannSolver;

  // What subgrid scale model to use.
  ENUM_SGS_MODEL mSGSModelType;

  // What type of fluctuations to add to the initial solution.
  ENUM_FLUCTUATIONS mFluctuationsInitSol;

  // Polynomial degree of the grid.
  int mNPolyGridDOFs;

  // Polynomial degree of the solution.
  int mNPolySolDOFs;

  // Polynomial degree of the integral rule.
  int mNPolyIntRule;

	// Mach number.
  su2double mMach;

  // Reynolds number.
  su2double mReynolds;

  // Length used for the Reynolds number.
  su2double mReynoldsLength;

  // Free stream temperature (only for external flows).
  su2double mTempFreeStream;

  // Free stream velocity direction (only for external flows).
  su2double mVelDirFreeStream[3];

  // Lift direction (only for external flows).
  su2double mLiftDirection[3];

  // Body force per unit volume.
  su2double mFBody[3];

  // Whether or not a body force has been specified.
  bool mFBodySpecified;

  // Lower coordinates of the box to which fluctuations may be added.
  su2double mLowerCoorBoxFluctuations[3];

  // Upper coordinates of the box to which fluctuations may be added.
  su2double mUpperCoorBoxFluctuations[3];

	// Pulse center.
	su2double mPulseCenter[3];
	// Pulse amplitude strength.
	su2double mPulseStrength;
	// Pulse radial width.
	su2double mPulseWidth;
	// Whether it is a 1D pulse.
	bool      mPulse1D;
	// Direction of the 1D (x-)pulse.
	su2double mPulseDirection1D;
	// Quiescent flow simulation (all outlets, M = 0).
	bool mQuiescentFlow;

  // Name of the grid file.
  std::string mGridFile;

  // Name of the solution file.
  std::string mSolFile;

  // Name of the restart file.
  std::string mRestartFile;

  // Name of the paraview visualization file.
  std::string mParaviewFile;

  // Name of the file for the spatially averaged data, if appropriate.
  std::string mSpaceAverageFile;

  // Name of the SU2 grid file.
  std::string mSU2GridFile;

  // Whether or not a restart must be carried out.
  bool mRestart;

  // Whether or not to write an SU2 grid file.
  bool mWriteSU2GridFile;

  // Whether or not the average solution must be computed.
  bool mComputeAverageSolution;

  // Whether or not the averaging must be continued from the restart.
  bool mContinueAveragingFromRestart;

  // Whether or not splitting (for MPI) in i-direction is allowed.
  bool mAllowSplitIDir;

  // Whether or not splitting (for MPI) in j-direction is allowed.
  bool mAllowSplitJDir;

  // Whether or not splitting (for MPI) in k-direction is allowed.
  bool mAllowSplitKDir;

  // Whether or not the forces must be monitored.
  bool mMonitorForces;

  // Reference area for the force coefficients.
  su2double mReferenceArea;

  // Maximum number of synchronization time steps.
  int mNTimeSynchrMax;

  // Saving frequency of the solution.
  int mNSave;

  // Writing frequency of the time step information.
  int mNWriteTimeHistory;

  // Synchronization time step. Note that multiple time steps can be carried
  // out within a synchronization time step.
  su2double mDeltaTSynchr;

  // CFL number to be used for global time stepping.
  su2double mCFL;

  // Parameter for the symmetrizing terms.
  su2double mThetaSym;

  // Relative penalty parameter.
  su2double mRelativePenaltyParameter;

  // Pointer to the SGS model object.
  CSGSModel *mSGSModel;

  // -----------Wall model info--------------

  // Boundary treatment for the viscous wall boundaries.
  ENUM_WALL_MODEL mWallModelType;

  // Pointer to Wall Model object
  CWallModel *mWallModel;

  // Number of points in the equilibrium wall model grid.
  int mNPointsGridWallModel;

  // Exchange location for the wall model.
  su2double mExchangeLocation;

  // Expansion ratio grid equilibrium wall model.
  su2double mExpansionRatioWallModel;

	// Option whether or not an NSCBC is specified.
	bool mNSCBC_Specified;

	// Type of transverse relaxation implementation: (1-beta) or beta.
	// Default is true, i.e. use: (1-beta) for both beta_l and beta_t.
	bool mNSCBC_Outlet_OneMinusBeta;
	// Type of averaged Mach number used: local(element) or global(boundary). 
	// Default is false, i.e. use: average over entire boundary per block.
	bool mNSCBC_Outlet_AverageLocal;
	// NSCBC outlet length scale.
	su2double mNSCBC_Outlet_len;	
	// NSCBC outlet relaxation coefficient for normal wave amplitude. 
	su2double mNSCBC_Outlet_sigma;
	// NSCBC outlet relaxation coefficient for coupled   transverse 
	// wave amplitude.
	su2double mNSCBC_Outlet_beta_l;	
	// NSCBC outlet relaxation coefficient for uncoupled transverse 
	// wave amplitude.
	su2double mNSCBC_Outlet_beta_t;

	// NSCBC inlet length scale.
	su2double mNSCBC_Inlet_len;	
	// NSCBC inlet relaxation coefficient for normal wave amplitude. 
	su2double mNSCBC_Inlet_sigma;
	// Type of reconstruction in the inlet NSCBC: whether or not to use
	// the outgoing wave amplitude in the reconstruction process.
	// Default is yes. Note, this means for a total inlet NSCBC, this is reflecting.
	bool mNSCBC_Inlet_incoming;

	// Whether or not to use a characteristic matching layer.
	bool mCharacteristicMatchingLayer;

	// Whether or not to use a sponge layer. Default is no.
	bool mSpongeLayer;
	// Whether or not a sponge layer is used, per block face.
	bool mSpongeLayerSpecified[6];
	// Number of elements in each layer, per block face.
	int  mNSpongeLayerElements[6];
	// Damping constant in each layer, per block face.
	su2double mDampingConstant[6];
	// Damping exponent in each layer, per block face.
	su2double mDampingExponent[6];

  // Name of the averaged solution file used in a sponge layer.
  std::string mAveragedSolutionFile;

private:
  //------------------------------------------------
  // Private member variables.
  //------------------------------------------------

  // Function, which extracts a possible parameter.
  void AnalyzeString(std::string &lineBuf);

  // Function to set a logical parameter depending on the input string.
  void CheckYesNoKeyword(const std::string &keyword,
                         std::string       &value,
                         bool              &param);

  // Function, which checks the input parameters.
  void CheckParameters(void);

  // Function to create a template parameter file.
  void CreateTemplateParamFile(const char *fileName);

  // Template function to read data from the given string.
  template <class T>
  void ReadDataFromString(const std::string &keyword,
                          const std::string &value,
                          T                 *param,
                          const int         nArg = 1);

  //------------------------------------------------
  // Disabled constructors and assignment operators.
  //------------------------------------------------

  // Copy constructor, disabled.
  InputParamClass(const InputParamClass&);

  // Assignment operator, disabled.
  InputParamClass& operator=(const InputParamClass&);
};

//------------------------------------------------------------------------------

// Implementation of the template function ReadDataFromString.
// This function attempts to read a certain amount of data from a given string.
template <class T>
void InputParamClass::ReadDataFromString(const std::string &keyword,
                                         const std::string &value,
                                         T                 *param,
                                         const int         nArg)
{
  // Couple the string value to the input stream istr and read the
  // data. Check if everything went okay. If not, create an error message
  // and exit.
  std::istringstream istr(value);

  for(int nn=0; nn<nArg; ++nn)
  {
    if(!(istr >> param[nn]))
    {
      std::ostringstream message;
      message << "\"" << keyword << "\" wrongly specified.";
      TerminateAll("InputParamClass::ReadDataFromString",
                   __FILE__, __LINE__, message.str());
    }
  }
}
