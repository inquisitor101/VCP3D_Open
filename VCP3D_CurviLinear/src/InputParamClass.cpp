//------------------------------------------------------------------------------
// File, which contains the implementation of the member functions of the
// class InputParamClass.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"

//------------------------------------------------------------------------------

// Constuctor. The default values are set.
InputParamClass::InputParamClass()
{
  mDOFLocation         = NO_LOCATION;
  mFEMVariables        = NO_VARIABLES;
  mRiemannSolver       = NO_RIEMANN;
  mFluctuationsInitSol = NO_FLUCTUATIONS;

  mNPolyGridDOFs = -1;
  mNPolySolDOFs  = -1;
  mNPolyIntRule  = -1;

  mMach           = -one;
  mReynolds       = -one;
  mReynoldsLength = -one;
  mTempFreeStream = -one;

  mVelDirFreeStream[0] = valLarge;
  mVelDirFreeStream[1] = valLarge;
  mVelDirFreeStream[2] = valLarge;

  mLiftDirection[0] = valLarge;
  mLiftDirection[1] = valLarge;
  mLiftDirection[2] = valLarge;

  mFBody[0] = zero;
  mFBody[1] = zero;
  mFBody[2] = zero;
  mFBodySpecified = false;

	mPulseCenter[0] = zero;
	mPulseCenter[1] = zero;
	mPulseCenter[2] = zero;
	mPulseStrength  = half;
	mPulseWidth     = half;
	mQuiescentFlow  = false;

  mLowerCoorBoxFluctuations[0] = -valLarge;
  mLowerCoorBoxFluctuations[1] = -valLarge;
  mLowerCoorBoxFluctuations[2] = -valLarge;

  mUpperCoorBoxFluctuations[0] = valLarge;
  mUpperCoorBoxFluctuations[1] = valLarge;
  mUpperCoorBoxFluctuations[2] = valLarge;

  mExchangeLocation        = -one;
  mNPointsGridWallModel    = -1;
  mExpansionRatioWallModel = -one;

  mGridFile         = "";
  mSolFile          = "Solution.dat";
  mRestartFile      = "";
  mParaviewFile     = "Solution.vtk";
  mSpaceAverageFile = "SpaceAverage.vtk";
  mSU2GridFile      = "Grid.su2";

  mRestart          = false;
  mWriteSU2GridFile = false;

  mComputeAverageSolution       = false;
  mContinueAveragingFromRestart = false;

  mAllowSplitIDir = true;
  mAllowSplitJDir = true;
  mAllowSplitKDir = true;

  mMonitorForces = true;
  mReferenceArea = one;

  mNTimeSynchrMax    = -1;
  mNSave             =  0;
  mNWriteTimeHistory =  1;

  mDeltaTSynchr = -one;
  mCFL          = -one;

  mThetaSym                 = one;
  mRelativePenaltyParameter = one;

  mSGSModelType = NO_SGS_MODEL;
  mSGSModel     = NULL;

  mWallModelType = NO_WALL_MODEL;
  mWallModel     = NULL;

	mNSCBC_Specified           = false;
	mNSCBC_Outlet_OneMinusBeta = true;
	mNSCBC_Outlet_AverageLocal = false;

	mNSCBC_Outlet_len    = one;
	mNSCBC_Outlet_sigma  = zero;
	mNSCBC_Outlet_beta_l = one;
	mNSCBC_Outlet_beta_t = one;

	mNSCBC_Inlet_len     = one;
	mNSCBC_Inlet_sigma   = one;
	mNSCBC_Inlet_beta    = zero;
}

//------------------------------------------------------------------------------

// Destructor. Release the memory, if needed.
InputParamClass::~InputParamClass()
{
  for(unsigned int i=0; i<mSubfaces.size(); ++i)
  {
    for(unsigned int j=0; j<mSubfaces[i].size(); ++j)
    {
      if( mSubfaces[i][j] ) delete mSubfaces[i][j];
      mSubfaces[i][j] = NULL;
    }
  }

  if( mSGSModel )
  {
    delete mSGSModel;
    mSGSModel = NULL;
  }

  if( mWallModel )
  {
    delete mWallModel;
    mWallModel = NULL;
  }
}

//------------------------------------------------------------------------------

// Function, which extracts a possible input parameter from the given string.
void InputParamClass::AnalyzeString(std::string &lineBuf)
{
  // Replace all tab and return characters by spaces and remove the comment
  // part from lineBuf.
  ReplaceTabsAndReturns(lineBuf);

  std::string::size_type pos = lineBuf.find("#");
  if(pos != std::string::npos) lineBuf.erase(pos, lineBuf.length()-pos);

  // Search for the delimiting character, which is :, in lineBuf. If not
  // present, return to incate that nothing went wrong. If present,
  // the strings keyword and value are created.
  pos = lineBuf.find(":");
  if(pos == std::string::npos) return;

  std::string keyword(lineBuf, 0, pos);
  std::string value(lineBuf, pos+1);

  // Remove the leading and trailing blanks in keyword and value. If either of
  // them has zero length after the removal no information is present and a
  // return is made, without doing anything.
  RemoveBlanks(keyword);
  RemoveBlanks(value);

  if( !keyword.length() || !value.length()) return;

  // Create a lower case version of keyword. The string value should
  // not be changed yet, because it may contain a file name.
  CreateLowerCase(keyword);

  //----------------------------------------------------
  // Extraction of the input parameters.
  //----------------------------------------------------

  if(keyword == "location of the dofs")
  {
    CreateLowerCase(value);
    if(     value == "equidistant") mDOFLocation = EQUIDISTANT;
    else if(value == "lgl")         mDOFLocation = LGL_POINTS;
    else
      TerminateAll("InputParamClass::AnalyzeString", __FILE__, __LINE__,
                   "Unknown \"Location of the DOFs\" specified");
    return;
  }

  if(keyword == "fem working variables")
  {
    CreateLowerCase(value);
    if(     value == "conservative") mFEMVariables = CONSERVATIVE_VARIABLES;
    else if(value == "entropy")      mFEMVariables = ENTROPY_VARIABLES;
    else
      TerminateAll("InputParamClass::AnalyzeString", __FILE__, __LINE__,
                   "Unknown \"FEM working variables\" specified");
    return;
  }

  if(keyword == "riemann solver")
  {
    CreateLowerCase(value);
    if(     value == "roe")       mRiemannSolver = ROE;
    else if(value == "ismailroe") mRiemannSolver = ISMAIL_ROE;
    else
      TerminateAll("InputParamClass::AnalyzeString", __FILE__, __LINE__,
                   "Unknown \"Riemann solver\" specified");
    return;
  }

  if(keyword == "wall modeling")
  {
    CreateLowerCase(value);
    if(     value == "no wall model"){
      mWallModelType = NO_WALL_MODEL;
    }          
    else if(value == "equilibrium wall model"){
      mWallModelType = EQUILIBRIUM_WALL_MODEL;
    }
    else if(value == "logarithmic wall model"){
      mWallModelType = LOGARITHMIC_WALL_MODEL;
    }
    else{
      TerminateAll("InputParamClass::AnalyzeString", __FILE__, __LINE__,
                   "Unknown \"Wall modeling\" specified");
    }
    return;
  }

  if(keyword == "subgrid scale model")
  {
    CreateLowerCase(value);
    if(     value == "none")   mSGSModelType = NO_SGS_MODEL;
    else if(value == "wale")   mSGSModelType = WALE_MODEL;
    else if(value == "vreman") mSGSModelType = VREMAN_MODEL;
    else
      TerminateAll("InputParamClass::AnalyzeString", __FILE__, __LINE__,
                   "Unknown \"Subgrid scale model\" specified");
    return;
  }

  if(keyword == "fluctuations to add to initial solution")
  {
    CreateLowerCase(value);
    if(     value == "no fluctuations")
      mFluctuationsInitSol = NO_FLUCTUATIONS;
    else if(value == "random fluctuations")
      mFluctuationsInitSol = RANDOM;
    else if(value == "sine_and_cosine")
      mFluctuationsInitSol = SINE_AND_COSINE;
    else if(value == "synthetic fluctuations")
      mFluctuationsInitSol = SYNTHETIC;
		else if(value == "pressure pulse")
			mFluctuationsInitSol = PRESSURE_PULSE;
    else
     TerminateAll("InputParamClass::AnalyzeString", __FILE__, __LINE__,
                   "Unknown \"Fluctuations to add to initial solution\" specified");
    return;
  }

	if(keyword == "pressure pulse center")
		return ReadDataFromString(keyword, value, mPulseCenter, 3);
	
	if(keyword == "pressure pulse strength")
		return ReadDataFromString(keyword, value, &mPulseStrength);

	if(keyword == "pressure pulse width")
		return ReadDataFromString(keyword, value, &mPulseWidth);

  if(keyword == "lower coordinates box for fluctuations")
   return ReadDataFromString(keyword, value, mLowerCoorBoxFluctuations, 3);

  if(keyword == "upper coordinates box for fluctuations")
   return ReadDataFromString(keyword, value, mUpperCoorBoxFluctuations, 3);

  if(keyword == "polynomial degree grid")
    return ReadDataFromString(keyword, value, &mNPolyGridDOFs);

  if(keyword == "polynomial degree solution")
    return ReadDataFromString(keyword, value, &mNPolySolDOFs);

  if(keyword == "polynomial degree integration rule")
    return ReadDataFromString(keyword, value, &mNPolyIntRule);

  if(keyword == "mach number")
    return ReadDataFromString(keyword, value, &mMach);

  if(keyword == "reynolds number")
    return ReadDataFromString(keyword, value, &mReynolds);

  if(keyword == "reynolds length")
    return ReadDataFromString(keyword, value, &mReynoldsLength);

  if(keyword == "free stream temperature")
    return ReadDataFromString(keyword, value, &mTempFreeStream);

  if(keyword == "free stream velocity direction")
    return ReadDataFromString(keyword, value, mVelDirFreeStream, 3);

  if(keyword == "lift direction")
    return ReadDataFromString(keyword, value, mLiftDirection, 3);

  if(keyword == "body force per unit volume (x,y,z)")
  {
    mFBodySpecified = true;
    return ReadDataFromString(keyword, value, mFBody, 3);
  }

  if(keyword == "exchange location wall model")
    return ReadDataFromString(keyword, value, &mExchangeLocation);

  if(keyword == "number of grid points equilibrium wall model")
    return ReadDataFromString(keyword, value, &mNPointsGridWallModel);

  if(keyword == "expansion ratio grid equilibrium wall model")
    return ReadDataFromString(keyword, value, &mExpansionRatioWallModel);

  if(keyword == "grid file name"){mGridFile = value; return;}

  if(keyword == "solution file name"){mSolFile = value; return;}

  if(keyword == "restart file name") {mRestartFile = value; return;}

  if(keyword == "paraview visualization file name") {mParaviewFile = value; return;}

  if(keyword == "space average data file name") {mSpaceAverageFile = value; return;}

  if(keyword == "su2 grid file name") {mSU2GridFile = value; return;}

  if(keyword == "restart")
    return CheckYesNoKeyword(keyword, value, mRestart);

  if(keyword == "write su2 grid file")
    return CheckYesNoKeyword(keyword, value, mWriteSU2GridFile);

  if(keyword == "compute average solution")
    return CheckYesNoKeyword(keyword, value, mComputeAverageSolution);

  if(keyword == "continue averaging from restart")
    return CheckYesNoKeyword(keyword, value, mContinueAveragingFromRestart);

  if(keyword == "allow block splitting in i-direction")
    return CheckYesNoKeyword(keyword, value, mAllowSplitIDir);

  if(keyword == "allow block splitting in j-direction")
    return CheckYesNoKeyword(keyword, value, mAllowSplitJDir);

  if(keyword == "allow block splitting in k-direction")
    return CheckYesNoKeyword(keyword, value, mAllowSplitKDir);

  if(keyword == "monitor forces")
    return CheckYesNoKeyword(keyword, value, mMonitorForces);

  if(keyword == "reference area force coefficients")
    return ReadDataFromString(keyword, value, &mReferenceArea);

  if(keyword == "number of synchronization time steps")
    return ReadDataFromString(keyword, value, &mNTimeSynchrMax);

  if(keyword == "save solution every")
    return ReadDataFromString(keyword, value, &mNSave);

  if(keyword == "write time history info every")
    return ReadDataFromString(keyword, value, &mNWriteTimeHistory);

  if(keyword == "synchronization time step")
    return ReadDataFromString(keyword, value, &mDeltaTSynchr);

  if(keyword == "cfl number")
    return ReadDataFromString(keyword, value, &mCFL);

  if(keyword == "theta parameter symmetrizing terms")
    return ReadDataFromString(keyword, value, &mThetaSym);

  if(keyword == "relative penalty parameter")
    return ReadDataFromString(keyword, value, &mRelativePenaltyParameter);

	if(keyword == "outlet use local (element) averaged mach number")
		return CheckYesNoKeyword(keyword, value, mNSCBC_Outlet_AverageLocal);
	if(keyword == "outlet use one minus beta transverse relaxation")
		return CheckYesNoKeyword(keyword, value, mNSCBC_Outlet_OneMinusBeta);
	if(keyword == "outlet characteristic length")
		return ReadDataFromString(keyword, value, &mNSCBC_Outlet_len);
	if(keyword == "outlet normal relaxation (sigma)")
		return ReadDataFromString(keyword, value, &mNSCBC_Outlet_sigma);
	if(keyword == "outlet coupled   transverse relaxation (beta_l)")
		return ReadDataFromString(keyword, value, &mNSCBC_Outlet_beta_l);
	if(keyword == "outlet uncoupled transverse relaxation (beta_t)")
		return ReadDataFromString(keyword, value, &mNSCBC_Outlet_beta_t);
	
	if(keyword == "inlet characteristic length")
		return ReadDataFromString(keyword, value, &mNSCBC_Inlet_len);
	if(keyword == "inlet normal relaxation (sigma)")
		return ReadDataFromString(keyword, value, &mNSCBC_Inlet_sigma);
	if(keyword == "inlet transverse relaxation (beta)")
		return ReadDataFromString(keyword, value, &mNSCBC_Inlet_beta);

  // Unknown keyword. Print an error message and terminate.
  std::ostringstream message;
  message << "Unknown keyword \"" << keyword << "\" encountered";
  TerminateAll("InputParamClass::AnalyzeString",
               __FILE__, __LINE__, message.str());
}

//------------------------------------------------------------------------------

// Function to set a logical parameter depending on the input string.
void InputParamClass::CheckYesNoKeyword(const std::string &keyword,
                                        std::string       &value,
                                        bool              &param)
{
  // Convert value to lower case, such that testing of its contents is
  // a lot easier. Determine the case afterwards.
  CreateLowerCase(value);

  if(value == "yes"){param = true;  return;}
  if(value == "no") {param = false; return;}

  // Neither yes or no specified. Write an error message and exit.

  std::ostringstream message;
  message << "\"" << keyword << "\" wrongly specified. "
          << "It should be yes or no.";
  TerminateAll("InputParamClass::CheckYesNoKeyword",
               __FILE__, __LINE__, message.str());
}

//------------------------------------------------------------------------------

// Function, which checks the input parameters.
void InputParamClass::CheckParameters(void)
{
  std::ostringstream message;

  // Check if the key parameters were specified correctly. If not,
  // write an error message and exit.
  if(mDOFLocation == NO_LOCATION)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"Location of the DOFs\" not specified in parameter file");

  if(mFEMVariables == NO_VARIABLES)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"FEM working variables\" not specified in parameter file");

  if(mRiemannSolver == NO_RIEMANN)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"Riemann solver\" not specified in parameter file");

  if(mNPolyGridDOFs <= 0)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"Polynomial degree grid\" not or wrongly "
                 "specified in parameter file");

  if(mNPolyGridDOFs > 9)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"Polynomial degree grid\" must be less than 10");

  if(mNPolySolDOFs <= 0)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"Polynomial degree solution\" not or wrongly "
                 "specified in parameter file");

  if(mNPolySolDOFs > 9)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"Polynomial degree solution\" must be less than 10");

  if(mNPolyIntRule <= 0)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"Polynomial degree integration rule\" not or wrongly "
                 "specified in parameter file");

  if(mNPolyIntRule < 2*mNPolySolDOFs)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"Polynomial degree integration rule\" not accurate enough");

  if(mWallModelType != NO_WALL_MODEL)
  {
    if(mExchangeLocation <= zero)
      TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                   "\"Exchange location wall model\" not or wrongly "
                   "specified in parameter file");

    if(mWallModelType == EQUILIBRIUM_WALL_MODEL)
    {
      if((mExpansionRatioWallModel < one) || (mExpansionRatioWallModel > two))
        TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                     "\"Expansion ratio grid equilibrium wall model\" not "
                     "or wrongly specified");

      if(mNPointsGridWallModel <= 1)
        TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                     "\"Number of grid points equilibrium wall model\" not "
                     "or wrongly specified");
    }
  }

  if(mGridFile == "")
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"Grid file name\" not specified in parameter file");

  if(mRestartFile == "" && mRestart)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"Restart file name\" not specified in parameter file");

  if(mNTimeSynchrMax <= 0)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"Number of synchronization time steps\" not or wrongly "
                 "specified in parameter file");

  if(mDeltaTSynchr <= 0.0)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"Synchronization time step\" not or wrongly "
                 "specified in parameter file");

  if(mCFL <= zero)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"CFL number\" not or wrongly specified "
                 "in parameter file");

  if(mRelativePenaltyParameter <= zero)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"Relative penalty parameter\" wrongly specified "
                 "in parameter file");

  if(mReferenceArea <= zero)
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "\"Reference area force coefficients\" wrongly specified "
                 "in parameter file");

  if(mLowerCoorBoxFluctuations[0] >= mUpperCoorBoxFluctuations[0])
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "Wrong x-coordinates of the box for fluctuations");

  if(mLowerCoorBoxFluctuations[1] >= mUpperCoorBoxFluctuations[1])
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "Wrong y-coordinates of the box for fluctuations");

  if(mLowerCoorBoxFluctuations[2] >= mUpperCoorBoxFluctuations[2])
    TerminateAll("InputParamClass::CheckParameters", __FILE__, __LINE__,
                 "Wrong z-coordinates of the box for fluctuations");

  // Overwrite some parameters, depending on the value of other parameters.
  if(mNSave < 0) mNSave = 0;

  if(mNWriteTimeHistory <= 0) mNWriteTimeHistory = 1;

  // If the last 4 characters of the paraview file are ".vtk", remove those.
  // During the writing these are added again.
  std::size_t found = mParaviewFile.find(".vtk");
  if(found!=std::string::npos)
    mParaviewFile.erase(mParaviewFile.end()-4, mParaviewFile.end());

  // If no restart is performed, set mContinueAveragingFromRestart to false.
  if( !mRestart ) mContinueAveragingFromRestart = false;

  // If no averaging needs to be performed, set mContinueAveragingFromRestart to false.
  if( !mComputeAverageSolution ) mContinueAveragingFromRestart = false;
}

//------------------------------------------------------------------------------

// Function, which computes the reference conditions.
void InputParamClass::ComputeReferenceConditions(void)
{
  // Compute the viscosity. Set it to a constant for now.
  mu = (su2double) 1.716e-5;

  // Attempt to extract pressure, density and temperature from the
  // in- and outflow boundaries.
  su2double p = -one, pTot = -one, rho = -one, rhoTot = -one, T = -one, TTot = -one;

  // Determine if this is expected to be a quiescent simulation.
  if( FABS(mMach) < epsThreshold ) mQuiescentFlow = true; 

  for(unsigned int i=0; i<mSubfaces.size(); ++i)
  {
    for(unsigned int j=0; j<mSubfaces[i].size(); ++j)
    {
      const std::vector<PrescribedDataClass> &prescribedData = mSubfaces[i][j]->GetPrescribedData();

      for(unsigned int k=0; k<prescribedData.size(); ++k)
      {
        if(     prescribedData[k].mNameVariable == "Pressure")         p      = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "TotalPressure")    pTot   = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "Density")          rho    = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "TotalDensity")     rhoTot = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "Temperature")      T      = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "TotalTemperature") TTot   = prescribedData[k].mData[0];
      }

      // Check if all BCs are outlet, in case of a quiescent flow.
      if( mQuiescentFlow ) 
				if( mSubfaces[i][j]->GetTypeBoundaryPrescribed() != BC_OUTFLOW_SUBSONIC
						&&
						mSubfaces[i][j]->GetTypeBoundaryPrescribed() != BC_OUTFLOW_CHARACTERISTIC )
					TerminateAll("InputParamClass::ComputeReferenceConditions", __FILE__, __LINE__,
											 "Must specify 6 outlet-type faces for a quiescent simulation.");
    }
  }

	// In case of quiescent flow configuration and no free-stream temperature specification, 
	// use below reference temperature of T = 300.0.
	if( mQuiescentFlow ) T = ( mTempFreeStream < zero ) ? 300.0 : mTempFreeStream;

  // Check if a reference pressure can be deducted.
  pRef = -one;
  if(     p    > zero) pRef = p;
  else if(pTot > zero) pRef = pTot;

  // Check if a reference density can be conducted directly.
  rhoRef = -one;
  if(     rho    > zero) rhoRef = rho;
  else if(rhoTot > zero) rhoRef = rhoTot;

  // If no reference density could be constructed, check if it possible via the
  // temperature.
  if((rhoRef < zero) && (pRef > zero))
  {
    su2double TRef = -one;
    if(     T    > zero) TRef = T;
    else if(TTot > zero) TRef = TTot;

    if(TRef > zero) rhoRef = pRef/(RGas*TRef);
  }

  //// DEBUG. Set the reference pressure and density, such that
  ////        it can be compared with the channel flow codes.
  // pRef   = 2.70736890373328;
  // rhoRef = 3.4528409751710832e-05;
  //// END DEBUG

  // Check if it was possible to obtain the reference values for the
  // pressure and density. If not, they must be conducted from the
  // free stream values.
  if((pRef < zero) || (rhoRef < zero))
  {
    // No inflow or outflow boundaries, so most likely an external flow
    // with a farfield. Check if the necessary free-stream data were specified.
    if(mMach <= zero)
      TerminateAll("InputParamClass::ComputeReferenceConditions", __FILE__, __LINE__,
                   "\"Mach number\" not or wrongly specified in parameter file");

    if(mReynolds <= zero)
      TerminateAll("InputParamClass::ComputeReferenceConditions", __FILE__, __LINE__,
                   "\"Reynolds number\" not or wrongly specified in parameter file");

    if(mReynoldsLength <= zero)
      TerminateAll("InputParamClass::ComputeReferenceConditions", __FILE__, __LINE__,
                   "\"Reynolds length\" not or wrongly specified in parameter file");

    if(mTempFreeStream < zero)
      TerminateAll("InputParamClass::ComputeReferenceConditions", __FILE__, __LINE__,
                   "\"Free stream temperature\" not or wrongly specified in parameter file");

    // Check if the free stream direction has been specified and create a unit vector.
    // If not, terminate.
    if((mVelDirFreeStream[0] < valLarge) && (mVelDirFreeStream[1] < valLarge) &&
       (mVelDirFreeStream[2] < valLarge))
    {
      const su2double invLen = one/SQRT(mVelDirFreeStream[0]*mVelDirFreeStream[0]
                             +          mVelDirFreeStream[1]*mVelDirFreeStream[1]
                             +          mVelDirFreeStream[2]*mVelDirFreeStream[2]);

      mVelDirFreeStream[0] *= invLen;
      mVelDirFreeStream[1] *= invLen;
      mVelDirFreeStream[2] *= invLen;
    }
    else
      TerminateAll("InputParamClass::ComputeReferenceConditions", __FILE__, __LINE__,
                   "\"Free stream velocity direction\" not or wrongly specified in parameter file");

    // Check if a lift direction has not been specified. If not, set the default value
    // which is normal to the free stream direction without a z-component.
    if( !((mLiftDirection[0] < valLarge) && (mLiftDirection[1] < valLarge) &&
          (mLiftDirection[2] < valLarge)) )
    {
      mLiftDirection[0] = -mVelDirFreeStream[1];
      mLiftDirection[1] =  mVelDirFreeStream[0];
      mLiftDirection[2] =  zero;
    }

    // Create a unit vector out of the specified lift direction.
    const su2double invLen = one/SQRT(mLiftDirection[0]*mLiftDirection[0]
                           +          mLiftDirection[1]*mLiftDirection[1]
                           +          mLiftDirection[2]*mLiftDirection[2]);

    mLiftDirection[0] *= invLen;
    mLiftDirection[1] *= invLen;
    mLiftDirection[2] *= invLen;

    // Check if the lift direction is normal to the free stream direction.
    const su2double dotProd = mVelDirFreeStream[0]*mLiftDirection[0]
                            + mVelDirFreeStream[1]*mLiftDirection[1]
                            + mVelDirFreeStream[2]*mLiftDirection[2];
    if(FABS(dotProd) > (su2double) 1.e-3)
      TerminateAll("InputParamClass::ComputeReferenceConditions", __FILE__, __LINE__,
                   "\"Lift direction\" is not normal to \"Free stream velocity direction\"");

    // Compute the free stream velocity from the temperature and Mach number.
    const su2double vInf = mMach*SQRT(GamConstant*RGas*mTempFreeStream);

    // Compute the reference density from the Reynolds number,
    // dynamic viscosity and velocity.
    rhoRef = mReynolds*mu/(mReynoldsLength*vInf);

    // Compute the reference pressure via the gas law.
    pRef = rhoRef*RGas*mTempFreeStream;
  }

  // Compute the reference velocity.
  uRef = SQRT(pRef/rhoRef);

  // Non-dimensionalize the viscosity. Note that the reference length is missing
  // in this non-dimensionalization, because it is taken to be 1.0
  mu *= uRef/pRef;
}

//------------------------------------------------------------------------------

// Function, which creates a template parameter file.
void InputParamClass::CreateTemplateParamFile(const char *fileName)
{
  std::ostringstream message;

  // Check for the master rank.
  if(rank == 0)
  {
    // Open the file for writing. If not successful, print an error message
    // and exit.
    std::ofstream paramFile(fileName);
    if( !paramFile )
    {
      message << "Template parameter file " << fileName << " could not be opened";
      Terminate("InputParamClass::CreateTemplateParamFile",
                __FILE__, __LINE__, message.str());
    }

    // Write the header.
    paramFile << "============================================================"
              << "=====================" << std::endl;
    paramFile << "This is a template parameter file for VCP3D_CurviLinear."
              << std::endl;
    paramFile << "It contains all possible parameters in several sections and "
              << "must be adjusted" << std::endl;
    paramFile << "by the user to suit his requirements."
              << std::endl << std::endl;

    paramFile << "Comments are indicated by #, i.e. all info following a # "
              << "is ignored." << std::endl;
    paramFile << "The sequence of the parameters is arbitrary and keywords "
              << "are case insensitive." << std::endl;
    paramFile << "When a keyword occurs multiple times, the last value "
              << "is taken." << std::endl;
    paramFile << "============================================================"
              << "=====================" << std::endl << std::endl;

    // Write the grid parameters.
    paramFile << "----------------------------------------------------------"
              << std::endl;
    paramFile << "               Grid parameters" << std::endl;
    paramFile << "----------------------------------------------------------"
              << std::endl << std::endl;

    paramFile << "             Polynomial degree grid: MISSING PARAMETER" << std::endl;
    paramFile << "         Polynomial degree solution: MISSING PARAMETER" << std::endl;
    paramFile << " Polynomial degree integration rule: MISSING PARAMETER" << std::endl;
    paramFile << "               Location of the DOFs: Equidistant" << std::endl;
    paramFile << "                # Other possibility: LGL"         << std::endl
              << std::endl;

    // Write the IO parameters.
    paramFile << "----------------------------------------------------------"
              << std::endl;
    paramFile << "                 IO parameters" << std::endl;
    paramFile << "----------------------------------------------------------"
              << std::endl << std::endl;

    paramFile << "                   Grid file name: MISSING FILE NAME" << std::endl
              << std::endl;

    paramFile << "                Restart file name: MISSING FILE NAME" << std::endl;
    paramFile << "                          Restart: no   # yes" << std::endl
              << std::endl;

    paramFile << "               Solution file name: Solution.dat" << std::endl
              << std::endl;

    paramFile << " Paraview visualization file name: Solution.vtk" << std::endl
              << std::endl;

    paramFile << "     Space average data file name: SpaceAverage.vtk" << std::endl
              << std::endl;

    paramFile << "               SU2 grid file name: GridSU2.su2" << std::endl;
    paramFile << "              Write SU2 grid file: no   # yes" << std::endl
              << std::endl;

    paramFile << "              Save solution every: MISSING PARAMETER" << std::endl;
    paramFile << "    Write time history info every: 1" << std::endl
              << std::endl;

    // Write the discretization parameters.
    paramFile << "----------------------------------------------------------"
              << std::endl;
    paramFile << "                 Discretization parameters" << std::endl;
    paramFile << "----------------------------------------------------------"
              << std::endl << std::endl;

    paramFile << "  FEM working variables: Conservative" << std::endl;
    paramFile << "    # Other possibility: Entropy" << std::endl << std::endl;

    paramFile << "         Riemann solver: Roe" << std::endl;
    paramFile << "    # Other possibility: IsmailRoe" << std::endl << std::endl;

    paramFile << "    Subgrid scale model: None" << std::endl;
    paramFile << "  # Other possibilities: WALE" << std::endl;
    paramFile << "  #                      Vreman" << std::endl
              << std::endl;

    paramFile << "          Wall modeling: No wall model" << std::endl;
    paramFile << "  # Other possibilities: Equilibrium wall model" << std::endl;
    paramFile << "  #                      Logarithmic wall model" << std::endl
              << std::endl;

    paramFile << "                 Exchange location wall model: MISSING PARAMETER" << std::endl;
    paramFile << " Number of grid points equilibrium wall model: MISSING PARAMETER" << std::endl;
    paramFile << "  Expansion ratio grid equilibrium wall model: MISSING PARAMETER" << std::endl
              << std::endl;

    paramFile << "    Theta parameter symmetrizing terms: 1.0" << std::endl;
    paramFile << "            Relative penalty parameter: 1.0" << std::endl
              << std::endl;

    // Write the physics parameters.
    paramFile << "----------------------------------------------------------"
              << std::endl;
    paramFile << "                       Physics parameters" << std::endl;
    paramFile << "----------------------------------------------------------"
              << std::endl << std::endl;

    paramFile << "                   Mach number: MISSING PARAMETER" << std::endl;
    paramFile << "               Reynolds number: MISSING PARAMETER" << std::endl;
    paramFile << "               Reynolds length: MISSING PARAMETER" << std::endl;
    paramFile << "       Free stream temperature: MISSING PARAMETER" << std::endl;
    paramFile << "Free stream velocity direction: 1.0 0.0 0.0"       << std::endl;
    paramFile << "                Lift direction: 0.0 1.0 0.0"       << std::endl
              << std::endl;

    paramFile << "      Body force per unit volume (x,y,z): 0.0 0.0 0.0" << std::endl
              << std::endl;

    paramFile << " Fluctuations to add to initial solution: No fluctuations" << std::endl;
    paramFile << "                   # Other possibilities: Random fluctuations" << std::endl;
    paramFile << "                   #                      Sine and cosine" << std::endl;
    paramFile << "                   #                      Synthetic fluctuations" << std::endl
              << std::endl;

    paramFile << "Lower coordinates box for fluctuations: -1.e+10 -1.e+10 -1.e+10" << std::endl;
    paramFile << "Upper coordinates box for fluctuations:  1.e+10  1.e+10  1.e+10" << std::endl
              << std::endl;

    paramFile << "                        Monitor forces: no # yes" << std::endl;
    paramFile << "     Reference area force coefficients: MISSING PARAMETER" << std::endl
              << std::endl;

    // Write the time stepping parameters.
    paramFile << "----------------------------------------------------------"
              << std::endl;
    paramFile << "                 Time stepping parameters" << std::endl;
    paramFile << "----------------------------------------------------------"
              << std::endl << std::endl;

    paramFile << " Number of synchronization time steps: MISSING PARAMETER"
              << std::endl;
    paramFile << "            Synchronization time step: MISSING PARAMETER"
              << std::endl;
    paramFile << "                           CFL number: MISSING PARAMETER"
              << std::endl;

    paramFile << "             Compute average solution: no    # yes" << std::endl;
    paramFile << "      Continue averaging from restart: no    # yes" << std::endl
              << std::endl;

    // Write the block splitting parameters.
    paramFile << "----------------------------------------------------------"
              << std::endl;
    paramFile << "              Block splitting parameters" << std::endl;
    paramFile << "----------------------------------------------------------"
              << std::endl << std::endl;

    paramFile << "  Allow block splitting in i-direction: yes    # no" << std::endl;
    paramFile << "  Allow block splitting in j-direction: yes    # no" << std::endl;
    paramFile << "  Allow block splitting in k-direction: yes    # no" << std::endl
              << std::endl;

    // Write a template for the subface information of the block.
    paramFile << "----------------------------------------------------------"
              << std::endl;
    paramFile << "                 Subface information" << std::endl;
    paramFile << "----------------------------------------------------------"
              << std::endl << std::endl;

    paramFile << "#Face nSfaces    BoundaryCondition     ib  ie  jb  je  kb  ke       idb ide jdb jde kdb kde     Orientation        Translation"
              << std::endl;
    paramFile << "    1    1" << std::endl;
    paramFile << "                    Internal1to1        1   1   1  17   1  11        49  49   1  17   1  11      1   2   3" << std::endl;
    paramFile << "    2    1" << std::endl;
    paramFile << "                    Internal1to1       49  49   1  17   1  11         1   1   1  17   1  11      1   2   3" << std::endl;
    paramFile << "    3    1" << std::endl;
    paramFile << "                  BCWallHeatFlux        1  49   1   1   1  11" << std::endl;
    paramFile << "    4    1" << std::endl;
    paramFile << "                      BCFarfield        1  49  17  17   1  11" << std::endl;
    paramFile << "    5    1" << std::endl;
    paramFile << "                 Periodic1to1Trans      1  49   1  17   1   1         1  49   1  17  11  11      1   2   3        0.0  0.0  0.2"
              << std::endl;
    paramFile << "    6    1" << std::endl;
    paramFile << "                 Periodic1to1Trans      1  49   1  17  11  11         1  49   1  17   1   1      1   2   3        0.0  0.0 -0.2"
              << std::endl;

    // Write a template for the prescribed boundary data.
    paramFile << "----------------------------------------------------------"
              << std::endl;
    paramFile << "                 Prescribed boundary data" << std::endl;
    paramFile << "----------------------------------------------------------"
              << std::endl << std::endl;

    paramFile << "# FaceID   SubfaceID   nPoints    Variable1  Variable2" << std::endl;
    paramFile << "    3          1         3            x      HeatFlux"  << std::endl;
    paramFile << "                                     0.0      8.0"      << std::endl;
    paramFile << "                                     0.5      4.0"      << std::endl;
    paramFile << "                                     1.0      6.0"      << std::endl;

    // Close the parameter file.
    paramFile.close();
  }

  // Create an error message that the parameter file was not found and a
  // template has been created.
  message << "Parameter file " << fileName << " not found. "
          << "A template has been created.";
  TerminateAll("InputParamClass::CreateTemplateParamFile",
               __FILE__, __LINE__, message.str());
}

//------------------------------------------------------------------------------

// Function, which determines the non-dimensional primitive variables
// of the free-stream.
void InputParamClass::DetermineNonDimPrimVarFreeStream(su2double *primVar) const
{
  // Attempt to extract the variables from the in- and outflow boundaries.
  su2double p = -one, pTot = -one, rho = -one, rhoTot = -one, T = -one, TTot = -one;
  su2double u = valLarge, v = valLarge, w = valLarge;
  su2double uDir = valLarge, vDir = valLarge, wDir = valLarge;

  for(unsigned int i=0; i<mSubfaces.size(); ++i)
  {
    for(unsigned int j=0; j<mSubfaces[i].size(); ++j)
    {
      const std::vector<PrescribedDataClass> &prescribedData = mSubfaces[i][j]->GetPrescribedData();

      for(unsigned int k=0; k<prescribedData.size(); ++k)
      {
        if(     prescribedData[k].mNameVariable == "Pressure")         p      = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "TotalPressure")    pTot   = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "Density")          rho    = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "TotalDensity")     rhoTot = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "Temperature")      T      = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "TotalTemperature") TTot   = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "VelocityX")        u      = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "VelocityY")        v      = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "VelocityZ")        w      = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "VelocityDirX")     uDir   = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "VelocityDirY")     vDir   = prescribedData[k].mData[0];
        else if(prescribedData[k].mNameVariable == "VelocityDirZ")     wDir   = prescribedData[k].mData[0];
      }
    }
  }

  // Check the situation.
  bool stateFromBoundaries = false;
  if((u < valLarge) && (v < valLarge) && (w < valLarge))
  {
    // Velocity is prescribed. Check if two out of three thermodynamic variables have
    // been specified. 
    int nVarSpecified = 0;
    if(p   > zero) ++nVarSpecified;
    if(rho > zero) ++nVarSpecified;
    if(T   > zero) ++nVarSpecified;

    // Check if at least two thermodynamic state variables have been specified.
    if(nVarSpecified >= 2)
    {
      // Reconstruct the pressure and density, if needed.
      if(p   <= zero) p   = rho*RGas*T;
      if(rho <= zero) rho = p/(RGas*T);

      // Set the dimensional primitive variables.
      primVar[0] = rho;
      primVar[1] = u;
      primVar[2] = v;
      primVar[3] = w;
      primVar[4] = p;

      // Set stateFromBoundaries to true.
      stateFromBoundaries = true;
    }
    else if(nVarSpecified == 1)
    {
      // Only one state variable specified. As long as this is not the temperature,
      // this may correspond to a subsonic inlet and a supersonic outlet, although
      // this is not adviceable.
      if(T <= zero)
      {
        // Compute the magnitude of the velocity squared.
        const su2double velMag2 = u*u + v*v + w*w;

        // Compute the pressure from the density and velocity magnitude (or vice versa).
        // Assume a Mach number of 0.95.
        if(p <= zero) p   = 0.95*rho*velMag2/GamConstant;
        else          rho = GamConstant*p/(0.95*velMag2);

        // Set the dimensional primitive variables.
        primVar[0] = rho;
        primVar[1] = u;
        primVar[2] = v;
        primVar[3] = w;
        primVar[4] = p;

        // Set stateFromBoundaries to true.
        stateFromBoundaries = true;
      }
    }
  }
  else if((uDir < valLarge) && (vDir < valLarge) && (wDir < valLarge))
  {
    // Velocity direction prescribed. Determine the corresponding unit vector.
    const su2double lenInv = one/SQRT(uDir*uDir + vDir*vDir + wDir*wDir);
    uDir *= lenInv;
    vDir *= lenInv;
    wDir *= lenInv;

    // Compute the value of gamma -1.
    const su2double gm1 = GamConstant - one;

    // Check if a Mach number can be extracted from a combination of a total and
    // a static quantity. If this is not the case, subsonic inflow with supersonic
    // outflow, assume a Mach number of 0.95.
    su2double M2 = (su2double) (0.95*0.95);

    if((pTot > zero) && (p > zero))
    {
      const su2double ratio = pTot/p;
      M2 = (POW(ratio,gm1/GamConstant)-one)*two/gm1;
    }
    else if((rhoTot > zero) && (rho > zero))
    {
      const su2double ratio = rhoTot/rho;
      M2 = (POW(ratio,one/GamConstant)-one)*two/gm1;
    }
    else if((TTot > zero) && (T > zero))
    {
      const su2double ratio = TTot/T;
      M2 = (ratio-one)*two/gm1;
    }

    // Check if at least two total conditions are prescribed.
    int nVarSpecified = 0;
    if(pTot   > zero) ++nVarSpecified;
    if(rhoTot > zero) ++nVarSpecified;
    if(TTot   > zero) ++nVarSpecified;

    if(nVarSpecified >= 2)
    {
      // Create the total pressure and total density, if not specified.
      if(pTot   <= zero) pTot   = rhoTot*RGas*TTot;
      if(rhoTot <= zero) rhoTot = pTot/(RGas*TTot);

      // Create the static pressure and density.
      const su2double tmp = one + half*gm1*M2;

      p   = pTot  /POW(tmp, GamConstant/gm1);
      rho = rhoTot/POW(tmp, one/gm1);

      // Compute the velocity magnitude.
      const su2double velMag = SQRT(M2*GamConstant*p/rho);

      // Set the dimensional primitive variables.
      primVar[0] = rho;
      primVar[1] = uDir*velMag;
      primVar[2] = vDir*velMag;
      primVar[3] = wDir*velMag;
      primVar[4] = p;

      // Set stateFromBoundaries to true.
      stateFromBoundaries = true;
    }
  }

	// Check if the primitive variables must be constructed from a quiescent flow condition.
	if( mQuiescentFlow )
	{
		// Set background temperature.
		su2double T; ( mTempFreeStream < zero ) ? T = 300.0 : T = mTempFreeStream;

		// Compute background density.
		const su2double rho = p/(RGas*T);

		// Set the dimensional primitive variables.
		primVar[0] = rho;
		primVar[1] = zero;
		primVar[2] = zero;
		primVar[3] = zero;
		primVar[4] = p;

		// Set stateFromBoundaries to true.
		stateFromBoundaries = true;
	}

  // Check if the primitive variables must be constructed from the specified
  // free stream conditions in the parameter file.
  if( !stateFromBoundaries )
  {
    // Compute the free stream velocity from the temperature and Mach number.
    const su2double vInf = mMach*SQRT(GamConstant*RGas*mTempFreeStream);

    // Compute the dimensional viscosity and the free stream density
    // from the Reynolds number.
    const su2double muDim  = mu*pRef/uRef;
    const su2double rhoInf = mReynolds*muDim/(mReynoldsLength*vInf);

    // Compute the free stream pressure from the perfect gas law.
    const su2double pInf = rhoInf*RGas*mTempFreeStream;

    // Set the dimensional primitive variables.
    primVar[0] = rhoInf;
    primVar[1] = vInf*mVelDirFreeStream[0];
    primVar[2] = vInf*mVelDirFreeStream[1];
    primVar[3] = vInf*mVelDirFreeStream[2];
    primVar[4] = pInf;
  }

  // Non-dimensionalize the primitive variables.
  primVar[0] /= rhoRef;
  primVar[1] /= uRef;
  primVar[2] /= uRef;
  primVar[3] /= uRef;
  primVar[4] /= pRef;
}

//------------------------------------------------------------------------------

// Function to read the prescribed boundary data.
void InputParamClass::ReadBoundaryData(const char *fileName)
{
  std::ostringstream message;

  // Open the file for reading. If not successful terminate.
  // file is created.
  std::ifstream boundaryDataFile(fileName);
  if( !boundaryDataFile )
  {
    message << "Boundary data file " << fileName << " not found.";
    TerminateAll("InputParamClass::ReadBoundaryData",
                 __FILE__, __LINE__, message.str());
  }

  // While loop to find the keyword "Prescribed boundary data".
  std::string lineBuf;
  bool stringFound = false;
  while( std::getline(boundaryDataFile, lineBuf) )
  {
    // Replace all tab and return characters, remove the leading
    // and trailing blanks and create a lower case version.
    ReplaceTabsAndReturns(lineBuf);
    RemoveBlanks(lineBuf);
    CreateLowerCase(lineBuf);

    // Check if the string equals "prescribed boundary data".
    // If so, set stringFound to true and break the loop.
    if(lineBuf == "prescribed boundary data")
    {
      stringFound = true;
      break;
    }
  }

  // Check if the subface information was found. If not return.
  if( !stringFound ) {boundaryDataFile.close(); return;}

  // The sequence of the prescribed data can be arbitrary.
  // Therefore use a while loop to extract the data.
  while( std::getline(boundaryDataFile, lineBuf) )
  {
    // Couple an istringstream to lineBuf.
    std::istringstream istr(lineBuf);

    // Attempt to read face ID, subface ID and number of points.
    int faceID, subfaceID, nPoints;
    if( (istr >> faceID >> subfaceID >> nPoints) )
    {
      // Subtract 1 from the face ID and subface ID to obtain the C-numbering.
      --faceID; --subfaceID;

      // Read the actual data.
      mSubfaces[faceID][subfaceID]->ReadBoundaryData(faceID, subfaceID, nPoints,
                                                     istr, boundaryDataFile);
    }
  }

  // Close the file again.
  boundaryDataFile.close();

  // Check if the expected boundary data is actually prescribed.
  for(unsigned int i=0; i<mSubfaces.size(); ++i)
    for(unsigned int j=0; j<mSubfaces[i].size(); ++j)
      mSubfaces[i][j]->CheckPrescribedData(i, j);
}

//------------------------------------------------------------------------------

// Function, which reads the input parameters from the given file.
void InputParamClass::ReadInputParameters(const char *fileName)
{
  // Open the file for reading. If not successful a template parameter
  // file is created.
  std::ifstream paramFile(fileName);
  if( !paramFile ) CreateTemplateParamFile(fileName);

  // Read the string as long as there is info in the parameter file.
  // Close the parameter file afterwards.
  std::string lineBuf;
  while( std::getline(paramFile, lineBuf) ) AnalyzeString(lineBuf);
  paramFile.close();

  // If a wall model used, initialize the wall model (set up y-values of points, etc)
  if(mWallModelType == EQUILIBRIUM_WALL_MODEL){
    mWallModel = new CWallModel1DEQ(this);
  }
  else if(mWallModelType == LOGARITHMIC_WALL_MODEL){
    mWallModel = new CWallModelLogLaw(this);
  }

  // If a subgrid scale model is used, initialize the model.
  if(     mSGSModelType == WALE_MODEL)   mSGSModel = new CWALEModel();
  else if(mSGSModelType == VREMAN_MODEL) mSGSModel = new CVremanModel();

  // Check the input parameters.
  CheckParameters();
}

//------------------------------------------------------------------------------

// Function to read the subface information.
void InputParamClass::ReadSubfaceInformation(const char *fileName)
{
  std::ostringstream message;

  // Open the file for reading. If not successful terminate.
  // file is created.
  std::ifstream subfaceFile(fileName);
  if( !subfaceFile ) 
  {
    message << "Subface file " << fileName << " not found.";
    TerminateAll("InputParamClass::ReadSubfaceInformation",
                 __FILE__, __LINE__, message.str());
  }

  // While loop to find the keyword "Subface information".
  std::string lineBuf;
  bool stringFound = false;
  while( std::getline(subfaceFile, lineBuf) )
  {
    // Replace all tab and return characters, remove the leading
    // and trailing blanks and create a lower case version.
    ReplaceTabsAndReturns(lineBuf);
    RemoveBlanks(lineBuf);
    CreateLowerCase(lineBuf);

    // Check if the string equals "subface information".
    // If so, set stringFound to true and break the loop.
    if(lineBuf == "subface information")
    {
      stringFound = true;
      break;
    }
  }

  // Check if the subface information was found. If not, terminate.
  if( !stringFound )
  {
    message << "Keyword \"Subface information\" not found in the file " << fileName << ".";
    TerminateAll("InputParamClass::ReadSubfaceInformation",
                 __FILE__, __LINE__, message.str());
  }

  // Set the boolean for wall subfaces present to false.
  bool wallSubfacesPresent = false;

  // Allocate the first index for mSubfaces, which is 6, the number of faces of a hexahedron.
  mSubfaces.resize(6);

  // Loop over the 6 faces of the block.
  for(int i=0; i<6; ++i)
  {
    // While loop until the information of the current face is found.
    stringFound = false;
    while( std::getline(subfaceFile, lineBuf) )
    {
      // Couple an istringstream to lineBuf.
      std::istringstream istr(lineBuf);

      // Check if this is the desired string. If so, read the number of subfaces, allocate
      // the memory mSubfaces[i] and set stringFound to true.
      int faceID;
      if( (istr >> faceID) )
      {
        if(faceID == (i+1))
        {
          int nSubfaces;
          istr >> nSubfaces;
          mSubfaces[i].resize(nSubfaces, NULL);
          stringFound = true;
        }
      }

      // Break the while loop if the string was found.
      if( stringFound ) break;
    }

    // Terminate if the subface information for this face was not found.
    if( !stringFound )
    {
      message << "Subface information for face " << i+1 << " not found in the file " << fileName << ".";
      TerminateAll("InputParamClass::ReadSubfaceInformation",
                   __FILE__, __LINE__, message.str());
    }

    // Loop over the subfaces of the current face.
    for(unsigned int j=0; j<mSubfaces[i].size(); ++j)
    {
      // Read the string with the subface information and couple it to an istringstream.
      std::getline(subfaceFile, lineBuf);
      std::istringstream istr(lineBuf);

      // Read the boundary condition and convert it to lower case.
      std::string BCSubface, BCSubfaceOr;
      istr >> BCSubfaceOr;
      BCSubface = BCSubfaceOr;
      CreateLowerCase(BCSubface);

      // Determine the boundary condition and allocate the correct class.
      if(     BCSubface == "internal1to1")            mSubfaces[i][j] = new Internal1to1SubfaceClass(istr);
      else if(BCSubface == "periodic1to1trans")       mSubfaces[i][j] = new Periodic1to1TransSubfaceClass(istr);
      else if(BCSubface == "bcfarfield")              mSubfaces[i][j] = new BCFarfieldSubfaceClass(istr);
      else if(BCSubface == "bcwallisothermal")        mSubfaces[i][j] = new BCIsothermalWallSubfaceClass(istr);
      else if(BCSubface == "bcwallisothermalwm")      mSubfaces[i][j] = new BCIsothermalWallSubfaceClass(istr, mWallModelType);
      else if(BCSubface == "bcwallheatflux")          mSubfaces[i][j] = new BCHeatFluxWallSubfaceClass(istr);
      else if(BCSubface == "bcwallheatfluxwm")        mSubfaces[i][j] = new BCHeatFluxWallSubfaceClass(istr, mWallModelType);
      else if(BCSubface == "bcwallinviscid")          mSubfaces[i][j] = new BCInviscidWallSubfaceClass(istr);
      else if(BCSubface == "bcsymmetryplane")         mSubfaces[i][j] = new BCSymmetrySubfaceClass(istr);
      else if(BCSubface == "bcinflowsubsonic")        mSubfaces[i][j] = new BCInflowSubsonicSubfaceClass(istr);
      else if(BCSubface == "bcinflowsupersonic")      mSubfaces[i][j] = new BCInflowSupersonicSubfaceClass(istr);
      else if(BCSubface == "bcoutflowsubsonic")       mSubfaces[i][j] = new BCOutflowSubsonicSubfaceClass(istr);
      else if(BCSubface == "bcoutflowsupersonic")     mSubfaces[i][j] = new BCOutflowSupersonicSubfaceClass(istr);
			else if(BCSubface == "bcinflowcharacteristic")  mSubfaces[i][j] = new BCInflowCharacteristicSubfaceClass(istr);
			else if(BCSubface == "bcoutflowcharacteristic") mSubfaces[i][j] = new BCOutflowCharacteristicSubfaceClass(istr);
			else
      {
        message << "Unknown boundary condition \"" << BCSubfaceOr << "\".";
        TerminateAll("InputParamClass::ReadSubfaceInformation",
                     __FILE__, __LINE__, message.str());
      }

      // Set wallSubfacesPresent if this is a wall.
      if( mSubfaces[i][j]->BCIsWall() ) wallSubfacesPresent = true;

			// Book-keep the iFace and jSubface indices.
			mSubfaces[i][j]->mBoundaryID = i;
			mSubfaces[i][j]->mSubfacesID = j;

			// Set NSCBC condition to true, if characteristic inlet/outlet are specified.
			if( mSubfaces[i][j]->GetTypeBoundaryPrescribed() == BC_OUTFLOW_CHARACTERISTIC ) mNSCBC_Specified = true;
			if( mSubfaces[i][j]->GetTypeBoundaryPrescribed() == BC_INFLOW_CHARACTERISTIC  ) mNSCBC_Specified = true;
    }
  }

  // Set mMonitorForces to false if no walls are present.
  if( !wallSubfacesPresent ) mMonitorForces = false;
}

//------------------------------------------------------------------------------

// Function, which writes the header for the time step history.
void InputParamClass::WriteHeaderTimeStepHistory(void)
{
  if(rank == 0)
  {
#ifdef HAVE_OPENMP
#pragma omp single nowait
#endif
    {
      // Write the header line.
      std::cout << "#-------------------------------------------------";
      if( mSGSModel      ) std::cout << "-------------";
      if( mMonitorForces ) std::cout << "----------------------------------------------------";
      std::cout << "|" << std::endl;

      // Write the first line of the monitoring names.
      std::cout << "#   Time   | Wall time  |  Physical  |  Maximum   |";
      if( mSGSModel      ) std::cout << "   Maximum  |";
      if( mMonitorForces ) std::cout << "     CL     |     CD     |     CL     |     CD     |";
      std::cout << std::endl;

      // Write the second line of the monitoring names.
      std::cout << "#   Step   |     (s)    |  time (s)  |    Mach    |";
      if( mSGSModel      ) std::cout << "muTurb/muLam|";
      if( mMonitorForces ) std::cout << "    Total   |    Total   |   Viscous  |   Viscous  |";
      std::cout << std::endl;

      // Write the trailing line, which is identical to the header line.
      std::cout << "#-------------------------------------------------";
      if( mSGSModel      ) std::cout << "-------------";
      if( mMonitorForces ) std::cout << "----------------------------------------------------";
      std::cout << "|" << std::endl;
    }
  }
}
