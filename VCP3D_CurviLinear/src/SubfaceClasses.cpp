//------------------------------------------------------------------------------
// File, which contains the implementation of the member functions of the
// subface classes.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"

//------------------------------------------------------------------------------

// Constructor. Nothing to be done.
PrescribedDataClass::PrescribedDataClass(){}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
PrescribedDataClass::~PrescribedDataClass(){}

//------------------------------------------------------------------------------

// Overloaded constructor of SubfaceBaseClass.
SubfaceBaseClass::SubfaceBaseClass(std::istringstream &istr)
{
  // Read the subface range. Substract one, because C numbering is used.
  istr >> mElemIBeg >> mElemIEnd >> mElemJBeg >> mElemJEnd >> mElemKBeg >> mElemKEnd;

  --mElemIBeg; --mElemIEnd;
  --mElemJBeg; --mElemJEnd;
  --mElemKBeg; --mElemKEnd;

	// Initialize the boundary and subface ID to -ve. This way, in case they 
	// are undefined, an error is thrown.
	mBoundaryID = mSubfacesID = -1;
}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
SubfaceBaseClass::~SubfaceBaseClass(){}

//------------------------------------------------------------------------------

// Virtual function, which sets the value of mMachAverageBoundary.
// It may be overwritten by the derived class.
void SubfaceBaseClass::SetMachAverageBoundary(const su2double Mavg) {}

//------------------------------------------------------------------------------

// Function, which configures the tuning parameters required for a NSCBC.
// It is overwritten, if NSCBC is used.
void SubfaceBaseClass::ConfigureParamNSCBC(const InputParamClass 			*inputParam, 
												 									 const StandardElementClass *standardHex,
												 									 su2double            		 **metric,
												 									 su2double                   factNorm) {}

//------------------------------------------------------------------------------

// Virtual function, which computes the weighted Mach number on this
// local element boundary that is used to compute the average.
// The average Mach number on this element surface is also computed.
su2double SubfaceBaseClass::WeightedMachElement(const InputParamClass *inputParam,
                                                const int              nInt,
					 					                            su2double            **sol,
					 					                            const su2double        factNorm,
					 					                            su2double            **metric,
					 					                            su2double             *intWeights) {return zero;}

//------------------------------------------------------------------------------

// Virtual function, which indicates whether or not a boundary condition
// is physical. It may be overwritten by the derived class.
bool SubfaceBaseClass::BCIsPhysical(void)
{
  // Return the default value, which is true.
  return true;
}

//------------------------------------------------------------------------------

// Virtual function, which indicates whether or not a boundary condition
// is a wall boundary condition. It should be overwritten by the
// wall boundary condition classes.
bool SubfaceBaseClass::BCIsWall(void)
{
  // Return the default value, which is false.
  return false;
}

//------------------------------------------------------------------------------

// Function, which checks if the expected data is prescribed.
void SubfaceBaseClass::CheckPrescribedData(const int faceID,
                                           const int subfaceID)
{
  // Get the possible number of sets of variables that can be expected.
  // If nothing needs to be prescribed, return immediately.
  const int nSets = this->GetNSetsPrescribedVariables();
  if(nSets == 0) return;

  // Loop over the number of possible sets.
  bool varSetFound = false;
  int  setFound    = -1;
  for(int i=0; i<nSets; ++i)
  {
    // Get the names of the required variables for this set.
    std::vector<std::string> varNames;
    this->GetNamesPrescibedVariables(i, varNames);

    // Loop over the required names.
    unsigned int nVarProvided = 0;
    for(unsigned int j=0; j<varNames.size(); ++j)
    {
      // Check if the name is present in the provided data.
      for(unsigned int k=0; k<mPrescribedData.size(); ++k)
        if(varNames[j] == mPrescribedData[k].mNameVariable) ++nVarProvided;
    }

    // Break the loop if the variable set is provided.
    if(nVarProvided == varNames.size())
    {
      varSetFound = true;
      setFound    = i;
      break;
    }
  }

  // Terminate if no variable set was found.
  if( !varSetFound )
  {
    std::ostringstream message;
    message << "Face " << faceID+1 << ", subface " << subfaceID+1
            << ": Provided variable set not sufficient.";
    TerminateAll("SubfaceBaseClass::CheckPrescribedData",
                 __FILE__, __LINE__, message.str());
  }

  // Set the indices for the prescribed data, if appropriate.
  this->SetIndicesPrescribedData(setFound);
}

//------------------------------------------------------------------------------

// Function, which converts the prescribed data to the required form.
// It should be overwritten be the derived classes, if this is needed.
// Therefore this function does not do anything.
void SubfaceBaseClass::ConvertPrescribedData(
                              const int                nIntegration,
                              const ENUM_FEM_VARIABLES FEMVariables,
                              std::vector<su2double *> &prescribedData){}

//------------------------------------------------------------------------------

// Function, which indicates whether prescribed data is expected.
// It may be overwritten be the derived classes.
bool SubfaceBaseClass::ExpectPrescribedData(void){return false;}

//------------------------------------------------------------------------------

// Function, which gets the index for in mPrescribedData of the given variable.
int SubfaceBaseClass::GetIndexPrescribedData(const char *varName)
{
  // Loop over the prescribed data and return the correct index.
  for(unsigned int k=0; k<mPrescribedData.size(); ++k)
    if(varName == mPrescribedData[k].mNameVariable) return ((int) k);

  // This point should not be reached.
  TerminateAll("SubfaceBaseClass::GetIndexPrescribedData", __FILE__, __LINE__,
               "This should not happen");

  // Add a return statement to avoid a compiler warning.
  return -1;
}

//------------------------------------------------------------------------------

// Function, which gets the names of the prescribed variables for the given set.
// It should be overwritten be the derived classes.
void SubfaceBaseClass::GetNamesPrescibedVariables(const int                set,
                                                  std::vector<std::string> &varNames)
{
  TerminateAll("SubfaceBaseClass::GetNamesPrescibedVariables", __FILE__, __LINE__,
               "This function must be overwritten by the derived class");
}

//------------------------------------------------------------------------------

// Function, which returns the number of data sets that can be prescribed.
// It may be overwritten be the derived classes.
int SubfaceBaseClass::GetNSetsPrescribedVariables(void){return 0;}

//------------------------------------------------------------------------------

// Virtual function, which makes available the beginning i-index
// of the element donor range. The default value is returned.
int SubfaceBaseClass::GetElemDonorIBeg(void) {return -10000;}

//------------------------------------------------------------------------------

// Virtual function, which makes available the beginning j-index
// of the element donor range. The default value is returned.
int SubfaceBaseClass::GetElemDonorJBeg(void) {return -10000;}

//------------------------------------------------------------------------------

// Virtual function, which makes available the beginning k-index
// of the element donor range. The default value is returned.
int SubfaceBaseClass::GetElemDonorKBeg(void) {return -10000;}

//------------------------------------------------------------------------------

// Virtual function, which makes available the number of prescribed
// variables. It must be overwritten by the physical boundary conditions.
int SubfaceBaseClass::GetNVarPrescribed(void)
{
  Terminate("SubfaceBaseClass::GetNVarPrescribed", __FILE__, __LINE__,
            "Function must be overloaded by the derived class");
  return 0;
}

//------------------------------------------------------------------------------

// Function, which makes the possible periodic transformation available.
// This base class function does not do anything.
void SubfaceBaseClass::GetPeriodicTransformation(su2double *trans){}

//------------------------------------------------------------------------------
 
// Function to make available the prescribed data as a const reference.
const std::vector<PrescribedDataClass> &SubfaceBaseClass::GetPrescribedData(void)
{
  return mPrescribedData;
}

//------------------------------------------------------------------------------

// Virtual function, which indicates whether or not this subface is 1 to 1 matching.
// It may be overwritten by the derived class, as this function sets the default.
bool SubfaceBaseClass::Is1To1Matching(void)
{
  return false;
}

//------------------------------------------------------------------------------

// Function, which determines the prescribed data from the free stream
// data. It should be overwritten be the derived classes.
void SubfaceBaseClass::PrescribedDataFromFreeStream(
                                 const int                nIntegration,
                                 const su2double          *nonDimPrimVarFreeStream,
                                 const ENUM_FEM_VARIABLES FEMVariables,
                                 std::vector<su2double *> &prescribedData)
{
  Terminate("SubfaceBaseClass::PrescribedDataFromFreeStream", __FILE__, __LINE__,
            "This function must be overwritten by the derived classes");
}

//------------------------------------------------------------------------------

// Function, which determines the prescribed data in the
// given number of integration points.
void SubfaceBaseClass::PrescribedDataIntegrationPoints(
                                 const int                nIntegration,
                                 const su2double          *x,
                                 const su2double          *y,
                                 const su2double          *z,
                                 const su2double          *nonDimPrimVarFreeStream,
                                 const ENUM_FEM_VARIABLES FEMVariables,
                                 std::vector<su2double *> &prescribedData)
{
  // Check if any data is prescribed at all. If not, the boundary state
  // must be determined from the free stream values.
  if(mPrescribedData.size() == 0)
    this->PrescribedDataFromFreeStream(nIntegration, nonDimPrimVarFreeStream,
                                       FEMVariables, prescribedData);
  else
  {
    // Data is prescribed. Check if constant data is prescribed.
    if(mPrescribedData[0].mData.size() == 1)
    {
      // Constant data. Loop over the elements of mIndicesPrescribedData
      // to store the data in the correct way.
      for(unsigned int l=0; l<mIndicesPrescribedData.size(); ++l)
      {
        const int ll = mIndicesPrescribedData[l];
        for(int i=0; i<nIntegration; ++i)
          prescribedData[l][i] = mPrescribedData[ll].mData[0];
      }
    }
    else
    {
      // Varying data prescribed.
      Terminate("SubfaceBaseClass::PrescribedDataIntegrationPoints", __FILE__, __LINE__,
                "Interpolation of varying data not implemented yet");
    }

    // Convert the prescribed data to the required form, if needed.
    this->ConvertPrescribedData(nIntegration, FEMVariables, prescribedData);
  }
}

//------------------------------------------------------------------------------

// Function, which reads the prescribed data for this subface.
void SubfaceBaseClass::ReadBoundaryData(const int          faceID,
                                        const int          subfaceID,
                                        const int          nPoints,
                                        std::istringstream &istrVarNames,
                                        std::ifstream      &boundaryDataFile)
{
  std::ostringstream message;

  // Check if data is actually expected for this subface.
  if( !(this->ExpectPrescribedData()) )
  {
    message << "Face " << faceID+1 << ", subface " << subfaceID+1
            << ": Boundary data is prescribed, but this is not expected.";
    TerminateAll("SubfaceBaseClass::ReadBoundaryData",
                 __FILE__, __LINE__, message.str());
  }

  // Read the variables names from istrVarNames.
  std::string lineBuf;
  while( istrVarNames >> lineBuf )
  {
    PrescribedDataClass newData;
    newData.mNameVariable = lineBuf;
    mPrescribedData.push_back(newData);
  }

  // Loop over the number of prescribed variables and allocate the memory
  // to store the data.
  for(unsigned int j=0; j<mPrescribedData.size(); ++j)
    mPrescribedData[j].mData.resize(nPoints);

  // Loop over the number of points and read the prescribed data.
  for(int i=0; i<nPoints; ++i)
  {
    // Read the string.
    if( !(std::getline(boundaryDataFile, lineBuf)) )
    {
      message << "Face " << faceID+1 << ", subface " << subfaceID+1
              << ": Unexpected end of the file while reading prescribed data";
      TerminateAll("SubfaceBaseClass::ReadBoundaryData",
                   __FILE__, __LINE__, message.str());
    }

    // Couple an istringstream to lineBuf.
    std::istringstream istr(lineBuf);

    // Loop over the number of variables and read the data.
    for(unsigned int j=0; j<mPrescribedData.size(); ++j)
    {
      if( !(istr >> mPrescribedData[j].mData[i]) )
      {
        message << "Face " << faceID+1 << ", subface " << subfaceID+1
                << ": Something wrong while reading prescribed data";
        TerminateAll("SubfaceBaseClass::ReadBoundaryData",
                     __FILE__, __LINE__, message.str());
      }
    }
  }
}

//------------------------------------------------------------------------------

// Virtual function, which sets the indices for the prescribed data.
// It should be overwritten be the derived classes, if this is needed.
// Therefore this function does not do anything.
void SubfaceBaseClass::SetIndicesPrescribedData(const int set){}

//------------------------------------------------------------------------------

// Virtual function, which indicates whether or not a wall model is used
// for the boundary. It may be overwritten by the derived class.
bool SubfaceBaseClass::UseWallModel(void)
{
  // Return the default value, which is false.
  return false;
}

//------------------------------------------------------------------------------

// Virtual function, which gets the type of boundary prescribed. 
// It should be overwritten by the derived class.
int SubfaceBaseClass::GetTypeBoundaryPrescribed(void) {return BC_UNDEFINED_TYPE;}

//------------------------------------------------------------------------------

// Virtual function, which determines whether or not the boundary is of 
// a characteristic nature (inlet or outlet).
bool SubfaceBaseClass::IsCharacteristicBC(void) {return false;}

//------------------------------------------------------------------------------

// Overloaded constructor of Internal1to1SubfaceClass.
Internal1to1SubfaceClass::Internal1to1SubfaceClass(std::istringstream &istr)
 : SubfaceBaseClass(istr)
{
  // Read the donor range. Substract one, because C numbering is used.
  istr >> mElemDonorIBeg >> mElemDonorIEnd >> mElemDonorJBeg >> mElemDonorJEnd
       >> mElemDonorKBeg >> mElemDonorKEnd;

  --mElemDonorIBeg; --mElemDonorIEnd;
  --mElemDonorJBeg; --mElemDonorJEnd;
  --mElemDonorKBeg; --mElemDonorKEnd;

  // Read the short hand of the transformation matrix.
  istr >> mL1 >> mL2 >> mL3;

  // Check if this corresponds to the identity transformation matrix.
  if((mL1 != 1) || (mL2 != 2) || (mL3 !=3))
    TerminateAll("Internal1to1SubfaceClass::Internal1to1SubfaceClass", __FILE__, __LINE__,
                 "Only identity transformation matrices are supported");

  // Check if the number of elements and number of donor elements match.
  const int nElemI = mElemIEnd - mElemIBeg;
  const int nElemJ = mElemJEnd - mElemJBeg;
  const int nElemK = mElemKEnd - mElemKBeg;

  const int nDonorElemI = mElemDonorIEnd - mElemDonorIBeg;
  const int nDonorElemJ = mElemDonorJEnd - mElemDonorJBeg;
  const int nDonorElemK = mElemDonorKEnd - mElemDonorKBeg;

  if((nElemI != nDonorElemI) || (nElemJ != nDonorElemJ) || (nElemK != nDonorElemK))
  {
    std::ostringstream message;
    message << "Subface " << mElemIBeg+1 << " " << mElemIEnd+1 << " "
                          << mElemJBeg+1 << " " << mElemJEnd+1 << " "
                          << mElemKBeg+1 << " " << mElemKEnd+1 << ": "
                          << "Donor range does not match subface range";
    TerminateAll("Internal1to1SubfaceClass::Internal1to1SubfaceClass",
                 __FILE__, __LINE__, message.str());
                          
  }
}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
Internal1to1SubfaceClass::~Internal1to1SubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which indicates whether or not a boundary condition is physical.
bool Internal1to1SubfaceClass::BCIsPhysical(void)
{
  // Return false, as an internal 1 to 1 boundary is not physical.
  return false;
}

//------------------------------------------------------------------------------

// Function, which makes available the beginning i-index
// of the element donor range.
int Internal1to1SubfaceClass::GetElemDonorIBeg(void) {return mElemDonorIBeg;}

//------------------------------------------------------------------------------

// Function, which makes available the beginning j-index
// of the element donor range.
int Internal1to1SubfaceClass::GetElemDonorJBeg(void) {return mElemDonorJBeg;}

//------------------------------------------------------------------------------

// Function, which makes available the beginning k-index
// of the element donor range.
int Internal1to1SubfaceClass::GetElemDonorKBeg(void) {return mElemDonorKBeg;}

//------------------------------------------------------------------------------

// Function, which indicates whether or not this subface is 1 to 1 matching.
// It is still virtual, because it is overwritten by Periodic1to1TransSubfaceClass.
bool Internal1to1SubfaceClass::Is1To1Matching(void)
{
  // Return true, as an internal 1 to 1 boundary is 1 to 1 matching.
  return true;
}

//------------------------------------------------------------------------------

// Function, which gets the type of boundary prescribed. 
int Internal1to1SubfaceClass::GetTypeBoundaryPrescribed(void) {return INTERNAL_1TO1;}

//------------------------------------------------------------------------------

// Overloaded constructor of Periodic1to1TransSubfaceClass.
Periodic1to1TransSubfaceClass::Periodic1to1TransSubfaceClass(std::istringstream &istr)
 :  Internal1to1SubfaceClass(istr)
{
  // Read the periodic translation.
  istr >> mPeriodicTranslation[0] >> mPeriodicTranslation[1] >> mPeriodicTranslation[2];
}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
Periodic1to1TransSubfaceClass::~Periodic1to1TransSubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which makes the periodic transformation available.
void Periodic1to1TransSubfaceClass::GetPeriodicTransformation(su2double *trans)
{
  trans[0] = mPeriodicTranslation[0];
  trans[1] = mPeriodicTranslation[1];
  trans[2] = mPeriodicTranslation[2];
}

//------------------------------------------------------------------------------

// Function, which indicates whether or not this subface is 1 to 1 matching
// and not periodic, which is false.
bool Periodic1to1TransSubfaceClass::Is1To1Matching(void)
{
  return false;
}

//------------------------------------------------------------------------------

// Function, which gets the type of boundary prescribed. 
int Periodic1to1TransSubfaceClass::GetTypeBoundaryPrescribed(void) {return PERIODIC_1TO1_TRANS;}

//------------------------------------------------------------------------------

// Overloaded constructor of BCFarfieldSubfaceClass.
BCFarfieldSubfaceClass::BCFarfieldSubfaceClass(std::istringstream &istr)
 : SubfaceBaseClass(istr){}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
BCFarfieldSubfaceClass::~BCFarfieldSubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which makes available the number of prescribed variables.
int BCFarfieldSubfaceClass::GetNVarPrescribed(void) {return nVar;}

//------------------------------------------------------------------------------

// Function, which determines the prescribed data from the free stream.
void BCFarfieldSubfaceClass::PrescribedDataFromFreeStream(
                             const int                nIntegration,
                             const su2double          *nonDimPrimVarFreeStream,
                             const ENUM_FEM_VARIABLES FEMVariables,
                             std::vector<su2double *> &prescribedData)
{
  // Easier storage of the free stream variables.
  const su2double rho = nonDimPrimVarFreeStream[0];
  const su2double u   = nonDimPrimVarFreeStream[1];
  const su2double v   = nonDimPrimVarFreeStream[2];
  const su2double w   = nonDimPrimVarFreeStream[3];
  const su2double p   = nonDimPrimVarFreeStream[4];

  // Make a distinction between the working variables.
  switch( FEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      // Conservative variables are used. Compute the free stream
      // conservative variables.
      const su2double ru = rho*u;
      const su2double rv = rho*v;
      const su2double rw = rho*w;
      const su2double rE = p/(GamConstant-one) + half*rho*(u*u + v*v + w*w);

      // Loop over the integration points and set the value.
#pragma omp simd
      for(int l=0; l<nIntegration; ++l)
      {
        prescribedData[0][l] = rho;
        prescribedData[1][l] = ru;
        prescribedData[2][l] = rv;
        prescribedData[3][l] = rw;
        prescribedData[4][l] = rE;
      }

      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      // Entropy variables are used. Compute the free stream entropy variables.
      const su2double ovgm1 = one/(GamConstant - one);
      const su2double s     = LOG(p/POW(rho,GamConstant));
      const su2double pInv  = one/p;

      const su2double V0 =  (GamConstant-s)*ovgm1 - half*pInv*rho*(u*u + v*v + w*w);
      const su2double V1 =  rho*u*pInv;
      const su2double V2 =  rho*v*pInv;
      const su2double V3 =  rho*w*pInv;
      const su2double V4 = -rho*pInv;

      // Loop over the integration points and set the value.
#pragma omp simd
      for(int l=0; l<nIntegration; ++l)
      {
        prescribedData[0][l] = V0;
        prescribedData[1][l] = V1;
        prescribedData[2][l] = V2;
        prescribedData[3][l] = V3;
        prescribedData[4][l] = V4;
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

// Function, which gets the type of boundary prescribed. 
int BCFarfieldSubfaceClass::GetTypeBoundaryPrescribed(void) {return BC_FARFIELD;}

//------------------------------------------------------------------------------

// Overloaded constructor of BCIsothermalWallSubfaceClass.
BCIsothermalWallSubfaceClass::BCIsothermalWallSubfaceClass(std::istringstream &istr)
 : SubfaceBaseClass(istr)
{
  mUseWallModel = false;
}

//------------------------------------------------------------------------------

// Overloaded constructor of BCIsothermalWallSubfaceClass.
BCIsothermalWallSubfaceClass::BCIsothermalWallSubfaceClass(std::istringstream &istr,
                                                           ENUM_WALL_MODEL    wallModelType)
 : SubfaceBaseClass(istr)
{
  if(wallModelType == NO_WALL_MODEL) mUseWallModel = false;
  else                               mUseWallModel = true;
}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
BCIsothermalWallSubfaceClass::~BCIsothermalWallSubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which indicates whether or not a boundary condition
// is a wall boundary condition. True for this BC.
bool BCIsothermalWallSubfaceClass::BCIsWall(void){return true;}

//------------------------------------------------------------------------------

// Function, which indicates whether prescribed data is expected.
bool BCIsothermalWallSubfaceClass::ExpectPrescribedData(void){return true;}

//------------------------------------------------------------------------------

// Function, which gets the names of the prescribed variables for the given set.
void BCIsothermalWallSubfaceClass::GetNamesPrescibedVariables(const int                set,
                                                              std::vector<std::string> &varNames)
{
  // For an isothermal wall, the temperature must be provided.
  varNames.push_back("Temperature");
}

//------------------------------------------------------------------------------

// Function, which returns the number of data sets that can be prescribed.
int BCIsothermalWallSubfaceClass::GetNSetsPrescribedVariables(void){return 1;}

//------------------------------------------------------------------------------

// Function, which makes available the number of prescribed variables.
int BCIsothermalWallSubfaceClass::GetNVarPrescribed(void) {return 1;}

//------------------------------------------------------------------------------

// Function, which determines the indices for the prescribed data.
void BCIsothermalWallSubfaceClass::SetIndicesPrescribedData(const int set)
{
  // Set the index for the temperature.
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("Temperature"));
}

//------------------------------------------------------------------------------

// Function, which indicates whether or not a wall model is used for the boundary.
bool BCIsothermalWallSubfaceClass::UseWallModel(void)
{
  // Return the value of mUseWallModel.
  return mUseWallModel;
}

//------------------------------------------------------------------------------

// Function, which gets the type of boundary prescribed. 
int BCIsothermalWallSubfaceClass::GetTypeBoundaryPrescribed(void) {return BC_ISOTHERMAL_WALL;}

//------------------------------------------------------------------------------

// Overloaded constructor of BCHeatFluxWallSubfaceClass.
BCHeatFluxWallSubfaceClass::BCHeatFluxWallSubfaceClass(std::istringstream &istr)
 : SubfaceBaseClass(istr)
{
  mUseWallModel = false;
}

//------------------------------------------------------------------------------

// Overloaded constructor of BCHeatFluxWallSubfaceClass.
BCHeatFluxWallSubfaceClass::BCHeatFluxWallSubfaceClass(std::istringstream &istr,
                                                       ENUM_WALL_MODEL    wallModelType)
 : SubfaceBaseClass(istr)
{
  if(wallModelType == NO_WALL_MODEL) mUseWallModel = false;
  else                               mUseWallModel = true;
}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
BCHeatFluxWallSubfaceClass::~BCHeatFluxWallSubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which indicates whether or not a boundary condition
// is a wall boundary condition. True for this BC.
bool BCHeatFluxWallSubfaceClass::BCIsWall(void){return true;}

//------------------------------------------------------------------------------

// Function, which converts the prescribed data to the required form,
// i.e. the heat flux is made non-dimensioal.
void BCHeatFluxWallSubfaceClass::ConvertPrescribedData(
                                        const int                nIntegration,
                                        const ENUM_FEM_VARIABLES FEMVariables,
                                        std::vector<su2double *> &prescribedData)
{
  // Determine the conversion factor to a non-dimensional heat flux.
  const su2double qRefInv = one/(pRef*uRef);

  // Convert the dimensional heat flux to a non-dimensional heat flux.
#pragma omp simd
  for(int l=0; l<nIntegration; ++l)
    prescribedData[0][l] *= qRefInv;
}

//------------------------------------------------------------------------------

// Function, which indicates whether prescribed data is expected.
bool BCHeatFluxWallSubfaceClass::ExpectPrescribedData(void){return true;}

//------------------------------------------------------------------------------

// Function, which gets the names of the prescribed variables for the given set.
void BCHeatFluxWallSubfaceClass::GetNamesPrescibedVariables(const int                set,
                                                            std::vector<std::string> &varNames)
{
  // For a heat flux wall, the heat flux must be provided.
  varNames.push_back("HeatFlux");
}

//------------------------------------------------------------------------------

// Function, which returns the number of data sets that can be prescribed.
int BCHeatFluxWallSubfaceClass::GetNSetsPrescribedVariables(void){return 1;}

//------------------------------------------------------------------------------

// Function, which makes available the number of prescribed variables.
int BCHeatFluxWallSubfaceClass::GetNVarPrescribed(void) {return 1;}

//------------------------------------------------------------------------------

// Function, which determines the indices for the prescribed data.
void BCHeatFluxWallSubfaceClass::SetIndicesPrescribedData(const int set)
{
  // Set the index for the heat flux.
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("HeatFlux"));
}

//------------------------------------------------------------------------------

// Function, which indicates whether or not a wall model is used for the boundary.
bool BCHeatFluxWallSubfaceClass::UseWallModel(void)
{
  // Return the value of mUseWallModel.
  return mUseWallModel;
}

//------------------------------------------------------------------------------

// Function, which gets the type of boundary prescribed. 
int BCHeatFluxWallSubfaceClass::GetTypeBoundaryPrescribed(void) {return BC_HEATFLUX_WALL;}

//------------------------------------------------------------------------------

// Overloaded constructor of BCInviscidWallSubfaceClass.
BCInviscidWallSubfaceClass::BCInviscidWallSubfaceClass(std::istringstream &istr)
 : SubfaceBaseClass(istr){}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
BCInviscidWallSubfaceClass::~BCInviscidWallSubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which indicates whether or not a boundary condition
// is a wall boundary condition. True for this BC.
bool BCInviscidWallSubfaceClass::BCIsWall(void){return true;}

//------------------------------------------------------------------------------

// Function, which makes available the number of prescribed variables.
int BCInviscidWallSubfaceClass::GetNVarPrescribed(void) {return 0;}

//------------------------------------------------------------------------------

// Function, which gets the type of boundary prescribed. 
int BCInviscidWallSubfaceClass::GetTypeBoundaryPrescribed(void) {return BC_INVISCID_WALL;}

//------------------------------------------------------------------------------

// Overloaded constructor of BCSymmetrySubfaceClass.
BCSymmetrySubfaceClass::BCSymmetrySubfaceClass(std::istringstream &istr)
 : SubfaceBaseClass(istr){}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
BCSymmetrySubfaceClass::~BCSymmetrySubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which makes available the number of prescribed variables.
int BCSymmetrySubfaceClass::GetNVarPrescribed(void) {return 0;}

//------------------------------------------------------------------------------

// Function, which gets the type of boundary prescribed. 
int BCSymmetrySubfaceClass::GetTypeBoundaryPrescribed(void) {return BC_SYMMETRY;}

//------------------------------------------------------------------------------

// Overloaded constructor of BCInflowSubsonicSubfaceClass.
BCInflowSubsonicSubfaceClass::BCInflowSubsonicSubfaceClass(std::istringstream &istr)
 : SubfaceBaseClass(istr){}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
BCInflowSubsonicSubfaceClass::~BCInflowSubsonicSubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which converts the prescribed data to the required form.
void BCInflowSubsonicSubfaceClass::ConvertPrescribedData(
                             const int                nIntegration,
                             const ENUM_FEM_VARIABLES FEMVariables,
                             std::vector<su2double *> &prescribedData)
{
  // Compute the inverse of the reference conditions.
  const su2double rhoRefInv = one/rhoRef;
  const su2double uRefInv   = one/uRef;
  const su2double pRefInv   = one/pRef;

  // Determine which set of variables have been prescribed.
  if( mTotalConditionsSpecified )
  {
    // Total conditions specified. Compute the conversion factor from the
    // dimension total temperature to the non-dimensional total enthalpy.
    const su2double Ttot2HtotNonDim = Cp*uRefInv*uRefInv;

    // Loop over the integration points.
#pragma omp simd
    for(int l=0; l<nIntegration; ++l)
    {
      // Compute the non-dimensional total pressure and total enthalpy.
      prescribedData[0][l] *= pRefInv;
      prescribedData[1][l] *= Ttot2HtotNonDim;

      // Make sure the unit vector really is a unit vector.
      const su2double len2   = prescribedData[2][l]*prescribedData[2][l]
                             + prescribedData[3][l]*prescribedData[3][l]
                             + prescribedData[4][l]*prescribedData[4][l];
      const su2double lenInv = one/SQRT(len2);

      prescribedData[2][l] *= lenInv;
      prescribedData[3][l] *= lenInv;
      prescribedData[4][l] *= lenInv;
    }
  }
  else
  {
    // Density and velocities are specified. Loop over the integration
    // points and non-dimensionalize them.
#pragma omp simd
    for(int l=0; l<nIntegration; ++l)
    {
      prescribedData[0][l] *= rhoRefInv;
      prescribedData[1][l] *= uRefInv;
      prescribedData[2][l] *= uRefInv;
      prescribedData[3][l] *= uRefInv;
    }
  }
}

//------------------------------------------------------------------------------

// Function, which indicates whether prescribed data is expected.
bool BCInflowSubsonicSubfaceClass::ExpectPrescribedData(void){return true;}

//------------------------------------------------------------------------------

// Function, which gets the names of the prescribed variables for the given set.
void BCInflowSubsonicSubfaceClass::GetNamesPrescibedVariables(const int                set,
                                                              std::vector<std::string> &varNames)
{
  // For a subsonic inflow boundary, there are two options for the prescribed
  // data. Make a distinction between the sets and set the data accordingly.
  switch( set )
  {
    case 0:
    {
      // The first set of variables. Total data and velocity direction are prescribed.
      varNames.push_back("TotalPressure");
      varNames.push_back("TotalTemperature");
      varNames.push_back("VelocityDirX");
      varNames.push_back("VelocityDirY");
      varNames.push_back("VelocityDirZ");

      break;
    }

    case 1:
    {
      // The second set of variables. Density and velocities are prescribed.
      varNames.push_back("Density");
      varNames.push_back("VelocityX");
      varNames.push_back("VelocityY");
      varNames.push_back("VelocityZ");

      break;
    }

    default:
      TerminateAll("BCInflowSubsonicSubfaceClass::GetNamesPrescibedVariables",
                   __FILE__, __LINE__, "This should not happen");
  }
}

//------------------------------------------------------------------------------

// Function, which returns the number of data sets that can be prescribed.
int BCInflowSubsonicSubfaceClass::GetNSetsPrescribedVariables(void){return 2;}

//------------------------------------------------------------------------------

// Function, which makes available the number of prescribed variables.
int BCInflowSubsonicSubfaceClass::GetNVarPrescribed(void)
{
  // If the total conditions are specified, nVar variables must be
  // prescribed, otherwise nVar-1.
  if( mTotalConditionsSpecified ) return nVar;
  else                            return nVar-1;
}

//------------------------------------------------------------------------------

// Function, which determines the indices for the prescribed data.
void BCInflowSubsonicSubfaceClass::SetIndicesPrescribedData(const int set)
{
  // Determine the set that has been prescribed.
  if(set == 0)
  {
    // Total conditions have been specified.
    // Set mTotalConditionsSpecified to true.
    mTotalConditionsSpecified = true;

    // Set the indices for the corresponding variables.
    mIndicesPrescribedData.push_back(GetIndexPrescribedData("TotalPressure"));
    mIndicesPrescribedData.push_back(GetIndexPrescribedData("TotalTemperature"));
    mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityDirX"));
    mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityDirY"));
    mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityDirZ"));
  }
  else
  {
    // Density and velocities are prescribed.
    // Set mTotalConditionsSpecified to false.
    mTotalConditionsSpecified = false;

    // Set the indices for the corresponding variables.
    mIndicesPrescribedData.push_back(GetIndexPrescribedData("Density"));
    mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityX"));
    mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityY"));
    mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityZ"));
  }
}

//------------------------------------------------------------------------------

// Function, which gets the type of boundary prescribed. 
int BCInflowSubsonicSubfaceClass::GetTypeBoundaryPrescribed(void) {return BC_INFLOW_SUBSONIC;}

//------------------------------------------------------------------------------

// Overloaded constructor of BCInflowSupersonicSubfaceClass.
BCInflowSupersonicSubfaceClass::BCInflowSupersonicSubfaceClass(std::istringstream &istr)
 : SubfaceBaseClass(istr){}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
BCInflowSupersonicSubfaceClass::~BCInflowSupersonicSubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which converts the prescribed data to the required form.
void BCInflowSupersonicSubfaceClass::ConvertPrescribedData(
                             const int                nIntegration,
                             const ENUM_FEM_VARIABLES FEMVariables,
                             std::vector<su2double *> &prescribedData)
{
  // Compute the inverse of the reference values
  const su2double rhoRefInv = one/rhoRef;
  const su2double uRefInv   = one/uRef;
  const su2double pRefInv   = one/pRef;

  // Compute 1/(gam-1), which is needed in the conversion
  // to the working variables.
  const su2double ovgm1 = one/(GamConstant - one);

  // Make a distinction between the working variables.
  switch( FEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      // The conservative variables are used as working variables.
      // Loop over the integration points.
#pragma omp simd
      for(int l=0; l<nIntegration; ++l)
      {
        // Compute the non-dimensional primitive variables.
        const su2double rho = prescribedData[0][l]*rhoRefInv;
        const su2double u   = prescribedData[1][l]*uRefInv;
        const su2double v   = prescribedData[2][l]*uRefInv;
        const su2double w   = prescribedData[3][l]*uRefInv;
        const su2double p   = prescribedData[4][l]*pRefInv;

        // Compute the conservative variables.
        prescribedData[0][l] = rho;
        prescribedData[1][l] = rho*u;
        prescribedData[2][l] = rho*v;
        prescribedData[3][l] = rho*w;
        prescribedData[4][l] = ovgm1*p + half*rho*(u*u + v*v + w*w);
      }

      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      // The entropy variables are used as working variables.
      // Loop over the integration points.
#pragma omp simd
      for(int l=0; l<nIntegration; ++l)
      {
        // Compute the non-dimensional primitive variables.
        const su2double rho = prescribedData[0][l]*rhoRefInv;
        const su2double u   = prescribedData[1][l]*uRefInv;
        const su2double v   = prescribedData[2][l]*uRefInv;
        const su2double w   = prescribedData[3][l]*uRefInv;
        const su2double p   = prescribedData[4][l]*pRefInv;

        // Compute the entropy variables.
        const su2double s     = LOG(p/POW(rho,GamConstant));
        const su2double pInv  = one/p;

        prescribedData[0][l] =  (GamConstant-s)*ovgm1 - half*pInv*rho*(u*u + v*v + w*w);
        prescribedData[1][l] =  rho*u*pInv;
        prescribedData[2][l] =  rho*v*pInv;
        prescribedData[3][l] =  rho*w*pInv;
        prescribedData[4][l] = -rho*pInv;
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

// Function, which indicates whether prescribed data is expected.
bool BCInflowSupersonicSubfaceClass::ExpectPrescribedData(void){return true;}

//------------------------------------------------------------------------------

// Function, which gets the names of the prescribed variables for the given set.
void BCInflowSupersonicSubfaceClass::GetNamesPrescibedVariables(const int                set,
                                                                std::vector<std::string> &varNames)
{
  // Supersonic inflow boundary. The primitive variables are prescribed.
  varNames.push_back("Density");
  varNames.push_back("VelocityX");
  varNames.push_back("VelocityY");
  varNames.push_back("VelocityZ");
  varNames.push_back("Pressure");
}

//------------------------------------------------------------------------------

// Function, which returns the number of data sets that can be prescribed.
int BCInflowSupersonicSubfaceClass::GetNSetsPrescribedVariables(void){return 1;}

//------------------------------------------------------------------------------

// Function, which makes available the number of prescribed variables.
int BCInflowSupersonicSubfaceClass::GetNVarPrescribed(void) {return nVar;}

//------------------------------------------------------------------------------

// Function, which determines the indices for the prescribed data.
void BCInflowSupersonicSubfaceClass::SetIndicesPrescribedData(const int set)
{
  // Set the indices for the corresponding variables.
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("Density"));
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityX"));
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityY"));
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityZ"));
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("Pressure"));
}

//------------------------------------------------------------------------------

// Function, which gets the type of boundary prescribed. 
int BCInflowSupersonicSubfaceClass::GetTypeBoundaryPrescribed(void) {return BC_INFLOW_SUPERSONIC;}

//------------------------------------------------------------------------------

// Overloaded constructor of BCOutflowSubsonicSubfaceClass.
BCOutflowSubsonicSubfaceClass::BCOutflowSubsonicSubfaceClass(std::istringstream &istr)
 : SubfaceBaseClass(istr){}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
BCOutflowSubsonicSubfaceClass::~BCOutflowSubsonicSubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which converts the prescribed data to the required form.
void BCOutflowSubsonicSubfaceClass::ConvertPrescribedData(
                             const int                nIntegration,
                             const ENUM_FEM_VARIABLES FEMVariables,
                             std::vector<su2double *> &prescribedData)
{
  // Convert the pressure from dimensional to non-dimensional.
  const su2double pRefInv = one/pRef;

#pragma omp simd
  for(int l=0; l<nIntegration; ++l)
    prescribedData[0][l] *= pRefInv;
}

//------------------------------------------------------------------------------

// Function, which indicates whether prescribed data is expected.
bool BCOutflowSubsonicSubfaceClass::ExpectPrescribedData(void){return true;}

//------------------------------------------------------------------------------

// Function, which gets the names of the prescribed variables for the given set.
void BCOutflowSubsonicSubfaceClass::GetNamesPrescibedVariables(const int                set,
                                                               std::vector<std::string> &varNames)
{
  // Subsonic outflow boundary. The pressure is prescribed.
  varNames.push_back("Pressure");
}

//------------------------------------------------------------------------------

// Function, which returns the number of data sets that can be prescribed.
int BCOutflowSubsonicSubfaceClass::GetNSetsPrescribedVariables(void){return 1;}

//------------------------------------------------------------------------------

// Function, which makes available the number of prescribed variables.
int BCOutflowSubsonicSubfaceClass::GetNVarPrescribed(void) {return 1;}

//------------------------------------------------------------------------------

// Function, which determines the indices for the prescribed data.
void BCOutflowSubsonicSubfaceClass::SetIndicesPrescribedData(const int set)
{
  // Set the index for the pressure.
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("Pressure"));
}

//------------------------------------------------------------------------------

// Function, which gets the type of boundary prescribed. 
int BCOutflowSubsonicSubfaceClass::GetTypeBoundaryPrescribed(void) {return BC_OUTFLOW_SUBSONIC;}

//------------------------------------------------------------------------------

// Overloaded constructor of BCOutflowSupersonicSubfaceClass.
BCOutflowSupersonicSubfaceClass::BCOutflowSupersonicSubfaceClass(std::istringstream &istr)
 : SubfaceBaseClass(istr){}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
BCOutflowSupersonicSubfaceClass::~BCOutflowSupersonicSubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which makes available the number of prescribed variables.
int BCOutflowSupersonicSubfaceClass::GetNVarPrescribed(void) {return 0;}

//------------------------------------------------------------------------------

// Function, which gets the type of boundary prescribed. 
int BCOutflowSupersonicSubfaceClass::GetTypeBoundaryPrescribed(void) {return BC_OUTFLOW_SUPERSONIC;}




//------------------------------------------------------------------------------
// Characteristic BC: Standard 
//------------------------------------------------------------------------------

// Overloaded constructor of BCStandardCharacteristicSubfaceClass.
BCStandardCharacteristicSubfaceClass::BCStandardCharacteristicSubfaceClass(std::istringstream &istr)
 : SubfaceBaseClass(istr){}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
BCStandardCharacteristicSubfaceClass::~BCStandardCharacteristicSubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which sets the value of mMachAverageBoundary.
void BCStandardCharacteristicSubfaceClass::SetMachAverageBoundary(const su2double Mavg) 
{
	mMachAverageBoundary = Mavg;
}

//------------------------------------------------------------------------------

// Function, which flags that this is a characteristic boundary.
bool BCStandardCharacteristicSubfaceClass::IsCharacteristicBC(void)
{
	return true;
}





//------------------------------------------------------------------------------
// Characteristic BC: Outlet 
//------------------------------------------------------------------------------

// Overloaded constructor of BCOutflowCharacteristicSubfaceClass.
BCOutflowCharacteristicSubfaceClass::BCOutflowCharacteristicSubfaceClass(std::istringstream &istr)
 : BCStandardCharacteristicSubfaceClass(istr){}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
BCOutflowCharacteristicSubfaceClass::~BCOutflowCharacteristicSubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which makes available the number of prescribed variables.
int BCOutflowCharacteristicSubfaceClass::GetNVarPrescribed(void) {return 1;}

//------------------------------------------------------------------------------

// Function, which gets the type of boundary prescribed. 
int BCOutflowCharacteristicSubfaceClass::GetTypeBoundaryPrescribed(void) {return BC_OUTFLOW_CHARACTERISTIC;}

//------------------------------------------------------------------------------

// Function, which converts the prescribed data to the required form.
void BCOutflowCharacteristicSubfaceClass::ConvertPrescribedData(
                             const int                nIntegration,
                             const ENUM_FEM_VARIABLES FEMVariables,
                             std::vector<su2double *> &prescribedData)
{
  // Convert the pressure from dimensional to non-dimensional.
  const su2double pRefInv = one/pRef;

#pragma omp simd
  for(int l=0; l<nIntegration; ++l)
    prescribedData[0][l] *= pRefInv;
}

//------------------------------------------------------------------------------

// Function, which gets the names of the prescribed variables for the given set.
void BCOutflowCharacteristicSubfaceClass::GetNamesPrescibedVariables(const int                set,
                                                               			 std::vector<std::string> &varNames)
{
  // Characteristic outflow boundary. The pressure is prescribed.
  varNames.push_back("Pressure");
}

//------------------------------------------------------------------------------

// Function, which returns the number of data sets that can be prescribed.
int BCOutflowCharacteristicSubfaceClass::GetNSetsPrescribedVariables(void){return 1;}

//------------------------------------------------------------------------------

// Function, which indicates whether prescribed data is expected.
bool BCOutflowCharacteristicSubfaceClass::ExpectPrescribedData(void){return true;}

//------------------------------------------------------------------------------

// Function, which determines the indices for the prescribed data.
void BCOutflowCharacteristicSubfaceClass::SetIndicesPrescribedData(const int set)
{
  // Set the index for the pressure.
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("Pressure"));
}

//------------------------------------------------------------------------------

// Function, which configures the tuning parameters required for a NSCBC.
void BCOutflowCharacteristicSubfaceClass::ConfigureParamNSCBC(
		                                         const InputParamClass 			*inputParam, 
													                   const StandardElementClass *standardHex,
													                   su2double            		 **metric,
													                   su2double                   factNorm)
{
	// If the working variables are not entropy based, terminate.
	if( inputParam->mFEMVariables != ENTROPY_VARIABLES )
		TerminateAll("BCOutflowCharacteristicSubfaceClass::ConfigureParamNSCBC", __FILE__, __LINE__,
				         "Only entropy variables are supported.");

	// Extract the max-face Lagrange derivative coefficient at the DOFs. 
	// Note, this is multiplied by the sign of the normal on this boundary
	// face, because it then gives either 
	// dL1D[0]: IMIN boundary or dL1D[end]: IMAX boundary. 
	// Note, dL1D[0] = -dL1D[end]. Also, since it is easier to simply use
	// dL1D[0] instead of dL1D[nNode1D*nNode1DPad+nNode1D], then 
	// dL1D[end] = -dL1D[0] is used for convenience only.
	mDerCoefficient = factNorm*(-standardHex->mDerLagrangeDOFs1D[0]);

	// Set characteristic length scale.
	mLengthScale = inputParam->mNSCBC_Outlet_len;

	// Set relaxation coefficients for outlet.
	mSigma       = inputParam->mNSCBC_Outlet_sigma;
	mBeta_l      = inputParam->mNSCBC_Outlet_beta_l;
	mBeta_t      = inputParam->mNSCBC_Outlet_beta_t;

	// Identify the incoming acoustic wave amplitude index.
	mIndexPhi    = (factNorm > zero) ? 0 : 4;   

	// Identify the relevant indices required for this characteristic boundary.
	IdentifyBoundaryIndices(metric);
}




//------------------------------------------------------------------------------
// Characteristic BC: Inlet (generic) 
//------------------------------------------------------------------------------

// Overloaded constructor of BCInflowCharacteristicSubfaceClass.
BCInflowCharacteristicSubfaceClass::BCInflowCharacteristicSubfaceClass(std::istringstream &istr)
 : BCStandardCharacteristicSubfaceClass(istr){}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
BCInflowCharacteristicSubfaceClass::~BCInflowCharacteristicSubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which configures the tuning parameters required for a NSCBC.
void BCInflowCharacteristicSubfaceClass::ConfigureParamNSCBC(
		                                         const InputParamClass 			*inputParam, 
													                   const StandardElementClass *standardHex,
													                   su2double            		 **metric,
													                   su2double                   factNorm)
{
	// If the working variables are not entropy based, terminate.
	if( inputParam->mFEMVariables != ENTROPY_VARIABLES )
		TerminateAll("BCInflowCharacteristicSubfaceClass::ConfigureParamNSCBC", __FILE__, __LINE__,
				         "Only entropy variables are supported.");

	// Extract the max-face Lagrange derivative coefficient at the DOFs. 
	// Note, this is multiplied by the sign of the normal on this boundary
	// face, because it then gives either 
	// dL1D[0]: IMIN boundary or dL1D[end]: IMAX boundary. 
	// Note, dL1D[0] = -dL1D[end]. Also, since it is easier to simply use
	// dL1D[0] instead of dL1D[nNode1D*nNode1DPad+nNode1D], then 
	// dL1D[end] = -dL1D[0] is used for convenience only.
	mDerCoefficient = factNorm*(-standardHex->mDerLagrangeDOFs1D[0]);

	// Set characteristic length scale.
	mLengthScale    = inputParam->mNSCBC_Inlet_len;

	// Set relaxation coefficients for inlet.
	mSigma          = inputParam->mNSCBC_Inlet_sigma;

	// Identify the incoming acoustic wave amplitude index.
	mIndexPhi       = (factNorm > zero) ? 0 : 4;   
	// Identify the outgoing acoustic wave amplitude index. 
	mIndexPsi       = (mIndexPhi == 0 ) ? 4 : 0;

	// Identify the relevant indices required for this characteristic boundary.
	IdentifyBoundaryIndices(metric);
}




//------------------------------------------------------------------------------
// Characteristic BC: Inlet (static) 
//------------------------------------------------------------------------------

// Overloaded constructor of BCInflowStaticCharacteristicSubfaceClass.
BCInflowStaticCharacteristicSubfaceClass::BCInflowStaticCharacteristicSubfaceClass(std::istringstream &istr)
 : BCInflowCharacteristicSubfaceClass(istr){}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
BCInflowStaticCharacteristicSubfaceClass::~BCInflowStaticCharacteristicSubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which makes available the number of prescribed variables.
int BCInflowStaticCharacteristicSubfaceClass::GetNVarPrescribed(void) 
{
	return nVar-1;
}

//------------------------------------------------------------------------------

// Function, which gets the type of boundary prescribed. 
int BCInflowStaticCharacteristicSubfaceClass::GetTypeBoundaryPrescribed(void) 
{
	return BC_INFLOW_STATIC_CHARACTERISTIC;
}

//------------------------------------------------------------------------------

// Function, which converts the prescribed data to the required form.
void BCInflowStaticCharacteristicSubfaceClass::ConvertPrescribedData(
                             const int                nIntegration,
                             const ENUM_FEM_VARIABLES FEMVariables,
                             std::vector<su2double *> &prescribedData)
{
  // Compute the inverse of the reference conditions.
  const su2double rhoRefInv = one/rhoRef;
  const su2double uRefInv   = one/uRef;
  //const su2double pRefInv   = one/pRef;

  // Density and velocities are specified. Loop over the integration
  // points and non-dimensionalize them.
#pragma omp simd
  for(int l=0; l<nIntegration; ++l)
  {
    prescribedData[0][l] *= rhoRefInv;
    prescribedData[1][l] *= uRefInv;
    prescribedData[2][l] *= uRefInv;
    prescribedData[3][l] *= uRefInv;
  }
}

//------------------------------------------------------------------------------

// Function, which gets the names of the prescribed variables for the given set.
void BCInflowStaticCharacteristicSubfaceClass::GetNamesPrescibedVariables(const int                set,
                                                                          std::vector<std::string> &varNames)
{
	// The static and total inflows, along with their respective data sets are separated in different classes.
	// Thus, make sure there exist only one set.
	if( set !=  0 )
		TerminateAll("BCInflowStaticCharacteristicSubfaceClass::GetNamesPrescibedVariables",
                  __FILE__, __LINE__, "This should not happen");

  // This is the set of static variables. Density and velocities are prescribed.
  varNames.push_back("Density");
  varNames.push_back("VelocityX");
  varNames.push_back("VelocityY");
  varNames.push_back("VelocityZ");
}

//------------------------------------------------------------------------------

// Function, which returns the number of data sets that can be prescribed.
int BCInflowStaticCharacteristicSubfaceClass::GetNSetsPrescribedVariables(void){return 1;}

//------------------------------------------------------------------------------

// Function, which indicates whether prescribed data is expected.
bool BCInflowStaticCharacteristicSubfaceClass::ExpectPrescribedData(void){return true;}

//------------------------------------------------------------------------------

// Function, which determines the indices for the prescribed data.
void BCInflowStaticCharacteristicSubfaceClass::SetIndicesPrescribedData(const int set)
{
	// The static and total inflows, along with their respective data sets are separated in different classes.
	// Thus, make sure there exist only one set.
	if( set !=  0 )
		TerminateAll("BCInflowStaticCharacteristicSubfaceClass::SetIndicesPrescribedData",
                  __FILE__, __LINE__, "This should not happen");

  // Set the indices for the corresponding variables.
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("Density"));
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityX"));
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityY"));
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityZ"));
}




//------------------------------------------------------------------------------
// Characteristic BC: Inlet (total) 
//------------------------------------------------------------------------------

// Overloaded constructor of BCInflowTotalCharacteristicSubfaceClass.
BCInflowTotalCharacteristicSubfaceClass::BCInflowTotalCharacteristicSubfaceClass(std::istringstream &istr)
 : BCInflowCharacteristicSubfaceClass(istr){}

//------------------------------------------------------------------------------

// Destructor. Nothing to be done.
BCInflowTotalCharacteristicSubfaceClass::~BCInflowTotalCharacteristicSubfaceClass(){}

//------------------------------------------------------------------------------

// Function, which makes available the number of prescribed variables.
int BCInflowTotalCharacteristicSubfaceClass::GetNVarPrescribed(void) 
{
	return nVar;
}

//------------------------------------------------------------------------------

// Function, which gets the type of boundary prescribed. 
int BCInflowTotalCharacteristicSubfaceClass::GetTypeBoundaryPrescribed(void) 
{
	return BC_INFLOW_TOTAL_CHARACTERISTIC;
}

//------------------------------------------------------------------------------

// Function, which converts the prescribed data to the required form.
void BCInflowTotalCharacteristicSubfaceClass::ConvertPrescribedData(
                             const int                nIntegration,
                             const ENUM_FEM_VARIABLES FEMVariables,
                             std::vector<su2double *> &prescribedData)
{
  // Compute the inverse of the reference values.
  const su2double pRefInv = one/pRef;

	// Compute the reference temperature and its inverse.
	// Note, division by RGas is skipped, since this produces the same scaling 
	// used when constructing the temperature in the NSCBC ComputeBoundary.
	const su2double TRef    = pRef/rhoRef;
	const su2double TRefInv = one/TRef;

  // Loop over the integration points.
#pragma omp simd
  for(int l=0; l<nIntegration; ++l)
  {
    // Compute the non-dimensional total pressure and total temperature.
    prescribedData[0][l] *= pRefInv;
    prescribedData[1][l] *= TRefInv;

    // Make sure the unit vector really is a unit vector.
    const su2double len2   = prescribedData[2][l]*prescribedData[2][l]
                           + prescribedData[3][l]*prescribedData[3][l]
                           + prescribedData[4][l]*prescribedData[4][l];
    const su2double lenInv = one/SQRT(len2);

    prescribedData[2][l] *= lenInv;
    prescribedData[3][l] *= lenInv;
    prescribedData[4][l] *= lenInv;
  }
}

//------------------------------------------------------------------------------

// Function, which gets the names of the prescribed variables for the given set.
void BCInflowTotalCharacteristicSubfaceClass::GetNamesPrescibedVariables(const int                set,
                                                                         std::vector<std::string> &varNames)
{
	// The static and total inflows, along with their respective data sets are separated in different classes.
	// Thus, make sure there exist only one set.
	if( set !=  0 )
		TerminateAll("BCInflowTotalCharacteristicSubfaceClass::GetNamesPrescibedVariables",
                  __FILE__, __LINE__, "This should not happen");

	// This is the set of total variables. Total data and velocity direction are prescribed.
  varNames.push_back("TotalPressure");
  varNames.push_back("TotalTemperature");
  varNames.push_back("VelocityDirX");
  varNames.push_back("VelocityDirY");
  varNames.push_back("VelocityDirZ");     
}

//------------------------------------------------------------------------------

// Function, which returns the number of data sets that can be prescribed.
int BCInflowTotalCharacteristicSubfaceClass::GetNSetsPrescribedVariables(void){return 1;}

//------------------------------------------------------------------------------

// Function, which indicates whether prescribed data is expected.
bool BCInflowTotalCharacteristicSubfaceClass::ExpectPrescribedData(void){return true;}

//------------------------------------------------------------------------------

// Function, which determines the indices for the prescribed data.
void BCInflowTotalCharacteristicSubfaceClass::SetIndicesPrescribedData(const int set)
{
	// The static and total inflows, along with their respective data sets are separated in different classes.
	// Thus, make sure there exist only one set.
	if( set !=  0 )
		TerminateAll("BCInflowTotalCharacteristicSubfaceClass::SetIndicesPrescribedData",
                  __FILE__, __LINE__, "This should not happen");

  // Set the indices for the corresponding variables.
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("TotalPressure"));
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("TotalTemperature"));
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityDirX"));
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityDirY"));
  mIndicesPrescribedData.push_back(GetIndexPrescribedData("VelocityDirZ"));
}




