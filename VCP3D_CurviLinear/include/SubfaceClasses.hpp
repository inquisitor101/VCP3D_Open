//------------------------------------------------------------------------------
// SubfaceClasses.hpp: Definition of the classes that contain the boundary
//                     condition information on the subfaces. This include file
//                     should not be included by any source file. These files
//                     should just include VCP3D_CurviLinear.hpp.
//------------------------------------------------------------------------------

#pragma once

// Forward definition of InputParamClass and StandardElementClass.
class InputParamClass;
class StandardElementClass;

// Definition of the class with prescribed data.
class PrescribedDataClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

  // Constructor. Nothing to be done.
  PrescribedDataClass();

  // Destructor. Nothing to be done.
  ~PrescribedDataClass();

  //------------------------------------------
  // Public member variables.
  //------------------------------------------

  // Name of the prescribed variable.
  std::string mNameVariable;

  // Vector with prescribed data.
  std::vector<su2double> mData;
};

//------------------------------------------------------------------------------

// Definition of the subface base class.
class SubfaceBaseClass
{
public:
  //------------------------------------------
  // Destructor.
  //------------------------------------------
  
  // Destructor. Nothing to be done.
  virtual ~SubfaceBaseClass();

  //------------------------------------------
  // Public member functions.
  //------------------------------------------

  // Function, which checks if the expected data is prescribed.
  void CheckPrescribedData(const int faceID,
                           const int subfaceID);

  // Function to make available the prescribed data as a const reference.
  const std::vector<PrescribedDataClass> &GetPrescribedData(void);

  // Function, which determines the prescribed data in the
  // given number of integration points.
  void PrescribedDataIntegrationPoints(const int                nIntegration,
                                       const su2double          *x,
                                       const su2double          *y,
                                       const su2double          *z,
                                       const su2double          *nonDimPrimVarFreeStream,
                                       const ENUM_FEM_VARIABLES FEMVariables,
                                       std::vector<su2double *> &prescribedData);

  // Function, which reads the prescribed data for this subface.
  void ReadBoundaryData(const int          faceID,
                        const int          subfaceID,
                        const int          nPoints,
                        std::istringstream &istrVarNames,
                        std::ifstream      &boundaryDataFile);

  //----------------------------------------------------------
  // Public virtual member functions. These may/must
  // be overloaded by the derived classes.
  //----------------------------------------------------------

	// Virtual function, which sets the value of mMachAverageBoundary.
	// It may be overwritten by the derived class.
	virtual void SetMachAverageBoundary(const su2double Mavg);

	// Virtual function, which configures the tuning parameters 
	// required for a NSCBC. It may be overwritten by the derived class.
	virtual void ConfigureParamNSCBC(const InputParamClass 			*inputParam, 
																	 const StandardElementClass *standardHex,
																	 su2double            		 **metric,
																	 su2double                   factNorm);

	// Virtual function, which computes the weighted Mach number on this
	// local element boundary that is used to compute the average.
	// The average Mach number on this element surface is also computed.
	// It may be overwritten by the derived class.
	virtual su2double WeightedMachElement(const InputParamClass *inputParam,
			                                  const int              nInt,
										 					          su2double            **sol,
										 					          const su2double        factNorm,
										 					          su2double            **metric,
										 					          su2double             *intWeights);

  // Virtual function, which indicates whether or not a boundary condition
  // is physical. It may be overwritten by the derived class.
  virtual bool BCIsPhysical(void);

  // Virtual function, which indicates whether or not a boundary condition
  // is a wall boundary condition. It should be overwritten by the
  // wall boundary condition classes.
  virtual bool BCIsWall(void);

  // Function, which computes the boundary state (the right state) from the
  // given left state and the prescribed boundary data. It must be overwritten
  // by the derived class.
  virtual void ComputeBoundaryState(const InputParamClass      *inputParam,
		                                const StandardElementClass *standardHex,
                                    const int                   nInt,
                                    su2double                 **solL,
                                    su2double                 **dSolDxL,
                                    su2double                 **dSolDyL,
                                    su2double                 **dSolDzL,
                                    const su2double             factNorm,
                                    su2double                 **metricL,
                                    su2double                 **prescribedData,
                                    su2double                 **solR,
                                    su2double                 **dSolDxR,
                                    su2double                 **dSolDyR,
                                    su2double                 **dSolDzR,
                                    bool                       &heatFluxPrescribed,
                                    su2double                 *&prescribedWallData,
                                    su2double                  &wallPerm);

  // Virtual function, which makes available the beginning i-index
  // of the element donor range. It may be overwritten by the derived class.
  virtual int GetElemDonorIBeg(void);

  // Virtual function, which makes available the beginning j-index
  // of the element donor range. It may be overwritten by the derived class.
  virtual int GetElemDonorJBeg(void);

  // Virtual function, which makes available the beginning k-index
  // of the element donor range. It may be overwritten by the derived class.
  virtual int GetElemDonorKBeg(void);
  
  // Virtual function, which makes available the number of prescribed
  // variables. It must be overwritten by the physical boundary conditions.
  virtual int GetNVarPrescribed(void);

  // Virtual function, which makes the possible periodic transformation available.
  // It may be overwritten by the derived class.
  virtual void GetPeriodicTransformation(su2double *trans);

  // Virtual function, which indicates whether or not this subface is 1 to 1 matching.
  // It may be overwritten by the derived class.
  virtual bool Is1To1Matching(void);

  // Virtual function, which indicates whether or not a wall model is used
  // for the boundary. It may be overwritten by the derived class.
  virtual bool UseWallModel(void);

	// Virtual function, which makes available the type of boundary 
	// condition prescribed. It is overwritten by the derived class.
	virtual int GetTypeBoundaryPrescribed(void);

  //------------------------------------------
  // Public member variables.
  //------------------------------------------

  // Begin i-index of the element subface range.
  int mElemIBeg;

  // End i-index of the element subface range.
  int mElemIEnd;

  // Begin j-index of the element subface range.
  int mElemJBeg;

  // End J-index of the element subface range.
  int mElemJEnd;

  // Begin k-index of the element subface range.
  int mElemKBeg;

  // End k-index of the element subface range.
  int mElemKEnd;

	// Boundary face ID.
	int mBoundaryID;
	// Subface index ID.
	int mSubfacesID;

protected:

  //------------------------------------------
  // Constructor.
  //------------------------------------------

  // Overloaded constructor. Read the subface range.
  SubfaceBaseClass(std::istringstream &istr);

  //------------------------------------------
  // Protected member functions.
  //------------------------------------------

  // Function, which applies the inviscid wall boundary conditions.
  void ApplyInviscidWallBC(const InputParamClass *inputParam,
                           const int             nInt,
                           su2double             **solL,
                           su2double             **metricL,
                           su2double             **solR);

  // Function, which gets the index for in mPrescribedData of the given variable.
  int GetIndexPrescribedData(const char *varName);

  //------------------------------------------
  // Protected member variables.
  //------------------------------------------

  // Indices in mPrescribedData that the actual boundary
  // conditions routines expect.
  std::vector<int> mIndicesPrescribedData;

private:
  //----------------------------------------------------------
  // Private virtual member functions. These may/must
  // be overloaded by the derived classes.
  //----------------------------------------------------------

  // Function, which converts the prescribed data to the required form.
  // It should be overwritten by the derived classes, if this is needed.
  virtual void ConvertPrescribedData(const int                nIntegration,
                                     const ENUM_FEM_VARIABLES FEMVariables,
                                     std::vector<su2double *> &prescribedData);

  // Function, which gets the names of the prescribed variables for the given set.
  // It should be overwritten by the derived classes.
  virtual void GetNamesPrescibedVariables(const int                set,
                                          std::vector<std::string> &varNames);

  // Function, which returns the number of data sets that can be prescribed.
  // It may be overwritten by the derived classes.
  virtual int GetNSetsPrescribedVariables(void);

  // Function, which indicates whether prescribed data is expected.
  // It may be overwritten by the derived classes.
  virtual bool ExpectPrescribedData(void);

  // Function, which determines the prescribed data from the free stream
  // data. It should be overwritten by the derived classes.
  virtual void PrescribedDataFromFreeStream(const int                nIntegration,
                                            const su2double          *nonDimPrimVarFreeStream,
                                            const ENUM_FEM_VARIABLES FEMVariables,
                                            std::vector<su2double *> &prescribedData);

  // Function, which determines the indices for the prescribed data.
  // It may be overwritten by the derived class.
  virtual void SetIndicesPrescribedData(const int set);

  //------------------------------------------
  // Private member variables.
  //------------------------------------------

  // Vector with the prescribed data.
  std::vector<PrescribedDataClass> mPrescribedData;

  //------------------------------------------
  // Private member functions.
  //------------------------------------------

  // Disabled default constructor.
  SubfaceBaseClass();
};

//------------------------------------------------------------------------------

// Definition of the internal 1 to 1 matching subface class.
class Internal1to1SubfaceClass : public SubfaceBaseClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

  // Overloaded constructor. Read the subface data.
  Internal1to1SubfaceClass(std::istringstream &istr);

  // Destructor. Nothing to be done.
  virtual ~Internal1to1SubfaceClass();

  //----------------------------------------------------------
  // Public member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which indicates whether or not a boundary condition is physical.
  bool BCIsPhysical(void);

  // Function, which makes available the beginning i-index
  // of the element donor range.
  int GetElemDonorIBeg(void);

  // Function, which makes available the beginning j-index
  // of the element donor range.
  int GetElemDonorJBeg(void);
  
  // Function, which makes available the beginning k-index
  // of the element donor range.
  int GetElemDonorKBeg(void);
 
  // Function, which indicates whether or not this subface is 1 to 1 matching.
  // It is still virtual, because it is overwritten by Periodic1to1TransSubfaceClass.
  virtual bool Is1To1Matching(void);

	// Virtual function, which makes available the type of boundary 
	// condition prescribed. It is still virtual, because it is 
	// overwritten by Periodic1to1TransSubfaceClass.
	virtual int GetTypeBoundaryPrescribed(void);

  //------------------------------------------
  // Public member variables.
  //------------------------------------------

  // Begin i-index of the element donor subface range.
  int mElemDonorIBeg;
  
  // End i-index of the element donor subface range.
  int mElemDonorIEnd;
  
  // Begin j-index of the element donor subface range.
  int mElemDonorJBeg;
  
  // End J-index of the element donor subface range.
  int mElemDonorJEnd;
  
  // Begin k-index of the element donor subface range.
  int mElemDonorKBeg;

  // End k-index of the element donor subface range.
  int mElemDonorKEnd;
  
  // First index of the transformation matrix between this subface and the donor subface.
  int mL1;

  // Second index of the transformation matrix between this subface and the donor subface.
  int mL2;

  // Third index of the transformation matrix between this subface and the donor subface.
  int mL3;

private:
  //------------------------------------------
  // Private member functions.
  //------------------------------------------

  // Disabled default constructor.
  Internal1to1SubfaceClass();
};

//------------------------------------------------------------------------------

// Definition of the translational periodic 1 to 1 matching subface class.
class Periodic1to1TransSubfaceClass : public Internal1to1SubfaceClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

  // Overloaded constructor. Read the subface data.
  Periodic1to1TransSubfaceClass(std::istringstream &istr);

  // Destructor. Nothing to be done.
  ~Periodic1to1TransSubfaceClass();

  //----------------------------------------------------------
  // Public member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which makes the periodic transformation available.
  void GetPeriodicTransformation(su2double *trans);

  // Function, which indicates whether or not this subface is 1 to 1 matching
  // and not periodic.
  bool Is1To1Matching(void);

	// Function, which makes available the type of boundary 
	// condition prescribed. 
	int GetTypeBoundaryPrescribed(void);

  //------------------------------------------
  // Public member variables.
  //------------------------------------------

  // Translation vector from the current subface to the donor subface.
  su2double mPeriodicTranslation[3];

private:
  //------------------------------------------
  // Private member functions.
  //------------------------------------------

  // Disabled default constructor.
  Periodic1to1TransSubfaceClass();
};

//------------------------------------------------------------------------------

// Definition of the farfield boundary subface class.
class BCFarfieldSubfaceClass : public SubfaceBaseClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

  // Overloaded constructor. Read the subface data.
  BCFarfieldSubfaceClass(std::istringstream &istr);

  // Destructor. Nothing to be done.
  ~BCFarfieldSubfaceClass();

  //----------------------------------------------------------
  // Public member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which computes the boundary state (the right state) from the
  // given left state and the prescribed boundary data.
  void ComputeBoundaryState(const InputParamClass      *inputParam,
		                        const StandardElementClass *standardHex,
                            const int                   nInt,
                            su2double                 **solL,
                            su2double                 **dSolDxL,
                            su2double                 **dSolDyL,
                            su2double                 **dSolDzL,
                            const su2double             factNorm,
                            su2double                 **metricL,
                            su2double                 **prescribedData,
                            su2double                 **solR,
                            su2double                 **dSolDxR,
                            su2double                 **dSolDyR,
                            su2double                 **dSolDzR,
                            bool                       &heatFluxPrescribed,
                            su2double                 *&prescribedWallData,
                            su2double                  &wallPerm);

  // Function, which makes available the number of prescribed variables.
  int GetNVarPrescribed(void);

	// Function, which makes available the type of boundary 
	// condition prescribed. 
	int GetTypeBoundaryPrescribed(void);

private:

  //----------------------------------------------------------
  // Private member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which determines the prescribed data from the free stream.
  void PrescribedDataFromFreeStream(const int                nIntegration,
                                    const su2double          *nonDimPrimVarFreeStream,
                                    const ENUM_FEM_VARIABLES FEMVariables,
                                    std::vector<su2double *> &prescribedData);

  //------------------------------------------
  // Private member functions.
  //------------------------------------------

  // Disabled default constructor.
  BCFarfieldSubfaceClass();
};

//------------------------------------------------------------------------------

// Definition of the isothermal wall boundary subface class.
class BCIsothermalWallSubfaceClass : public SubfaceBaseClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

  // Overloaded constructor. Read the subface data.
  BCIsothermalWallSubfaceClass(std::istringstream &istr);

  // Overloaded constructor. Read the subface data and set the wall model information.
  BCIsothermalWallSubfaceClass(std::istringstream &istr,
                               ENUM_WALL_MODEL    wallModelType);

  // Destructor. Nothing to be done.
  ~BCIsothermalWallSubfaceClass();

  //----------------------------------------------------------
  // Public member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which indicates whether or not a boundary condition
  // is a wall boundary condition. True for this BC.
  bool BCIsWall(void);

  // Function, which computes the boundary state (the right state) from the
  // given left state and the prescribed boundary data.
  void ComputeBoundaryState(const InputParamClass      *inputParam,
		                        const StandardElementClass *standardHex,
                            const int                   nInt,
                            su2double                 **solL,
                            su2double                 **dSolDxL,
                            su2double                 **dSolDyL,
                            su2double                 **dSolDzL,
                            const su2double             factNorm,
                            su2double                 **metricL,
                            su2double                 **prescribedData,
                            su2double                 **solR,
                            su2double                 **dSolDxR,
                            su2double                 **dSolDyR,
                            su2double                 **dSolDzR,
                            bool                       &heatFluxPrescribed,
                            su2double                 *&prescribedWallData,
                            su2double                  &wallPerm);

  // Function, which makes available the number of prescribed variables.
  int GetNVarPrescribed(void);

  // Function, which indicates whether or not a wall model is used
  // for the boundary.
  bool UseWallModel(void);

	// Function, which makes available the type of boundary 
	// condition prescribed. 
	int GetTypeBoundaryPrescribed(void);

private:
  //----------------------------------------------------------
  // Private member variables.
  //----------------------------------------------------------

  // Whether or not wall modeling must be used.
  bool mUseWallModel;

  //----------------------------------------------------------
  // Private member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which gets the names of the prescribed variables for the given set.
  void GetNamesPrescibedVariables(const int                set,
                                  std::vector<std::string> &varNames);

  // Function, which returns the number of data sets that can be prescribed.
  int GetNSetsPrescribedVariables(void);
  
  // Function, which indicates whether prescribed data is expected.
  bool ExpectPrescribedData(void);

  // Function, which determines the indices for the prescribed data.
  void SetIndicesPrescribedData(const int set);

  //------------------------------------------
  // Private member functions.
  //------------------------------------------

  // Disabled default constructor.
  BCIsothermalWallSubfaceClass();
};

//------------------------------------------------------------------------------

// Definition of the heat flux wall boundary subface class.
class BCHeatFluxWallSubfaceClass : public SubfaceBaseClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

  // Overloaded constructor. Read the subface data.
  BCHeatFluxWallSubfaceClass(std::istringstream &istr);

  // Overloaded constructor. Read the subface data and set the wall model information.
  BCHeatFluxWallSubfaceClass(std::istringstream &istr,
                             ENUM_WALL_MODEL    wallModelType);

  // Destructor. Nothing to be done.
  ~BCHeatFluxWallSubfaceClass();

  //----------------------------------------------------------
  // Public member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which indicates whether or not a boundary condition
  // is a wall boundary condition. True for this BC.
  bool BCIsWall(void);

  // Function, which computes the boundary state (the right state) from the
  // given left state and the prescribed boundary data.
  void ComputeBoundaryState(const InputParamClass      *inputParam,
		                        const StandardElementClass *standardHex,
                            const int                   nInt,
                            su2double                 **solL,
                            su2double                 **dSolDxL,
                            su2double                 **dSolDyL,
                            su2double                 **dSolDzL,
                            const su2double             factNorm,
                            su2double                 **metricL,
                            su2double                 **prescribedData,
                            su2double                 **solR,
                            su2double                 **dSolDxR,
                            su2double                 **dSolDyR,
                            su2double                 **dSolDzR,
                            bool                       &heatFluxPrescribed,
                            su2double                 *&prescribedWallData,
                            su2double                  &wallPerm);

  // Function, which makes available the number of prescribed variables.
  int GetNVarPrescribed(void);

  // Function, which indicates whether or not a wall model is used
  // for the boundary.
  bool UseWallModel(void);
 
	// Function, which makes available the type of boundary 
	// condition prescribed. 
	int GetTypeBoundaryPrescribed(void);
 
private:
  //----------------------------------------------------------
  // Private member variables.
  //----------------------------------------------------------

  // Whether or not wall modeling must be used.
  bool mUseWallModel;

  //----------------------------------------------------------
  // Private member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which converts the prescribed data to the required form,
  // i.e. the heat flux is made non-dimensioal.
  void ConvertPrescribedData(const int                nIntegration,
                             const ENUM_FEM_VARIABLES FEMVariables,
                             std::vector<su2double *> &prescribedData);

  // Function, which gets the names of the prescribed variables for the given set.
  void GetNamesPrescibedVariables(const int                set,
                                  std::vector<std::string> &varNames);

  // Function, which returns the number of data sets that can be prescribed.
  int GetNSetsPrescribedVariables(void);

  // Function, which indicates whether prescribed data is expected.
  bool ExpectPrescribedData(void);

  // Function, which determines the indices for the prescribed data.
  void SetIndicesPrescribedData(const int set);

  //------------------------------------------
  // Private member functions.
  //------------------------------------------

  // Disabled default constructor.
  BCHeatFluxWallSubfaceClass();
};

//------------------------------------------------------------------------------

// Definition of the inviscid wall boundary subface class.
class BCInviscidWallSubfaceClass : public SubfaceBaseClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

  // Overloaded constructor. Read the subface data.
  BCInviscidWallSubfaceClass(std::istringstream &istr);

  // Destructor. Nothing to be done.
  ~BCInviscidWallSubfaceClass();

  //----------------------------------------------------------
  // Public member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which indicates whether or not a boundary condition
  // is a wall boundary condition. True for this BC.
  bool BCIsWall(void);

  // Function, which computes the boundary state (the right state) from the
  // given left state and the prescribed boundary data.
  void ComputeBoundaryState(const InputParamClass      *inputParam,
		                        const StandardElementClass *standardHex,
                            const int                   nInt,
                            su2double                 **solL,
                            su2double                 **dSolDxL,
                            su2double                 **dSolDyL,
                            su2double                 **dSolDzL,
                            const su2double             factNorm,
                            su2double                 **metricL,
                            su2double                 **prescribedData,
                            su2double                 **solR,
                            su2double                 **dSolDxR,
                            su2double                 **dSolDyR,
                            su2double                 **dSolDzR,
                            bool                       &heatFluxPrescribed,
                            su2double                 *&prescribedWallData,
                            su2double                  &wallPerm);

  // Function, which makes available the number of prescribed variables.
  int GetNVarPrescribed(void);
 
	// Function, which makes available the type of boundary 
	// condition prescribed. 
	int GetTypeBoundaryPrescribed(void);
 
private:
  //------------------------------------------
  // Private member functions.
  //------------------------------------------

  // Disabled default constructor.
  BCInviscidWallSubfaceClass();
};

//------------------------------------------------------------------------------

// Definition of the symmetry plane boundary subface class.
class BCSymmetrySubfaceClass : public SubfaceBaseClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

  // Overloaded constructor. Read the subface data.
  BCSymmetrySubfaceClass(std::istringstream &istr);

  // Destructor. Nothing to be done.
  ~BCSymmetrySubfaceClass();

  //----------------------------------------------------------
  // Public member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which computes the boundary state (the right state) from the
  // given left state and the prescribed boundary data.
  void ComputeBoundaryState(const InputParamClass      *inputParam,
		                        const StandardElementClass *standardHex,
                            const int                   nInt,
                            su2double                 **solL,
                            su2double                 **dSolDxL,
                            su2double                 **dSolDyL,
                            su2double                 **dSolDzL,
                            const su2double             factNorm,
                            su2double                 **metricL,
                            su2double                 **prescribedData,
                            su2double                 **solR,
                            su2double                 **dSolDxR,
                            su2double                 **dSolDyR,
                            su2double                 **dSolDzR,
                            bool                       &heatFluxPrescribed,
                            su2double                 *&prescribedWallData,
                            su2double                  &wallPerm);

  // Function, which makes available the number of prescribed variables.
  int GetNVarPrescribed(void);
 
	// Function, which makes available the type of boundary 
	// condition prescribed. 
	int GetTypeBoundaryPrescribed(void);
 
private:
  //------------------------------------------
  // Private member functions.
  //------------------------------------------

  // Disabled default constructor.
  BCSymmetrySubfaceClass();
};

//------------------------------------------------------------------------------

// Definition of the subsonic inflow boundary subface class.
class BCInflowSubsonicSubfaceClass : public SubfaceBaseClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

  // Overloaded constructor. Read the subface data.
  BCInflowSubsonicSubfaceClass(std::istringstream &istr);

  // Destructor. Nothing to be done.
  ~BCInflowSubsonicSubfaceClass();

  //----------------------------------------------------------
  // Public member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which computes the boundary state (the right state) from the
  // given left state and the prescribed boundary data.
  void ComputeBoundaryState(const InputParamClass      *inputParam,
		                        const StandardElementClass *standardHex,
                            const int                   nInt,
                            su2double                 **solL,
                            su2double                 **dSolDxL,
                            su2double                 **dSolDyL,
                            su2double                 **dSolDzL,
                            const su2double             factNorm,
                            su2double                 **metricL,
                            su2double                 **prescribedData,
                            su2double                 **solR,
                            su2double                 **dSolDxR,
                            su2double                 **dSolDyR,
                            su2double                 **dSolDzR,
                            bool                       &heatFluxPrescribed,
                            su2double                 *&prescribedWallData,
                            su2double                  &wallPerm);

  // Function, which makes available the number of prescribed variables.
  int GetNVarPrescribed(void);
  
	// Function, which makes available the type of boundary 
	// condition prescribed. 
	int GetTypeBoundaryPrescribed(void);
  
private:
  //----------------------------------------------------------
  // Private member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which converts the prescribed data to the required form.
  void ConvertPrescribedData(const int                nIntegration,
                             const ENUM_FEM_VARIABLES FEMVariables,
                             std::vector<su2double *> &prescribedData);

  // Function, which gets the names of the prescribed variables for the given set.
  void GetNamesPrescibedVariables(const int                set,
                                  std::vector<std::string> &varNames);

  // Function, which returns the number of data sets that can be prescribed.
  int GetNSetsPrescribedVariables(void);

  // Function, which indicates whether prescribed data is expected.
  bool ExpectPrescribedData(void);

  // Function, which determines the indices for the prescribed data.
  void SetIndicesPrescribedData(const int set);

  //------------------------------------------
  // Private member variables.
  //------------------------------------------

  // Boolean, which indicates whether or not total
  // conditions have been specified.
  bool mTotalConditionsSpecified;

  //------------------------------------------
  // Private member functions.
  //------------------------------------------

  // Disabled default constructor.
  BCInflowSubsonicSubfaceClass();
};

//------------------------------------------------------------------------------

// Definition of the supersonic inflow boundary subface class.
class BCInflowSupersonicSubfaceClass : public SubfaceBaseClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

  // Overloaded constructor. Read the subface data.
  BCInflowSupersonicSubfaceClass(std::istringstream &istr);

  // Destructor. Nothing to be done.
  ~BCInflowSupersonicSubfaceClass();

  //----------------------------------------------------------
  // Public member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which computes the boundary state (the right state) from the
  // given left state and the prescribed boundary data.
  void ComputeBoundaryState(const InputParamClass      *inputParam,
		                        const StandardElementClass *standardHex,
                            const int                   nInt,
                            su2double                 **solL,
                            su2double                 **dSolDxL,
                            su2double                 **dSolDyL,
                            su2double                 **dSolDzL,
                            const su2double             factNorm,
                            su2double                 **metricL,
                            su2double                 **prescribedData,
                            su2double                 **solR,
                            su2double                 **dSolDxR,
                            su2double                 **dSolDyR,
                            su2double                 **dSolDzR,
                            bool                       &heatFluxPrescribed,
                            su2double                 *&prescribedWallData,
                            su2double                  &wallPerm);

  // Function, which makes available the number of prescribed variables.
  int GetNVarPrescribed(void);
  
	// Function, which makes available the type of boundary 
	// condition prescribed. 
	int GetTypeBoundaryPrescribed(void);
  
private:
  //----------------------------------------------------------
  // Private member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which converts the prescribed data to the required form,
  // i.e. the heat flux is made non-dimensioal.
  void ConvertPrescribedData(const int                nIntegration,
                             const ENUM_FEM_VARIABLES FEMVariables,
                             std::vector<su2double *> &prescribedData);

  // Function, which gets the names of the prescribed variables for the given set.
  void GetNamesPrescibedVariables(const int                set,
                                  std::vector<std::string> &varNames);

  // Function, which returns the number of data sets that can be prescribed.
  int GetNSetsPrescribedVariables(void);

  // Function, which indicates whether prescribed data is expected.
  bool ExpectPrescribedData(void);

  // Function, which determines the indices for the prescribed data.
  void SetIndicesPrescribedData(const int set);

  //------------------------------------------
  // Private member functions.
  //------------------------------------------

  // Disabled default constructor.
  BCInflowSupersonicSubfaceClass();
};

//------------------------------------------------------------------------------

// Definition of the subsonic outflow boundary subface class.
class BCOutflowSubsonicSubfaceClass : public SubfaceBaseClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

  // Overloaded constructor. Read the subface data.
  BCOutflowSubsonicSubfaceClass(std::istringstream &istr);

  // Destructor. Nothing to be done.
  ~BCOutflowSubsonicSubfaceClass();

  //----------------------------------------------------------
  // Public member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which computes the boundary state (the right state) from the
  // given left state and the prescribed boundary data.
  void ComputeBoundaryState(const InputParamClass      *inputParam,
		                        const StandardElementClass *standardHex,
                            const int                   nInt,
                            su2double                 **solL,
                            su2double                 **dSolDxL,
                            su2double                 **dSolDyL,
                            su2double                 **dSolDzL,
                            const su2double             factNorm,
                            su2double                 **metricL,
                            su2double                 **prescribedData,
                            su2double                 **solR,
                            su2double                 **dSolDxR,
                            su2double                 **dSolDyR,
                            su2double                 **dSolDzR,
                            bool                       &heatFluxPrescribed,
                            su2double                 *&prescribedWallData,
                            su2double                  &wallPerm);

  // Function, which makes available the number of prescribed variables.
  int GetNVarPrescribed(void);
 
	// Function, which makes available the type of boundary 
	// condition prescribed. 
	int GetTypeBoundaryPrescribed(void);
  
private:
  //----------------------------------------------------------
  // Private member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which converts the prescribed data to the required form,
  // i.e. the heat flux is made non-dimensioal.
  void ConvertPrescribedData(const int                nIntegration,
                             const ENUM_FEM_VARIABLES FEMVariables,
                             std::vector<su2double *> &prescribedData);

  // Function, which gets the names of the prescribed variables for the given set.
  void GetNamesPrescibedVariables(const int                set,
                                  std::vector<std::string> &varNames);

  // Function, which returns the number of data sets that can be prescribed.
  int GetNSetsPrescribedVariables(void);

  // Function, which indicates whether prescribed data is expected.
  bool ExpectPrescribedData(void);

  // Function, which determines the indices for the prescribed data.
  void SetIndicesPrescribedData(const int set);

  //------------------------------------------
  // Private member functions.
  //------------------------------------------

  // Disabled default constructor.
  BCOutflowSubsonicSubfaceClass();
};

//------------------------------------------------------------------------------

// Definition of the supersonic outflow boundary subface class.
class BCOutflowSupersonicSubfaceClass : public SubfaceBaseClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

  // Overloaded constructor. Read the subface data.
  BCOutflowSupersonicSubfaceClass(std::istringstream &istr);

  // Destructor. Nothing to be done.
  ~BCOutflowSupersonicSubfaceClass();

  //----------------------------------------------------------
  // Public member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which computes the boundary state (the right state) from the
  // given left state and the prescribed boundary data.
  void ComputeBoundaryState(const InputParamClass      *inputParam,
		                        const StandardElementClass *standardHex,
                            const int                   nInt,
                            su2double                 **solL,
                            su2double                 **dSolDxL,
                            su2double                 **dSolDyL,
                            su2double                 **dSolDzL,
                            const su2double             factNorm,
                            su2double                 **metricL,
                            su2double                 **prescribedData,
                            su2double                 **solR,
                            su2double                 **dSolDxR,
                            su2double                 **dSolDyR,
                            su2double                 **dSolDzR,
                            bool                       &heatFluxPrescribed,
                            su2double                 *&prescribedWallData,
                            su2double                  &wallPerm);

  // Function, which makes available the number of prescribed variables.
  int GetNVarPrescribed(void);
 
	// Function, which makes available the type of boundary 
	// condition prescribed. 
	int GetTypeBoundaryPrescribed(void);
 
private:
  //------------------------------------------
  // Private member functions.
  //------------------------------------------

  // Disabled default constructor.
  BCOutflowSupersonicSubfaceClass();
};




//------------------------------------------------------------------------------
// Characteristic BC: Standard
//------------------------------------------------------------------------------

class BCStandardCharacteristicSubfaceClass : public SubfaceBaseClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

	// Overloaded constructor. Read the subface data.
	BCStandardCharacteristicSubfaceClass(std::istringstream &istr);

	// Destructor. Nothing to be done.
	virtual ~BCStandardCharacteristicSubfaceClass();


  //------------------------------------------
  // Public member variables.
  //------------------------------------------

  //----------------------------------------------------------
  // Public member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

	// Function, which sets the value of mMachAverageBoundary.
	void SetMachAverageBoundary(const su2double Mavg);

  // Pure virtual function, which makes available the number 
	// of prescribed variables.
  virtual int GetNVarPrescribed(void) = 0;

	// Pure virtual function, which makes available the type of 
	// boundary condition prescribed. 
	virtual int GetTypeBoundaryPrescribed(void) = 0;

	// Pure virtual function, which configures the tuning parameters 
	// required for a NSCBC.
	virtual void ConfigureParamNSCBC(const InputParamClass 			*inputParam, 
																	 const StandardElementClass *standardHex,
																	 su2double            		 **metric,
																	 su2double                   factNorm) = 0;

	// Function, which computes the weighted Mach number on this
	// local element boundary that is used to compute the average.
	// The average Mach number on this element surface is also computed.
	su2double WeightedMachElement(const InputParamClass *inputParam,
			                          const int              nInt,
										 					  su2double            **sol,
										 					  const su2double        factNorm,
										 					  su2double            **metric,
										 					  su2double             *intWeights);

  // Pure virtual function, which computes the boundary state 
	// (the right state) from the given left state and the 
	// prescribed boundary data.
  virtual void ComputeBoundaryState(const InputParamClass      *inputParam,
		                                const StandardElementClass *standardHex,
                                    const int                   nInt,
                                    su2double                 **solL,
                                    su2double                 **dSolDxL,
                                    su2double                 **dSolDyL,
                                    su2double                 **dSolDzL,
                                    const su2double             factNorm,
                                    su2double                 **metricL,
                                    su2double                 **prescribedData,
                                    su2double                 **solR,
                                    su2double                 **dSolDxR,
                                    su2double                 **dSolDyR,
                                    su2double                 **dSolDzR,
                                    bool                       &heatFluxPrescribed,
                                    su2double                 *&prescribedWallData,
                                    su2double                  &wallPerm) = 0;

protected:
	//------------------------------------------
  // Protected member variables.
  //------------------------------------------

	// Average Mach number over this element.
	su2double mMachAverageElement;
	// Average Mach number over entire boundary.
	su2double mMachAverageBoundary;

	// Index of the normal and transverse metrics.
	// [0]: starting index for normal metric,
	// [1]: starting index of one of the transverse metrics,
	// [2]: starting index of the remaining transverse metric.
	int mIndexMetric[3];

	// Index of the normal and transverse derivatives.
	// [0]: normal derivatie index,
	// [1]: transverse index for one of the derivatives,
	// [2]: transverse index for the remaining derivative.
	int mIndexGradient[3];

	//------------------------------------------
  // Protected member functions.
  //------------------------------------------

	// Function, which identifies the normal and transverse metrics
	// as well as the index of the normal wave amplitude.
  void IdentifyBoundaryIndices(su2double **metric);

	// Function, which converts the entropy variables and their gradient
	// into the primitive variables and their gradient. Note, the left 
	// state is the interior state and the right state is the external 
	// which is overwritten in this function by the primitive data.
	void ConvertEntropyToPrimitive(const int   nInt,
			                           su2double **solL,
			                           su2double **dSolDrL,
																 su2double **dSolDsL,
																 su2double **dSolDtL,
																 su2double **solR, 
																 su2double **dSolDrR,
																 su2double **dSolDsR,
																 su2double **dSolDtR);

	// Function, which reconstructs a boundary state from the 
	// boundary-conforming normal derivative of the primitive variables. 
	// This also reconstructs the gradient of that state. The result is 
	// a boundary state and gradient in entropy variables.
	void ReconstructBoundaryStateAndGradient(const StandardElementClass *standardHex, 
			                                     const su2double             factNorm,
			                                     su2double                 **solR,
			                                     su2double                 **dSolDrR,
																					 su2double                 **dSolDsR,
																					 su2double                 **dSolDtR);

private: 
	//------------------------------------------
  // Private member functions.
  //------------------------------------------

	// Pure virtual function, which converts the prescribed data 
	// to the required form, i.e. the pressure is made non-dimensional.
  virtual void ConvertPrescribedData(const int                nIntegration,
                                     const ENUM_FEM_VARIABLES FEMVariables,
                                     std::vector<su2double *> &prescribedData) = 0;

  // Pure virtual function, which gets the names of the prescribed 
	// variables for the given set.
  virtual void GetNamesPrescibedVariables(const int                set,
                                          std::vector<std::string> &varNames) = 0;

  // Pure virtual function, which returns the number of data sets 
	// that can be prescribed.
  virtual int GetNSetsPrescribedVariables(void) = 0;

  // Pure virtual function, which indicates whether prescribed data 
	// is expected.
  virtual bool ExpectPrescribedData(void) = 0;

  // Pure virtual function, which determines the indices for the 
	// prescribed data.
  virtual void SetIndicesPrescribedData(const int set) = 0;
};




//------------------------------------------------------------------------------
// Characteristic BC: Outlet
//------------------------------------------------------------------------------

class BCOutflowCharacteristicSubfaceClass : public BCStandardCharacteristicSubfaceClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

	// Overloaded constructor. Read the subface data.
	BCOutflowCharacteristicSubfaceClass(std::istringstream &istr);

	// Destructor. Nothing to be done.
	~BCOutflowCharacteristicSubfaceClass();


  //----------------------------------------------------------
  // Public member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which makes available the number of prescribed variables.
  int GetNVarPrescribed(void);

	// Function, which makes available the type of boundary 
	// condition prescribed. 
	int GetTypeBoundaryPrescribed(void);
 
	// Function, which configures the tuning parameters required for a NSCBC.
	void ConfigureParamNSCBC(const InputParamClass 			*inputParam, 
													 const StandardElementClass *standardHex,
													 su2double            		 **metric,
													 su2double                   factNorm);

  // Function, which computes the boundary state (the right state) from the
  // given left state and the prescribed boundary data.
  void ComputeBoundaryState(const InputParamClass      *inputParam,
		                        const StandardElementClass *standardHex,
                            const int                   nInt,
                            su2double                 **solL,
                            su2double                 **dSolDxL,
                            su2double                 **dSolDyL,
                            su2double                 **dSolDzL,
                            const su2double             factNorm,
                            su2double                 **metricL,
                            su2double                 **prescribedData,
                            su2double                 **solR,
                            su2double                 **dSolDxR,
                            su2double                 **dSolDyR,
                            su2double                 **dSolDzR,
                            bool                       &heatFluxPrescribed,
                            su2double                 *&prescribedWallData,
                            su2double                  &wallPerm);

private:
  //------------------------------------------
  // Private member variables.
  //------------------------------------------

	// Index of the incoming(+) normal wave amplitude.
  int mIndexPhi;
	// Lagrange derivative coefficient at boundary w.r.t 
	// normal (rr-)dimension.
	su2double mDerCoefficient;

	// Characteristic length scale.
	su2double mLengthScale;
	// Normal     relaxation coefficient.
	su2double mSigma;
	// Transverse relaxation coefficient (coupled).
	su2double mBeta_l;
	// Transverse relaxation coefficient (uncoupled).
	su2double mBeta_t;



	//------------------------------------------
  // Private member functions.
  //------------------------------------------

  // Function, which converts the prescribed data to the required form,
  // i.e. the pressure is made non-dimensional.
  void ConvertPrescribedData(const int                nIntegration,
                             const ENUM_FEM_VARIABLES FEMVariables,
                             std::vector<su2double *> &prescribedData);

  // Function, which gets the names of the prescribed variables for the given set.
  void GetNamesPrescibedVariables(const int                set,
                                  std::vector<std::string> &varNames);

  // Function, which returns the number of data sets that can be prescribed.
  int GetNSetsPrescribedVariables(void);

  // Function, which indicates whether prescribed data is expected.
  bool ExpectPrescribedData(void);

  // Function, which determines the indices for the prescribed data.
  void SetIndicesPrescribedData(const int set);

};




//------------------------------------------------------------------------------
// Characteristic BC: Inlet
//------------------------------------------------------------------------------

class BCInflowCharacteristicSubfaceClass : public BCStandardCharacteristicSubfaceClass
{
public:
  //------------------------------------------
  // Constructors and destructors.
  //------------------------------------------

	// Overloaded constructor. Read the subface data.
	BCInflowCharacteristicSubfaceClass(std::istringstream &istr);

	// Destructor. Nothing to be done.
	~BCInflowCharacteristicSubfaceClass();

  //----------------------------------------------------------
  // Public member functions, which overload the virtual
  // functions of the base class.
  //----------------------------------------------------------

  // Function, which makes available the number of prescribed variables.
  int GetNVarPrescribed(void);

	// Function, which makes available the type of boundary 
	// condition prescribed. 
	int GetTypeBoundaryPrescribed(void);
 
	// Function, which configures the tuning parameters required for a NSCBC.
	void ConfigureParamNSCBC(const InputParamClass 			*inputParam, 
													 const StandardElementClass *standardHex,
													 su2double            		 **metric,
													 su2double                   factNorm);

  // Function, which computes the boundary state (the right state) from the
  // given left state and the prescribed boundary data.
  void ComputeBoundaryState(const InputParamClass      *inputParam,
		                        const StandardElementClass *standardHex,
                            const int                   nInt,
                            su2double                 **solL,
                            su2double                 **dSolDxL,
                            su2double                 **dSolDyL,
                            su2double                 **dSolDzL,
                            const su2double             factNorm,
                            su2double                 **metricL,
                            su2double                 **prescribedData,
                            su2double                 **solR,
                            su2double                 **dSolDxR,
                            su2double                 **dSolDyR,
                            su2double                 **dSolDzR,
                            bool                       &heatFluxPrescribed,
                            su2double                 *&prescribedWallData,
                            su2double                  &wallPerm);

private:
  //------------------------------------------
  // Private member variables.
  //------------------------------------------

	// Index of the incoming(+) acoustic wave amplitude.
  int mIndexPhi;
	// Lagrange derivative coefficient at boundary w.r.t 
	// normal (rr-)dimension.
	su2double mDerCoefficient;

	// Characteristic length scale.
	su2double mLengthScale;
	// Normal     relaxation coefficient.
	su2double mSigma;
	// Transverse relaxation coefficient.
	su2double mBeta;



	//------------------------------------------
  // Private member functions.
  //------------------------------------------

	// Boolean, which indicates whether or not total conditions are specified.
	bool mTotalConditionsSpecified;

  // Function, which converts the prescribed data to the required form,
  // i.e. the pressure is made non-dimensional.
  void ConvertPrescribedData(const int                nIntegration,
                             const ENUM_FEM_VARIABLES FEMVariables,
                             std::vector<su2double *> &prescribedData);

  // Function, which gets the names of the prescribed variables for the given set.
  void GetNamesPrescibedVariables(const int                set,
                                  std::vector<std::string> &varNames);

  // Function, which returns the number of data sets that can be prescribed.
  int GetNSetsPrescribedVariables(void);

  // Function, which indicates whether prescribed data is expected.
  bool ExpectPrescribedData(void);

  // Function, which determines the indices for the prescribed data.
  void SetIndicesPrescribedData(const int set);

};



