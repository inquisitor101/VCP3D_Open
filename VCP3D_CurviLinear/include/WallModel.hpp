//------------------------------------------------------------------------------
// WallModel.hpp: Definition of the classes, which contain the wall models.
//                This include file should not be included by any source
//                code file directly. These files should just include
//                VCP3D_CurviLinear.hpp.
//------------------------------------------------------------------------------


#pragma once

// Forward declaration of the input parameter class.
class InputParamClass;

// Base class for the wall models.
class CWallModel {

public:

  // Constructor of the class.
  CWallModel();

  //Destructor of the class.
  virtual ~CWallModel(void);

  // Virtual function, which makes available the number of wall points.
  virtual int GetNWallModelPoints(void);
  
  // Virtual function, which computes the wall shear stress and heat flux
  virtual void WallShearStressAndHeatFlux(const su2double tExchange,
                                          const su2double velExchange,
                                          const su2double muExchange,
                                          const su2double pExchange,
                                          const bool      heatFluxPrescribed,
                                          const su2double wallData,
                                          su2double       &tauWall,
                                          su2double       &qWall,
                                          su2double       &ViscosityWall,
                                          su2double       &kOverCvWall,
                                          su2double       **workArray1,
                                          su2double       **workArray2);
};

// Class for the equilibrium wall model.
class CWallModel1DEQ : public CWallModel {

public:

  // Overloaded constructor of the class.
  CWallModel1DEQ(const InputParamClass *inputParam);

  // Destructor of the class.
  ~CWallModel1DEQ(void);

  // Function, which makes available the number of wall points.
  int GetNWallModelPoints(void);

  // Function, which computes the wall shear stress and heat flux from the data at the exchange location.
  void WallShearStressAndHeatFlux(const su2double tExchange,
                                  const su2double velExchange,
                                  const su2double muExchange,
                                  const su2double pExchange,
                                  const bool      heatFluxPrescribed,
                                  const su2double wallData,
                                  su2double       &tauWall,
                                  su2double       &qWall,
                                  su2double       &ViscosityWall,
                                  su2double       &kOverCvWall,
                                  su2double       **workArray1,
                                  su2double       **workArray2);
private:

  su2double h_wm;
  su2double expansionRatio;  // Stretching factor used for the wall model grid..
  int       numPoints;       // Number of points used in the wall model grid.
  su2double *y_cv;           // The coordinates in normal direction of the wall model grid (control volumes).
  su2double *y_fa;           // The coordinates in normal direction of the wall model grid (faces of CV).

  // Default constructor, disabled.
  CWallModel1DEQ();
};

// Class for the logarithmic wall model.
class CWallModelLogLaw : public CWallModel {
  
public:
  
  // Overloaded constructor.
  CWallModelLogLaw(const InputParamClass *inputParam);
  
  // Destructor of the class.
  ~CWallModelLogLaw(void);
  
  // Function, which computes the wall shear stress and heat flux from the data at the exchange location.
  void WallShearStressAndHeatFlux(const su2double tExchange,
                                  const su2double velExchange,
                                  const su2double muExchange,
                                  const su2double pExchange,
                                  const bool      heatFluxPrescribed,
                                  const su2double wallData,
                                  su2double       &tauWall,
                                  su2double       &qWall,
                                  su2double       &ViscosityWall,
                                  su2double       &kOverCvWall,
                                  su2double       **workArray1,
                                  su2double       **workArray2);
private:
  su2double h_wm;
  su2double C;  //Constant to match the Reichardt BL profile

  // Default constructor, disabled.
  CWallModelLogLaw();
};
