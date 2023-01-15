//------------------------------------------------------------------------------
// SGSModel.hpp: Definition of the classes, which contain the subgrid scale
//               models for LES. This include file should not be included by
//               any source code file directly. These files should just include
//               VCP3D_CurviLinear.hpp.
//------------------------------------------------------------------------------


#pragma once

// Base class for the subgrid scale models.
class CSGSModel
{
public:

  // Constructor of the class.
  CSGSModel();

  //Destructor of the class.
  virtual ~CSGSModel(void);

  // Virtual function, which computes the eddy viscosity for the
  // given items.
  virtual void EddyViscosity(const int       nItems,
                             const su2double lenScale,
                             su2double       **primVar,
                             su2double       **gradVel,
                             su2double       *muSGS);
};

// Class for the WALE subgrid scale model.
class CWALEModel : public CSGSModel
{
public:

  // Constructor of the class.
  CWALEModel();

  // Destructor of the class.
  ~CWALEModel(void);

  // Function, which computes the eddy viscosity for the given items.
  void EddyViscosity(const int       nItems,
                     const su2double lenScale,
                     su2double       **primVar,
                     su2double       **gradVel,
                     su2double       *muSGS);
};

// Class for the Vreman subgrid scale model.
class CVremanModel : public CSGSModel
{
public:
  
  // Constructor of the class.
  CVremanModel();
  
  // Destructor of the class.
  ~CVremanModel(void);

  // Function, which computes the eddy viscosity for the given items.
  void EddyViscosity(const int       nItems,
                     const su2double lenScale,
                     su2double       **primVar,
                     su2double       **gradVel,
                     su2double       *muSGS);
};
