//------------------------------------------------------------------------------
// File, which contains the implementation for the member functions of the
// subgrid scale classes.
//------------------------------------------------------------------------------

#include "VCP3D_CurviLinear.hpp"

//------------------------------------------------------------------------------
// Implementation of the member functions of the CSGSModel base class.
//------------------------------------------------------------------------------

CSGSModel::CSGSModel(){}

CSGSModel::~CSGSModel(){}

void CSGSModel::EddyViscosity(const int       nItems,
                              const su2double lenScale,
                              su2double       **primVar,
                              su2double       **gradVel,
                              su2double       *muSGS){}

//------------------------------------------------------------------------------
// Implementation of the member functions of the CWALEModel class.
//------------------------------------------------------------------------------

// Constructor of the class.
CWALEModel::CWALEModel()
  :  CSGSModel() {}

// Destructor of the class.
CWALEModel::~CWALEModel(){}

// Function, which computes the eddy viscosity for the given items.
void CWALEModel::EddyViscosity(const int       nItems,
                               const su2double lenScale,
                               su2double       **primVar,
                               su2double       **gradVel,
                               su2double       *muSGS)
{
  // Abbreviate the constant 1/3.
  const su2double third = one/three;

  // Compute the square of the length scale of the WALE model.
  const su2double Cw = (su2double) 0.325;

  const su2double len  = Cw*lenScale;
  const su2double len2 = len*len;

  // Loop over the number of items.
#pragma omp simd
  for(int l=0; l<nItems; ++l)
  {
    // Easier storage of density.
    const su2double rho = primVar[0][l];

    // Easier storage of the velocity gradients.
    const su2double dudx = gradVel[0][l];
    const su2double dvdx = gradVel[1][l];
    const su2double dwdx = gradVel[2][l];

    const su2double dudy = gradVel[3][l];
    const su2double dvdy = gradVel[4][l];
    const su2double dwdy = gradVel[5][l];

    const su2double dudz = gradVel[6][l];
    const su2double dvdz = gradVel[7][l];
    const su2double dwdz = gradVel[8][l];

    // Compute the strain rate tensor, which is symmetric.
    const su2double S11 = dudx, S22 = dvdy, S33 = dwdz;
    const su2double S12 = half*(dudy + dvdx);
    const su2double S13 = half*(dudz + dwdx);
    const su2double S23 = half*(dvdz + dwdy);

    // Compute the values of the Sd tensor. First without the trace
    // correction of the diagonal terms.
    su2double Sd11 = dudx*dudx + dudy*dvdx + dudz*dwdx;
    su2double Sd22 = dvdx*dudy + dvdy*dvdy + dvdz*dwdy;
    su2double Sd33 = dwdx*dudz + dwdy*dvdz + dwdz*dwdz;

    const su2double Sd12 = half*(dudx*dudy + dudy*dvdy + dudz*dwdy
                         +       dvdx*dudx + dvdy*dvdx + dvdz*dwdx);
    const su2double Sd13 = half*(dudx*dudz + dudy*dvdz + dudz*dwdz
                         +       dwdx*dudx + dwdy*dvdx + dwdz*dwdx);
    const su2double Sd23 = half*(dvdx*dudz + dvdy*dvdz + dvdz*dwdz
                         +       dwdx*dudy + dwdy*dvdy + dwdz*dwdy);

    // Correct the diagonal elements, such that the trace of the Sd tensor is zero.
    const su2double thirdTrace = third*(Sd11 + Sd22 + Sd33);

    Sd11 -= thirdTrace;
    Sd22 -= thirdTrace;
    Sd33 -= thirdTrace;

    // Compute the summation of both tensors.
    const su2double sumS  = S11*S11 + S22*S22 + S33*S33
                          + two*(S12*S12 + S13*S13 + S23*S23);
    const su2double sumSd = Sd11*Sd11 + Sd22*Sd22 + Sd33*Sd33
                          + two*(Sd12*Sd12 + Sd13*Sd13 + Sd23*Sd23);
 
    // Compute the kinematic eddy viscosity.
    const su2double sumSdPow3_2 = sumSd*SQRT(sumSd);
    const su2double sumSdPow5_4 = SQRT(sumSdPow3_2*sumSd);
    const su2double sumSPow5_2  = sumS*sumS*SQRT(sumS);
    const su2double denom       = sumSPow5_2 + sumSdPow5_4;

    const su2double nuEddy = len2*sumSdPow3_2/std::max(denom, epsSmall);

    // Compute the SGS dynamic viscosity.
    muSGS[l] = rho*nuEddy;
  }
}

//------------------------------------------------------------------------------
// Implementation of the member functions of the CVremanModel class.
//------------------------------------------------------------------------------

// Constructor of the class.
CVremanModel::CVremanModel()
  :  CSGSModel() {}

// Destructor of the class.
CVremanModel::~CVremanModel(){}

// Function, which computes the eddy viscosity for the given items.
void CVremanModel::EddyViscosity(const int       nItems,
                                 const su2double lenScale,
                                 su2double       **primVar,
                                 su2double       **gradVel,
                                 su2double       *muSGS)
{
  // Abbreviate the constant 1/3.
  const su2double third = one/three;

  // Constant of the Vreman model = 2.5*Cs*Cs where Cs is the Smagorinsky constant.
  const su2double cVreman = (su2double) 0.07;

  // Compute the square of the length scale.
  const su2double len2 = lenScale*lenScale;

  // Loop over the number of items.
#pragma omp simd
  for(int l=0; l<nItems; ++l)
  {
    // Easier storage of density.
    const su2double rho = primVar[0][l];

    // Easier storage of the velocity gradients.
    const su2double dudx = gradVel[0][l];
    const su2double dvdx = gradVel[1][l];
    const su2double dwdx = gradVel[2][l];

    const su2double dudy = gradVel[3][l];
    const su2double dvdy = gradVel[4][l];
    const su2double dwdy = gradVel[5][l];

    const su2double dudz = gradVel[6][l];
    const su2double dvdz = gradVel[7][l];
    const su2double dwdz = gradVel[8][l];

    // Compute the diagonal elements of the velocity gradient tensor.
    su2double alpha11 = dudx;
    su2double alpha22 = dvdy;
    su2double alpha33 = dwdz;

    // Remove the trace (is this really necessary?).
    const su2double tmp = third*(alpha11 + alpha22 + alpha33);
    alpha11 -= tmp;
    alpha22 -= tmp;
    alpha33 -= tmp;

    // Compute the off-diagonal terms of the velocity gradient tensor.
    const su2double alpha12 = dudy;
    const su2double alpha13 = dudz;
    const su2double alpha23 = dvdz;

    const su2double alpha21 = dvdx;
    const su2double alpha31 = dwdx;
    const su2double alpha32 = dwdy;

    // Compute the square of the alpha tensor, multiplied by the
    // square of the length scale.
    const su2double beta11  = len2*(alpha11*alpha11 + alpha12*alpha12 + alpha13*alpha13);
    const su2double beta12  = len2*(alpha11*alpha21 + alpha12*alpha22 + alpha13*alpha23);
    const su2double beta13  = len2*(alpha11*alpha31 + alpha12*alpha32 + alpha13*alpha33);
    const su2double beta22  = len2*(alpha21*alpha21 + alpha22*alpha22 + alpha23*alpha23);
    const su2double beta23  = len2*(alpha21*alpha31 + alpha22*alpha32 + alpha23*alpha33);
    const su2double beta33  = len2*(alpha31*alpha31 + alpha32*alpha32 + alpha33*alpha33);

    // Compute the kinematic eddy viscosity.
    su2double B = beta11*beta22-beta12*beta12+beta11*beta33-beta13*beta13+beta22*beta33-beta23*beta23;
    B = (B + FABS(B))*half;

    const su2double denom = alpha11*alpha11+alpha22*alpha22+alpha33*alpha33 +
                            alpha12*alpha12+alpha13*alpha13+alpha23*alpha23 +
                            alpha21*alpha21+alpha31*alpha31+alpha32*alpha32;

    const su2double nuEddy = SQRT(B/std::max(denom, epsSmall));

    // Compute the SGS dynamic viscosity.
    muSGS[l] = rho*cVreman*nuEddy;
  }
}
