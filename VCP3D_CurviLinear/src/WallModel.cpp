//------------------------------------------------------------------------------
// File, which contains the implementation for the wall model functions
// for large eddy simulations.
//------------------------------------------------------------------------------

#include "VCP3D_CurviLinear.hpp"

// Prototypes for Lapack functions.
// Not needed when the MKL is used.
#ifndef HAVE_MKL
extern "C"
{
#ifdef USE_SINGLE_PRECISION
void sgtsv_(int *n, int *nrhs, float *dl, float *d,
            float *du, float *b, int *ldb, int *info);
#else
void dgtsv_(int *n, int *nrhs, double *dl, double *d,
            double *du, double *b, int *ldb, int *info);
#endif
}
#endif

//------------------------------------------------------------------------------
// Implementation of the member functions of the CWallModel base class.
//------------------------------------------------------------------------------

CWallModel::CWallModel(){}

CWallModel::~CWallModel(void){}

int CWallModel::GetNWallModelPoints(void){return 0;}

void CWallModel::WallShearStressAndHeatFlux(const su2double tExchange,
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
                                            su2double       **workArray2) {}

//------------------------------------------------------------------------------
// Implementation of the member functions of the CWallModel1DEQ class.
//------------------------------------------------------------------------------

// Overloaded constructor of the class.
CWallModel1DEQ::CWallModel1DEQ(const InputParamClass *inputParam)
  :  CWallModel() {

  /* Copy the data into the member variables. */
  numPoints      = inputParam->mNPointsGridWallModel;
  h_wm           = inputParam->mExchangeLocation;
  expansionRatio = inputParam->mExpansionRatioWallModel;

  /* Modify numPoints, such that the number of faces is a multiple of
     the vector length. */
  int nfa   = numPoints - 1;
  nfa       = ((nfa+vecLen1D-1)/vecLen1D)*vecLen1D;
  numPoints = nfa + 1;

  /* Allocate the memory for the coordinates of the grid points used
     in the 1D equilibrium wall model. */
  y_cv = (su2double *) AllocateMemory(numPoints*sizeof(su2double));
  y_fa = (su2double *) AllocateMemory(nfa*sizeof(su2double));

  if(!y_cv || !y_fa)
    Terminate("CWallModel1DEQ::CWallModel1DEQ", __FILE__, __LINE__,
              "Memory allocation failure for y_cv and y_fa");

  /* Determine the scaled version of the normal coordinates, where the
    first normal coordinate is simply 1.0. */
  su2double currentHeight = one;

  y_cv[0] = zero;
  for(int i=1; i<numPoints; ++i) {
    y_cv[i]   = y_cv[i-1] + currentHeight;
    y_fa[i-1] = half * (y_cv[i] + y_cv[i-1]);
    currentHeight *= expansionRatio;
  }

  /* Scale the coordinates, such that the exchange point coincides
     with the last cell center. */
  const su2double scaleFact = h_wm/y_cv[numPoints-1];
  for(int i=0; i<nfa; ++i)
    y_fa[i] *= scaleFact;

  for(int i=0; i<numPoints; ++i)
    y_cv[i] *= scaleFact;
}

// Destructor of the class. Release the allocated memory.
CWallModel1DEQ::~CWallModel1DEQ(void){

  FreeMemory((void **) &y_cv);
  FreeMemory((void **) &y_fa);
}

// Function, which makes available the number of wall points.
int CWallModel1DEQ::GetNWallModelPoints(void){return numPoints;}

// Function, which computes the wall shear stress and heat flux from the data at the exchange location.
void CWallModel1DEQ::WallShearStressAndHeatFlux(const su2double tExchange,
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
                                                su2double       **workArray2) {

  /* Set the maximum number of iterations of the iterative algorithm.
     and the tolerance for this algorithm.  */
  const int max_iter  = 25;
  const su2double tol = (su2double) 1e-3;

  /* Model constant in wall units where the laminar and eddy viscosity are
     of the same order of magnitude. */
  const su2double A = (su2double) 17.0;

  /* Compute the dimensional laminar viscosity and thermal conductity. */
  const su2double mu_lam = mu*pRef/uRef;
  const su2double k_lam  = Cp*mu_lam/Prandtl_Lam;

  /* Set the wall temperature. If the heat flux is prescribed use the temperature
     in the exchange point, otherwise the wall temperature (stored in wallData). */
  su2double TWall = heatFluxPrescribed ? tExchange : wallData;

  /* Compute the temperature difference between the wall and first internal
     point when the heatflux is prescribed. Note that the heat flux value
     must be made dimensional. */
  const su2double deltaTHeatflux = k_lam*(y_cv[1]-y_cv[0])*wallData*pRef*uRef;

  /* Set the return value for the viscosity of the wall and kOverCv.
     As mu is constant, these values can be set here and there is
     no need to change them. */
  ViscosityWall = mu_lam;
  kOverCvWall   = k_lam/Cv;

  /* Set the pointers for the help variables. */
  su2double *u        = workArray1[0];
  su2double *T        = workArray1[1];
  su2double *jacFace  = workArray1[2];
  su2double *fluxFace = workArray1[3];
  su2double *lower    = workArray2[0];
  su2double *diagonal = workArray2[1];
  su2double *upper    = workArray2[2];
  su2double *rhs      = workArray2[3];

  /* Set tau wall to initial guess. */
  tauWall = mu_lam*velExchange/h_wm;
  qWall   = -k_lam*(tExchange-TWall)/h_wm;

  /* Store the number of faces. */
  const int nfa = numPoints - 1;

  /* Initialize the velocity and temperature to a linear distribution.
     The last value is set outside the loop to improve vector performance
     as nfa is a multiple of the vector length. */
#pragma omp simd
  for(int i=0; i<nfa; ++i)
  {
    const su2double fact = y_cv[i]*h_wm;
    u[i] = fact*velExchange;
    T[i] = TWall + fact*(tExchange - TWall);
  }

  u[nfa] = velExchange;
  T[nfa] = tExchange;

  /* Set parameters for control. */
  bool converged = false;
  int iter = 0;

  /* Start of the iterative algorithm. */
  while (converged == false){
    
    /* Update the iteration counter and terminate if too many iterations
       are needed. */
    iter += 1;
    if (iter == max_iter)
    {
      std::cout.precision(15);
      std::cout << std::endl;
      std::cout << "tExchange:   " << tExchange   << std::endl;
      std::cout << "velExchange: " << velExchange << std::endl;
      std::cout << "muExchange:  " << muExchange  << std::endl;
      std::cout << "pExchange:   " << pExchange   << std::endl;
      std::cout << "wallData:    " << wallData    << std::endl;
      std::cout << std::endl;

      Terminate("CWallModel1DEQ::WallShearStressAndHeatFlux", __FILE__, __LINE__,
                "Iterative algorithm did not converge");
    }

    /* Store the previous values of the wall shear stress, heat flux
       and wall temperature. */
    const su2double tauWall_prev = tauWall;
    const su2double qWall_prev   = qWall;
    const su2double TWall_prev   = TWall;

    /* Loop over the faces to compute the face Jacobian of the momentum equation. */
#pragma omp simd
    for(int i=0; i<nfa; ++i)
    {
      /* Compute the density at the face. */
      const su2double TFace = half*(T[i+1] + T[i]);
      const su2double rho   = pExchange/(RGas*TFace);

      /* Compute the y+ value of the face. */
      const su2double utau = SQRT(tauWall/rho);
      const su2double yplus = rho*y_fa[i]*utau/mu_lam;

      /* Compute the eddy viscosity. */
      const su2double tmp = one - EXP((-yplus/A));
      const su2double D   = tmp*tmp;
      const su2double mut = rho*Karman*y_fa[i]*utau*D;

      /* Compute the value of the Jacobian at this face. */
      jacFace[i] = (mu_lam + mut)/(y_cv[i+1] - y_cv[i]);
    }

    /* Create the tri-diagonal matrix for the momentum equation.
       Note that only the interior points must be solved.
       First internal point. */
    rhs[0]      = -jacFace[0]*u[0];
    diagonal[0] = -(jacFace[0]+jacFace[1]);
    upper[0]    =  jacFace[1];

    /* The regular internal points. */
#pragma omp simd
    for(int i=1; i<(numPoints-3); ++i)
    {
      rhs[i]      =  zero;
      lower[i-1]  =  jacFace[i];
      diagonal[i] = -(jacFace[i] + jacFace[i+1]);
      upper[i]    =  jacFace[i+1];
    }

    /* The last internal point. */
    rhs[numPoints-3]      = -jacFace[numPoints-2]*u[numPoints-1];
    diagonal[numPoints-3] = -(jacFace[numPoints-3]+jacFace[numPoints-2]);
    lower[numPoints-4]    =  jacFace[numPoints-3];

    /* Solve the tri-diagonal matrix for the velocities. */
    int info, nrhs = 1, size = numPoints-2;

#ifdef USE_SINGLE_PRECISION
    sgtsv_(&size, &nrhs, lower, diagonal, upper, rhs, &size, &info);
#else
    dgtsv_(&size, &nrhs, lower, diagonal, upper, rhs, &size, &info);
#endif
    if (info != 0)
      Terminate("CWallModel1DEQ::WallShearStressAndHeatFlux", __FILE__, __LINE__,
                "\"Unsuccessful call to dgtsv_ for the velocity field\"");

    /* Copy the data of the internal points to u. */
    for(int i=0; i<size; ++i)
      u[i+1] = rhs[i];

    /* Compute the new value of tauWall. */
    tauWall = mu_lam*(u[1] - u[0])/(y_cv[1] - y_cv[0]);

    /* Loop over the faces to compute the face flux and Jacobian
       of the energy equation. */
#pragma omp simd
    for(int i=0; i<nfa; ++i)
    {
      /* Compute the density at the face. */
      const su2double TFace = half*(T[i+1] + T[i]);
      const su2double rho   = pExchange/(RGas*TFace);

      /* Compute the y+ value of the face. */
      const su2double utau = SQRT(tauWall/rho);
      const su2double yplus = rho*y_fa[i]*utau/mu_lam;

      /* Compute the eddy viscosity. */
      const su2double tmp = one - EXP((-yplus/A));
      const su2double D   = tmp*tmp;
      const su2double mut = rho*Karman*y_fa[i]*utau*D;

      /* Compute the value of the velocity flux and Jacobian for the
         temperature at this face. */
      const su2double dyInv = one/(y_cv[i+1] - y_cv[i]);

      fluxFace[i] = half*dyInv*(mu_lam + mut)*(u[i+1]+u[i])*(u[i+1]-u[i]);
      jacFace[i]  = Cp*dyInv*(mu_lam/Prandtl_Lam + mut/Prandtl_Turb);
    }

    /* Create the tri-diagonal matrix for the energy equation.
       Note that only the interior points must be solved.
       First internal point. Note that for an adiabatic
       wall the wall temperature depends on the first interior
       point, which must be taken into account in the right
       hand side and diagonal. */
    rhs[0]      =  fluxFace[0] - fluxFace[1];
    diagonal[0] = -jacFace[1];
    upper[0]    =  jacFace[1];

    if( heatFluxPrescribed ) rhs[0] -= jacFace[0]*deltaTHeatflux;
    else {rhs[0] -= jacFace[0]*T[0]; diagonal[0] -= jacFace[0];}

    /* The regular internal points. */
#pragma omp simd
    for(int i=1; i<(numPoints-3); ++i)
    {
      rhs[i]      =  fluxFace[i] - fluxFace[i+1];
      lower[i-1]  =  jacFace[i];
      diagonal[i] = -(jacFace[i] + jacFace[i+1]);
      upper[i]    =  jacFace[i+1];
    }

    /* The last internal point. */
    rhs[numPoints-3]      = -jacFace[numPoints-2]*T[numPoints-1]
                          +  fluxFace[numPoints-3] - fluxFace[numPoints-2];
    diagonal[numPoints-3] = -(jacFace[numPoints-3]+jacFace[numPoints-2]);
    lower[numPoints-4]    =  jacFace[numPoints-3];

    /* Solve the tri-diagonal matrix for the temperatures. */
#ifdef USE_SINGLE_PRECISION
    sgtsv_(&size, &nrhs, lower, diagonal, upper, rhs, &size, &info);
#else
    dgtsv_(&size, &nrhs, lower, diagonal, upper, rhs, &size, &info);
#endif
    if (info != 0)
      Terminate("CWallModel1DEQ::WallShearStressAndHeatFlux", __FILE__, __LINE__,
                "\"Unsuccessful call to dgtsv_ for the temperature field\"");

    /* Copy the data of the internal points to T. */
    for(int i=0; i<size; ++i)
      T[i+1] = rhs[i];

    /* Compute the new wall temperature for a prescribed heat flux. */
    if( heatFluxPrescribed ) TWall = T[0] = T[1] + deltaTHeatflux;

    /* Compute the new value of qWall. */
    qWall = -k_lam * (T[1] - T[0])/(y_cv[1]-y_cv[0]);

    /* Compute the deltas for the wall shear stress, heat flux and
       wall temperature and the corresponding tolerances. */
    const su2double dTauWall = FABS(tauWall - tauWall_prev);
    const su2double dQWall   = FABS(qWall - qWall_prev);
    const su2double dTWall   = FABS(TWall - TWall_prev);

    const su2double tauWallAbs = FABS(tauWall);
    const su2double tauWallTol = tol*std::max(tol,tauWallAbs);

    const su2double qWallAbs = FABS(qWall);
    const su2double qWallTol = tol*std::max(tol,qWallAbs);

    const su2double TWallAbs = FABS(TWall);
    const su2double TWallTol = tol*std::max(tol,TWallAbs);

    /* Check for convergence. */
    if((dTauWall <= tauWallTol) && (dQWall <= qWallTol) &&
       (dTWall   <= TWallTol)) converged = true;
  }

  /* Final check of the Y+. */
  const su2double rho = pExchange / (RGas * T[0]);
  if (rho*y_cv[1]*SQRT(tauWall/rho)/mu > one)
    Terminate("CWallModel1DEQ::WallShearStressAndHeatFlux", __FILE__, __LINE__,
              "\"Y+ greater than one: Increase the number of points or growth ratio.");
}

//------------------------------------------------------------------------------
// Implementation of the member functions of the CWallModelLogLaw class.
//------------------------------------------------------------------------------

// Overloaded constructor.
CWallModelLogLaw::CWallModelLogLaw(const InputParamClass *inputParam)
  :  CWallModel() {

  C = 5.25; /* Constant to match the Reichardt BL profile ->  C = 4.1;  or 5.25. */

  /* Add the term with the von Karman constant. */
  C -=  LOG(Karman)/Karman;

  /* Copy the exchange height from the input parameters. */
  h_wm = inputParam->mExchangeLocation;
}

// Destructor of the class.
CWallModelLogLaw::~CWallModelLogLaw(void){}

// Function, which computes the wall shear stress and heat flux from the data at the exchange location.
void CWallModelLogLaw::WallShearStressAndHeatFlux(const su2double tExchange,
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
                                                  su2double       **workArray2) {

  /* Some constants that appear in Reichardt's profile and Kader's law. */
  const su2double oneOv11  = one/eleven;
  const su2double third    = one/three;
  const su2double oneOv100 = (su2double) 0.01;
  const su2double C1Kader  = (su2double) 3.85;
  const su2double C2Kader  = (su2double) 1.3;
  const su2double C3Kader  = (su2double) 2.12;

  /* Compute the dimensional laminar viscosity. */
  const su2double mu_lam = mu*pRef/uRef;

  /* Set the wall temperature. If the heat flux is prescribed use the temperature
     in the exchange point, otherwise the wall temperature (stored in wallData). */
  const su2double TWall = heatFluxPrescribed ? tExchange : wallData;

  /* Get the required data from the fluid model. */
  const su2double rho_wall = pExchange / (TWall * RGas);
  const su2double nu_wall  = mu_lam / rho_wall;

  /* Initial guess of the friction velocity. */ 
  su2double u_tau = std::max(0.01*velExchange, 1.e-5);
  
  /* Set parameters for control of the Newton iteration. */
  bool converged = false;
  unsigned short iter = 0, max_iter = 50;
  const su2double tol = (su2double) 1e-3;
  
  while (converged == false){
    
    /* Increase the number of iterations and terminate if too
       many are needed. */
    iter += 1;
    if (iter == max_iter)
      Terminate("CWallModelLogLaw::WallShearStressAndHeatFlux", __FILE__,
                __LINE__, "Newton method for u_tau did not converge");

    /* Compute the yPlus value and the derivative w.r.t. u_tau. */
    const su2double u_tau0        = u_tau;
    const su2double dy_plusdu_tau = h_wm/nu_wall;
    const su2double y_plus        = u_tau0*dy_plusdu_tau;

    /* Determine the different parts of the Reichardt boundary layer
       analytical law and its derivative w.r.t. y_plus. */
    const su2double val1  = LOG(Karman*y_plus + one)/Karman;
    const su2double dVal1 = one/(Karman*y_plus + one);

    const su2double tmp2  = C*EXP(-y_plus*oneOv11);
    const su2double val2  = C - tmp2;
    const su2double dVal2 = tmp2*oneOv11;

    const su2double tmp3  = -C*EXP(-y_plus*third)*oneOv11;
    const su2double val3  =  tmp3*y_plus;
    const su2double dVal3 =  tmp3 - third*val3;

    /* Determine the function for which the root must be found, which is
       velExchange/u_tau - profile = 0, where profile is the Reichardt
       function, which is val1 + val2 + val3. */
    const su2double fval = velExchange/u_tau0 - val1 - val2 - val3;

    /* Determine the derivative of fval w.r.t. u_tau0. */
    const su2double fprime = -velExchange/(u_tau0*u_tau0)
                           - dy_plusdu_tau*(dVal1 + dVal2 + dVal3);
  
    /* Newton method.  */
    const su2double newton_step = fval/fprime;
    u_tau = u_tau0 - newton_step;
    
    /* Define a norm */
    if (FABS(one - u_tau/u_tau0) < tol) converged = true;
  }
  
  /* Compute the wall shear stress. */
  tauWall = rho_wall * u_tau*u_tau;

  /* Check if the heat flux is prescribed. */
  if( heatFluxPrescribed )
  {
    /* Heat flux is prescribed and stored in wallData. Set the value. */
    qWall = wallData;
  }
  else
  {
    /* Isothermal wall. Kader's law will be used to approximate the
       variations of the temperature inside the boundary layer.
       First compute y+, y+^4 and expressions involving Prandtl.  */
    const su2double y_plus   = u_tau*h_wm/nu_wall;
    const su2double abvyPPr  = y_plus*Prandtl_Lam;
    const su2double abvyPPr4 = abvyPPr*abvyPPr*abvyPPr*abvyPPr;
    const su2double Pr3      = Prandtl_Lam*Prandtl_Lam*Prandtl_Lam;
    const su2double abvPr    = C1Kader*POW(Prandtl_Lam,third) - C2Kader;

    /* Compute the actual heat flux. */
    const su2double lhs = - ((tExchange - TWall) * rho_wall * Cp * u_tau);
    const su2double Gamma = - (oneOv100 * abvyPPr4)/(one + five*y_plus*Pr3);
    const su2double rhs_1 = Prandtl_Lam * y_plus * EXP(Gamma);
    const su2double rhs_2 = (C3Kader*LOG(one+y_plus) + abvPr*abvPr + C3Kader*LOG(Prandtl_Lam)) * EXP(one/Gamma);
    qWall = lhs/(rhs_1 + rhs_2);
  }

  /* Set the values of the viscosity and thermal conductivity divided by Cv
     for the wall. */
  ViscosityWall = mu_lam;
  kOverCvWall   = (mu_lam*Cp/Prandtl_Lam)/Cv;
}
