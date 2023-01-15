//------------------------------------------------------------------------------
// File, which contains the implementation of the member function
// ComputeBoundaryState of subface classes, which correspond to physical
// boundary conditions.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"

//------------------------------------------------------------------------------

// Function, which computes the boundary state (the right state) from the
// given left state and the prescribed boundary data. It must be overwritten
// by the derived class.
void SubfaceBaseClass::ComputeBoundaryState(const InputParamClass      *inputParam,
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
                                            su2double                  &wallPerm)
{
  Terminate("SubfaceBaseClass::ComputeBoundaryState", __FILE__, __LINE__,
            "Function must be overwritten by the derived class");
}

//------------------------------------------------------------------------------

// Function, which applies the inviscid wall boundary conditions.
void SubfaceBaseClass::ApplyInviscidWallBC(const InputParamClass *inputParam,
                                           const int             nInt,
                                           su2double             **solL,
                                           su2double             **metricL,
                                           su2double             **solR)
{
  // Make a distinction between the working variables.
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      // The conservative variables are used as working variables.
      // Abbreviate gamma-1 and its inverse.
      const su2double gm1   = GamConstant - one;
      const su2double ovgm1 = one/gm1;

      // Loop over the number of integration points.
#pragma omp simd
      for(int l=0; l<nInt; ++l)
      {
        // Compute the velocities and pressure of the left state.
        su2double rhoInv = one/solL[0][l];
        su2double u      = rhoInv*solL[1][l];
        su2double v      = rhoInv*solL[2][l];
        su2double w      = rhoInv*solL[3][l];
        su2double p      = gm1*(solL[4][l] - half*solL[0][l]*(u*u + v*v + w*w));

        // Compute the normal velocity.
        su2double un = u*metricL[0][l] + v*metricL[1][l] + w*metricL[2][l];

        // Subtract the normal velocity to obtain the tangential velocity.
        u -= un*metricL[0][l];
        v -= un*metricL[1][l];
        w -= un*metricL[2][l];

        // Compute the right state.
        solR[0][l] = solL[0][l];
        solR[1][l] = solL[0][l]*u;
        solR[2][l] = solL[0][l]*v;
        solR[3][l] = solL[0][l]*w;
        solR[4][l] = p*ovgm1 + half*solL[0][l]*(u*u + v*v + w*w);
      }

      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      // The entropy variables are used as working variables.
      // Loop over the number of integration points.
#pragma omp simd
      for(int l=0; l<nInt; ++l)
      {
        // Compute the normal component of the variables involving the velocity.
        const su2double Vn = solL[1][l]*metricL[0][l] + solL[2][l]*metricL[1][l]
                           + solL[3][l]*metricL[2][l];

        // Compute the right state.
        solR[0][l] = solL[0][l] - half*Vn*Vn/solL[4][l];
        solR[1][l] = solL[1][l] - Vn*metricL[0][l];
        solR[2][l] = solL[2][l] - Vn*metricL[1][l];
        solR[3][l] = solL[3][l] - Vn*metricL[2][l];
        solR[4][l] = solL[4][l];
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

// Function, which computes the boundary state (the right state) from the
// given left state and the prescribed boundary data.
void BCFarfieldSubfaceClass::ComputeBoundaryState(const InputParamClass      *inputParam,
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
                                                  su2double                  &wallPerm)
{
  // Farfield boundary condition. For now use a simple approach and just set
  // the right state to the prescribed state.
  for(int l=0; l<nVar; ++l)
  {
#pragma omp simd
    for(int m=0; m<nInt; ++m)
      solR[l][m] = prescribedData[l][m];
  }
}

//------------------------------------------------------------------------------

// Function, which computes the boundary state (the right state) from the
// given left state and the prescribed boundary data.
void BCIsothermalWallSubfaceClass::ComputeBoundaryState(const InputParamClass      *inputParam,
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
                                                        su2double                  &wallPerm)
{
  // Set the boolean for heatFluxPrescribed, set wallPerm to zero and
  // set the pointer for the prescribed data for the wall temperature.
  heatFluxPrescribed = false;
  wallPerm           = zero;
  prescribedWallData = prescribedData[0];

  // Check for a wall model.
  if( mUseWallModel )
  {
    // A a wall model is used, so inviscid boundary conditions are applied.
    ApplyInviscidWallBC(inputParam, nInt, solL, metricL, solR);
  }
  else
  {
    // No wall model. Make a distinction between the working variables.
    switch( inputParam->mFEMVariables )
    {
      case CONSERVATIVE_VARIABLES:
      {
        // The conservative variables are used. Determine the conversion factor from the
        // dimensional temperature to a non-dimensional internal energy per unit mass.
        const su2double eRefInv = Cv/(uRef*uRef);

        // Loop over the integration points, set the velocities to zero, extrapolate the
        // density and use the prescribed temperature.
#pragma omp simd
        for(int l=0; l<nInt; ++l)
        {
          solR[0][l] = solL[0][l];
          solR[1][l] = zero;
          solR[2][l] = zero;
          solR[3][l] = zero;
          solR[4][l] = solL[0][l]*prescribedData[0][l]*eRefInv;
        }

        break;
      }

      //------------------------------------------------------------------------

      case ENTROPY_VARIABLES:
      {
        // The entropy variables are used as working variables.
        // Determine the factor multiplying the inverse of the dimensional
        // temperature to obain the non-dimensional last entropy variable.
        const su2double fact = -pRef/(RGas*rhoRef);

        // Loop over the integration points, set the velocities to zero, extrapolate the
        // entropy and use the prescribed temperature.
#pragma omp simd
        for(int l=0; l<nInt; ++l)
        {
          solR[0][l] = solL[0][l] - half*(solL[1][l]*solL[1][l] + solL[2][l]*solL[2][l]
                     +                    solL[3][l]*solL[3][l])/solL[4][l];
          solR[1][l] = zero;
          solR[2][l] = zero;
          solR[3][l] = zero;
          solR[4][l] = fact/prescribedData[0][l];
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
  }
}

//------------------------------------------------------------------------------

// Function, which computes the boundary state (the right state) from the
// given left state and the prescribed boundary data.
void BCHeatFluxWallSubfaceClass::ComputeBoundaryState(const InputParamClass      *inputParam,
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
                                                      su2double                  &wallPerm)
{
  // Set the boolean for heatFluxPrescribed, set wallPerm to zero and set
  // the pointer for the prescribed data for the heat flux.
  heatFluxPrescribed = true;
  wallPerm           = zero;
  prescribedWallData = prescribedData[0];

  // Check for a wall model.
  if( mUseWallModel )
  {
    // A a wall model is used, so inviscid boundary conditions are applied.
    ApplyInviscidWallBC(inputParam, nInt, solL, metricL, solR);
  }
  else
  {
    // No wall model. Make a distinction between the working variables.
    switch( inputParam->mFEMVariables )
    {
      case CONSERVATIVE_VARIABLES:
      {
        // The conservative variables are used. Set the velocities to zero
        // and extrapolate the density and pressure for the given number
        // of points.
#pragma omp simd
        for(int l=0; l<nInt; ++l)
        {
          solR[0][l] = solL[0][l];
          solR[1][l] = zero;
          solR[2][l] = zero;
          solR[3][l] = zero;
          solR[4][l] = solL[4][l] - half*(solL[1][l]*solL[1][l]
                     +                    solL[2][l]*solL[2][l]
                     +                    solL[3][l]*solL[3][l])/solL[0][l];
        }

        break;
      }

      //-------------------------------------------------------------------------

      case ENTROPY_VARIABLES:
      {
        // The entropy variables are used as working variables. Set Set the velocities
        // to zero and extrapolate the density and pressure for the given number
        // of points. Extrapolating density and pressure means that also the entropy
        // is extrapolated.
#pragma omp simd
        for(int l=0; l<nInt; ++l)
        {
          solR[0][l] = solL[0][l] - half*(solL[1][l]*solL[1][l] + solL[2][l]*solL[2][l]
                     +                    solL[3][l]*solL[3][l])/solL[4][l];
          solR[1][l] = zero;
          solR[2][l] = zero;
          solR[3][l] = zero;
          solR[4][l] = solL[4][l];
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
  }
}

//------------------------------------------------------------------------------

// Function, which computes the boundary state (the right state) from the
// given left state and the prescribed boundary data.
void BCInviscidWallSubfaceClass::ComputeBoundaryState(const InputParamClass      *inputParam,
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
                                                      su2double                  &wallPerm)
{
  // Set wallPerm to zero.
  wallPerm = zero;

  // Call ApplyInviscidWallBC to compute the right state.
  ApplyInviscidWallBC(inputParam, nInt, solL, metricL, solR);
}

//------------------------------------------------------------------------------

// Function, which computes the boundary state (the right state) from the
// given left state and the prescribed boundary data.
void BCSymmetrySubfaceClass::ComputeBoundaryState(const InputParamClass      *inputParam,
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
                                                  su2double                  &wallPerm)
{
  Terminate("BCSymmetrySubfaceClass::ComputeBoundaryState",
             __FILE__, __LINE__, "Not implemented yet");
}

//------------------------------------------------------------------------------

// Function, which computes the boundary state (the right state) from the
// given left state and the prescribed boundary data.
void BCInflowSubsonicSubfaceClass::ComputeBoundaryState(const InputParamClass      *inputParam,
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
                                                        su2double                  &wallPerm)
{
  // Some abbreviations involving gamma.
  const su2double gm1    = GamConstant - one;
  const su2double ovgm1  = one/gm1;
  const su2double govgm1 = GamConstant*ovgm1;

  // Make a distinction between the working variables.
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      // Conservative variables are used.
      // Loop over the number of integration points.
#pragma omp simd
      for(int l=0; l<nInt; ++l)
      {
        // Compute the primitive variables from the conservative ones
        // for the left state.
        const su2double rho    = solL[0][l];
        const su2double rhoInv = one/rho;
        const su2double u      = rhoInv*solL[1][l];
        const su2double v      = rhoInv*solL[2][l];
        const su2double w      = rhoInv*solL[3][l];
        const su2double p      = gm1*(solL[4][l] - half*rho*(u*u + v*v + w*w));

        // Compute the components of the unit outward normal.
        const su2double nx = factNorm*metricL[0][l];
        const su2double ny = factNorm*metricL[1][l];
        const su2double nz = factNorm*metricL[2][l];

        // Compute the normal velocity and the speed of sound.
        const su2double un = u*nx + v*ny + w*nz;
        const su2double a2 = GamConstant*p*rhoInv;
        const su2double a  = SQRT(a2);

        // Compute the Riemann invariant to be extrapolated.
        su2double riemann = two*a*ovgm1 + un;

        // Apply total enthalpy scaling to increase stability.
        // If this is not desired, comment the following lines.
        const su2double H = rhoInv*(solL[4][l] + p);
        riemann          *= SQRT(prescribedData[1][l]/H);

        // Compute the dot product between the normal and the
        // velocity direction. This value should be negative.
        const su2double alpha = nx*prescribedData[2][l]
                              + ny*prescribedData[3][l]
                              + nz*prescribedData[4][l];

        // Coefficients in the quadratic equation for
        // the magnitude of the velocity.
        const su2double aa =  one + half*gm1*alpha*alpha;
        const su2double bb = -gm1*alpha*riemann;
        const su2double cc =  half*gm1*riemann*riemann - two*prescribedData[1][l];

        // Solve the equation for the magnitude of the velocity. As this value
        // must be positive and both aa and bb are positive (alpha is negative and
        // riemann is positive up till Mach = 5.0 or so, which is not really subsonic
        // anymore), it is clear which of the two possible solutions must be taken.
        // Some clipping is present, but this is normally not active.
        su2double dd = bb*bb - four*aa*cc;  dd = SQRT(std::max(zero, dd));
        su2double qR = (-bb + dd)/(two*aa); qR = std::max(zero, qR);

        // Compute the square of the Mach number and clip it between 0 and 1,
        // because this is a subsonic inflow boundary.
        su2double qR2 = qR*qR;
        su2double aR2 = gm1*(prescribedData[1][l] - half*qR2);
        su2double MR2 = qR2/aR2; MR2 = std::min(one, MR2);

        // Compute the final value of the magnitude of the velocity.
        const su2double tmp = one/(one + half*gm1*MR2);
        aR2 = gm1*prescribedData[1][l]*tmp;
        qR2 = MR2*aR2;
        qR  = SQRT(qR2);

        // Compute the pressure from the prescribed total pressure
        // and the density from the speed of sound and pressure.
        const su2double pR   = prescribedData[2][l]*POW(tmp, govgm1);
        const su2double rhoR = GamConstant*pR/aR2;

        // Compute the conservative variables of the right state.
        solR[0][l] = rhoR;
        solR[1][l] = rhoR*qR*prescribedData[2][l];
        solR[2][l] = rhoR*qR*prescribedData[3][l];
        solR[3][l] = rhoR*qR*prescribedData[4][l];
        solR[4][l] = pR*ovgm1 + half*rhoR*qR2;
      }

      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      // Entropy variables are used.
      // Loop over the number of integration points.
#pragma omp simd
      for(int l=0; l<nInt; ++l)
      {
        // Compute the velocities and speed of sound for the left state.
        const su2double V4Inv =  one/solL[4][l];
        const su2double u     = -V4Inv*solL[1][l];
        const su2double v     = -V4Inv*solL[2][l];
        const su2double w     = -V4Inv*solL[3][l];
        const su2double a2    = -GamConstant*V4Inv;
        const su2double a     =  SQRT(a2);

        // Compute the components of the unit outward normal.
        const su2double nx = factNorm*metricL[0][l];
        const su2double ny = factNorm*metricL[1][l];
        const su2double nz = factNorm*metricL[2][l];

        // Compute the normal velocity.
        const su2double un = u*nx + v*ny + w*nz;

        // Compute the Riemann invariant to be extrapolated.
        su2double riemann = two*a*ovgm1 + un;

        // Apply total enthalpy scaling to increase stability.
        // If this is not desired, comment the following lines.
        const su2double H = a2*ovgm1 + half*(u*u + v*v + w*w);
        riemann          *= SQRT(prescribedData[1][l]/H);

        // Compute the dot product between the normal and the
        // velocity direction. This value should be negative.
        const su2double alpha = nx*prescribedData[2][l]
                              + ny*prescribedData[3][l]
                              + nz*prescribedData[4][l];

        // Coefficients in the quadratic equation for
        // the magnitude of the velocity.
        const su2double aa =  one + half*gm1*alpha*alpha;
        const su2double bb = -gm1*alpha*riemann;
        const su2double cc =  half*gm1*riemann*riemann - two*prescribedData[1][l];

        // Solve the equation for the magnitude of the velocity. As this value
        // must be positive and both aa and bb are positive (alpha is negative and
        // riemann is positive up till Mach = 5.0 or so, which is not really subsonic
        // anymore), it is clear which of the two possible solutions must be taken.
        // Some clipping is present, but this is normally not active.
        su2double dd = bb*bb - four*aa*cc;  dd = SQRT(std::max(zero, dd));
        su2double qR = (-bb + dd)/(two*aa); qR = std::max(zero, qR);

        // Compute the square of the Mach number and clip it between 0 and 1,
        // because this is a subsonic inflow boundary.
        su2double qR2 = qR*qR;
        su2double aR2 = gm1*(prescribedData[1][l] - half*qR2);
        su2double MR2 = qR2/aR2; MR2 = std::min(one, MR2);

        // Compute the final value of the magnitude of the velocity.
        const su2double tmp = one/(one + half*gm1*MR2);
        aR2 = gm1*prescribedData[1][l]*tmp;
        qR2 = MR2*aR2;
        qR  = SQRT(qR2);

        // Compute the pressure from the prescribed total pressure
        // and the density from the speed of sound and pressure.
        const su2double pR   = prescribedData[0][l]*POW(tmp, govgm1);
        const su2double rhoR = GamConstant*pR/aR2;

        // Compute the entropy variables from the primitive ones.
        const su2double pRInv = one/pR;
        const su2double sR    = LOG(pR/POW(rhoR,GamConstant));

        solR[0][l] =  (GamConstant-sR)*ovgm1 - half*pRInv*rhoR*qR2;
        solR[1][l] =  rhoR*pRInv*qR*prescribedData[2][l];
        solR[2][l] =  rhoR*pRInv*qR*prescribedData[3][l];
        solR[3][l] =  rhoR*pRInv*qR*prescribedData[4][l];
        solR[4][l] = -rhoR*pRInv;
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

// Function, which computes the boundary state (the right state) from the
// given left state and the prescribed boundary data.
void BCInflowSupersonicSubfaceClass::ComputeBoundaryState(const InputParamClass      *inputParam,
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
                                                          su2double                  &wallPerm)
{
  // Supersonic inflow. Set the right state to the prescribed state.
  for(int l=0; l<nVar; ++l)
  {
#pragma omp simd
    for(int m=0; m<nInt; ++m)
      solR[l][m] = prescribedData[l][m];
  }
}

//------------------------------------------------------------------------------

// Function, which computes the boundary state (the right state) from the
// given left state and the prescribed boundary data.
void BCOutflowSubsonicSubfaceClass::ComputeBoundaryState(const InputParamClass      *inputParam,
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
                                                         su2double                  &wallPerm)
{
  // Some abbreviations involving gamma.
  const su2double gm1   = GamConstant - one;
  const su2double ovgm1 = one/gm1;
  const su2double ovg   = one/GamConstant;

  // Make a distinction between the working variables.
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      // Conservative variables are used.
      // Loop over the number of integration points.
#pragma omp simd
      for(int l=0; l<nInt; ++l)
      {
        // Compute the primitive variables from the conservative ones
        // for the left state.
        const su2double rho    = solL[0][l];
        const su2double rhoInv = one/rho;
        const su2double u      = rhoInv*solL[1][l];
        const su2double v      = rhoInv*solL[2][l];
        const su2double w      = rhoInv*solL[3][l];
        const su2double p      = gm1*(solL[4][l] - half*rho*(u*u + v*v + w*w));

        // Compute the components of the unit outward normal.
        const su2double nx = factNorm*metricL[0][l];
        const su2double ny = factNorm*metricL[1][l];
        const su2double nz = factNorm*metricL[2][l];

        // Compute the normal velocity and the speed of sound.
        const su2double un = u*nx + v*ny + w*nz;
        const su2double a2 = GamConstant*p*rhoInv;
        const su2double a  = SQRT(a2);

        // Subsonic exit flow: there is one incoming characteristic,
        // therefore one variable can be specified (back pressure) and is used
        // to compute the variables of the right state. Compute the entropy and
        // the acoustic Riemann variable. These invariants, as well as the
        // tangential velocity components, are extrapolated. Note that s
        // actually contains EXP(entropy). This is to avoid an EXP now and
        // a LOG evaluation later (when the density is computed).
        const su2double s       = p*POW(rhoInv, GamConstant);
        const su2double riemann = two*a*ovgm1 + un;

        // Compute the density and normal velocity of the right state.
        const su2double pR   = prescribedData[0][l];
        const su2double rhoR = POW(pR/s, ovg);
        const su2double aR   = SQRT(GamConstant*pR/rhoR);
        const su2double unR  = riemann - two*aR*ovgm1;

        // Compute the velocity of the right state from the known normal
        // velocity and the extrapolated tangential velocity components.
        const su2double uR = u + (unR - un)*nx;
        const su2double vR = v + (unR - un)*ny;
        const su2double wR = w + (unR - un)*nz;

        // Compute the conservative variables of the right state.
        solR[0][l] = rhoR;
        solR[1][l] = rhoR*uR;
        solR[2][l] = rhoR*vR;
        solR[3][l] = rhoR*wR;
        solR[4][l] = pR*ovgm1 + half*rhoR*(uR*uR + vR*vR + wR*wR);
      }

      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      // Entropy variables are used.
      // Loop over the number of integration points.
#pragma omp simd
      for(int l=0; l<nInt; ++l)
      {
        // Compute the velocities, entropy and speed of sound for the left state.
        const su2double V4Inv =  one/solL[4][l];
        const su2double u     = -V4Inv*solL[1][l];
        const su2double v     = -V4Inv*solL[2][l];
        const su2double w     = -V4Inv*solL[3][l];
        const su2double eKin  =  half*(u*u + v*v + w*w);
        const su2double s     =  GamConstant - gm1*(solL[0][l] - solL[4][l]*eKin);
        const su2double a     =  SQRT(-GamConstant*V4Inv);

        // Compute the components of the unit outward normal.
        const su2double nx = factNorm*metricL[0][l];
        const su2double ny = factNorm*metricL[1][l];
        const su2double nz = factNorm*metricL[2][l];

        // Compute the normal velocity.
        const su2double un = u*nx + v*ny + w*nz;

        // Subsonic exit flow: there is one incoming characteristic,
        // therefore one variable can be specified (back pressure) and is used
        // to compute the variables of the right state. Compute the acoustic
        // Riemann variable, as this variable as well as the entropy and the
        // tangential velocity components, are extrapolated.
        const su2double riemann = two*a*ovgm1 + un;

        // Compute the density and normal velocity of the right state.
        const su2double pR   = prescribedData[0][l];
        const su2double rhoR = POW((pR/EXP(s)), ovg);
        const su2double aR   = SQRT(GamConstant*pR/rhoR);
        const su2double unR  = riemann - two*aR*ovgm1;

        // Compute the velocity of the right state from the known normal
        // velocity and the extrapolated tangential velocity components.
        const su2double uR = u + (unR - un)*nx;
        const su2double vR = v + (unR - un)*ny;
        const su2double wR = w + (unR - un)*nz;

        // Compute the entropy variables from the primitive ones.
        const su2double pRInv = one/pR;

        solR[0][l] =  (GamConstant-s)*ovgm1 - half*pRInv*rhoR*(uR*uR + vR*vR + wR*wR);
        solR[1][l] =  rhoR*pRInv*uR;
        solR[2][l] =  rhoR*pRInv*vR;
        solR[3][l] =  rhoR*pRInv*wR;
        solR[4][l] = -rhoR*pRInv;
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

// Function, which computes the boundary state (the right state) from the
// given left state and the prescribed boundary data.
void BCOutflowSupersonicSubfaceClass::ComputeBoundaryState(const InputParamClass      *inputParam,
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
                                                           su2double                  &wallPerm)
{
  // Supersonic outflow. Set the right state to the left state.
  for(int l=0; l<nVar; ++l)
  {
#pragma omp simd
    for(int m=0; m<nInt; ++m)
      solR[l][m] = solL[l][m];
  }
}




//------------------------------------------------------------------------------
// Characteristic BC: Standard 
//------------------------------------------------------------------------------

// Function, which identifies the normal and transverse metrics
// as well as the index of the normal wave amplitude.
void BCStandardCharacteristicSubfaceClass::IdentifyBoundaryIndices(su2double **metric)
{
	// Reference node.
	const unsigned short iRef = 0;

	// Extract unit normal metrics.
	const su2double nx = metric[0][iRef];
	const su2double ny = metric[1][iRef];
	const su2double nz = metric[2][iRef];

	// Metric term: drdX
	const su2double drdx = metric[4][iRef]; 
	const su2double drdy = metric[5][iRef]; 
	const su2double drdz = metric[6][iRef]; 
	// Metric term: dsdX
	const su2double dsdx = metric[7][iRef]; 
	const su2double dsdy = metric[8][iRef]; 
	const su2double dsdz = metric[9][iRef]; 
	// Metric term: dtdX
	const su2double dtdx = metric[10][iRef];
	const su2double dtdy = metric[11][iRef];
	const su2double dtdz = metric[12][iRef];

	// Dot-product: drdX*n
	const su2double dotprod_nr = nx*drdx + ny*drdy + nz*drdz;
	// Dot-product: dsdX*n
	const su2double dotprod_ns = nx*dsdx + ny*dsdy + nz*dsdz;
	// Dot-product: dtdX*n
	const su2double dotprod_nt = nx*dtdx + ny*dtdy + nz*dtdz;

	// Length of each metric.
	const su2double len_drdX = SQRT(drdx*drdx + drdy*drdy + drdz*drdz);
	const su2double len_dsdX = SQRT(dsdx*dsdx + dsdy*dsdy + dsdz*dsdz);
	const su2double len_dtdX = SQRT(dtdx*dtdx + dtdy*dtdy + dtdz*dtdz);
	
	// Cosine of the angle between the normal and each of the metrics.
	su2double cosTheta[3];
	cosTheta[0] = dotprod_nr/len_drdX; // cos(theta(r))
	cosTheta[1] = dotprod_ns/len_dsdX; // cos(theta(s))
	cosTheta[2] = dotprod_nt/len_dtdX; // cos(theta(t))

	// Offset for metric index.
	unsigned short offset[3];
	offset[0] = 4;  // starting index for r
	offset[1] = 7;  // starting index for s
	offset[2] = 10; // starting index for t

	// Indices found are initialized to false.
	bool idxFound = false;

	// Loop over all three directions and determine the orientation of
	// this surface.
	for(unsigned short i=0; i<3; i++){
	
		// Check if the i-index corresponds to the normal direction.
		if( FABS(cosTheta[i]-one) < epsThreshold ){
	
			// Make sure this is executed once.
			assert(idxFound == false);

			// Specify map of normal index of gradient (i.e d/dr).
			mIndexGradient[0] = i;

			// Make sure other directions are transverse.
			for(unsigned short j=0; j<3; j++){
				if(j != i){
					// Consistency check.
					assert( FABS(cosTheta[j]-one) > epsThreshold );
				}
			}

			// Assign the transverse directions appropriately to ensure that the global
			// reference framework is the same. There are 3 choices, regardless of min/max
			// faces. If normal direction is:  
			//   xi: [0],  then  eta: [1]  and  zeta: [2], 
			//   xi: [1],  then  eta: [0]  and  zeta: [2], 
			//   xi: [2],  then  eta: [0]  and  zeta: [1]. 
			switch( i )
			{
				case( 0 ): { mIndexGradient[1] = 1; mIndexGradient[2] = 2; break; } 
				case( 1 ): { mIndexGradient[1] = 0; mIndexGradient[2] = 2; break; } 
				case( 2 ): { mIndexGradient[1] = 0; mIndexGradient[2] = 1; break; } 
				default: break;
			}

			// Assign indices.
			mIndexMetric[0] = offset[mIndexGradient[0]];
			mIndexMetric[1] = offset[mIndexGradient[1]];
			mIndexMetric[2] = offset[mIndexGradient[2]];

			// Indices have been found.
			idxFound = true;
		}
	}

	// Check whether the indices are all identified.
	if(!idxFound)
		Terminate("BCStandardCharacteristicSubfaceClass::IdentifyBoundaryIndices", __FILE__, __LINE__,
							"Could not identify relevant boundary indices!");
}

//------------------------------------------------------------------------------

// Function, which computes the weighted Mach number on this
// local element boundary that is used to compute the average.
// The average Mach number on this element surface is also computed.
su2double BCStandardCharacteristicSubfaceClass::WeightedMachElement(const InputParamClass *inputParam,
		                                                                const int              nInt,
															                                      su2double            **sol,
															                                      const su2double        factNorm,
															                                      su2double            **metric,
															                                      su2double             *intWeights)
{
	// Initialize weighted Mach number on this element boundary.
	su2double Mwtd = zero;
	// Initialize local surface area.
	su2double area = zero;

	// Extract the starting index of the normal metric.
	const int I0 = mIndexMetric[0], I1 = I0+1, I2 = I1+1;

	// Loop over the number of integration points.
#pragma omp simd
	for(int l=0; l<nInt; ++l)
	{
		// Extract normal length.
		const su2double magn  = metric[3][l];
		// Standard integration weights.
		const su2double wts   = intWeights[l];

		// Extract normal metric components.
		const su2double rrx   = metric[I0][l];
		const su2double rry   = metric[I1][l];
		const su2double rrz   = metric[I2][l];

		// Compute the magnitude of the normal metric.
		const su2double magrr = SQRT( rrx*rrx + rry*rry + rrz*rrz );

		// Compute the velocities.
		const su2double V4Inv =  one/sol[4][l];
		const su2double u     = -V4Inv*sol[1][l];
		const su2double v     = -V4Inv*sol[2][l];
		const su2double w     = -V4Inv*sol[3][l];

		// Compute the local speed of sound.
		const su2double a     = SQRT( -GamConstant*V4Inv );
		// Compute the projected normal velocity.
		const su2double urr   = ( u*rrx + v*rry + w*rrz )/magrr;

		// Compute the local projected Mach number.
		const su2double Mrr   = urr/a;

		// Temporary storage of ||n||*weights.
		const su2double wmagn = wts*magn;

		// Average the local projected Mach number on this element boundary.
		Mwtd += wmagn*Mrr;
		area += wmagn;
	}

	// Compute the current averaged projected Mach number on this element.
	mMachAverageElement = Mwtd/area;

	// Return the weighted projected Mach number on this element boundary.
	return Mwtd;
}

//------------------------------------------------------------------------------

// Function, which converts the entropy variables and their gradient
// into the primitive variables and their gradient. Note, the left 
// state is the interior state and the right state is the external 
// which is overwritten in this function by the primitive data.
void BCStandardCharacteristicSubfaceClass::ConvertEntropyToPrimitive(const int   nInt,
		                                                                 su2double **solL,
			                                                               su2double **dSolDrL,
																                                     su2double **dSolDsL,
																                                     su2double **dSolDtL,
																                                     su2double **solR, 
																                                     su2double **dSolDrR,
																                                     su2double **dSolDsR,
																                                     su2double **dSolDtR)
{
	// Some abbreviations involving gamma.
  const su2double gm1   = GamConstant - one;
  const su2double ovgm1 = one/gm1;

	// Loop over the number of integration points.
#pragma omp simd
	for(int l=0; l<nInt; ++l)
	{
		// Initialize the Jacobian of the transformation.
		su2double dQdV[16];

		// Compute the velocities and speed of sound for the left state.
    const su2double V4Inv  =  one/solL[4][l];
    const su2double u      = -V4Inv*solL[1][l];
    const su2double v      = -V4Inv*solL[2][l];
    const su2double w      = -V4Inv*solL[3][l];
		const su2double ek     =  half*(u*u + v*v + w*w);
		const su2double s      =  GamConstant - gm1*(solL[0][l] - solL[4][l]*ek);
    const su2double a2     = -GamConstant*V4Inv;
		const su2double rho    =  POW(-solL[4][l]*EXP(s), -ovgm1);
		const su2double p      =  a2*rho/GamConstant;
		const su2double povrho = -V4Inv;
		const su2double E      =  ovgm1*povrho + ek;
    const su2double H      =  E + povrho; 

		// Assemble non-zero entries of the Jacobian of transformation: entropy-to-primitive.
		// Note, the commented indices start from index 1, not 0 -- similar to MATLAB.
		// The 1st row.
		dQdV[ 0] = rho;          // dQdV(1,1) 
		dQdV[ 1] = rho*u;        // dQdV(1,2)
		dQdV[ 2] = rho*v;        // dQdV(1,3)
		dQdV[ 3] = rho*w;        // dQdV(1,4)
		dQdV[ 4] = rho*E;        // dQdV(1,5)

		// Then 2nd row.
		dQdV[ 5] = povrho;       // dQdV(2,2)
		dQdV[ 6] = povrho*u;     // dQdV(2,5)

		// The 3rd row.
		dQdV[ 7] = povrho;       // dQdV(3,3)
		dQdV[ 8] = povrho*v;     // dQdV(3,5)

		// The 4th row.
		dQdV[ 9] = povrho;       // dQdV(4,4)
		dQdV[10] = povrho*w;     // dQdV(4,5)

		// The 5th row.
		dQdV[11] = p;            // dQdV(5,1)
		dQdV[12] = p*u;          // dQdV(5,2)
		dQdV[13] = p*v;          // dQdV(5,3)
		dQdV[14] = p*w;          // dQdV(5,4)
		dQdV[15] = p*H;          // dQdV(5,5)

		// Obtain the derivative of the primitive variables from 
		// the derivative of the entropy variables in the three directions.
		// Note, the right state is overwritten by the gradient of the 
		// primitive variables. For readability purposes, loop over all 
		// three directions and assign the relevant pointer depending on
		// the index: i.
		for(int i=0; i<3; i++)
		{
			// Initialize the gradient pointers for the entropy and primitive 
			// variables. Note, the left state is always the internal state 
			// given in entropy variables. The right state is the unknown 
			// state, which for now will be overwritten by the primitive 
			// gradient.
			su2double **dV = NULL, **dQ = NULL;
			switch( i )
			{ 
				case(0): { dV = dSolDrL; dQ = dSolDrR; break; } // r-direction.
				case(1): { dV = dSolDsL; dQ = dSolDsR; break; } // s-direction.
				case(2): { dV = dSolDtL; dQ = dSolDtR; break; } // t-direction.
				default: break;
			}

			// Compute: drho.
  		dQ[0][l] = dQdV[ 0]*dV[0][l] + dQdV[ 1]*dV[1][l] + dQdV[ 2]*dV[2][l] + dQdV[ 3]*dV[3][l] + dQdV[ 4]*dV[4][l];
  		// Compute: du. 
  		dQ[1][l] =                     dQdV[ 5]*dV[1][l]                                         + dQdV[ 6]*dV[4][l];
  		// Compute: dv.
  		dQ[2][l] =                                         dQdV[ 7]*dV[2][l]                     + dQdV[ 8]*dV[4][l];
  		// Compute: dw.
  		dQ[3][l] =                                                             dQdV[ 9]*dV[3][l] + dQdV[10]*dV[4][l];
  		// Compute: dp.
  		dQ[4][l] = dQdV[11]*dV[0][l] + dQdV[12]*dV[1][l] + dQdV[13]*dV[2][l] + dQdV[14]*dV[3][l] + dQdV[15]*dV[4][l];
		}

		// Copy the primitive variables variables.
		solR[0][l] = rho;
		solR[1][l] = u;
		solR[2][l] = v;
		solR[3][l] = w;
		solR[4][l] = p;
	}
}

//------------------------------------------------------------------------------

// Function, which reconstructs a boundary state from the boundary-conforming 
// normal derivative of the primitive variables. This also reconstructs the 
// gradient of that state. The result is a boundary state and gradient in 
// entropy variables. These overwrite the entire right state variables.
void BCStandardCharacteristicSubfaceClass::ReconstructBoundaryStateAndGradient(
		                                       const StandardElementClass *standardHex, 
																					 const su2double             factNorm,
			                                     su2double                 **solR,
			                                     su2double                 **dSolDrR,
																					 su2double                 **dSolDsR,
																					 su2double                 **dSolDtR)
{
	// Some abbreviations involving gamma.
  const su2double gm1    = GamConstant - one;
  const su2double ovgm1  = one/gm1;
  const su2double govgm1 = GamConstant*ovgm1;

	// Extract the relevant number of integration points in 1D and 2D.
  const int nInt1D    = standardHex->mNIntegration1D;
	const int nInt1DPad = standardHex->mNIntegration1DPad;
  const int nInt2D    = standardHex->mNIntegration2D;
  const int nInt2DPad = standardHex->mNIntegration2DPad;

	// Extract the transpose of the 1D Lagrange differentiation function 
	// based on the integration points as basis nodes and applied on the 
	// same integration points as well.
	const su2double *dL1DqT = standardHex->mDerLagrangeInt1D_BaseIntTranspose;

	// Extract the polynomial-correction matrix: A.
	const su2double *A      = standardHex->mPolynomialCorrectionMatrix;

	// Compute the new boundary-conforming primitive variable state 
	// from the B matrix, currently stored in dSolDsR. The result
	// is then appropriately stored in solR.
	for(int k=0; k<nVar; k++)
	{
		// Make sure the padded data in dSolDsR is set to zero to
		// avoid any floating point errors during the multiplication.
		for(int ii=nInt2D; ii<nInt2DPad; ii++) dSolDsR[k][ii] = zero;

		// Perform a (padded) matrix-vector multiplication.
		for(int i=0; i<nInt2D; i++)
		{
			su2double tmp = zero;
#pragma omp simd reduction(+:tmp)
			for(int j=0; j<nInt2DPad; j++)
			{
				tmp += A[i*nInt2DPad+j]*dSolDsR[k][j];
			}
			solR[k][i] = tmp;
		}
	}

	// To avoid numerical issues, set the padded values to the first entry.
	for(int i=0; i<nVar; i++) 
	{
		for(int l=nInt2D; l<nInt2DPad; l++)
		{
			solR[i][l]    = solR[i][0];
			dSolDrR[i][l] = dSolDrR[i][0];
			dSolDsR[i][l] = dSolDsR[i][0];
			dSolDtR[i][l] = dSolDtR[i][0];
		}
	}

	// Multiply the resulting solution by -ve, in case this is a min face. 
	// The reason is because the PC-matrix is based on the max face. Thus,
	// multiplying by -ve gives the correct value.
	if( factNorm < zero )
	{
		for(int i=0; i<nVar; i++)
		{
#pragma omp simd
			for(int l=0; l<nInt2DPad; l++) solR[i][l] *= -one;
		}
	}


	// Convert the variables and their derivative in the normal 
	// direction of the primitive variables into entropy ones.
#pragma omp simd
	for(int l=0; l<nInt2DPad; l++)
	{
		// Initialize the Jacobian of the transformation.
		su2double dVdQ[16];

		// Extract the primitive variables.
		const su2double rho   = solR[0][l];
		const su2double u     = solR[1][l];
		const su2double v     = solR[2][l];
		const su2double w     = solR[3][l];
		const su2double p     = solR[4][l];
		const su2double s     = LOG(p/POW(rho, GamConstant)); 

		// Extract the normal-based derivative of the primitive variables.
		const su2double dr    = dSolDrR[0][l];
		const su2double du    = dSolDrR[1][l];
		const su2double dv    = dSolDrR[2][l];
		const su2double dw    = dSolDrR[3][l];
		const su2double dp    = dSolDrR[4][l];

		// Abbreviations.
		const su2double ovrho = one/rho;
		const su2double ovp   = one/p;
		const su2double rovp  = rho*ovp;
		const su2double ek    = half*( u*u + v*v + w*w );	

		// Assemble the entropy variables.
		solR[0][l] = -rovp*ek + (GamConstant-s)*ovgm1;
		solR[1][l] =  rovp*u;
		solR[2][l] =  rovp*v;
		solR[3][l] =  rovp*w;
		solR[4][l] = -rovp;

		// Assemble non-zero entries of the Jacobian of transformation: primitive-to-entropy.
		// Note, the commented indices start from index 1, not 0 -- similar to MATLAB.
		// The 1st row.
		dVdQ[ 0] =  govgm1*ovrho - ovp*ek;      // dVdQ(1,1)
		dVdQ[ 1] = -solR[1][l];                 // dVdQ(1,2) 
		dVdQ[ 2] = -solR[2][l];                 // dVdQ(1,3)
		dVdQ[ 3] = -solR[3][l];                 // dVdQ(1,4)
		dVdQ[ 4] =  ovp*( rovp*ek - ovgm1 );    // dVdQ(1,5)

		// The 2nd tow.
		dVdQ[ 5] = u*ovp;                       // dVdQ(2,1)
		dVdQ[ 6] =  rovp;                       // dVdQ(2,2)
		dVdQ[ 7] =  -ovp*solR[1][l];            // dVdQ(2,5)

		// The 3rd row.
		dVdQ[ 8] = v*ovp;                       // dVdQ(3,1)
		dVdQ[ 9] =  rovp;                       // dVdQ(3,3)
		dVdQ[10] =  -ovp*solR[2][l];            // dVdQ(3,5)

		// The 4th row.
		dVdQ[11] = w*ovp;                       // dVdQ(4,1)
		dVdQ[12] =  rovp;                       // dVdQ(4,4)
		dVdQ[13] =  -ovp*solR[3][l];            // dVdQ(4,5)

		// The 5th row.
		dVdQ[14] = -ovp;                        // dVdQ(5,1)
		dVdQ[15] =  ovp*rovp;                   // dVdQ(5,5)

		// Overwrite the entries in dSolDrR by the derivative of the entropy variables in
		// the normal (rr-)direction.
		dSolDrR[0][l] = dVdQ[ 0]*dr + dVdQ[ 1]*du + dVdQ[ 2]*dv + dVdQ[ 3]*dw + dVdQ[ 4]*dp;
		dSolDrR[1][l] = dVdQ[ 5]*dr + dVdQ[ 6]*du                             + dVdQ[ 7]*dp;
		dSolDrR[2][l] = dVdQ[ 8]*dr               + dVdQ[ 9]*dv               + dVdQ[10]*dp;
		dSolDrR[3][l] = dVdQ[11]*dr                             + dVdQ[12]*dw + dVdQ[13]*dp;
		dSolDrR[4][l] = dVdQ[14]*dr                                           + dVdQ[15]*dp;
	}

	// Compute the transverse derivate components as applied on the boundary-conforming 
	// entropy variables located in solR.
	for(int iVar=0; iVar<nVar; iVar++)
	{
		int ij = 0;
		for(int i=0; i<nInt1D; i++)
		{
			for(int j=0; j<nInt1D; j++)
			{
				// Reset data.
				dSolDsR[iVar][ij] = dSolDtR[iVar][ij] = zero;

				for(int k=0; k<nInt1D; k++)
				{
					// Select solution indices.
					const int s = k*nInt1D + j;
					const int t = i*nInt1D + k;	

					// Compute the derivative in the s-direction.
					dSolDsR[iVar][ij] += dL1DqT[i*nInt1DPad+k]*solR[iVar][s];
					// Compute the derivative in the t-direction.
					dSolDtR[iVar][ij] += dL1DqT[j*nInt1DPad+k]*solR[iVar][t];
				}
				// Increment global nodal index.
				ij++;
			}
		}
	}
}




//------------------------------------------------------------------------------
// Characteristic BC: Outlet 
//------------------------------------------------------------------------------

// Function, which computes the boundary state (the right state) from the
// given left state and the prescribed boundary data.
void BCOutflowCharacteristicSubfaceClass::ComputeBoundaryState(const InputParamClass      *inputParam,
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
                                                               su2double                  &wallPerm)
{
	// Abbreviation for the incoming wave amplitude index.
	const int iPhi       = mIndexPhi;
	// Abbreviation for the nodal Lagrange coefficient used in the PC-reconstruction.
	const su2double dell = mDerCoefficient;

	// Select type of averaging used in the Mach number calculations.
	const su2double Mach  = ( inputParam->mNSCBC_Outlet_AverageLocal ) ? mMachAverageElement : mMachAverageBoundary; 
	// Compute the square of the boundary-averaged Mach number.
	const su2double Mavg2 = Mach*Mach;
	
	// Compute the normal relaxation coefficient.
	const su2double kk = half*mSigma*(one - Mavg2)/mLengthScale;
	
	// Assign the transverse relaxation coefficients. If either of the terms 
	// is negative, then the relevant coefficient should be overwritten by 
	// the averaged Mach number instead.
	const su2double bl = ( mBeta_l < zero ) ? Mach : mBeta_l;
	const su2double bt = ( mBeta_t < zero ) ? Mach : mBeta_t;

	// Abbreviation involving the transverse terms. Check if the transverse 
	// relaxation coefficients are directly used or if (1-beta) is how they 
	// are factored in. These are clipped at an upper limit of one, although
	// should not be needed.
	const su2double cbl = ( inputParam->mNSCBC_Outlet_OneMinusBeta ) ? std::min(one, one-bl) : std::min(one, bl);
	const su2double cbt = ( inputParam->mNSCBC_Outlet_OneMinusBeta ) ? std::min(one, one-bt) : std::min(one, bt);

	// In the following, the directions (r, s, t) correspond to:
	//  r: is the normal     (xi  -)direction,
	//  s: first  transverse (eta -)direction,
	//  t: second transverse (zeta-)direction.
	//  ... these depend on the normal (r-)direction, specifically:
	//      r = x, then: s = y, t = z (original framework),
	//      r = y, then: s = x, t = z (rotation around  t),
	//      r = z, then: s = y, t = x (rotation around  s).

	// Gradients w.r.t. boundary orientation (internal/known   left state).
	su2double **dSolDrL, **dSolDsL, **dSolDtL;
	// Gradients w.r.t. boundary orientation (external/unknown right state).
	su2double **dSolDrR, **dSolDsR, **dSolDtR;

	// Assign the gradient and metric directions according to the normal 
	// convention on this boundary face.
	for(int i=0; i<3; i++)
	{
		// Temporary gradient pointers.
		su2double **tmpL = NULL, **tmpR = NULL;
		// Step 1: Select old boundary orientation.
    switch( mIndexGradient[i] )
		{
			case(0): { tmpL = dSolDxL; tmpR = dSolDxR; break; }
			case(1): { tmpL = dSolDyL; tmpR = dSolDyR; break; }
			case(2): { tmpL = dSolDzL; tmpR = dSolDzR; break; }
			default: break;
		}

		// Step 2: Map the old orientation to the new boundary-normal orientation.
		switch( i )
		{
			case(0): { dSolDrL = tmpL; dSolDrR = tmpR; break; }
			case(1): { dSolDsL = tmpL; dSolDsR = tmpR; break; }
			case(2): { dSolDtL = tmpL; dSolDtR = tmpR; break; }
			default: break;
		}
	}

	// Starting metric indices on this face for the normal, and the two 
	// transverse directions, respectively.
	const int I0 = mIndexMetric[0], I1 = I0+1, I2 = I1+1; // Normal
	const int J0 = mIndexMetric[1], J1 = J0+1, J2 = J1+1; // Transverse first
	const int K0 = mIndexMetric[2], K1 = K0+1, K2 = K1+1; // Transverse second

	// Convert the entropy variables and their entropy-based gradient 
	// into primitive variables and primitive-based gradient. These 
	// primitive-based data overwrite the right-state variables for now.
  ConvertEntropyToPrimitive(nInt,
			                      solL, dSolDrL, dSolDsL, dSolDtL,
			                      solR, dSolDrR, dSolDsR, dSolDtR);

	// Loop over the number of integration points.
#pragma omp simd
	for(int l=0; l<nInt; ++l)
	{
		// Initialize the wave amplitudes vector.
		su2double LL[5] = {zero, zero, zero, zero, zero}; // normal waves.
		su2double LT[5] = {zero, zero, zero, zero, zero}; // coupled   transverse.
		su2double TT[5] = {zero, zero, zero, zero, zero}; // uncoupled transverse.

		// Extract prescribed pressure.
		const su2double pInf = prescribedData[0][l];

		// Extract the metrics according to the boundary-normal convention. Note, 
		// according to Lodato et al. these are: 
		// [rr]: \nu, [ss]: \tau, [tt]: \kappa.
		const su2double rr[3] = { metricL[I0][l], metricL[I1][l], metricL[I2][l] };
		const su2double ss[3] = { metricL[J0][l], metricL[J1][l], metricL[J2][l] };
		const su2double tt[3] = { metricL[K0][l], metricL[K1][l], metricL[K2][l] };

		// Magnitude of the normal metric vector.
		const su2double magrr   = SQRT( rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2] );
		// Inverse of the magnitude of the normal metric vector.
		const su2double ovmagrr = one/magrr;
		// Normalize the metric in the normal (rr-)direction. This is \bar{nu} in Lodato et al.
		const su2double rn[3]   = { rr[0]*ovmagrr, rr[1]*ovmagrr, rr[2]*ovmagrr };

		// Extract the primitive variables.
		const su2double rho    = solR[0][l];
		const su2double u      = solR[1][l];
		const su2double v      = solR[2][l];
		const su2double w      = solR[3][l];
		const su2double p      = solR[4][l];

		// Extract the derivative of the normal derivative, needed in forming the B (PC-)matrix.
		const su2double drhodr = dSolDrR[0][l];
		const su2double dudr   = dSolDrR[1][l];
		const su2double dvdr   = dSolDrR[2][l];
		const su2double dwdr   = dSolDrR[3][l];
		const su2double dpdr   = dSolDrR[4][l];

		// Inverse of rho.
		const su2double ovrho  = one/rho;

		// Compute the speed of sound and its square.
		const su2double a2     = GamConstant*p*ovrho;
		const su2double a      = SQRT(a2);

		// Compute some abbreviations required.
		const su2double ova2   = one/a2;
		const su2double rhoa   = rho*a;
		const su2double ovrhoa = one/rhoa;
		const su2double rhova  = rho/a;

		// Projected the velocity vector.
		const su2double ur = u*rr[0] + v*rr[1] + w*rr[2]; // w.r.t. rr-direction.
		const su2double us = u*ss[0] + v*ss[1] + w*ss[2]; // w.r.t. ss-direction.
		const su2double ut = u*tt[0] + v*tt[1] + w*tt[2]; // w.r.t. tt-direction.

		// Projected speed of sound w.r.t. rr-direction.
		const su2double arr = a*magrr;

		// Projected speed of sound w.r.t. the coupled rs-directions.
		const su2double ars = a*( rn[0]*ss[0] + rn[1]*ss[1] + rn[2]*ss[2] );
		// Projected speed of sound w.r.t. the coupled rt-directions.
		const su2double art = a*( rn[0]*tt[0] + rn[1]*tt[1] + rn[2]*tt[2] );
		
		// Define the projected (acoustic) eigenvalues. Note, the convective 
		// eigenvalue (vorticity/entropy) is not needed, since it avoids a
		// numerical singularity by skipping division by lmb2 later in the 
		// reconstruction process altogether.
		const su2double lmb1 = ur - arr;
		const su2double lmb3 = ur + arr;

		// Inverse of the acoustic eigenvalues.
		const su2double ovlmb1 = one/lmb1;
		const su2double ovlmb3 = one/lmb3;

		// In the below, the variable indices in dSolDrR, dSolDsR, dSolDtR are:
		//  [0]: rho,    [1]: u,    [2]: v,    [3]: w,    [4]: p.

		// Compute the differential of the entropy: ds = drho - ova2*dp.
		const su2double dsdr = drhodr - ova2*dpdr; 

		// Compute the differential of the projected vorticity w.r.t. rr-direction: 
		//   dwi = e_{ijk}*rn[j]*du[k],    where e_{ijk} is the levi-civita symbol.
		const su2double dw1dr = rn[1]*dwdr - rn[2]*dvdr;
		const su2double dw2dr = rn[2]*dudr - rn[0]*dwdr;
		const su2double dw3dr = rn[0]*dvdr - rn[1]*dudr;

		// The derivative of the rr-projected velocity in the normal rr-direction.
		const su2double durdr = rn[0]*dSolDrR[1][l] 
			                    + rn[1]*dSolDrR[2][l] 
													+ rn[2]*dSolDrR[3][l]; 
		
		// The derivative of the rr-projected velocity in the transverse ss-direction.
		const su2double durds = rn[0]*dSolDsR[1][l] 
			                    + rn[1]*dSolDsR[2][l] 
													+ rn[2]*dSolDsR[3][l]; 
		
		// The derivative of the rr-projected velocity in the transverse tt-direction.
		const su2double durdt = rn[0]*dSolDtR[1][l] 
			                    + rn[1]*dSolDtR[2][l] 
													+ rn[2]*dSolDtR[3][l];

		// Assemble the internal-based (-)normal wave amplitude. 
		// Note, the entropy/vorticity wave amplitudes intentionally
		// skip the multiplication with their eigenvalues in this step.	
		LL[0] = half*lmb1*( ovrhoa*dpdr - durdr );
		LL[1] =           (  rn[0]*dsdr - dw1dr );
		LL[2] =           (  rn[1]*dsdr - dw2dr );
		LL[3] =           (  rn[2]*dsdr - dw3dr );
		LL[4] = half*lmb3*( ovrhoa*dpdr + durdr );

		// Assemble the coupled-transverse wave amplitudes for the acoustic waves only.
		// The acoustic wave relative to (u-a).
		LT[0] = half*ovrhoa*(us - ars)*dSolDsR[4][l] - half*us*durds
			    + half*ovrhoa*(ut - art)*dSolDtR[4][l] - half*ut*durdt;
		// The acoustic wave relative to (u+a).
		LT[4] = half*ovrhoa*(us + ars)*dSolDsR[4][l] + half*us*durds
			    + half*ovrhoa*(ut + art)*dSolDtR[4][l] + half*ut*durdt;
		
		// Assemble the uncoupled-transverse wave amplitudes for the acoustic waves only.
		TT[0] = -half*a*( ss[0]*dSolDsR[1][l] + ss[1]*dSolDsR[2][l] + ss[2]*dSolDsR[3][l]
				  +           tt[0]*dSolDtR[1][l] + tt[1]*dSolDtR[2][l] + tt[2]*dSolDtR[3][l] );
		// Both uncoupled-transverse acoustic wave amplitudes share the same expression.
		TT[4] = TT[0];

		// Correct for the incoming waves. Note, Lt - TT = transverse terms.
		LL[iPhi] = kk*ovrho*(p - pInf) - cbl*LT[iPhi] + cbt*TT[iPhi]; 

		// Reconstruct the boundary-conforming normal derivative of the
		// primitive variables. Note, these are stored in dSolDrR.
		
		// Construct: drhodr.
		dSolDrR[0][l] =  ovlmb1*rhova*LL[0]
			            +         rn[0]*LL[1]
									+         rn[1]*LL[2]
									+         rn[2]*LL[3]
									+  ovlmb3*rhova*LL[4];
		
		// Construct: dudr.
		dSolDrR[1][l] = -ovlmb1*rn[0]*LL[0]
			            -         rn[2]*LL[2]
									+         rn[1]*LL[3]
									+  ovlmb3*rn[0]*LL[4];

		// Construct: dvdr.
		dSolDrR[2][l] = -ovlmb1*rn[1]*LL[0]
			            +         rn[2]*LL[1]
									-         rn[0]*LL[3]
									+  ovlmb3*rn[1]*LL[4];

		// Construct: dwdr.
		dSolDrR[3][l] = -ovlmb1*rn[2]*LL[0]
			            -         rn[1]*LL[1]
									+         rn[0]*LL[2]
									+  ovlmb3*rn[2]*LL[4];

		// Construct: dpdr.
		dSolDrR[4][l] =  ovlmb1*rhoa*LL[0]
			            +  ovlmb3*rhoa*LL[4];

		// Form the B matrix, which is used in the PC-reconstruction process.
		// Note, use dSolDsR as temporary storage for the matrix B. Recall, 
		// dSolDrR contains the boundary-conforming primitive derivative. 
		// Obviously, both derivatives are taken in the normal (rr-)direction.
		dSolDsR[0][l] = dSolDrR[0][l] - ( drhodr - rho*dell ); // density.  
		dSolDsR[1][l] = dSolDrR[1][l] - ( dudr   -   u*dell ); // u-velocity.
		dSolDsR[2][l] = dSolDrR[2][l] - ( dvdr   -   v*dell ); // v-velocity.
		dSolDsR[3][l] = dSolDrR[3][l] - ( dwdr   -   w*dell ); // w-velocity.
		dSolDsR[4][l] = dSolDrR[4][l] - ( dpdr   -   p*dell ); // pressure.
	}

	// Approximate the primitive variable state from the boundary-conforming
	// normal wave amplitude via a PC-approach. Then, convert the primitive 
	// state into the entropy state which are the working variables.
	ReconstructBoundaryStateAndGradient(standardHex, factNorm, solR, dSolDrR, dSolDsR, dSolDtR);
}




//------------------------------------------------------------------------------
// Characteristic BC: Inlet 
//------------------------------------------------------------------------------

// Function, which computes the boundary state (the right state) from the
// given left state and the prescribed boundary data.
void BCInflowCharacteristicSubfaceClass::ComputeBoundaryState(const InputParamClass      *inputParam,
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
                                                              su2double                  &wallPerm)
{
	// Abbreviation for the incoming wave amplitude index.
	const int iPhi       = mIndexPhi;
	// Abbreviation for the nodal Lagrange coefficient used in the PC-reconstruction.
	const su2double dell = mDerCoefficient;
	// Deduce whether the sign used in the incoming acoustic wave coefficient.
	const su2double sign = ( factNorm > zero ) ? -one : one;

	// Select the Mach number based on the global boundary always.
	const su2double Mach  = mMachAverageBoundary; 
	// Compute the square of the boundary-averaged Mach number.
	const su2double Mavg2 = Mach*Mach;

	// Compute the normal relaxation coefficient for the acoustic   wave.
	const su2double ka = sign*half*mSigma*(one - Mavg2)/mLengthScale;
	// Compute the normal relaxation coefficient for the convective waves.
	const su2double ku = mSigma/mLengthScale;

	// In the following, the directions (r, s, t) correspond to:
	//  r: is the normal     (xi  -)direction,
	//  s: first  transverse (eta -)direction,
	//  t: second transverse (zeta-)direction.
	//  ... these depend on the normal (r-)direction, specifically:
	//      r = x, then: s = y, t = z (original framework),
	//      r = y, then: s = x, t = z (rotation around  t),
	//      r = z, then: s = y, t = x (rotation around  s).

	// Gradients w.r.t. boundary orientation (internal/known   left state).
	su2double **dSolDrL, **dSolDsL, **dSolDtL;
	// Gradients w.r.t. boundary orientation (external/unknown right state).
	su2double **dSolDrR, **dSolDsR, **dSolDtR;

	// Assign the gradient and metric directions according to the normal 
	// convention on this boundary face.
	for(int i=0; i<3; i++)
	{
		// Temporary gradient pointers.
		su2double **tmpL = NULL, **tmpR = NULL;
		// Step 1: Select old boundary orientation.
    switch( mIndexGradient[i] )
		{
			case(0): { tmpL = dSolDxL; tmpR = dSolDxR; break; }
			case(1): { tmpL = dSolDyL; tmpR = dSolDyR; break; }
			case(2): { tmpL = dSolDzL; tmpR = dSolDzR; break; }
			default: break;
		}

		// Step 2: Map the old orientation to the new boundary-normal orientation.
		switch( i )
		{
			case(0): { dSolDrL = tmpL; dSolDrR = tmpR; break; }
			case(1): { dSolDsL = tmpL; dSolDsR = tmpR; break; }
			case(2): { dSolDtL = tmpL; dSolDtR = tmpR; break; }
			default: break;
		}
	}

	// Starting metric indices on this face for the normal, and the two 
	// transverse directions, respectively.
	const int I0 = mIndexMetric[0], I1 = I0+1, I2 = I1+1; // Normal
	const int J0 = mIndexMetric[1], J1 = J0+1, J2 = J1+1; // Transverse first
	const int K0 = mIndexMetric[2], K1 = K0+1, K2 = K1+1; // Transverse second

	// Convert the entropy variables and their entropy-based gradient 
	// into primitive variables and primitive-based gradient. These 
	// primitive-based data overwrite the right-state variables for now.
  ConvertEntropyToPrimitive(nInt,
			                      solL, dSolDrL, dSolDsL, dSolDtL,
			                      solR, dSolDrR, dSolDsR, dSolDtR);

	// Loop over the number of integration points.
#pragma omp simd
	for(int l=0; l<nInt; ++l)
	{
		// Initialize the normal wave amplitudes vector.
		su2double LL[5] = {zero, zero, zero, zero, zero};

		// Extract prescribed static conditions: rho, u, v, w.
		const su2double rhoInf = prescribedData[0][l];
		const su2double uInf   = prescribedData[1][l];
		const su2double vInf   = prescribedData[2][l];
		const su2double wInf   = prescribedData[3][l];

		// Extract the metrics according to the boundary-normal convention. Note, 
		// according to Lodato et al. these are: 
		// [rr]: \nu, [ss]: \tau, [tt]: \kappa.
		const su2double rr[3] = { metricL[I0][l], metricL[I1][l], metricL[I2][l] };
		const su2double ss[3] = { metricL[J0][l], metricL[J1][l], metricL[J2][l] };
		const su2double tt[3] = { metricL[K0][l], metricL[K1][l], metricL[K2][l] };

		// Magnitude of the normal metric vector.
		const su2double magrr   = SQRT( rr[0]*rr[0] + rr[1]*rr[1] + rr[2]*rr[2] );
		// Inverse of the magnitude of the normal metric vector.
		const su2double ovmagrr = one/magrr;
		// Normalize the metric in the normal (rr-)direction. This is \bar{nu} in Lodato et al.
		const su2double rn[3]   = { rr[0]*ovmagrr, rr[1]*ovmagrr, rr[2]*ovmagrr };

		// Extract the primitive variables.
		const su2double rho    = solR[0][l];
		const su2double u      = solR[1][l];
		const su2double v      = solR[2][l];
		const su2double w      = solR[3][l];
		const su2double p      = solR[4][l];

		// Extract the derivative of the normal derivative, needed in forming the B (PC-)matrix.
		const su2double drhodr = dSolDrR[0][l];
		const su2double dudr   = dSolDrR[1][l];
		const su2double dvdr   = dSolDrR[2][l];
		const su2double dwdr   = dSolDrR[3][l];
		const su2double dpdr   = dSolDrR[4][l];

		// Inverse of rho.
		const su2double ovrho  = one/rho;

		// Compute the speed of sound and its square.
		const su2double a2     = GamConstant*p*ovrho;
		const su2double a      = SQRT(a2);

		// Compute some abbreviations required.
		const su2double ova2   = one/a2;
		const su2double rhoa   = rho*a;
		const su2double ovrhoa = one/rhoa;
		const su2double rhova  = rho/a;

		// Projected the velocity vector.
		const su2double ur = u*rr[0] + v*rr[1] + w*rr[2]; // w.r.t. rr-direction.
		const su2double us = u*ss[0] + v*ss[1] + w*ss[2]; // w.r.t. ss-direction.
		const su2double ut = u*tt[0] + v*tt[1] + w*tt[2]; // w.r.t. tt-direction.

		// Projected speed of sound w.r.t. rr-direction.
		const su2double arr = a*magrr;

		// Projected speed of sound w.r.t. the coupled rs-directions.
		const su2double ars = a*( rn[0]*ss[0] + rn[1]*ss[1] + rn[2]*ss[2] );
		// Projected speed of sound w.r.t. the coupled rt-directions.
		const su2double art = a*( rn[0]*tt[0] + rn[1]*tt[1] + rn[2]*tt[2] );
		
		// Define the projected eigenvalues: acoustic and convective. 
		const su2double lmb1 = ur - arr;
		const su2double lmb2 = ur;
		const su2double lmb3 = ur + arr;

		// Inverse of the eigenvalues.
		const su2double ovlmb1 = one/lmb1;
		const su2double ovlmb2 = one/lmb2; 
		const su2double ovlmb3 = one/lmb3;

		// In the below, the variable indices in dSolDrR, dSolDsR, dSolDtR are:
		//  [0]: rho,    [1]: u,    [2]: v,    [3]: w,    [4]: p.

		// Compute the differential of the entropy: ds = drho - ova2*dp.
		const su2double dsdr = drhodr - ova2*dpdr; 

		// Compute the differential of the projected vorticity w.r.t. rr-direction: 
		//   dwi = e_{ijk}*rn[j]*du[k],    where e_{ijk} is the levi-civita symbol.
		const su2double dw1dr = rn[1]*dwdr - rn[2]*dvdr;
		const su2double dw2dr = rn[2]*dudr - rn[0]*dwdr;
		const su2double dw3dr = rn[0]*dvdr - rn[1]*dudr;

		// The derivative of the rr-projected velocity in the normal rr-direction.
		const su2double durdr = rn[0]*dSolDrR[1][l] 
			                    + rn[1]*dSolDrR[2][l] 
													+ rn[2]*dSolDrR[3][l]; 
		
		// The derivative of the rr-projected velocity in the transverse ss-direction.
		const su2double durds = rn[0]*dSolDsR[1][l] 
			                    + rn[1]*dSolDsR[2][l] 
													+ rn[2]*dSolDsR[3][l]; 
		
		// The derivative of the rr-projected velocity in the transverse tt-direction.
		const su2double durdt = rn[0]*dSolDtR[1][l] 
			                    + rn[1]*dSolDtR[2][l] 
													+ rn[2]*dSolDtR[3][l];

		// Assemble the internal-based (-)normal wave amplitude. 
		// Note, the entropy/vorticity wave amplitudes need not be constructed
		// since they will always be overwritten by the external part.
		LL[0] = half*lmb1*( ovrhoa*dpdr - durdr );
		LL[4] = half*lmb3*( ovrhoa*dpdr + durdr );

		// Difference between the internal and external/reference data.
		const su2double drho = rho - rhoInf;
		const su2double du   = u   - uInf;
		const su2double dv   = v   - vInf;
		const su2double dw   = w   - wInf;

		// Abbreviation used in the convective relaxation coefficient.
		// Note, for efficiency reasons, include the lmb2 eigenvalue 
		// division in the coefficient.
		const su2double kaol = ku*a*ovlmb2;

		// Correct for the incoming acoustic   wave.
		LL[iPhi] = ka*a*( rn[0]*du + rn[1]*dv + rn[2]*dw   ); 
		// Correct for the incoming convective waves.
		LL[1]    = kaol*( rn[2]*dv - rn[1]*dw - rn[0]*drho );
		LL[2]    = kaol*( rn[0]*dw - rn[2]*du - rn[1]*drho );
		LL[3]    = kaol*( rn[1]*du - rn[0]*dv - rn[2]*drho );


		// Reconstruct the boundary-conforming normal derivative of the
		// primitive variables. Note, these are stored in dSolDrR. Note,
		// there is no need to multiply the convective entries in LL by 
		// the inverse of lmb2, since these are already done so during
		// the assembly of the external LL state.
		
		// Construct: drhodr.
		dSolDrR[0][l] =  ovlmb1*rhova*LL[0]
			            +         rn[0]*LL[1]
									+         rn[1]*LL[2]
									+         rn[2]*LL[3]
									+  ovlmb3*rhova*LL[4];
		
		// Construct: dudr.
		dSolDrR[1][l] = -ovlmb1*rn[0]*LL[0]
			            -         rn[2]*LL[2]
									+         rn[1]*LL[3]
									+  ovlmb3*rn[0]*LL[4];

		// Construct: dvdr.
		dSolDrR[2][l] = -ovlmb1*rn[1]*LL[0]
			            +         rn[2]*LL[1]
									-         rn[0]*LL[3]
									+  ovlmb3*rn[1]*LL[4];

		// Construct: dwdr.
		dSolDrR[3][l] = -ovlmb1*rn[2]*LL[0]
			            -         rn[1]*LL[1]
									+         rn[0]*LL[2]
									+  ovlmb3*rn[2]*LL[4];

		// Construct: dpdr.
		dSolDrR[4][l] =  ovlmb1*rhoa*LL[0]
			            +  ovlmb3*rhoa*LL[4];

		// Form the B matrix, which is used in the PC-reconstruction process.
		// Note, use dSolDsR as temporary storage for the matrix B. Recall, 
		// dSolDrR contains the boundary-conforming primitive derivative. 
		// Obviously, both derivatives are taken in the normal (rr-)direction.
		dSolDsR[0][l] = dSolDrR[0][l] - ( drhodr - rho*dell ); // density.  
		dSolDsR[1][l] = dSolDrR[1][l] - ( dudr   -   u*dell ); // u-velocity.
		dSolDsR[2][l] = dSolDrR[2][l] - ( dvdr   -   v*dell ); // v-velocity.
		dSolDsR[3][l] = dSolDrR[3][l] - ( dwdr   -   w*dell ); // w-velocity.
		dSolDsR[4][l] = dSolDrR[4][l] - ( dpdr   -   p*dell ); // pressure.
	}

	// Approximate the primitive variable state from the boundary-conforming
	// normal wave amplitude via a PC-approach. Then, convert the primitive 
	// state into the entropy state which are the working variables.
	ReconstructBoundaryStateAndGradient(standardHex, factNorm, solR, dSolDrR, dSolDsR, dSolDtR);
}


