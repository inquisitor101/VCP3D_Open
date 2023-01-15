//------------------------------------------------------------------------------
// File, which contains the implementation of the standalone functions to
// compute the fluxes in the integration points of a face. These functions
// are a highly likely candidates for performance optimization and are
// therefore put in a separate file.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"

//------------------------------------------------------------------------------

// Function, which computes the sum of the inviscid and viscous fluxes in the
// integration points of a boundary face. Note, in this convention, the left
// state is always the interior state.
void FluxesBoundaryFace(const InputParamClass                                    *inputParam,
		                    const StandardElementClass                               *standardHex,
                        std::vector<std::vector< std::vector<ElementClass *> > > &elements,
                        const int                                                 nInt,
                        const int                                                 nIntPad,
                        SubfaceBaseClass                                         *BC,
                        const su2double                                           lenScale,
                        const su2double                                           lenScaleLES,
                        const su2double                                          *intWeights,
                        su2double                                               **solL,
                        su2double                                               **dSolDxL,
                        su2double                                               **dSolDyL,
                        su2double                                               **dSolDzL,
                        const su2double                                           factNorm,
                        su2double                                               **metricL,
                        ExchangeDataWallModelClass                               *exchangeDataWM,
                        su2double                                               **prescribedData,
                        su2double                                               **solR,
                        su2double                                               **dSolDxR,
                        su2double                                               **dSolDyR,
                        su2double                                               **dSolDzR,
                        su2double                                                *eddyVis,
                        su2double                                               **fluxTot,
												su2double                                               **dFluxSymDxL,
												su2double                                               **dFluxSymDyL,
												su2double                                               **dFluxSymDzL,
												su2double                                               **dFluxSymDxR,
												su2double                                               **dFluxSymDyR,
												su2double                                               **dFluxSymDzR,
                        const bool                                                ComputeMonitoringData,
                        su2double                                                &EddyVisMax,
                        su2double                                                *forceCoef)
{
  // Determine whether or not a wall model is used.
  const bool useWallModel = BC->UseWallModel();
 	// Determine whether or not a NSCBC boundary is used.
	const bool useNSCBC     =  ( BC->GetTypeBoundaryPrescribed() == BC_OUTFLOW_CHARACTERISTIC
												    || BC->GetTypeBoundaryPrescribed() == BC_INFLOW_CHARACTERISTIC );

  // Initialize the padded values to avoid problems.
  for(int l=0; l<nVar; ++l)
  {
    for(int m=nInt; m<nIntPad; ++m)
    {
      solL[l][m] = solL[l][0];
      dSolDxL[l][m] = dSolDyL[l][m] = dSolDzL[l][m] = zero;
    }
  }

  // At entry the gradients w.r.t. parametric coordinates are stored.
  // Convert them to gradients w.r.t. Cartesian coordinates if no wall
  // model or NSCBC is used.
  if( !useWallModel && !useNSCBC )
  {
    for(int l=0; l<nVar; ++l)
    {
#pragma omp simd
      for(int m=0; m<nIntPad; ++m)
      {
        // Left state. Store the gradients w.r.t. the parametric coordinates.
        su2double dvdr = dSolDxL[l][m], dvds = dSolDyL[l][m], dvdt = dSolDzL[l][m];

        // Compute the Cartesian gradients of the left state.
        dSolDxL[l][m] = dvdr*metricL[4][l] + dvds*metricL[7][l] + dvdt*metricL[10][l];
        dSolDyL[l][m] = dvdr*metricL[5][l] + dvds*metricL[8][l] + dvdt*metricL[11][l];
        dSolDzL[l][m] = dvdr*metricL[6][l] + dvds*metricL[9][l] + dvdt*metricL[12][l]; 
      }
    }
  }

  // Compute the right state, which corresponds to the boundary state.
  bool heatFluxPrescribed = false;
  su2double wallPerm = one, *prescribedWallData = NULL;
	// Note, in case an NSCBC is used, the right-state gradient is computed
	// from the estimated wave amplitude and its PC-reconstruction state. 
	// Recall, the left-state is always the internal state and hence is not 
	// modified.
  BC->ComputeBoundaryState(inputParam, standardHex, nIntPad,  
			                     solL, dSolDxL, dSolDyL, dSolDzL,
                           factNorm, metricL, prescribedData, 
													 solR, dSolDxR, dSolDyR, dSolDzR, 
													 heatFluxPrescribed, prescribedWallData, wallPerm);

  // Determine whether or not the forces must be computed for this face.
  bool ComputeForces = false;
  if(ComputeMonitoringData && inputParam->mMonitorForces && BC->BCIsWall())
    ComputeForces = true;

  // Compute the inviscid fluxes in the integration points of the face.
  // Use dFluxSymDxR as temporary storage.
  InviscidFluxesFace(inputParam, nIntPad, factNorm, metricL, intWeights, solL,
                     solR, wallPerm, fluxTot, ComputeForces, forceCoef, dFluxSymDxR);


	// Set up the gradient of the solution which is needed to compute the viscous numerical fluxes.
	// Note, in case no NSCBC is prescribed, then we already use the left-state gradient which is 
	// always based on the internal state. This also computes the gradients in Cartesian coordinates.
	if( useNSCBC )
	{
    for(int l=0; l<nVar; ++l)
    {
#pragma omp simd
      for(int m=0; m<nIntPad; ++m)
      {
				// Extract metrics, based on the left (internal) state.
				const su2double drdx = metricL[4][l], dsdx = metricL[7][l], dtdx = metricL[10][l],
												drdy = metricL[5][l], dsdy = metricL[8][l], dtdy = metricL[11][l],
												drdz = metricL[6][l], dsdz = metricL[9][l], dtdz = metricL[12][l];

        // Right state. Store the gradients w.r.t. the parametric coordinates.
        const su2double dvdrR = dSolDxR[l][m], dvdsR = dSolDyR[l][m], dvdtR = dSolDzR[l][m];
				// Left state. Store the gradients w.r.t. the parametric coordinates.
        const su2double dvdrL = dSolDxL[l][m], dvdsL = dSolDyL[l][m], dvdtL = dSolDzL[l][m];

        // Compute the Cartesian gradients of the right state.
        dSolDxR[l][m] = dvdrR*drdx + dvdsR*dsdx + dvdtR*dtdx;
        dSolDyR[l][m] = dvdrR*drdy + dvdsR*dsdy + dvdtR*dtdy;
        dSolDzR[l][m] = dvdrR*drdz + dvdsR*dsdz + dvdtR*dtdz; 
        
				// Compute the Cartesian gradients of the left state.
        dSolDxL[l][m] = dvdrL*drdx + dvdsL*dsdx + dvdtL*dtdx;
        dSolDyL[l][m] = dvdrL*drdy + dvdsL*dsdy + dvdtL*dtdy;
        dSolDzL[l][m] = dvdrL*drdz + dvdsL*dsdz + dvdtL*dtdz; 
      
				// Compute the average of the internal and external gradients. Recall, the left-state is 
				// always the gradient of the internal solution and the right-state is the reconstructed
				// NSCBC-compliant gradient. Note, these are in Cartesian coordinates, as they should be.
        dSolDxL[l][m] = half*( dSolDxL[l][m] + dSolDxR[l][m] );
        dSolDyL[l][m] = half*( dSolDyL[l][m] + dSolDyR[l][m] );
        dSolDzL[l][m] = half*( dSolDzL[l][m] + dSolDzR[l][m] ); 
			}
    }
	}

  // Check if a wall model is used to compute the viscous fluxes.
  if( useWallModel )
  {
    // Wall model is used. Compute the viscous shear stress and heat flux
    // from the wall model and compute all the viscous fluxes, i.e. also
    // penalty and symmetrizing fluxes, in the integration points of the face.
    // The symmetrizing fluxes are stored in dFluxSymDxR, dFluxSymDyR and dFluxSymDzR.
    ViscousFluxesFaceWallModel(exchangeDataWM, elements, inputParam, nIntPad, factNorm,
                               metricL, intWeights, lenScale, lenScaleLES, 
															 solL, solR, 
															 heatFluxPrescribed, prescribedWallData, fluxTot,
                               dFluxSymDxR, dFluxSymDyR, dFluxSymDzR, 
															 eddyVis, ComputeForces, forceCoef);
  }
  else
  {
    // No wall model is used. Compute all the viscous fluxes, i.e. also
    // penalty and symmetrizing fluxes, in the integration points of the face.
    // The symmetrizing fluxes are stored in dFluxSymDxR, dFluxSymDyR and dFluxSymDzR.
    ViscousFluxesFace(inputParam, nIntPad, factNorm, metricL, 
				              intWeights, lenScale, lenScaleLES, 
											solL, solR, 
											dSolDxL, dSolDyL, dSolDzL,
                      heatFluxPrescribed, prescribedWallData, fluxTot, 
											dFluxSymDxR, dFluxSymDyR, dFluxSymDzR, 
											eddyVis, ComputeForces, forceCoef);
  }

  // The symmetrizing fluxes are multiplied by the gradients of the basis
  // functions for the residual. These must be the Cartesian gradients,
  // but for the standard element the gradients w.r.t. the parametric
  // coordinates are stored. To account for this, the Cartesian symmetrizing
  // fluxes are converted parametric fluxes. Note that this function is called
  // for both min and max boundaries, hence store symmetrizing fluxes for both
  // the left and right state, although they are the same. Therefore the final
  // symmetrizing fluxes are stored in 
	// dFluxSymDxL, dFluxSymDyL and dFluxSymDzL for the left state, and in
	// dFluxSymDxR, dFluxSymDyR and dFluxSymDzR for the right state.
  for(int l=0; l<nVar; ++l)
  {
#pragma omp simd
    for(int m=0; m<nIntPad; ++m)
    {
      // Store the Cartesian fluxes.
      const su2double fx = dFluxSymDxR[l][m], fy = dFluxSymDyR[l][m], fz = dFluxSymDzR[l][m];

      // Compute the symmetrizing fluxes for the left element.
      dFluxSymDxL[l][m] = fx*metricL[ 4][l]  + fy*metricL[ 5][l]  + fz*metricL[ 6][l];
      dFluxSymDyL[l][m] = fx*metricL[ 7][l]  + fy*metricL[ 8][l]  + fz*metricL[ 9][l];
      dFluxSymDzL[l][m] = fx*metricL[10][l]  + fy*metricL[11][l]  + fz*metricL[12][l];

      // Compute the symmetrizing fluxes for the right element.
      dFluxSymDxR[l][m] = dFluxSymDxL[l][m];
      dFluxSymDyR[l][m] = dFluxSymDyL[l][m];
      dFluxSymDzR[l][m] = dFluxSymDzL[l][m];
    }
  }

  // Update the maximum of the eddy viscosity, if the monitoring data
  // must be computed.
  if( ComputeMonitoringData )
  {
    for(int m=0; m<nIntPad; ++m)
      EddyVisMax = std::max(EddyVisMax, eddyVis[m]);
  }
}

//------------------------------------------------------------------------------

// Function, which computes the sum of the inviscid and viscous fluxes in the
// integration points of an internal face.
void FluxesInternalFace(const InputParamClass *inputParam,
                        const int              nInt,
                        const int              nIntPad,
                        const su2double        lenScale,
                        const su2double        lenScaleLES,
                        const su2double       *intWeights,
                        su2double            **solL,
                        su2double            **dSolDxL,
                        su2double            **dSolDyL,
                        su2double            **dSolDzL,
                        su2double            **metricL,
                        su2double            **solR,
                        su2double            **dSolDxR,
                        su2double            **dSolDyR,
                        su2double            **dSolDzR,
                        su2double            **metricR,
                        su2double             *eddyVis,
                        su2double            **fluxTot,
												su2double            **dFluxSymDxL,
												su2double            **dFluxSymDyL,
												su2double            **dFluxSymDzL,
												su2double            **dFluxSymDxR,
												su2double            **dFluxSymDyR,
												su2double            **dFluxSymDzR,
                        const bool             ComputeMonitoringData,
                        su2double             &EddyVisMax)
{
  // Initialize the padded values to avoid problems.
  for(int l=0; l<nVar; ++l)
  {
    for(int m=nInt; m<nIntPad; ++m)
    {
      solL[l][m] = solL[l][0];
      solR[l][m] = solR[l][0];

      dSolDxL[l][m] = dSolDyL[l][m] = dSolDzL[l][m] = zero;
      dSolDxR[l][m] = dSolDyR[l][m] = dSolDzR[l][m] = zero;
    }
  }

  // At entry the gradients w.r.t. parametric coordinates are stored.
  // Convert them to gradients w.r.t. Cartesian coordinates and average
  // them afterwards.
  for(int l=0; l<nVar; ++l)
  {
#pragma omp simd
    for(int m=0; m<nIntPad; ++m)
    {
      // Left state. Store the gradients w.r.t. the parametric coordinates.
      su2double dvdr = dSolDxL[l][m], dvds = dSolDyL[l][m], dvdt = dSolDzL[l][m];

      // Compute the Cartesian gradients of the left state.
      dSolDxL[l][m] = dvdr*metricL[4][l] + dvds*metricL[7][l] + dvdt*metricL[10][l];
      dSolDyL[l][m] = dvdr*metricL[5][l] + dvds*metricL[8][l] + dvdt*metricL[11][l];
      dSolDzL[l][m] = dvdr*metricL[6][l] + dvds*metricL[9][l] + dvdt*metricL[12][l]; 

      // Right state. Store the gradients w.r.t. the parametric coordinates.
      dvdr = dSolDxR[l][m]; dvds = dSolDyR[l][m]; dvdt = dSolDzR[l][m];

      // Compute the Cartesian gradients of the right state.
      dSolDxR[l][m] = dvdr*metricR[4][l] + dvds*metricR[7][l] + dvdt*metricR[10][l];
      dSolDyR[l][m] = dvdr*metricR[5][l] + dvds*metricR[8][l] + dvdt*metricR[11][l];
      dSolDzR[l][m] = dvdr*metricR[6][l] + dvds*metricR[9][l] + dvdt*metricR[12][l];

      // Store the average gradients in dSolDxL, dSolDyL and dSolDzL.
      dSolDxL[l][m] = half*(dSolDxL[l][m] + dSolDxR[l][m]);
      dSolDyL[l][m] = half*(dSolDyL[l][m] + dSolDyR[l][m]);
      dSolDzL[l][m] = half*(dSolDzL[l][m] + dSolDzR[l][m]);
    }
  }

  // Compute the inviscid fluxes in the integration points of the face.
  // Use dFluxSymDxR as temporary storage.
  InviscidFluxesFace(inputParam, nIntPad, one, metricL, intWeights, solL,
                     solR, one, fluxTot, false, NULL, dFluxSymDxR);

  // Compute all the viscous fluxes, i.e. also penalty and symmetrizing
  // fluxes, in the integration points of the face. The symmetrizing
  // fluxes are stored in dFluxSymDxR, dFluxSymDyR and dFluxSymDzR.
  ViscousFluxesFace(inputParam, nIntPad, one, metricL, 
			              intWeights, lenScale, lenScaleLES, 
										solL, solR, 
										dSolDxL, dSolDyL, dSolDzL,
                    false, NULL, fluxTot, 
										dFluxSymDxR, dFluxSymDyR, dFluxSymDzR,
                    eddyVis, false, NULL);

  // The symmetrizing fluxes are multiplied by the gradients of the basis
  // functions for the residual. These must be the Cartesian gradients,
  // but for the standard element the gradients w.r.t. the parametric
  // coordinates are stored. To account for this, the Cartesian symmetrizing
  // fluxes are converted parametric fluxes. Note that this conversion is
  // usually different for the left and the right element. Therefore the final
  // symmetrizing fluxes for the left element are stored in 
	// dFluxSymDxL, dFluxSymDyL and dFluxSymDzL and in
  // dFluxSymDxR, dFluxSymDyR and dFluxSymDzR in the right element.
  for(int l=0; l<nVar; ++l)
  {
#pragma omp simd
    for(int m=0; m<nIntPad; ++m)
    {
      // Store the Cartesian fluxes.
      const su2double fx = dFluxSymDxR[l][m], fy = dFluxSymDyR[l][m], fz = dFluxSymDzR[l][m];

      // Compute the symmetrizing fluxes for the left element.
      dFluxSymDxL[l][m] = fx*metricL[4][l]  + fy*metricL[5][l]  + fz*metricL[6][l];
      dFluxSymDyL[l][m] = fx*metricL[7][l]  + fy*metricL[8][l]  + fz*metricL[9][l];
      dFluxSymDzL[l][m] = fx*metricL[10][l] + fy*metricL[11][l] + fz*metricL[12][l];

      // Compute the symmetrizing fluxes for the right element.
      dFluxSymDxR[l][m] = fx*metricR[4][l]  + fy*metricR[5][l]  + fz*metricR[6][l];
      dFluxSymDyR[l][m] = fx*metricR[7][l]  + fy*metricR[8][l]  + fz*metricR[9][l];
      dFluxSymDzR[l][m] = fx*metricR[10][l] + fy*metricR[11][l] + fz*metricR[12][l];
    }
  }

  // Update the maximum of the eddy viscosity, if the monitoring data
  // must be computed.
  if( ComputeMonitoringData )
  {
    for(int m=0; m<nIntPad; ++m)
      EddyVisMax = std::max(EddyVisMax, eddyVis[m]);
  }
}

//------------------------------------------------------------------------------

// Function, which computes the inviscid fluxes in the integration points
// of a face for the given left and right solution.
void InviscidFluxesFace(const InputParamClass *inputParam,
                        const int              nIntPad,
                        const su2double        factNorm,
                        su2double            **metric,
                        const su2double       *intWeights,
                        su2double            **solL,
                        su2double            **solR,
                        const su2double        factUNorm,
                        su2double            **invFlux,
                        const bool             ComputeForces,
                        su2double             *forceCoef,
                        su2double            **pInt)
{
  // Values to scale the acoustic eigenvalues to obtain an adequate amount
  // of dissipation to be entropy satisfying in the Ismail_Roe flux.
  // Note that only alphaMax is taken, assuming that the jump in Mach number
  // over the interface is less than 0.5. For alphaMax = 0 the EC1 flux is obtained.
  const su2double beta     = one/six;
//const su2double alphaMax = two;
  const su2double alphaMax = zero;

  // Value of the entropy fix for standard Roe scheme.
  const su2double Delta = (su2double) 0.001;

  // Expressions involving gamma.
  const su2double gam    =  GamConstant;
  const su2double gm1    =  gam-one;
  const su2double ovgm1  =  one/(gam-one);
  const su2double ov1mg  = -ovgm1;
  const su2double gp1Ovg =  (gam+one)/gam;
  const su2double gm1Ovg =  gm1/gam;

  // Make a distinction between the working variables.
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      // Conservative variables are used. Determine the Riemann solver.
      switch( inputParam->mRiemannSolver )
      {
        case ROE:
        {
          // Standard Roe scheme.
          // Loop over the (padded) number of integration points.
#pragma omp simd
          for(int l=0; l<nIntPad; ++l)
          {
            // Determine the coefficient multiplying the fluxes and compute
            // the components of the unit normal that points from the left
            // to the right state.
            const su2double halfArea = half*intWeights[l]*metric[3][l];
            const su2double nx = factNorm*metric[0][l],
                            ny = factNorm*metric[1][l],
                            nz = factNorm*metric[2][l];

            // Compute the primitive variables of the left and right state.
            su2double tmp       = one/solL[0][l];
            const su2double vxL = tmp*solL[1][l];
            const su2double vyL = tmp*solL[2][l];
            const su2double vzL = tmp*solL[3][l];
            const su2double pL  = gm1*(solL[4][l]
                                - half*(vxL*solL[1][l] + vyL*solL[2][l] + vzL*solL[3][l]));

            tmp                 = one/solR[0][l];
            const su2double vxR = tmp*solR[1][l];
            const su2double vyR = tmp*solR[2][l];
            const su2double vzR = tmp*solR[3][l];
            const su2double pR  = gm1*(solR[4][l]
                                - half*(vxR*solR[1][l] + vyR*solR[2][l] + vzR*solR[3][l]));

            // Store the left and right pressure in pInt. It will be used
            // later on for computing the forces, if desired.
            pInt[0][l] = pL;
            pInt[1][l] = pR;

            // Compute the difference of the conservative mean flow variables.
            const su2double dr  = solR[0][l] - solL[0][l];
            const su2double dru = solR[1][l] - solL[1][l];
            const su2double drv = solR[2][l] - solL[2][l];
            const su2double drw = solR[3][l] - solL[3][l];
            const su2double drE = solR[4][l] - solL[4][l];

            // Compute the Roe average state.
            const su2double zL = SQRT(solL[0][l]);
            const su2double zR = SQRT(solR[0][l]);
            tmp                = one/(zL + zR);

            const su2double rHL = solL[4][l] + pL;
            const su2double rHR = solR[4][l] + pR;

            const su2double uAvg = tmp*(zL*vxL + zR*vxR);
            const su2double vAvg = tmp*(zL*vyL + zR*vyR);
            const su2double wAvg = tmp*(zL*vzL + zR*vzR);
            const su2double HAvg = tmp*(rHL/zL + rHR/zR);

            // Compute from the Roe average state some variables, which occur
            // quite often in the matrix vector product to be computed.
            const su2double alphaAvg = half*(uAvg*uAvg + vAvg*vAvg + wAvg*wAvg);
            tmp                      = gm1*(HAvg - alphaAvg);
            const su2double a2Avg    = FABS(tmp);
            const su2double aAvg     = SQRT(a2Avg);
            const su2double vnAvg    = uAvg*nx + vAvg*ny + wAvg*nz;
            const su2double ovaAvg   = one/aAvg;
            const su2double ova2Avg  = one/a2Avg;

            // Compute the absolute values of the three eigenvalues and
            // apply the entropy correction.
            su2double lam1 = FABS(vnAvg + aAvg);
            su2double lam2 = FABS(vnAvg - aAvg);
            su2double lam3 = FABS(vnAvg);

            tmp  = Delta*std::max(lam1, lam2);
            lam1 = std::max(lam1, tmp);
            lam2 = std::max(lam2, tmp);
            lam3 = std::max(lam3, tmp);

            // Some abbreviations, which occur quite often in the dissipation terms.
            const su2double abv1 = half*(lam1 + lam2);
            const su2double abv2 = half*(lam1 - lam2);
            const su2double abv3 = abv1 - lam3;

            const su2double abv4 = gm1*(alphaAvg*dr - uAvg*dru - vAvg*drv -wAvg*drw + drE);
            const su2double abv5 = nx*dru + ny*drv + nz*drw - vnAvg*dr;
            const su2double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
            const su2double abv7 = abv2*abv4*ovaAvg  + abv3*abv5;

            // Compute the Roe flux vector, which is 0.5*(FL + FR - |A|(UR-UL)).
            const su2double vnL = vxL*nx + vyL*ny + vzL*nz;
            const su2double vnR = vxR*nx + vyR*ny + vzR*nz;
            const su2double pa  = halfArea*(pL + pR);

            const su2double fact = factUNorm*halfArea;

            invFlux[0][l] = fact*(solL[0][l]*vnL + solR[0][l]*vnR - (lam3*dr + abv6));
            invFlux[1][l] = fact*(solL[1][l]*vnL + solR[1][l]*vnR
                          -       (lam3*dru + uAvg*abv6 + nx*abv7)) + pa*nx;
            invFlux[2][l] = fact*(solL[2][l]*vnL + solR[2][l]*vnR
                          -       (lam3*drv + vAvg*abv6 + ny*abv7)) + pa*ny;
            invFlux[3][l] = fact*(solL[3][l]*vnL + solR[3][l]*vnR
                          -       (lam3*drw + wAvg*abv6 + nz*abv7)) + pa*nz;
            invFlux[4][l] = fact*(solL[4][l]*vnL + solR[4][l]*vnR + pL*vnL + pR*vnR
                          -       (lam3*drE + HAvg*abv6 + vnAvg*abv7));
          }
       
          break;
        }

        //----------------------------------------------------------------------

        case ISMAIL_ROE:
        {
          // Entropy satisfying Riemann flux of Ismail and Roe.
          // Loop over the (padded) number of integration points.
#pragma omp simd
          for(int l=0; l<nIntPad; ++l)
          {
            // Determine the coefficient multiplying the fluxes and compute
            // the components of the unit normal that points from the left
            // to the right state.
            const su2double Area = intWeights[l]*metric[3][l];
            const su2double nx = factNorm*metric[0][l],
                            ny = factNorm*metric[1][l],
                            nz = factNorm*metric[2][l];

            // Compute the primitive variables of the left and right state.
            su2double tmp       = one/solL[0][l];
            const su2double vxL = tmp*solL[1][l];
            const su2double vyL = tmp*solL[2][l];
            const su2double vzL = tmp*solL[3][l];
            const su2double pL  = gm1*(solL[4][l]
                                 - half*(vxL*solL[1][l] + vyL*solL[2][l] + vzL*solL[3][l]));

            tmp                 = one/solR[0][l];
            const su2double vxR = tmp*solR[1][l];
            const su2double vyR = tmp*solR[2][l];
            const su2double vzR = tmp*solR[3][l];
            const su2double pR  = gm1*(solR[4][l]
                                - half*(vxR*solR[1][l] + vyR*solR[2][l] + vzR*solR[3][l]));

            const su2double rhoPInvL = solL[0][l]/pL;
            const su2double rhoPInvR = solR[0][l]/pR;

            // Store the left and right pressure in pInt. It will be used
            // later on for computing the forces, if desired.
            pInt[0][l] = pL;
            pInt[1][l] = pR;

            // Compute the entropy variables of the left and right state.
            tmp = LOG(pL/POW(solL[0][l],GamConstant));

            const su2double V0L =  (GamConstant-tmp)*ovgm1
                                -  half*rhoPInvL*(vxL*vxL + vyL*vyL + vzL*vzL);
            const su2double V1L =  rhoPInvL*vxL;
            const su2double V2L =  rhoPInvL*vyL;
            const su2double V3L =  rhoPInvL*vzL;
            const su2double V4L = -rhoPInvL;

            tmp = LOG(pR/POW(solR[0][l],GamConstant));

            const su2double V0R =  (GamConstant-tmp)*ovgm1
                                -  half*rhoPInvR*(vxR*vxR + vyR*vyR + vzR*vzR);
            const su2double V1R =  rhoPInvR*vxR;
            const su2double V2R =  rhoPInvR*vyR;
            const su2double V3R =  rhoPInvR*vzR;
            const su2double V4R = -rhoPInvR;

            // Compute the difference in entropy variables.
            const su2double dV0 = V0R - V0L, dV1 = V1R - V1L, dV2 = V2R - V2L,
                            dV3 = V3R - V3L, dV4 = V4R - V4L;

            // Compute the z-variables of the left and right states.
            const su2double z0L = SQRT(rhoPInvL),      z0R = SQRT(rhoPInvR);
            const su2double z1L = vxL*z0L,             z1R = vxR*z0R;
            const su2double z2L = vyL*z0L,             z2R = vyR*z0R;
            const su2double z3L = vzL*z0L,             z3R = vzR*z0R;
            const su2double z4L = SQRT(solL[0][l]*pL), z4R = SQRT(solR[0][l]*pR);

            // Compute the arithmetic average of the z-variables.
            const su2double z0Avg = half*(z0L + z0R);
            const su2double z1Avg = half*(z1L + z1R);
            const su2double z2Avg = half*(z2L + z2R);
            const su2double z3Avg = half*(z3L + z3R);
            const su2double z4Avg = half*(z4L + z4R);

            // Compute the logarithmic mean of z0.
            su2double zeta, f, u, F;

            zeta = z0L/z0R;
            f    = (zeta-one)/(zeta+one);
            u    = f*f;
            if(u < (su2double) 0.01)
              F = one + u/three + u*u/five + u*u*u/seven + u*u*u*u/nine;
            else
              F = LOG(zeta)/(two*f);

            const su2double z0LogAvg = z0Avg/F;

            // Compute the logarithmic mean of z4.
            zeta = z4L/z4R;
            f    = (zeta-one)/(zeta+one);
            u    = f*f;
            if(u < (su2double) 0.01)
              F = one + u/three + u*u/five + u*u*u/seven + u*u*u*u/nine;
            else
              F = LOG(zeta)/(two*f);

            const su2double z4LogAvg = z4Avg/F;

            // Compute the other averaged quantities that are necessary.
            const su2double oneOvz0Avg = one/z0Avg;
            const su2double rhoAvg = z0Avg*z4LogAvg;
            const su2double p1Avg  = oneOvz0Avg*z4Avg;
            const su2double p2Avg  = half*(gp1Ovg*z4LogAvg/z0LogAvg + gm1Ovg*p1Avg);
            const su2double vxAvg  = oneOvz0Avg*z1Avg;
            const su2double vyAvg  = oneOvz0Avg*z2Avg;
            const su2double vzAvg  = oneOvz0Avg*z3Avg;

            const su2double vnAvg  = vxAvg*nx + vyAvg*ny + vzAvg*nz;
            const su2double kinAvg = half*(vxAvg*vxAvg + vyAvg*vyAvg + vzAvg*vzAvg);
            const su2double a2Avg  = gam*p2Avg/rhoAvg;
            const su2double aAvg   = SQRT(a2Avg);
            const su2double HAvg   = a2Avg*ovgm1 + kinAvg;
            const su2double EAvg   = HAvg - p2Avg/rhoAvg;

            const su2double ovaAvg  = one/aAvg;
            const su2double ova2Avg = one/a2Avg;

            // Define the difference in conservative variables as dU/dV deltaV, where
            // the transformation matrix dU/dV must be evaluated at the averaged state.
            su2double abv1 = rhoAvg*(vxAvg*dV1 + vyAvg*dV2 + vzAvg*dV3);
            su2double abv2 = abv1 + rhoAvg*(dV0 + HAvg*dV4);

            const su2double dr  = abv1 + rhoAvg*(dV0 + EAvg*dV4);
            const su2double dru = vxAvg*abv2 + p2Avg*dV1;
            const su2double drv = vyAvg*abv2 + p2Avg*dV2;
            const su2double drw = vzAvg*abv2 + p2Avg*dV3;
            const su2double drE = HAvg*abv1 + rhoAvg*EAvg*(dV0 + HAvg*dV4)
                                + p2Avg*kinAvg*dV4;

            // Compute the absolute values of the eigenvalues of the flux Jacobian.
            su2double lam1 = FABS(vnAvg + aAvg);
            su2double lam2 = FABS(vnAvg - aAvg);
            su2double lam3 = FABS(vnAvg);

            // Scale the acoustic eigenvalue, such that the EC2 (or EC1) flux of Ismail
            // and Roe is obtained. Also multiply all eigenvalues by half to obtain
            // the correct scaling.
            lam1 *= half*(one + beta + alphaMax);
            lam2 *= half*(one + beta + alphaMax);
            lam3 *= half;

            // Some abbreviations, which occur quite often in the dissipation terms.
            abv1 = half*(lam1 + lam2);
            abv2 = half*(lam1 - lam2);

            const su2double abv3 = abv1 - lam3;
            const su2double abv4 = gm1*(kinAvg*dr - vxAvg*dru - vyAvg*drv - vzAvg*drw + drE);
            const su2double abv5 = nx*dru + ny*drv + nz*drw - vnAvg*dr;
            const su2double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
            const su2double abv7 = abv2*abv4*ovaAvg  + abv3*abv5;

            // Compute the central flux.
            invFlux[0][l] = factUNorm*rhoAvg*vnAvg;
            invFlux[1][l] = factUNorm*rhoAvg*vnAvg*vxAvg + p1Avg*nx;
            invFlux[2][l] = factUNorm*rhoAvg*vnAvg*vyAvg + p1Avg*ny;
            invFlux[3][l] = factUNorm*rhoAvg*vnAvg*vzAvg + p1Avg*nz;
            invFlux[4][l] = factUNorm*rhoAvg*vnAvg*HAvg;

            // Subtract the dissipation terms.
            invFlux[0][l] -= factUNorm*(lam3*dr + abv6);
            invFlux[1][l] -= factUNorm*(lam3*dru + vxAvg*abv6 + nx*abv7);
            invFlux[2][l] -= factUNorm*(lam3*drv + vyAvg*abv6 + ny*abv7);
            invFlux[3][l] -= factUNorm*(lam3*drw + vzAvg*abv6 + nz*abv7);
            invFlux[4][l] -= factUNorm*(lam3*drE + HAvg*abv6 + vnAvg*abv7);

            // Scale the fluxes with Area.
            invFlux[0][l] *= Area;
            invFlux[1][l] *= Area;
            invFlux[2][l] *= Area;
            invFlux[3][l] *= Area;
            invFlux[4][l] *= Area;
          }

          break;
        }

        //----------------------------------------------------------------------

        default:
        {
          // This is just to avoid a compiler warning.
          break;
        }
      }

      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      // Entropy variables are used. Determine the Riemann solver.
      switch( inputParam->mRiemannSolver )
      {
        case ROE:
        {
          // Standard Roe scheme.
          // Loop over the (padded) number of integration points.
#pragma omp simd
          for(int l=0; l<nIntPad; ++l)
          {
            // Determine the coefficient multiplying the fluxes and compute
            // the components of the unit normal that points from the left
            // to the right state.
            const su2double halfArea = half*intWeights[l]*metric[3][l];
            const su2double nx = factNorm*metric[0][l],
                            ny = factNorm*metric[1][l],
                            nz = factNorm*metric[2][l];

            // Compute the primitive variables of the left and right state.
            const su2double V4InvL =  one/solL[4][l];
            const su2double vxL    = -V4InvL*solL[1][l];
            const su2double vyL    = -V4InvL*solL[2][l];
            const su2double vzL    = -V4InvL*solL[3][l]; 
            const su2double vTot2L =  vxL*vxL + vyL*vyL + vzL*vzL;
            const su2double sL     =  gam - gm1*(solL[0][l] - half*solL[4][l]*vTot2L);
            const su2double tmpL   = -solL[4][l]*EXP(sL);
            const su2double rhoL   =  POW(tmpL, ov1mg);
            const su2double pL     = -V4InvL*rhoL;

            const su2double V4InvR =  one/solR[4][l];
            const su2double vxR    = -V4InvR*solR[1][l];
            const su2double vyR    = -V4InvR*solR[2][l];
            const su2double vzR    = -V4InvR*solR[3][l]; 
            const su2double vTot2R =  vxR*vxR + vyR*vyR + vzR*vzR;
            const su2double sR     =  gam - gm1*(solR[0][l] - half*solR[4][l]*vTot2R);
            const su2double tmpR   = -solR[4][l]*EXP(sR);
            const su2double rhoR   =  POW(tmpR, ov1mg);
            const su2double pR     = -V4InvR*rhoR;

            // Store the left and right pressure in pInt. It will be used
            // later on for computing the forces, if desired.
            pInt[0][l] = pL;
            pInt[1][l] = pR;

            // Compute the total energy for the left and right state.
            const su2double rEL = pL*ovgm1 + half*rhoL*vTot2L;
            const su2double rER = pR*ovgm1 + half*rhoR*vTot2R;

            // Compute the difference of the conservative mean flow variables.
            const su2double dr  = rhoR     - rhoL;
            const su2double dru = rhoR*vxR - rhoL*vxL;
            const su2double drv = rhoR*vyR - rhoL*vyL;
            const su2double drw = rhoR*vzR - rhoL*vzL;
            const su2double drE = rER      - rEL;

            // Compute the Roe average state.
            const su2double zL = SQRT(rhoL);
            const su2double zR = SQRT(rhoR);
            su2double tmp      = one/(zL + zR);

            const su2double rHL = rEL + pL;
            const su2double rHR = rER + pR;

            const su2double uAvg = tmp*(zL*vxL + zR*vxR);
            const su2double vAvg = tmp*(zL*vyL + zR*vyR);
            const su2double wAvg = tmp*(zL*vzL + zR*vzR);
            const su2double HAvg = tmp*(rHL/zL + rHR/zR);

            // Compute from the Roe average state some variables, which occur
            // quite often in the matrix vector product to be computed.
            const su2double alphaAvg = half*(uAvg*uAvg + vAvg*vAvg + wAvg*wAvg);
            tmp                      = gm1*(HAvg - alphaAvg);
            const su2double a2Avg    = FABS(tmp);
            const su2double aAvg     = SQRT(a2Avg);
            const su2double vnAvg    = uAvg*nx + vAvg*ny + wAvg*nz;
            const su2double ovaAvg   = one/aAvg;
            const su2double ova2Avg  = one/a2Avg;

            // Compute the absolute values of the three eigenvalues and
            // apply the entropy correction.
            su2double lam1 = FABS(vnAvg + aAvg);
            su2double lam2 = FABS(vnAvg - aAvg);
            su2double lam3 = FABS(vnAvg);

            tmp  = Delta*std::max(lam1, lam2);
            lam1 = std::max(lam1, tmp);
            lam2 = std::max(lam2, tmp);
            lam3 = std::max(lam3, tmp);

            // Some abbreviations, which occur quite often in the dissipation terms.
            const su2double abv1 = half*(lam1 + lam2);
            const su2double abv2 = half*(lam1 - lam2);
            const su2double abv3 = abv1 - lam3;

            const su2double abv4 = gm1*(alphaAvg*dr - uAvg*dru - vAvg*drv -wAvg*drw + drE);
            const su2double abv5 = nx*dru + ny*drv + nz*drw - vnAvg*dr;
            const su2double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
            const su2double abv7 = abv2*abv4*ovaAvg  + abv3*abv5;

            // Compute the Roe flux vector, which is 0.5*(FL + FR - |A|(UR-UL)).
            const su2double vnL = vxL*nx + vyL*ny + vzL*nz;
            const su2double vnR = vxR*nx + vyR*ny + vzR*nz;
            const su2double pa  = halfArea*(pL + pR);

            const su2double fact = factUNorm*halfArea;

            invFlux[0][l] = fact*(rhoL*vnL + rhoR*vnR - (lam3*dr + abv6));
            invFlux[1][l] = fact*(rhoL*vxL*vnL + rhoR*vxR*vnR
                          -       (lam3*dru + uAvg*abv6 + nx*abv7)) + pa*nx;
            invFlux[2][l] = fact*(rhoL*vyL*vnL + rhoR*vyR*vnR
                          -       (lam3*drv + vAvg*abv6 + ny*abv7)) + pa*ny;
            invFlux[3][l] = fact*(rhoL*vzL*vnL + rhoR*vzR*vnR
                          -       (lam3*drw + wAvg*abv6 + nz*abv7)) + pa*nz;
            invFlux[4][l] = fact*(rEL*vnL + rER*vnR + pL*vnL + pR*vnR
                          -       (lam3*drE + HAvg*abv6 + vnAvg*abv7));
          }
       
          break;
        }

        //----------------------------------------------------------------------

        case ISMAIL_ROE:
        {
          // Entropy satisfying Riemann flux of Ismail and Roe.
          // Loop over the (padded) number of integration points.
#pragma omp simd
          for(int l=0; l<nIntPad; ++l)
          {
            // Determine the coefficient multiplying the fluxes and compute
            // the components of the unit normal that points from the left
            // to the right state.
            const su2double Area = intWeights[l]*metric[3][l];
            const su2double nx = factNorm*metric[0][l],
                            ny = factNorm*metric[1][l],
                            nz = factNorm*metric[2][l];

            // Compute the primitive variables of the left and right state.
            const su2double V4InvL =  one/solL[4][l];
            const su2double vxL    = -V4InvL*solL[1][l];
            const su2double vyL    = -V4InvL*solL[2][l];
            const su2double vzL    = -V4InvL*solL[3][l]; 
            const su2double vTot2L =  vxL*vxL + vyL*vyL + vzL*vzL;
            const su2double sL     =  gam - gm1*(solL[0][l] - half*solL[4][l]*vTot2L);
            const su2double tmpL   = -solL[4][l]*EXP(sL);
            const su2double rhoL   =  POW(tmpL, ov1mg);
            const su2double pL     = -V4InvL*rhoL;

            const su2double V4InvR =  one/solR[4][l];
            const su2double vxR    = -V4InvR*solR[1][l];
            const su2double vyR    = -V4InvR*solR[2][l];
            const su2double vzR    = -V4InvR*solR[3][l]; 
            const su2double vTot2R =  vxR*vxR + vyR*vyR + vzR*vzR;
            const su2double sR     =  gam - gm1*(solR[0][l] - half*solR[4][l]*vTot2R);
            const su2double tmpR   = -solR[4][l]*EXP(sR);
            const su2double rhoR   =  POW(tmpR, ov1mg);
            const su2double pR     = -V4InvR*rhoR;

            // Store the left and right pressure in pInt. It will be used
            // later on for computing the forces, if desired.
            pInt[0][l] = pL;
            pInt[1][l] = pR;

            // Compute the difference in entropy variables.
            const su2double dV0 = solR[0][l] - solL[0][l];
            const su2double dV1 = solR[1][l] - solL[1][l];
            const su2double dV2 = solR[2][l] - solL[2][l];
            const su2double dV3 = solR[3][l] - solL[3][l];
            const su2double dV4 = solR[4][l] - solL[4][l];

            // Compute the z-variables of the left and right states.
            const su2double z0L = SQRT(rhoL/pL), z0R = SQRT(rhoR/pR);
            const su2double z1L = vxL*z0L,       z1R = vxR*z0R;
            const su2double z2L = vyL*z0L,       z2R = vyR*z0R;
            const su2double z3L = vzL*z0L,       z3R = vzR*z0R;
            const su2double z4L = SQRT(rhoL*pL), z4R = SQRT(rhoR*pR);

            // Compute the arithmetic average of the z-variables.
            const su2double z0Avg = half*(z0L + z0R);
            const su2double z1Avg = half*(z1L + z1R);
            const su2double z2Avg = half*(z2L + z2R);
            const su2double z3Avg = half*(z3L + z3R);
            const su2double z4Avg = half*(z4L + z4R);

            // Compute the logarithmic mean of z0.
            su2double zeta, f, u, F;

            zeta = z0L/z0R;
            f    = (zeta-one)/(zeta+one);
            u    = f*f;
            if(u < (su2double) 0.01)
              F = one + u/three + u*u/five + u*u*u/seven + u*u*u*u/nine;
            else
              F = LOG(zeta)/(two*f);

            const su2double z0LogAvg = z0Avg/F;

            // Compute the logarithmic mean of z4.
            zeta = z4L/z4R;
            f    = (zeta-one)/(zeta+one);
            u    = f*f;
            if(u < (su2double) 0.01)
              F = one + u/three + u*u/five + u*u*u/seven + u*u*u*u/nine;
            else
              F = LOG(zeta)/(two*f);

            const su2double z4LogAvg = z4Avg/F;

            // Compute the other averaged quantities that are necessary.
            const su2double oneOvz0Avg = one/z0Avg;
            const su2double rhoAvg = z0Avg*z4LogAvg;
            const su2double p1Avg  = oneOvz0Avg*z4Avg;
            const su2double p2Avg  = half*(gp1Ovg*z4LogAvg/z0LogAvg + gm1Ovg*p1Avg);
            const su2double vxAvg  = oneOvz0Avg*z1Avg;
            const su2double vyAvg  = oneOvz0Avg*z2Avg;
            const su2double vzAvg  = oneOvz0Avg*z3Avg;

            const su2double vnAvg  = vxAvg*nx + vyAvg*ny + vzAvg*nz;
            const su2double kinAvg = half*(vxAvg*vxAvg + vyAvg*vyAvg + vzAvg*vzAvg);
            const su2double a2Avg  = gam*p2Avg/rhoAvg;
            const su2double aAvg   = SQRT(a2Avg);
            const su2double HAvg   = a2Avg*ovgm1 + kinAvg;
            const su2double EAvg   = HAvg - p2Avg/rhoAvg;

            const su2double ovaAvg  = one/aAvg;
            const su2double ova2Avg = one/a2Avg;

            // Define the difference in conservative variables as dU/dV deltaV, where
            // the transformation matrix dU/dV must be evaluated at the averaged state.
            su2double abv1 = rhoAvg*(vxAvg*dV1 + vyAvg*dV2 + vzAvg*dV3);
            su2double abv2 = abv1 + rhoAvg*(dV0 + HAvg*dV4);

            const su2double dr  = abv1 + rhoAvg*(dV0 + EAvg*dV4);
            const su2double dru = vxAvg*abv2 + p2Avg*dV1;
            const su2double drv = vyAvg*abv2 + p2Avg*dV2;
            const su2double drw = vzAvg*abv2 + p2Avg*dV3;
            const su2double drE = HAvg*abv1 + rhoAvg*EAvg*(dV0 + HAvg*dV4)
                                + p2Avg*kinAvg*dV4;

            // Compute the absolute values of the eigenvalues of the flux Jacobian.
            su2double lam1 = FABS(vnAvg + aAvg);
            su2double lam2 = FABS(vnAvg - aAvg);
            su2double lam3 = FABS(vnAvg);

            // Scale the acoustic eigenvalue, such that the EC2 (or EC1) flux of Ismail
            // and Roe is obtained. Also multiply all eigenvalues by half to obtain
            // the correct scaling.
            lam1 *= half*(one + beta + alphaMax);
            lam2 *= half*(one + beta + alphaMax);
            lam3 *= half;

            // Some abbreviations, which occur quite often in the dissipation terms.
            abv1 = half*(lam1 + lam2);
            abv2 = half*(lam1 - lam2);

            const su2double abv3 = abv1 - lam3;
            const su2double abv4 = gm1*(kinAvg*dr - vxAvg*dru - vyAvg*drv - vzAvg*drw + drE);
            const su2double abv5 = nx*dru + ny*drv + nz*drw - vnAvg*dr;
            const su2double abv6 = abv3*abv4*ova2Avg + abv2*abv5*ovaAvg;
            const su2double abv7 = abv2*abv4*ovaAvg  + abv3*abv5;

            // Compute the central flux.
            invFlux[0][l] = factUNorm*rhoAvg*vnAvg;
            invFlux[1][l] = factUNorm*rhoAvg*vnAvg*vxAvg + p1Avg*nx;
            invFlux[2][l] = factUNorm*rhoAvg*vnAvg*vyAvg + p1Avg*ny;
            invFlux[3][l] = factUNorm*rhoAvg*vnAvg*vzAvg + p1Avg*nz;
            invFlux[4][l] = factUNorm*rhoAvg*vnAvg*HAvg;

            // Subtract the dissipation terms.
            invFlux[0][l] -= factUNorm*(lam3*dr + abv6);
            invFlux[1][l] -= factUNorm*(lam3*dru + vxAvg*abv6 + nx*abv7);
            invFlux[2][l] -= factUNorm*(lam3*drv + vyAvg*abv6 + ny*abv7);
            invFlux[3][l] -= factUNorm*(lam3*drw + vzAvg*abv6 + nz*abv7);
            invFlux[4][l] -= factUNorm*(lam3*drE + HAvg*abv6 + vnAvg*abv7);

            // Scale the fluxes with Area.
            invFlux[0][l] *= Area;
            invFlux[1][l] *= Area;
            invFlux[2][l] *= Area;
            invFlux[3][l] *= Area;
            invFlux[4][l] *= Area;
          }

          break;
        }

        //----------------------------------------------------------------------

        default:
        {
          // This is just to avoid a compiler warning.
          break;
        }
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

  // Update the inviscid forces, if needed.
  if( ComputeForces )
  {
    // Loop over the padded number of integration points.
    for(int l=0; l<nIntPad; ++l)
    {
      // Compute the averaged pressure times the area, integration weights
      // and the factor for the normal to make it pointing into the boundary.
      const su2double abv = half*factNorm*intWeights[l]*metric[3][l]
                          * (pInt[0][l] + pInt[1][l]);

      // Update the inviscid forces in the three directions.
      forceCoef[0] += abv*metric[0][l];
      forceCoef[1] += abv*metric[1][l];
      forceCoef[2] += abv*metric[2][l];
    }
  }
}

//------------------------------------------------------------------------------

// Function, which computes the penalty and symmetrizing fluxes when the
// conservative variables are used as working variables.
void PenAndSymFluxesConsVar(const InputParamClass *inputParam,
                            const int             nIntPad,
                            const su2double       factNorm,
                            su2double             **metric,
                            const su2double       *intWeights,
                            const su2double       lenScale,
                            su2double             **solAvg,
                            su2double             **dSol,
                            const su2double       *eddyVis,
                            const su2double       factHF_Lam,
                            const su2double       factHF_Turb,
                            su2double             **flux,
                            su2double             **fluxSymX,
                            su2double             **fluxSymY,
                            su2double             **fluxSymZ)
{
  // Compute the penalty parameter for the penalty terms.
  const su2double tau = (inputParam->mNPolySolDOFs+1)*(inputParam->mNPolySolDOFs+1)
                      *  inputParam->mRelativePenaltyParameter/lenScale;

  // Determine the viscous spectral radius over nu for both the laminar and the turbulent part.
  const su2double radVisOverNuLam  = std::max(std::max(one, two+lambdaOverMu), factHF_Lam);
  const su2double radVisOverNuTurb = std::max(std::max(one, two+lambdaOverMu), factHF_Turb);

  // Set the factor alpha for the two possibilities for the symmetrizing fluxes.
  const su2double alpha = lambdaOverMu;
  // const su2double alpha = one;

  // Other constants, which appear in the symmetrizing fluxes.
  const su2double beta     = lambdaOverMu + one - alpha;
  const su2double alphaP1  = alpha + one;
  const su2double lambdaP1 = lambdaOverMu + one;

  // Loop over the (padded) number of integration points.
#pragma omp simd
  for(int l=0; l<nIntPad; ++l)
  {
    // Determine the coefficient multiplying the fluxes and compute
    // the components of the scaled normal that points from the left
    // to the right state.
    const su2double Area = intWeights[l]*metric[3][l];
    su2double nx = Area*factNorm*metric[0][l],
              ny = Area*factNorm*metric[1][l],
              nz = Area*factNorm*metric[2][l];

    // Easier storage of the primitive variables.
    const su2double rho  = solAvg[0][l];
    const su2double u    = solAvg[1][l];
    const su2double v    = solAvg[2][l];
    const su2double w    = solAvg[3][l];
    const su2double Etot = solAvg[4][l];

    // Compute the inverse of the density.
    const su2double rhoInv = one/rho;

    // Easier storage of the difference of the conservative variables.
    const su2double dU0 = dSol[0][l];
    const su2double dU1 = dSol[1][l];
    const su2double dU2 = dSol[2][l];
    const su2double dU3 = dSol[3][l];
    const su2double dU4 = dSol[4][l];

    // Compute the penalty fluxes. Add them to the normal fluxes.
    // Use a scalar penalty for now, which means that also flux0
    // receives a contribution.
    const su2double tauScalar = tau*Area*rhoInv*(mu*radVisOverNuLam
                              +                  eddyVis[l]*radVisOverNuTurb);

    flux[0][l] += tauScalar*dU0;
    flux[1][l] += tauScalar*dU1;
    flux[2][l] += tauScalar*dU2;
    flux[3][l] += tauScalar*dU3;
    flux[4][l] += tauScalar*dU4;

    // Multiply the scaled normals by theta/2 to obtain the correct scaling
    // for the symmetrizing fluxes.
    nx *= half*inputParam->mThetaSym;
    ny *= half*inputParam->mThetaSym;
    nz *= half*inputParam->mThetaSym;

    // Abbreviations that occur in the symmetrizing fluxes.
    const su2double nu                  = (mu+eddyVis[l])*rhoInv;
    const su2double kScal               = (mu*factHF_Lam + eddyVis[l]*factHF_Turb)*rhoInv;
    const su2double velSquared          = u*u + v*v + w*w;
    const su2double velNorm             = u*nx + v*ny + w*nz;
    const su2double nuVelSquared        = nu*velSquared;
    const su2double nuVelNorm           = nu*velNorm;
    const su2double kScalEminVelSquared = kScal*(Etot - velSquared);

    const su2double nuVel[]    = {nu*u, nu*v, nu*w};
    const su2double kScalVel[] = {kScal*u, kScal*v, kScal*w};
    const su2double nuVelVel[] = {nu*u*velNorm, nu*v*velNorm, nu*w*velNorm};

    const su2double abv1 = nx*dU1 + ny*dU2 + nz*dU3;
    const su2double abv2 = nuVel[0]*dU1 + nuVel[1]*dU2 + nuVel[2]*dU3;

    const su2double abv2kScal = kScalVel[0]*dU1 + kScalVel[1]*dU2 + kScalVel[2]*dU3;

    const su2double abv3 = beta*(nu*abv1 - nuVelNorm*dU0);
    const su2double abv4 = kScal*dU4 - abv2kScal - kScalEminVelSquared*dU0 + abv2;

    // Compute the symmetrizing flux in x-direction.
    fluxSymX[0][l] =  zero;
    fluxSymX[1][l] = -(abv3 + alphaP1*nx*(nu*dU1 - nuVel[0]*dU0));
    fluxSymX[2][l] = -(nu*(nx*dU2 + alpha*ny*dU1) - (nx*nuVel[1] + alpha*ny*nuVel[0])*dU0);
    fluxSymX[3][l] = -(nu*(nx*dU3 + alpha*nz*dU1) - (nx*nuVel[2] + alpha*nz*nuVel[0])*dU0);
    fluxSymX[4][l] = -(nx*abv4 - (lambdaP1*nuVelVel[0] + nuVelSquared*nx)*dU0
                   +   alpha*nuVelNorm*dU1 + beta*nuVel[0]*abv1);

    // Compute the symmetrizing flux in y-direction.
    fluxSymY[0][l] =  zero;
    fluxSymY[1][l] = -(nu*(ny*dU1 + alpha*nx*dU2) - (ny*nuVel[0] + alpha*nx*nuVel[1])*dU0);
    fluxSymY[2][l] = -(abv3 + alphaP1*ny*(nu*dU2 - nuVel[1]*dU0));
    fluxSymY[3][l] = -(nu*(ny*dU3 + alpha*nz*dU2) - (ny*nuVel[2] + alpha*nz*nuVel[1])*dU0);
    fluxSymY[4][l] = -(ny*abv4 - (lambdaP1*nuVelVel[1] + nuVelSquared*ny)*dU0
                   +   alpha*nuVelNorm*dU2 + beta*nuVel[1]*abv1);

    // Compute the symmetrizing flux in z-direction.
    fluxSymZ[0][l] =  zero;
    fluxSymZ[1][l] = -(nu*(nz*dU1 + alpha*nx*dU3) - (nz*nuVel[0] + alpha*nx*nuVel[2])*dU0);
    fluxSymZ[2][l] = -(nu*(nz*dU2 + alpha*ny*dU3) - (nz*nuVel[1] + alpha*ny*nuVel[2])*dU0);
    fluxSymZ[3][l] = -(abv3 + alphaP1*nz*(nu*dU3 - nuVel[2]*dU0));
    fluxSymZ[4][l] = -(nz*abv4 - (lambdaP1*nuVelVel[2] + nuVelSquared*nz)*dU0
                   +   alpha*nuVelNorm*dU3 + beta*nuVel[2]*abv1);
  }
}

//------------------------------------------------------------------------------

// Function, which computes the penalty and symmetrizing fluxes when the
// entropy variables are used as working variables.
void PenAndSymFluxesSymmVar(const InputParamClass *inputParam,
                            const int              nIntPad,
                            const su2double        factNorm,
                            su2double            **metric,
                            const su2double       *intWeights,
                            const su2double        lenScale,
                            su2double            **solAvg,
                            su2double            **dSol,
                            const su2double       *eddyVis,
                            const su2double        factHF_Lam,
                            const su2double        factHF_Turb,
                            su2double            **flux,
                            su2double            **fluxSymX,
                            su2double            **fluxSymY,
                            su2double            **fluxSymZ)
{
  // Compute the penalty parameter for the penalty terms.
  const su2double tau = (inputParam->mNPolySolDOFs+1)*(inputParam->mNPolySolDOFs+1)
                      *  inputParam->mRelativePenaltyParameter/lenScale;

  // Compute some abbreviations involving Gamma.
  const su2double gm1   =  GamConstant - one;
  const su2double ovgm1 =  one/gm1;

  // Abbreviations involving lambda.
  const su2double lam   = lambdaOverMu;
  const su2double lamp1 = lambdaOverMu + one;
  const su2double lamp2 = lambdaOverMu + two;

  // Loop over the (padded) number of integration points.
#pragma omp simd
  for(int l=0; l<nIntPad; ++l)
  {
    // Compute the unit normal, which is the metric term multiplied by factNorm
    // to make sure that the normal points from the left state to the right state.
    const su2double nnx = factNorm*metric[0][l];
    const su2double nny = factNorm*metric[1][l];
    const su2double nnz = factNorm*metric[2][l];

    // Easier storage of the velocities and p/rho.
    const su2double u      = solAvg[1][l];
    const su2double v      = solAvg[2][l];
    const su2double w      = solAvg[3][l];
    const su2double pOvRho = solAvg[4][l];
    const su2double e      = pOvRho*ovgm1;

    // Easier storage of the difference in entropy variables.
    const su2double dV1 = dSol[1][l];
    const su2double dV2 = dSol[2][l];
    const su2double dV3 = dSol[3][l];
    const su2double dV4 = dSol[4][l];

    // Determine the coefficient multiplying the fluxes. This includes the
    // determinant, the integration weight of this point and p/rho.
    const su2double factFlux = pOvRho*intWeights[l]*metric[3][l]*(mu + eddyVis[l]);

    // Compute the normals scaled by factFlux.
    su2double nx = factFlux*nnx, ny = factFlux*nny, nz = factFlux*nnz;

    // Abbreviate an expression that appears in the viscous Jacobians.
    const su2double factHeatFlux = (mu*factHF_Lam + eddyVis[l]*factHF_Turb)
                                 / (mu + eddyVis[l]);
    const su2double abvE = u*u + v*v + w*w + factHeatFlux*e;

    // Compute the non-zero entries of the viscous Jacobians.
    const su2double K11_11 = lamp2, K11_14 = u*lamp2, K11_22 = one, K11_24 = v,
                    K11_33 = one, K11_34 = w, K11_44 = abvE + u*u*lamp1;
    const su2double K22_11 = one, K22_14 = u, K22_22 = lamp2, K22_24 = v*lamp2,
                    K22_33 = one, K22_34 = w, K22_44 = abvE + v*v*lamp1;
    const su2double K33_11 = one, K33_14 = u, K33_22 = one, K33_24 = v,
                    K33_33 = lamp2, K33_34 = w*lamp2, K33_44 = abvE + w*w*lamp1;

    const su2double K12_12 = lam, K12_14 = v*lam, K12_21 = one, K12_24 = u,
                    K12_41 = v,   K12_42 = u*lam, K12_44 = u*v*lamp1;
    const su2double K13_13 = lam, K13_14 = w*lam, K13_31 = one, K13_34 = u,
                    K13_41 = w,   K13_43 = u*lam, K13_44 = u*w*lamp1;
    const su2double K23_23 = lam, K23_24 = w*lam, K23_32 = one, K23_34 = v,
                    K23_42 = w,   K23_43 = v*lam, K23_44 = v*w*lamp1;

    // Compute the normal velocity and the normal difference of the symmetrizing variables.
    const su2double un  =   u*nnx +   v*nny +   w*nnz;
    const su2double dVn = dV1*nnx + dV2*nny + dV3*nnz;

    // Compute the penalty terms, which are added to flux1 to flux4. Note that flux0
    // does not get a penalty term when a matrix penalty is used.
    const su2double ttau   = tau*factFlux;
    const su2double abvMom = lamp1*(dVn + un*dV4);

    flux[1][l] += ttau*(abvMom*nnx + dV1 + u*dV4);
    flux[2][l] += ttau*(abvMom*nny + dV2 + v*dV4);
    flux[3][l] += ttau*(abvMom*nnz + dV3 + w*dV4);
    flux[4][l] += ttau*(lamp1*un*dVn + u*dV1 + v*dV2 + w*dV3 + (abvE + lamp1*un*un)*dV4);

    // Multiply the scaled normals by theta/2 to obtain the correct scaling
    // for the symmetrizing fluxes.
    nx *= half*inputParam->mThetaSym;
    ny *= half*inputParam->mThetaSym;
    nz *= half*inputParam->mThetaSym;

    // Compute the symmetrizing fluxes in x-direction, which are
    // -(K11 nx + K12 ny + K13 nz) dV.
    fluxSymX[0][l] =  zero;
    fluxSymX[1][l] = -(K11_11*dV1 + K11_14*dV4)*nx - (K12_12*dV2 + K12_14*dV4)*ny
                   -  (K13_13*dV3 + K13_14*dV4)*nz;
    fluxSymX[2][l] = -(K11_22*dV2 + K11_24*dV4)*nx - (K12_21*dV1 + K12_24*dV4)*ny;
    fluxSymX[3][l] = -(K11_33*dV3 + K11_34*dV4)*nx - (K13_31*dV1 + K13_34*dV4)*nz;
    fluxSymX[4][l] = -(K11_14*dV1 + K11_24*dV2 + K11_34*dV3 + K11_44*dV4)*nx
                   -  (K12_41*dV1 + K12_42*dV2 + K12_44*dV4)*ny
                   -  (K13_41*dV1 + K13_43*dV3 + K13_44*dV4)*nz;

    // Compute the symmetrizing fluxes in y-direction, which are
    // -(K21 nx + K22 ny + K23 nz) dV, where K21 = K12^T.
    fluxSymY[0][l] =  zero;
    fluxSymY[1][l] = -(K12_21*dV2 + K12_41*dV4)*nx - (K22_11*dV1 + K22_14*dV4)*ny;
    fluxSymY[2][l] = -(K12_12*dV1 + K12_42*dV4)*nx - (K22_22*dV2 + K22_24*dV4)*ny
                   -  (K23_23*dV3 + K23_24*dV4)*nz;
    fluxSymY[3][l] = -(K22_33*dV3 + K22_34*dV4)*ny - (K23_32*dV2 + K23_34*dV4)*nz;
    fluxSymY[4][l] = -(K12_14*dV1 + K12_24*dV2 + K12_44*dV4)*nx
                   -  (K22_14*dV1 + K22_24*dV2 + K22_34*dV3 + K22_44*dV4)*ny
                   -  (K23_42*dV2 + K23_43*dV3 + K23_44*dV4)*nz;

    // Compute the symmetrizing fluxes in z-direction, which are
    // -(K31 nx + K32 ny + K33 nz) dV, where K31 = K13^T and K32 = K23^T.
    fluxSymZ[0][l] =  zero;
    fluxSymZ[1][l] = -(K13_31*dV3 + K13_41*dV4)*nx - (K33_11*dV1 + K33_14*dV4)*nz;
    fluxSymZ[2][l] = -(K23_32*dV3 + K23_42*dV4)*ny - (K33_22*dV2 + K33_24*dV4)*nz;
    fluxSymZ[3][l] = -(K13_13*dV1 + K13_43*dV4)*nx - (K23_23*dV2 + K23_43*dV4)*ny
                   -  (K33_33*dV3 + K33_34*dV4)*nz;
    fluxSymZ[4][l] = -(K13_14*dV1 + K13_34*dV3 + K13_44*dV4)*nx
                   -  (K23_24*dV2 + K23_34*dV3 + K23_44*dV4)*ny
                   -  (K33_14*dV1 + K33_24*dV2 + K33_34*dV3 + K33_44*dV4)*nz;
  }
}

//------------------------------------------------------------------------------

// Function, which computes the viscous fluxes, including the penalty and
// symmetrizing fluxes on a face for the given solution and gradients.
void ViscousFluxesFace(const InputParamClass *inputParam,
                       const int              nIntPad,
                       const su2double        factNorm,
                       su2double            **metric,
                       const su2double       *intWeights,
                       const su2double        lenScale,
                       const su2double        lenScaleLES,
                       su2double            **solL,
                       su2double            **solR,
                       su2double            **dSolDx,
                       su2double            **dSolDy,
                       su2double            **dSolDz,
                       const bool             heatFluxPrescribed,
                       const su2double       *heatFlux,
                       su2double            **flux,
                       su2double            **fluxSymX,
                       su2double            **fluxSymY,
                       su2double            **fluxSymZ,
                       su2double             *eddyVis,
                       const bool             ComputeForces,
                       su2double             *forceCoef)
{
  // Set the multiplication factor of the heat flux, depending on whether
  // the heat flux is prescribed or not.
  const su2double multFactHeatFlux = heatFluxPrescribed ? zero : one;

  // Determine the corrected versions of factHeatFlux_Lam and factHeatFlux_Turb.
  const su2double factHF_Lam  = multFactHeatFlux*factHeatFlux_Lam;
  const su2double factHF_Turb = multFactHeatFlux*factHeatFlux_Turb;

  //----------------------------------------------------------------------------
  // Step 1: Compute the average solution, the difference of the solution
  //         and the gradients of velocity and internal energy. Use fluxSymX,
  //         solL and solR as storage. It must be used carefully in order not
  //         to overwrite data, but in this way no additional storage arrays
  //         are needed.
  //----------------------------------------------------------------------------

  // Set the pointers for temporary storage.
  su2double **solAvg  = solL;
  su2double **dSol    = solR;
  su2double **gradVel = fluxSymX;
  su2double **gradE   = gradVel + 9;

  // Make a distinction between the working variables.
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      // Loop over the padded number of integration points.
#pragma omp simd
      for(int l=0; l<nIntPad; ++l)
      {
        // Compute the average of the conservative variables.
        const su2double rho = half*(solL[0][l] + solR[0][l]);
        const su2double ru  = half*(solL[1][l] + solR[1][l]);
        const su2double rv  = half*(solL[2][l] + solR[2][l]);
        const su2double rw  = half*(solL[3][l] + solR[3][l]);
        const su2double rE  = half*(solL[4][l] + solR[4][l]);

        // Compute the velocities and the total energy per unit mass.
        const su2double rhoInv = one/rho;
        const su2double u      = rhoInv*ru;
        const su2double v      = rhoInv*rv;
        const su2double w      = rhoInv*rw;
        const su2double Etot   = rhoInv*rE;

        // Compute the gradients of the velocities.
        gradVel[0][l] = rhoInv*(dSolDx[1][l] - u*dSolDx[0][l]);   // dudx
        gradVel[1][l] = rhoInv*(dSolDx[2][l] - v*dSolDx[0][l]);   // dvdx.
        gradVel[2][l] = rhoInv*(dSolDx[3][l] - w*dSolDx[0][l]);   // dwdx.

        gradVel[3][l] = rhoInv*(dSolDy[1][l] - u*dSolDy[0][l]);   // dudy.
        gradVel[4][l] = rhoInv*(dSolDy[2][l] - v*dSolDy[0][l]);   // dvdy.
        gradVel[5][l] = rhoInv*(dSolDy[3][l] - w*dSolDy[0][l]);   // dwdy.

        gradVel[6][l] = rhoInv*(dSolDz[1][l] - u*dSolDz[0][l]);   // dudz.
        gradVel[7][l] = rhoInv*(dSolDz[2][l] - v*dSolDz[0][l]);   // dvdz.
        gradVel[8][l] = rhoInv*(dSolDz[3][l] - w*dSolDz[0][l]);   // dwdz.

        // Compute the gradients of the internal energy.
        gradE[0][l] = rhoInv*(dSolDx[4][l] - Etot*dSolDx[0][l])            // dedx.
                    - u*gradVel[0][l] - v*gradVel[1][l] - w*gradVel[2][l];
        gradE[1][l] = rhoInv*(dSolDy[4][l] - Etot*dSolDy[0][l])            // dedy.
                    - u*gradVel[3][l] - v*gradVel[4][l] - w*gradVel[5][l];
        gradE[2][l] = rhoInv*(dSolDz[4][l] - Etot*dSolDz[0][l])            // dedz.
                    - u*gradVel[6][l] - v*gradVel[7][l] - w*gradVel[8][l];

        // Store the difference of the conservative variables. Note that these must
        // be stored before the average solution, because solAvg overwrites solL.
        dSol[0][l] = solL[0][l] - solR[0][l];
        dSol[1][l] = solL[1][l] - solR[1][l];
        dSol[2][l] = solL[2][l] - solR[2][l];
        dSol[3][l] = solL[3][l] - solR[3][l];
        dSol[4][l] = solL[4][l] - solR[4][l];

        // Store the average primitive solution.
        solAvg[0][l] = rho;
        solAvg[1][l] = u;
        solAvg[2][l] = v;
        solAvg[3][l] = w;
        solAvg[4][l] = Etot;
      }

      break;
    }

    //------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      // Entropy variables are used. Compute some abbreviations involving Gamma.
      const su2double gm1   =  GamConstant - one;
      const su2double ovgm1 =  one/gm1;
      const su2double ov1mg = -ovgm1;

      // Loop over the padded number of integration points.
#pragma omp simd
      for(int l=0; l<nIntPad; ++l)
      {
        // Compute the average variables.
        const su2double V0 = half*(solL[0][l] + solR[0][l]);
        const su2double V1 = half*(solL[1][l] + solR[1][l]);
        const su2double V2 = half*(solL[2][l] + solR[2][l]);
        const su2double V3 = half*(solL[3][l] + solR[3][l]);
        const su2double V4 = half*(solL[4][l] + solR[4][l]);

        // Compute the primitive variables.
        const su2double V4Inv =  one/V4;
        const su2double u     = -V4Inv*V1;
        const su2double v     = -V4Inv*V2;
        const su2double w     = -V4Inv*V3;
        const su2double eKin  =  half*(u*u + v*v + w*w);
        const su2double s     =  GamConstant - gm1*(V0 - V4*eKin);
        const su2double tmp   = -V4*EXP(s);
        const su2double rho   =  POW(tmp, ov1mg);

        // Compute the velocity gradients.
        gradVel[0][l] = -V4Inv*(dSolDx[1][l] + u*dSolDx[4][l]);   // dudx.
        gradVel[1][l] = -V4Inv*(dSolDx[2][l] + v*dSolDx[4][l]);   // dvdx.
        gradVel[2][l] = -V4Inv*(dSolDx[3][l] + w*dSolDx[4][l]);   // dwdx.

        gradVel[3][l] = -V4Inv*(dSolDy[1][l] + u*dSolDy[4][l]);   // dudy.
        gradVel[4][l] = -V4Inv*(dSolDy[2][l] + v*dSolDy[4][l]);   // dvdy.
        gradVel[5][l] = -V4Inv*(dSolDy[3][l] + w*dSolDy[4][l]);   // dwdy.

        gradVel[6][l] = -V4Inv*(dSolDz[1][l] + u*dSolDz[4][l]);   // dudz.
        gradVel[7][l] = -V4Inv*(dSolDz[2][l] + v*dSolDz[4][l]);   // dvdz.
        gradVel[8][l] = -V4Inv*(dSolDz[3][l] + w*dSolDz[4][l]);   // dwdz.

        // Compute the gradients of the internal energy.
        const su2double fact = -V4Inv*V4Inv*ov1mg;

        gradE[0][l] = fact*dSolDx[4][l];   // dedx.
        gradE[1][l] = fact*dSolDy[4][l];   // dedy.
        gradE[2][l] = fact*dSolDz[4][l];   // dedz.

        // Store the difference of the entropy variables. Note that these must
        // be stored before the average solution, because solAvg overwrites solL.
        dSol[0][l] = solL[0][l] - solR[0][l];
        dSol[1][l] = solL[1][l] - solR[1][l];
        dSol[2][l] = solL[2][l] - solR[2][l];
        dSol[3][l] = solL[3][l] - solR[3][l];
        dSol[4][l] = solL[4][l] - solR[4][l];

        // Store the average primitive solution.
        solAvg[0][l] =  rho;
        solAvg[1][l] =  u;
        solAvg[2][l] =  v;
        solAvg[3][l] =  w;
        solAvg[4][l] = -V4Inv;  // p/rho.
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

  //----------------------------------------------------------------------------
  // Step 2: Compute the eddy viscosity if a subgrid scale model is used.
  //         If no SGS model is used, the eddy viscosity is set to zero.
  //----------------------------------------------------------------------------

  // Test for a subgrid scale model.
  if(inputParam->mSGSModelType != NO_SGS_MODEL)
  {
    // Compute the eddy viscosities.
    inputParam->mSGSModel->EddyViscosity(nIntPad, lenScaleLES, solAvg,
                                         gradVel, eddyVis);
  }
  else
  {
    // No subgrid scale model is used. Set the eddy viscosities to zero.
#pragma omp simd
    for(int l=0; l<nIntPad; ++l)
      eddyVis[l] = zero;
  }

  //----------------------------------------------------------------------------
  // Step 3: Compute the viscous fluxes in the integration points. As in step 1
  //         the gradients of velocity and internal energy have been computed,
  //         there is no need to make a distinction between the working variables.
  //----------------------------------------------------------------------------

  // Set the pointer to store the stress tensor, which may be needed
  // for computing the forces. Use dSolDx (and dSolDy) as storage.
  su2double **tauVis = dSolDx;

  // Loop over the padded number of integration points.
#pragma omp simd
  for(int l=0; l<nIntPad; ++l)
  {
    // Determine the components of the scaled (scaled by both area and
    // integration weight) normal that points from the left to the right state.
    const su2double Area = factNorm*intWeights[l]*metric[3][l];
    const su2double nx   = Area*metric[0][l],
                    ny   = Area*metric[1][l],
                    nz   = Area*metric[2][l];

    // More readable storage of the averaged velocities.
    const su2double u = solAvg[1][l];
    const su2double v = solAvg[2][l];
    const su2double w = solAvg[3][l];

    // More readable storage of the velocity gradients.
    const su2double dudx = gradVel[0][l];
    const su2double dvdx = gradVel[1][l];
    const su2double dwdx = gradVel[2][l];

    const su2double dudy = gradVel[3][l];
    const su2double dvdy = gradVel[4][l];
    const su2double dwdy = gradVel[5][l];

    const su2double dudz = gradVel[6][l];
    const su2double dvdz = gradVel[7][l];
    const su2double dwdz = gradVel[8][l];

    // More readable storage of the gradients of internal energy.
    const su2double dedx = gradE[0][l];
    const su2double dedy = gradE[1][l];
    const su2double dedz = gradE[2][l];
 
    // Compute the divergence of the velocity, multiplied by lambdaOverMu.
    const su2double lamDivVel = lambdaOverMu*(dudx + dvdy + dwdz);

    // Compute the elements of the viscous stress tensor.
    const su2double muTot = mu + eddyVis[l];

    const su2double tauxx = muTot*(two*dudx + lamDivVel);
    const su2double tauyy = muTot*(two*dvdy + lamDivVel);
    const su2double tauzz = muTot*(two*dwdz + lamDivVel);

    const su2double tauxy = muTot*(dudy + dvdx);
    const su2double tauxz = muTot*(dudz + dwdx);
    const su2double tauyz = muTot*(dvdz + dwdy);

    // Compute the elements of the heat flux vector.
    const su2double factHeatTot = mu*factHF_Lam + eddyVis[l]*factHF_Turb;

    const su2double qx = -factHeatTot*dedx;
    const su2double qy = -factHeatTot*dedy;
    const su2double qz = -factHeatTot*dedz;

    // Subtract the viscous normal flux from the currently stored values.
    // It is subtracted, because the viscous fluxes are defined with a
    // minus sign compared to the already stored inviscid fluxes.
    flux[1][l] -= tauxx*nx + tauxy*ny + tauxz*nz;
    flux[2][l] -= tauxy*nx + tauyy*ny + tauyz*nz;
    flux[3][l] -= tauxz*nx + tauyz*ny + tauzz*nz;
    flux[4][l] -= (u*tauxx + v*tauxy + w*tauxz - qx)*nx
                + (u*tauxy + v*tauyy + w*tauyz - qy)*ny
                + (u*tauxz + v*tauyz + w*tauzz - qz)*nz;

    // Store the viscous stress tensor.
    tauVis[0][l] = tauxx; tauVis[1][l] = tauyy; tauVis[2][l] = tauzz;
    tauVis[3][l] = tauxy; tauVis[4][l] = tauxz; tauVis[5][l] = tauyz;
  }

  // If the heat flux is prescribed, it should be added to the energy flux,
  // multiplied by the appropriate weighting factor.
  if( heatFluxPrescribed )
  {
    // Loop over the (padded) number of integration points.
#pragma omp simd
    for(int l=0; l<nIntPad; ++l)
      flux[4][l] += intWeights[l]*metric[3][l]*heatFlux[l];
  }

  //----------------------------------------------------------------------------
  // Step 4: Update the viscous forces, if forces must be computed.
  //----------------------------------------------------------------------------

  // Check if forces must be computed.
  if( ComputeForces )
  {
    // Loop over the padded number of integration points.
    for(int l=0; l<nIntPad; ++l)
    {
      // Determine the components of the scaled (scaled by both area and
      // integration weight) normal that points into the boundary.
      const su2double Area = factNorm*intWeights[l]*metric[3][l];
      const su2double nx   = Area*metric[0][l],
                      ny   = Area*metric[1][l],
                      nz   = Area*metric[2][l];

      // Easier storage of the elements of the stress tensor.
      const su2double tauxx = tauVis[0][l];
      const su2double tauyy = tauVis[1][l];
      const su2double tauzz = tauVis[2][l];
      const su2double tauxy = tauVis[3][l];
      const su2double tauxz = tauVis[4][l];
      const su2double tauyz = tauVis[5][l];

      // Update the viscous forces in the three directions.
      forceCoef[3] -= tauxx*nx + tauxy*ny + tauxz*nz;
      forceCoef[4] -= tauxy*nx + tauyy*ny + tauyz*nz;
      forceCoef[5] -= tauxz*nx + tauyz*ny + tauzz*nz;
    }
  }

  //----------------------------------------------------------------------------
  // Step 5: Compute the penalty and symmetrizing fluxes.
  //----------------------------------------------------------------------------

  // Make a distinction between the working variables and call the corresponding
  // function to do the actual work.
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      PenAndSymFluxesConsVar(inputParam, nIntPad, factNorm, metric, intWeights,
                             lenScale, solAvg, dSol, eddyVis, factHF_Lam,
                             factHF_Turb, flux, fluxSymX, fluxSymY, fluxSymZ);
      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      PenAndSymFluxesSymmVar(inputParam, nIntPad, factNorm, metric, intWeights,
                             lenScale, solAvg, dSol, eddyVis, factHF_Lam,
                             factHF_Turb, flux, fluxSymX, fluxSymY, fluxSymZ);
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

void ViscousFluxesFaceWallModel(ExchangeDataWallModelClass                               *exchangeDataWM,
                                std::vector<std::vector< std::vector<ElementClass *> > > &elements,
                                const InputParamClass                                    *inputParam,
                                const int                                                 nIntPad,
                                const su2double                                           factNorm,
                                su2double                                               **metric,
                                const su2double                                          *intWeights,
                                const su2double                                           lenScale,
                                const su2double                                           lenScaleLES,
                                su2double                                               **solL,
                                su2double                                               **solR,
                                const bool                                                heatFluxPrescribed,
                                const su2double                                          *prescribedWallData,
                                su2double                                               **flux,
                                su2double                                               **fluxSymX,
                                su2double                                               **fluxSymY,
                                su2double                                               **fluxSymZ,
                                su2double                                                *eddyVis,
                                const bool                                                ComputeForces,
                                su2double                                                *forceCoef)
{
  // Set the multiplication factor of the heat flux, depending on whether
  // the heat flux is prescribed or not.
  const su2double multFactHeatFlux = heatFluxPrescribed ? zero : one;

  // Determine the corrected versions of factHeatFlux_Lam and factHeatFlux_Turb.
  const su2double factHF_Lam  = multFactHeatFlux*factHeatFlux_Lam;
  const su2double factHF_Turb = multFactHeatFlux*factHeatFlux_Turb;

  //----------------------------------------------------------------------------
  // Step 1: Compute the data in the exchange points.
  //----------------------------------------------------------------------------

  // Use fluxSymY as storage for the solution in the exchange points.
  // Initialize the solution to zero.
  su2double **solExchange = fluxSymY;

  for(int m=0; m<nVar; ++m)
  {
#pragma omp simd
    for(int l=0; l<nIntPad; ++l)
      solExchange[m][l] = zero;
  }

  // Loop over the number of donor elements.
  for(unsigned int m=0; m<(exchangeDataWM->mNExchangePointsPerDonorElement.size()-1); ++m)
  {
    // Easier storage of the local indices of the donor element.
    const int i = exchangeDataWM->mIndicesDonorElements[3*m];
    const int j = exchangeDataWM->mIndicesDonorElements[3*m+1];
    const int k = exchangeDataWM->mIndicesDonorElements[3*m+2];

    // Easier storage of the solution of the exchange element.
    su2double **solElem = elements[i][j][k]->mSol.data();

    // Loop over the solution DOFs.
    for(unsigned int l=0; l<exchangeDataWM->mWeightsExchangePoints.size(); ++l)
    {
      // Loop over the integration points for this donor. This loop can be vectorized.
#pragma omp simd
      for(int n=exchangeDataWM->mNExchangePointsPerDonorElement[m];
              n<exchangeDataWM->mNExchangePointsPerDonorElement[m+1]; ++n)
      {
        // Easier storage of the interpolation weight.
        const su2double w = exchangeDataWM->mWeightsExchangePoints[l][n];

        // Update the interpolated solution.
        solExchange[0][n] += w*solElem[0][l];
        solExchange[1][n] += w*solElem[1][l];
        solExchange[2][n] += w*solElem[2][l];
        solExchange[3][n] += w*solElem[3][l];
        solExchange[4][n] += w*solElem[4][l];
      }
    }
  }

  // Store the solution of the exchange points in the right position and
  // convert it to dimensional primitive variables. Use fluxSymX as storage
  // for the primitive solution in the exchange points.
  su2double **solEx = fluxSymX;

  // Make a distinction between the working variables.
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      // The conservative variables are used as working variables.
      // Abbreviate gamma-1.
      const su2double gm1 = GamConstant - one;

      // Loop over the number of points for which the exchange solution
      // has been computed.
#pragma omp simd
      for(unsigned int l=0; l<exchangeDataWM->mExchangePointsPerDonorElement.size(); ++l)
      {
        // Easier storage of the corresponding integration point for this
        // exchange point.
        const int ll = exchangeDataWM->mExchangePointsPerDonorElement[l];

        // Compute the non-dimensional primitive variables.
        const su2double rho    = solExchange[0][l];
        const su2double rhoInv = one/rho;
        const su2double u      = rhoInv*solExchange[1][l];
        const su2double v      = rhoInv*solExchange[2][l];
        const su2double w      = rhoInv*solExchange[3][l];
        const su2double p      = gm1*(solExchange[4][l] - half*rho*(u*u + v*v + w*w));

        // Store the dimensional primitive variables.
        solEx[0][ll] = rhoRef*rho;
        solEx[1][ll] = uRef  *u;
        solEx[2][ll] = uRef  *v;
        solEx[3][ll] = uRef  *w;
        solEx[4][ll] = pRef  *p;
      }

      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      // The entropy variables are used as working variables.
      // Easier storage of some expressions involving gamma.
      const su2double gm1   =  GamConstant - one;
      const su2double ov1mg = -one/gm1;

      // Loop over the number of points for which the exchange solution
      // has been computed.
#pragma omp simd
      for(unsigned int l=0; l<exchangeDataWM->mExchangePointsPerDonorElement.size(); ++l)
      {
        // Easier storage of the corresponding integration point for this
        // exchange point.
        const int ll = exchangeDataWM->mExchangePointsPerDonorElement[l];
 
        // Compute the non-dimensional primitive variables.
        const su2double V4Inv =  one/solExchange[4][l];
        const su2double u     = -V4Inv*solExchange[1][l];
        const su2double v     = -V4Inv*solExchange[2][l];
        const su2double w     = -V4Inv*solExchange[3][l];
        const su2double eKin  =  half*(u*u + v*v + w*w);
        const su2double s     =  GamConstant - gm1*(solExchange[0][l] - solExchange[4][l]*eKin);
        const su2double tmp   = -solExchange[4][l]*EXP(s);
        const su2double rho   =  POW(tmp, ov1mg);
        const su2double p     = -rho*V4Inv; 

        // Store the dimensional primitive variables.
        solEx[0][ll] = rhoRef*rho;
        solEx[1][ll] = uRef  *u;
        solEx[2][ll] = uRef  *v;
        solEx[3][ll] = uRef  *w;
        solEx[4][ll] = pRef  *p;
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

  //----------------------------------------------------------------------------
  // Step 2: Compute the wall shear stress and heat flux, if needed,
  //         from the wall model.
  //----------------------------------------------------------------------------

  // Conversion factors from dimensional to non-dimensional for the wall
  // shear stress and heat flux.
  const su2double tauWallConversion = one/pRef;
  const su2double qWallConversion   = one/(pRef*uRef);

  // Compute the dimensional laminar viscosity.
  const su2double muDim = mu*pRef/uRef;

  // Again use fluxSymX as storage for the wall shear stress, wall heat flux
  // and the direction of the wall shear stress, which is direction of wall
  // parallel velocity.
  su2double  *tauWall = fluxSymX[0];
  su2double  *qWall   = fluxSymX[1];
  su2double **velDir  = fluxSymX + 2;

  // Loop over the number of points for which the wall shear stress must be computed.
  for(unsigned int l=0; l<exchangeDataWM->mExchangePointsPerDonorElement.size(); ++l)
  {
    // Easier storage of the primitive variables in the exchange point.
    const su2double rho = solEx[0][l];
    const su2double u   = solEx[1][l];
    const su2double v   = solEx[2][l];
    const su2double w   = solEx[3][l];
    const su2double p   = solEx[4][l];

    // Easier storage of the prescribed wall data for this integration point.
    // For an isothermal wall this is the wall temperature, while for a heat
    // flux wall this is the prescribed heat flux.
    const su2double wallData = prescribedWallData[l];

    // Compute the temperature in the exchange point.
    const su2double T = p/(RGas*rho);

    // Compute the normal velocity.
    const su2double un = u*metric[0][l] + v*metric[1][l] + w*metric[2][l];

    // Compute the components of the wall parallel velocity and its magnitude.
    const su2double uPar = u - un*metric[0][l];
    const su2double vPar = v - un*metric[1][l];
    const su2double wPar = w - un*metric[2][l];

    const su2double vParMag = SQRT(uPar*uPar + vPar*vPar + wPar*wPar);

    // Compute the wall shear stress and heat flux using the wall model.
    // The arrays fluxSymY and fluxSymZ are used as temporary storage.
    su2double tauW, qW, muW, kOverCvW;
    inputParam->mWallModel->WallShearStressAndHeatFlux(T, vParMag, muDim,
                                                       p, heatFluxPrescribed,
                                                       wallData, tauW, qW, muW,
                                                       kOverCvW, fluxSymY,
                                                       fluxSymZ);

    // Store the wall shear stress, heat flux and the direction of the
    // wall parallel velocity. Note that the wall shear stress and
    // heat flux must be made non-dimensional again.
    const su2double vParMagInv = one/std::max(vParMag, epsSmall);

    tauWall[l]   = tauW*tauWallConversion;
    qWall[l]     = qW*qWallConversion;
    velDir[0][l] = vParMagInv*uPar;
    velDir[1][l] = vParMagInv*vPar;
    velDir[2][l] = vParMagInv*wPar;
  }

  // Copy the data for the padded values to avoid problems.
  for(int l=(int) exchangeDataWM->mExchangePointsPerDonorElement.size(); l<nIntPad; ++l)
  {
    tauWall[l]   = tauWall[0];
    qWall[l]     = qWall[0];
    velDir[0][l] = velDir[0][0];
    velDir[1][l] = velDir[1][0];
    velDir[2][l] = velDir[2][0];
  }

  //----------------------------------------------------------------------------
  // Step 3: Update the viscous forces, if forces must be computed.
  //----------------------------------------------------------------------------

  // Check if forces must be computed.
  if( ComputeForces )
  {
    // Loop over the padded number of integration points.
    for(int l=0; l<nIntPad; ++l)
    {
      // Determine the Jacobian multiplied by the integration weight
      // and wall shear stress.
      const su2double tauScaled = intWeights[l]*metric[3][l]*tauWall[l];

      // Update the viscous forces in the three directions. Note that the
      // force direction is the direction of the  wall parallel velocity.
      forceCoef[3] += tauScaled*velDir[0][l];
      forceCoef[4] += tauScaled*velDir[1][l];
      forceCoef[5] += tauScaled*velDir[2][l];
    }
  }

  //----------------------------------------------------------------------------
  // Step 4: Compute the average solution and the difference of the solution.
  //         Use solL and solR as storage. It must be used carefully in order
  //         not to overwrite data, but in this way no additional storage arrays
  //         are needed.
  //----------------------------------------------------------------------------

  // Set the pointers for temporary storage.
  su2double **solAvg  = solL;
  su2double **dSol    = solR;

  // Make a distinction between the working variables.
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      // Loop over the padded number of integration points.
#pragma omp simd
      for(int l=0; l<nIntPad; ++l)
      {
        // Compute the average of the conservative variables.
        const su2double rho = half*(solL[0][l] + solR[0][l]);
        const su2double ru  = half*(solL[1][l] + solR[1][l]);
        const su2double rv  = half*(solL[2][l] + solR[2][l]);
        const su2double rw  = half*(solL[3][l] + solR[3][l]);
        const su2double rE  = half*(solL[4][l] + solR[4][l]);

        // Compute the velocities and the total energy per unit mass.
        const su2double rhoInv = one/rho;
        const su2double u      = rhoInv*ru;
        const su2double v      = rhoInv*rv;
        const su2double w      = rhoInv*rw;
        const su2double Etot   = rhoInv*rE;

        // Store the difference of the conservative variables. Note that these must
        // be stored before the average solution, because solAvg overwrites solL.
        dSol[0][l] = solL[0][l] - solR[0][l];
        dSol[1][l] = solL[1][l] - solR[1][l];
        dSol[2][l] = solL[2][l] - solR[2][l];
        dSol[3][l] = solL[3][l] - solR[3][l];
        dSol[4][l] = solL[4][l] - solR[4][l];

        // Store the average primitive solution.
        solAvg[0][l] = rho;
        solAvg[1][l] = u;
        solAvg[2][l] = v;
        solAvg[3][l] = w;
        solAvg[4][l] = Etot;
      }

      break;
    }

    //------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      // Entropy variables are used. Compute some abbreviations involving Gamma.
      const su2double gm1   =  GamConstant - one;
      const su2double ovgm1 =  one/gm1;
      const su2double ov1mg = -ovgm1;

      // Loop over the padded number of integration points.
#pragma omp simd
      for(int l=0; l<nIntPad; ++l)
      {
        // Compute the average variables.
        const su2double V0 = half*(solL[0][l] + solR[0][l]);
        const su2double V1 = half*(solL[1][l] + solR[1][l]);
        const su2double V2 = half*(solL[2][l] + solR[2][l]);
        const su2double V3 = half*(solL[3][l] + solR[3][l]);
        const su2double V4 = half*(solL[4][l] + solR[4][l]);

        // Compute the primitive variables.
        const su2double V4Inv =  one/V4;
        const su2double u     = -V4Inv*V1;
        const su2double v     = -V4Inv*V2;
        const su2double w     = -V4Inv*V3;
        const su2double eKin  =  half*(u*u + v*v + w*w);
        const su2double s     =  GamConstant - gm1*(V0 - V4*eKin);
        const su2double tmp   = -V4*EXP(s);
        const su2double rho   =  POW(tmp, ov1mg);

        // Store the difference of the entropy variables. Note that these must
        // be stored before the average solution, because solAvg overwrites solL.
        dSol[0][l] = solL[0][l] - solR[0][l];
        dSol[1][l] = solL[1][l] - solR[1][l];
        dSol[2][l] = solL[2][l] - solR[2][l];
        dSol[3][l] = solL[3][l] - solR[3][l];
        dSol[4][l] = solL[4][l] - solR[4][l];

        // Store the average primitive solution.
        solAvg[0][l] =  rho;
        solAvg[1][l] =  u;
        solAvg[2][l] =  v;
        solAvg[3][l] =  w;
        solAvg[4][l] = -V4Inv;  // p/rho.
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

  //----------------------------------------------------------------------------
  // Step 5: Compute the viscous fluxes in the integration points using the
  //         wall model. As the wall data has already been computed, there is
  //         no need to make a distinction between the working variables.
  //----------------------------------------------------------------------------

  // Loop over the padded number of integration points.
  for(int l=0; l<nIntPad; ++l)
  {
    // Easier storage of the velocity components.
    const su2double u = solAvg[1][l];
    const su2double v = solAvg[2][l];
    const su2double w = solAvg[3][l];

    // Compute the wall parallel velocity.
    const su2double uPar = u*velDir[0][l] + v*velDir[1][l] + w*velDir[2][l];

    // Determine the Jacobian multiplied by the integration weight
    // and wall shear stress.
    const su2double tauScaled = intWeights[l]*metric[3][l]*tauWall[l];

    // Update the fluxes. Note that the viscous fluxes are added instead of
    // subtracted, which is the convention in this code. The reason is that
    // the direction from the integration point to the exchange point is
    // opposite to the direction of the outward normal.
    flux[1][l] += tauScaled*velDir[0][l];
    flux[2][l] += tauScaled*velDir[1][l];
    flux[3][l] += tauScaled*velDir[2][l];
    flux[4][l] += tauScaled*uPar - intWeights[l]*metric[3][l]*qWall[l];
  }

  //----------------------------------------------------------------------------
  // Step 6: Compute the penalty and symmetrizing fluxes.
  //----------------------------------------------------------------------------

  // Set the eddy viscosities to zero.
#pragma omp simd
  for(int l=0; l<nIntPad; ++l)
    eddyVis[l] = zero;

  // Make a distinction between the working variables and call the corresponding
  // function to do the actual work.
  switch( inputParam->mFEMVariables )
  {
    case CONSERVATIVE_VARIABLES:
    {
      PenAndSymFluxesConsVar(inputParam, nIntPad, factNorm, metric, intWeights,
                             lenScale, solAvg, dSol, eddyVis, factHF_Lam,
                             factHF_Turb, flux, fluxSymX, fluxSymY, fluxSymZ);
      break;
    }

    //--------------------------------------------------------------------------

    case ENTROPY_VARIABLES:
    {
      PenAndSymFluxesSymmVar(inputParam, nIntPad, factNorm, metric, intWeights,
                             lenScale, solAvg, dSol, eddyVis, factHF_Lam,
                             factHF_Turb, flux, fluxSymX, fluxSymY, fluxSymZ);
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
