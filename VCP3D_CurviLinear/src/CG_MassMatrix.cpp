//------------------------------------------------------------------------------
// File, which contains the implementation of the standalone functions
// CG_ConservativeVarMassMatrix and CG_EntropyVarMassMatrix. These functions
// are a highly likely candidates for performance optimization and are
// therefore put in a separate file.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"

// Function prototypes of some BLAS routines.
extern "C"
{
#ifdef USE_SINGLE_PRECISION
  float sdot_(int *n, float *x, int *incx, float *y, int *incy);
  void saxpy_(int *n, float *a, float *x, int *incx, float *y, int *incy);
#else
  double ddot_(int *n, double *x, int *incx, double *y, int *incy);
  void  daxpy_(int *n, double *a, double *x, int *incx, double *y, int *incy);
#endif
}

//------------------------------------------------------------------------------

// Local function, which performs the axpy operation on the two dimensional
// vectors x and y.
void axpy(int       NVAR,
          int       n,
          su2double a,
          su2double **x,
          int       incx,
          su2double **y,
          int       incy)
{
  // Loop over the first dimension of the arrays.
  for(int l=0; l<NVAR; ++l)
  {
    // Call the appropriate function, depending on the precision used.
#ifdef USE_SINGLE_PRECISION
    saxpy_(&n, &a, x[l], &incx, y[l], &incy);
#else
    daxpy_(&n, &a, x[l], &incx, y[l], &incy);
#endif
  }
}

//------------------------------------------------------------------------------

// Local function, which computes the dot product of the two dimensional
// vectors x and y.
su2double dot(int       NVAR,
              int       n,
              su2double **x,
              int       incx,
              su2double **y,
              int       incy)
{
  // Initialize the result to zero.
  su2double result = zero;

  // Loop over NVAR to accumulate the value of the dot product.
  for(int l=0; l<NVAR; ++l)
  {
    // Call the appropriate function, depending on the precision used.
#ifdef USE_SINGLE_PRECISION
    result += sdot_(&n, x[l], &incx, y[l], &incy);
#else
    result += ddot_(&n, x[l], &incx, y[l], &incy);
#endif
  }

  // Return the result.
  return result;
}

//------------------------------------------------------------------------------

// Function, that solves the linear system of equations to obtain the residual
// for the conservative variables. A matrix free preconditioned flexible
// conjugate gradient method is used to do this. This function is a modified
// version of CG_EntropyVarMassMatrix, where the transformations dVdU and dUdV
// are left out. Hence it may not be fully optimal.
void CG_ConservativeVarMassMatrix(const StandardElementClass *standardHex,
                                  su2double                  **sol,
                                  const su2double            *Jac,
                                  const su2double            volElem,
                                  su2double                  **res,
                                  su2double                  **workArray)
{
  // Easier storage of the number of DOFs and integration points of the element.
  const int nDOFs1D  = standardHex->mNDOFs1D;
  const int nDOFsPad = standardHex->mNDOFsPad;
  const int nInt1D   = standardHex->mNIntegration1D;

  // Compute the inverse of the average Jacobian of the element.
  const su2double JacInv = eight/volElem;

  //----------------------------------------------------------------------------
  // Set the pointers for all the variables to the required memory locations.
  //----------------------------------------------------------------------------

  // Set the pointer for rVec. As the initialization is such that rVec == res,
  // the memory of res is used to store rVec.
  su2double **rVec = res;

  // Set the pointer to store the vectors array z, p, Ap and sol.
  su2double **zVec   = workArray;
  su2double **pVec   = zVec  + nVar;
  su2double **ApVec  = pVec  + nVar;
  su2double **solVec = ApVec + nVar;

  // Set the pointer for the interpolated p-vector in the integration points.
  su2double **pInt = solVec + nVar;
 
  //----------------------------------------------------------------------------
  // Initialization phase of the preconditioned conjugate gradient method.
  //----------------------------------------------------------------------------

  // Determine the maximum value of rVec. This serves as a scaling value
  // to determine the convergence.
  su2double rVecInitMax = zero;
  for(int l=0; l<nVar; ++l)
    for(int m=0; m<nDOFsPad; ++m)
      rVecInitMax = std::max(rVecInitMax, FABS(rVec[l][m]));

  // Initialize the solution to zero.
  for(int l=0; l<nVar; ++l)
#pragma omp simd
    for(int m=0; m<nDOFsPad; ++m)
      solVec[l][m] = zero;

  // If the value of rVecInitMax is very small, it does not make sense to
  // solve the linear system, because the solution is zero.
  // Test if this is not the case.
  if(rVecInitMax*JacInv > ((su2double) 1.e-10))
  {
    // Carry out the preconditioning step (multiply rVec with the inverse
    // of the Jacobian) and copy zVec into pVec.
    for(int l=0; l<nVar; ++l)
    {
#pragma omp simd
      for(int m=0; m<nDOFsPad; ++m)
      {
        zVec[l][m] = JacInv*rVec[l][m];
        pVec[l][m] = zVec[l][m];
      }
    }

    //--------------------------------------------------------------------------
    // Iterative phase of the preconditioned conjugate gradient method.
    //--------------------------------------------------------------------------

    // Loop over the number of iterations.
    for(int iter=0;;++iter)
    {
      // Safeguard to avoid an infinite loop.
      if(iter == 100)
        Terminate("CG_ConservativeVarMassMatrix", __FILE__, __LINE__,
                  "Convergence not reached for CG algorithm");

      // Compute the dot product of rVec and zVec.
      const su2double dotrz = dot(nVar, nDOFsPad, rVec, 1, zVec, 1);

      // The matrix vector product of the modified mass matrix and the vector pVec
      // must be determined. This is done in three steps. Step 1 is to interpolate
      // the data of pVec to the integration points. The result will be stored in pInt.
      TensorProductSolAndGradVolume(nInt1D, nVar, nDOFs1D, standardHex->mLegendreInt1D,
                                    NULL, pVec, pInt, NULL, NULL, NULL);

      // Step 2 of the matrix vector product. Multiply pInt with the integration
      // weight and the Jacobian. The loop is carried out over the padded number of
      // integration points for performance reasons.
      for(int l=0; l<nVar; ++l)
#pragma omp simd
        for(int m=0; m<nDOFsPad; ++m)
          pInt[l][m] *= standardHex->mIntWeights[l]*Jac[l];

      // Step 3 of the matrix vector product. Scatter the results of pInt back to the
      // DOFs. This is the final result, which is stored in ApVec. This array must be
      // initialized to zero, as the function TensorProductVolumeResidual accumulates
      // the values.
      for(int l=0; l<nVar; ++l)
#pragma omp simd
        for(int m=0; m<nDOFsPad; ++m)
          ApVec[l][m] = zero;

      TensorProductVolumeResidual(nInt1D, nVar, nDOFs1D,
                                  standardHex->mLegendreInt1DTranspose,
                                  standardHex->mLegendreInt1DTranspose,
                                  standardHex->mLegendreInt1DTranspose,
                                  pInt, ApVec);

      // Determine the dot product between pVec and ApVec.
      const su2double dotpAp = dot(nVar, nDOFsPad, pVec, 1, ApVec, 1);

      // Determine the coefficient alpha, which is the contribution of
      // pVec to the solution and -ApVec to the right hand side.
      const su2double alpha = dotrz/dotpAp;

      // Compute the new solution and right hand side.
      axpy(nVar, nDOFsPad,  alpha,  pVec, 1, solVec, 1);
      axpy(nVar, nDOFsPad, -alpha, ApVec, 1, rVec,   1);

      // Determine the Linf norm of rVec. Needed to check the convergence of
      // the iterative algorithm.
      su2double rVecMax = zero;
      for(int l=0; l<nVar; ++l)
        for(int m=0; m<nDOFsPad; ++m)
          rVecMax = std::max(rVecMax, FABS(rVec[l][m]));

      // Convergence criterion. At the moment the convergence criterion is
      // hard coded, but could be replaced by a user defined parameter.
      if(rVecMax < ((su2double) 1.e-10*rVecInitMax)) break;

      // Preconditioning step. Multiply rVec with the inverse of the Jacobian and
      // store the results in zVec. Loop over the padded value of nDOFs for
      // performance reasons.
      for(int l=0; l<nVar; ++l)
#pragma omp simd
        for(int m=0; m<nDOFsPad; ++m)
          zVec[l][m] = JacInv*rVec[l][m];

      // Compute the value of the coefficient beta. This is a flexible CG, hence the
      // Polak-Ribiere formula must be used.
      const su2double beta = -alpha*dot(nVar, nDOFsPad, zVec, 1, ApVec, 1)/dotrz;

      // Compute the new p vector.
      for(int l=0; l<nVar; ++l)
#pragma omp simd
        for(int m=0; m<nDOFsPad; ++m)
          pVec[l][m] = beta*pVec[l][m] + zVec[l][m];
    }
  }

  //----------------------------------------------------------------------------
  // Copy the data from solVec back into res.
  //----------------------------------------------------------------------------

  for(int l=0; l<nVar; ++l)
#pragma omp simd
    for(int m=0; m<nDOFsPad; ++m)
      res[l][m] = solVec[l][m];
}

//------------------------------------------------------------------------------

// Function, that solves the linear system of equations to obtain the residual
// for the entropy variables. A matrix free preconditioned flexible conjugate
// gradient method is used to do this.
void CG_EntropyVarMassMatrix(const StandardElementClass *standardHex,
                             su2double                  **sol,
                             const su2double            *Jac,
                             const su2double            volElem,
                             su2double                  **dUdV,
                             su2double                  **res,
                             su2double                  **workArray)
{
  // Easier storage of the number of DOFs and integration points of the element.
  const int nDOFs1D  = standardHex->mNDOFs1D;
  const int nDOFs    = standardHex->mNDOFs;
  const int nDOFsPad = standardHex->mNDOFsPad;
  const int nInt1D   = standardHex->mNIntegration1D;
  const int nIntPad  = standardHex->mNIntegrationPad;

  // Compute the inverse of the average Jacobian of the element.
  const su2double JacInv = eight/volElem;

  //----------------------------------------------------------------------------
  // Set the pointers for all the variables to the required memory locations.
  //----------------------------------------------------------------------------

  // Set the pointer for rVec. As the initialization is such that rVec == res,
  // the memory of res is used to store rVec.
  su2double **rVec = res;

  // Set the pointer to store the vectors array z, p, Ap and sol.
  su2double **zVec   = workArray;
  su2double **pVec   = zVec  + nVar;
  su2double **ApVec  = pVec  + nVar;
  su2double **solVec = ApVec + nVar;

  // Set the pointer for the interpolated p-vector in the integration points.
  su2double **pInt = solVec + nVar;

  // Determine the maximum value of rVec. This serves as a scaling value
  // to determine the convergence.
  su2double rVecInitMax = zero;
  for(int l=0; l<nVar; ++l)
    for(int m=0; m<nDOFsPad; ++m)
      rVecInitMax = std::max(rVecInitMax, FABS(rVec[l][m]));

  // Initialize the solution to zero.
  for(int l=0; l<nVar; ++l)
#pragma omp simd
    for(int m=0; m<nDOFsPad; ++m)
      solVec[l][m] = zero;

  // If the value of rVecInitMax is very small, it does not make sense to
  // solve the linear system, because the solution is zero.
  // Test if this is not the case.
  if(rVecInitMax*JacInv > ((su2double) 1.e-10))
  {
    //--------------------------------------------------------------------------
    // Compute the elements of the transformation matrix dVdU in the center of
    // the element.
    //--------------------------------------------------------------------------

    // Abbreviations involving gamma.
    const su2double gm1    = GamConstant - one;
    const su2double ov1mg  = one/(one-GamConstant);
    const su2double govgm1 = GamConstant/gm1;

    // Compute the entropy variables in the center of the element.
    su2double V0 = zero, V1 = zero, V2 = zero, V3 = zero, V4 = zero;
    for(int l=0; l<nDOFs; ++l)
    {
      V0 += standardHex->mBasisCenter[l]*sol[0][l];
      V1 += standardHex->mBasisCenter[l]*sol[1][l];
      V2 += standardHex->mBasisCenter[l]*sol[2][l];
      V3 += standardHex->mBasisCenter[l]*sol[3][l];
      V4 += standardHex->mBasisCenter[l]*sol[4][l];
    }

    // Compute the primitive variables from the entropy ones.
    const su2double V4Inv =  one/V4;
    const su2double u     = -V4Inv*V1;
    const su2double v     = -V4Inv*V2;
    const su2double w     = -V4Inv*V3;
    const su2double eKin  =  half*(u*u + v*v + w*w);
    const su2double s     =  GamConstant - gm1*(V0 - V4*eKin);
    const su2double tmp   = -V4*EXP(s);
    const su2double rho   =  POW(tmp, ov1mg);
    const su2double p     = -rho*V4Inv;
    const su2double pInv  =  one/p;

    // Two abbreviations that appear in the elements of dVdU.
    const su2double abv1 = gm1*rho*pInv*pInv;
    const su2double abv2 = eKin*abv1;

    // Store the elements of the transformation matrix dVdU in this DOF.
    // Note that this matrix is symmetric and hence only the upper-diagonal
    // part (or lower diagonal part) is stored. Multiply the transformation
    // matrix by the inverse of the Jacobian to account for the transformation
    // to the standard element.
    su2double dVdU[15];

    dVdU[0]  =  JacInv*(govgm1/rho + abv2*eKin);      // dVdU(0,0)
    dVdU[1]  = -JacInv*u*abv2;                        // dVdU(0,1) = dVdU(1,0)
    dVdU[2]  = -JacInv*v*abv2;                        // dVdU(0,2) = dVdU(2,0)
    dVdU[3]  = -JacInv*w*abv2;                        // dVdU(0,3) = dVdU(3,0)
    dVdU[4]  =  JacInv*(abv2 - pInv);                 // dVdU(0,4) = dVdU(4,0)
    dVdU[5]  =  JacInv*(abv1*u*u + pInv);             // dVdU(1,1)
    dVdU[6]  =  JacInv*abv1*u*v;                      // dVdU(1,2) = dVdU(2,1)
    dVdU[7]  =  JacInv*abv1*u*w;                      // dVdU(1,3) = dVdU(3,1)
    dVdU[8]  = -JacInv*abv1*u;                        // dVdU(1,4) = dVdU(4,1)
    dVdU[9]  =  JacInv*(abv1*v*v + pInv);             // dVdU(2,2)
    dVdU[10] =  JacInv*abv1*v*w;                      // dVdU(2,3) = dVdU(3,2)
    dVdU[11] = -JacInv*abv1*v;                        // dVdU(2,4) = dVdU(4,2)
    dVdU[12] =  JacInv*(abv1*w*w + pInv);             // dVdU(3,3)
    dVdU[13] = -JacInv*abv1*w;                        // dVdU(3,4) = dVdU(4,3)
    dVdU[14] =  JacInv*abv1;                          // dVdU(4,4)

    //--------------------------------------------------------------------------
    // Initialization phase of the preconditioned conjugate gradient method.
    //--------------------------------------------------------------------------

    // Preconditioning step. Multiply rVec with the dVdU and store the result in zVec.
    // Loop over the padded value of nDOFs for performance reasons.
#pragma omp simd
    for(int l=0; l<nDOFsPad; ++l)
    {
      zVec[0][l] = dVdU[0] *rVec[0][l] + dVdU[1] *rVec[1][l] + dVdU[2] *rVec[2][l]
                 + dVdU[3] *rVec[3][l] + dVdU[4] *rVec[4][l];
      zVec[1][l] = dVdU[1] *rVec[0][l] + dVdU[5] *rVec[1][l] + dVdU[6] *rVec[2][l]
                 + dVdU[7] *rVec[3][l] + dVdU[8] *rVec[4][l];
      zVec[2][l] = dVdU[2] *rVec[0][l] + dVdU[6] *rVec[1][l] + dVdU[9] *rVec[2][l]
                 + dVdU[10]*rVec[3][l] + dVdU[11]*rVec[4][l];
      zVec[3][l] = dVdU[3] *rVec[0][l] + dVdU[7] *rVec[1][l] + dVdU[10]*rVec[2][l]
                 + dVdU[12]*rVec[3][l] + dVdU[13]*rVec[4][l];
      zVec[4][l] = dVdU[4] *rVec[0][l] + dVdU[8] *rVec[1][l] + dVdU[11]*rVec[2][l]
                 + dVdU[13]*rVec[3][l] + dVdU[14]*rVec[4][l];
    }

    // Copy zVec into pVec.
    for(int l=0; l<nVar; ++l)
#pragma omp simd
      for(int m=0; m<nDOFsPad; ++m)
        pVec[l][m] = zVec[l][m];

    //--------------------------------------------------------------------------
    // Iterative phase of the preconditioned conjugate gradient method.
    //--------------------------------------------------------------------------

    // Loop over the number of iterations.
    for(int iter=0;;++iter)
    {
      // Safeguard to avoid an infinite loop.
      if(iter == 100)
        Terminate("CG_EntropyVarMassMatrix", __FILE__, __LINE__,
                  "Convergence not reached for CG algorithm");

      // Compute the dot product of rVec and zVec.
      const su2double dotrz = dot(nVar, nDOFsPad, rVec, 1, zVec, 1);

      // The matrix vector product of the modified mass matrix and the vector pVec
      // must be determined. This is done in three steps. Step 1 is to interpolate
      // the data of pVec to the integration points. The result will be stored in pInt.
      TensorProductSolAndGradVolume(nInt1D, nVar, nDOFs1D, standardHex->mLegendreInt1D,
                                    NULL, pVec, pInt, NULL, NULL, NULL);

      // Step 2 of the matrix vector product. Multiply pInt with dUdV, the integration
      // weight and the Jacobian. The loop is carried out over the padded number of
      // integration points for performance reasons.
#pragma omp simd
      for(int l=0; l<nIntPad; ++l)
      {
        const su2double intWeight = standardHex->mIntWeights[l]*Jac[l];

        const su2double p0 = intWeight*pInt[0][l], p1 = intWeight*pInt[1][l],
                        p2 = intWeight*pInt[2][l], p3 = intWeight*pInt[3][l],
                        p4 = intWeight*pInt[4][l];

        pInt[0][l] = dUdV[0][l]*p0 + dUdV[1][l]*p1 + dUdV[2][l] *p2 + dUdV[3][l] *p3 + dUdV[4][l] *p4;
        pInt[1][l] = dUdV[1][l]*p0 + dUdV[5][l]*p1 + dUdV[6][l] *p2 + dUdV[7][l] *p3 + dUdV[8][l] *p4;
        pInt[2][l] = dUdV[2][l]*p0 + dUdV[6][l]*p1 + dUdV[9][l] *p2 + dUdV[10][l]*p3 + dUdV[11][l]*p4;
        pInt[3][l] = dUdV[3][l]*p0 + dUdV[7][l]*p1 + dUdV[10][l]*p2 + dUdV[12][l]*p3 + dUdV[13][l]*p4;
        pInt[4][l] = dUdV[4][l]*p0 + dUdV[8][l]*p1 + dUdV[11][l]*p2 + dUdV[13][l]*p3 + dUdV[14][l]*p4;
      }

      // Step 3 of the matrix vector product. Scatter the results of pInt back to the
      // DOFs. This is the final result, which is stored in ApVec. This array must be
      // initialized to zero, as the function TensorProductVolumeResidual accumulates
      // the values.
      for(int l=0; l<nVar; ++l)
#pragma omp simd
        for(int m=0; m<nDOFsPad; ++m)
          ApVec[l][m] = zero;

      TensorProductVolumeResidual(nInt1D, nVar, nDOFs1D,
                                  standardHex->mLegendreInt1DTranspose,
                                  standardHex->mLegendreInt1DTranspose,
                                  standardHex->mLegendreInt1DTranspose,
                                  pInt, ApVec);

      // Determine the dot product between pVec and ApVec.
      const su2double dotpAp = dot(nVar, nDOFsPad, pVec, 1, ApVec, 1);

      // Determine the coefficient alpha, which is the contribution of
      // pVec to the solution and -ApVec to the right hand side.
      const su2double alpha = dotrz/dotpAp;

      // Compute the new solution and right hand side.
      axpy(nVar, nDOFsPad,  alpha,  pVec, 1, solVec, 1);
      axpy(nVar, nDOFsPad, -alpha, ApVec, 1, rVec,   1);

      // Determine the Linf norm of rVec. Needed to check the convergence of
      // the iterative algorithm.
      su2double rVecMax = zero;
      for(int l=0; l<nVar; ++l)
        for(int m=0; m<nDOFsPad; ++m)
          rVecMax = std::max(rVecMax, FABS(rVec[l][m]));

      // Convergence criterion. At the moment the convergence criterion is
      // hard coded, but could be replaced by a user defined parameter.
      if(rVecMax < ((su2double) 1.e-10*rVecInitMax)) break;

      // Preconditioning step. Multiply rVec with the dVdU and store the results in zVec.
      // Loop over the padded value of nDOFs for performance reasons.
#pragma omp simd
      for(int l=0; l<nDOFsPad; ++l)
      {
        zVec[0][l] = dVdU[0] *rVec[0][l] + dVdU[1] *rVec[1][l] + dVdU[2] *rVec[2][l]
                   + dVdU[3] *rVec[3][l] + dVdU[4] *rVec[4][l];
        zVec[1][l] = dVdU[1] *rVec[0][l] + dVdU[5] *rVec[1][l] + dVdU[6] *rVec[2][l]
                   + dVdU[7] *rVec[3][l] + dVdU[8] *rVec[4][l];
        zVec[2][l] = dVdU[2] *rVec[0][l] + dVdU[6] *rVec[1][l] + dVdU[9] *rVec[2][l]
                   + dVdU[10]*rVec[3][l] + dVdU[11]*rVec[4][l];
        zVec[3][l] = dVdU[3] *rVec[0][l] + dVdU[7] *rVec[1][l] + dVdU[10]*rVec[2][l]
                   + dVdU[12]*rVec[3][l] + dVdU[13]*rVec[4][l];
        zVec[4][l] = dVdU[4] *rVec[0][l] + dVdU[8] *rVec[1][l] + dVdU[11]*rVec[2][l]
                   + dVdU[13]*rVec[3][l] + dVdU[14]*rVec[4][l];
      }

      // Compute the value of the coefficient beta. This is a flexible CG, hence the
      // Polak-Ribiere formula must be used.
      const su2double beta = -alpha*dot(nVar, nDOFsPad, zVec, 1, ApVec, 1)/dotrz;

      // Compute the new p vector.
      for(int l=0; l<nVar; ++l)
#pragma omp simd
        for(int m=0; m<nDOFsPad; ++m)
          pVec[l][m] = beta*pVec[l][m] + zVec[l][m];
    }
  }

  //----------------------------------------------------------------------------
  // Copy the data from solVec back into res.
  //----------------------------------------------------------------------------

  for(int l=0; l<nVar; ++l)
#pragma omp simd
    for(int m=0; m<nDOFsPad; ++m)
      res[l][m] = solVec[l][m];
}
