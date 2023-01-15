//------------------------------------------------------------------------------
// TensorProductResKFaceK7.cpp: Implementation of the function that carries out
//                              the tensor product multiplications to compute
//                              the contribution from a K face to the residual.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"
#include "TensorProductResKFace.hpp"

//------------------------------------------------------------------------------

// Function, which performs the tensor product multiplication to accumulate the
// residual with a contribution from a face in k-direction. The C tensor is a
// rank 4 tensor of dimension C[K][K][K][N] stored in column major order, B is a
// rank 3 tensor of dimension B[M][M][N] stored in column major order, ATx and
// ATy are rank 2 tensors (i.e. matrices) of dimension AT[K][M] in column major
// order where the first dimension is padded to a multiple of the required vector
// length, and AFacez is a rank 1 tensor of dimension AFacez[K].
// The tensor product carried out is
// C += ATx*(ATy*(AFacez*B))
void TensorProductResKFaceK7(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *ATy,
                             su2double *AFacez,
                             su2double **B,
                             su2double **C)
{
  // In order to obtain the best performance the sizes of M and K must be known at
  // compile time. Therefore determine the size of these variables and call the
  // corresponding template instantiation of the function that carries out the
  // actual multiplication.
  switch( M )
  {
    case  9: TPRKF<7, 9>(N,ATx,ATy,AFacez,B,C); break;
    default: TerminateTensorProduct("TensorProductResKFace", 7, M);
  }
}
