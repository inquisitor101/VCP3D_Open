//------------------------------------------------------------------------------
// TensorProductResIFace.cpp: Implementation of the functions that carry out the
//                            tensor product multiplications to compute the
//                            contribution from an I face to the residual.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"

//------------------------------------------------------------------------------

// Function, which performs the tensor product multiplication to accumulate the
// residual with a contribution from a face in i-direction. The C tensor is a
// rank 4 tensor of dimension C[K][K][K][N] stored in column major order, B is a
// rank 3 tensor of dimension B[M][M][N] stored in column major order, ATy and
// ATz are rank 2 tensors (i.e. matrices) of dimension AT[K][M] in column major
// order where the first dimension is padded to a multiple of the required vector
// length, and AFacex is a rank 1 tensor of dimension AFacex[K].
// The tensor product carried out is
// C += AFacex*(ATy*(ATz*B))
void TensorProductResIFace(const int M,
                           const int N,
                           const int K,
                           su2double *AFacex,
                           su2double *ATy,
                           su2double *ATz,
                           su2double **B,
                           su2double **C)
{
  // In order to obtain the best performance the sizes of M and K must be known at
  // compile time. Therefore determine the size of these variables and call the
  // corresponding function that carries out the actual multiplication.
  switch( K )
  {
    case  2: TensorProductResIFaceK2(M,N,AFacex,ATy,ATz,B,C); break;
    case  3: TensorProductResIFaceK3(M,N,AFacex,ATy,ATz,B,C); break;
    case  4: TensorProductResIFaceK4(M,N,AFacex,ATy,ATz,B,C); break;
    case  5: TensorProductResIFaceK5(M,N,AFacex,ATy,ATz,B,C); break;
    case  6: TensorProductResIFaceK6(M,N,AFacex,ATy,ATz,B,C); break;
    case  7: TensorProductResIFaceK7(M,N,AFacex,ATy,ATz,B,C); break;
    case  8: TensorProductResIFaceK8(M,N,AFacex,ATy,ATz,B,C); break;
    case  9: TensorProductResIFaceK9(M,N,AFacex,ATy,ATz,B,C); break;
    case 10: TensorProductResIFaceK10(M,N,AFacex,ATy,ATz,B,C); break;
    default:
      TerminateAll("TensorProductResIFace", __FILE__, __LINE__,
                   "Value of K not considered");
  }
}
