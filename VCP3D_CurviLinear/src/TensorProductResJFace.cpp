//------------------------------------------------------------------------------
// TensorProductResJFace.cpp: Implementation of the functions that carry out the
//                            tensor product multiplications to compute the
//                            contribution from a J face to the residual.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"

//------------------------------------------------------------------------------

// Function, which performs the tensor product multiplication to accumulate the
// residual with a contribution from a face in j-direction. The C tensor is a
// rank 4 tensor of dimension C[K][K][K][N] stored in column major order, B is a
// rank 3 tensor of dimension B[M][M][N] stored in column major order, ATx and
// ATz are rank 2 tensors (i.e. matrices) of dimension AT[K][M] in column major
// order where the first dimension is padded to a multiple of the required vector
// length, and AFacey is a rank 1 tensor of dimension AFacey[K].
// The tensor product carried out is
// C += ATx*(AFacey*(ATz*B))
void TensorProductResJFace(const int M,
                           const int N,
                           const int K,
                           su2double *ATx,
                           su2double *AFacey,
                           su2double *ATz,
                           su2double **B,
                           su2double **C)
{
  // In order to obtain the best performance the sizes of M and K must be known at
  // compile time. Therefore determine the size of these variables and call the
  // corresponding template instantiation of the function that carries out the
  // actual multiplication.
  switch( K )
  {
    case  2: TensorProductResJFaceK2(M,N,ATx,AFacey,ATz,B,C); break;
    case  3: TensorProductResJFaceK3(M,N,ATx,AFacey,ATz,B,C); break;
    case  4: TensorProductResJFaceK4(M,N,ATx,AFacey,ATz,B,C); break;
    case  5: TensorProductResJFaceK5(M,N,ATx,AFacey,ATz,B,C); break;
    case  6: TensorProductResJFaceK6(M,N,ATx,AFacey,ATz,B,C); break;
    case  7: TensorProductResJFaceK7(M,N,ATx,AFacey,ATz,B,C); break;
    case  8: TensorProductResJFaceK8(M,N,ATx,AFacey,ATz,B,C); break;
    case  9: TensorProductResJFaceK9(M,N,ATx,AFacey,ATz,B,C); break;
    case 10: TensorProductResJFaceK10(M,N,ATx,AFacey,ATz,B,C); break;
    default:
      TerminateAll("TensorProductResJFace", __FILE__, __LINE__,
                   "Value of K not considered");
  }
}
