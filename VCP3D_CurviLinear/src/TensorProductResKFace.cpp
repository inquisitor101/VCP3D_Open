//------------------------------------------------------------------------------
// TensorProductResKFace.cpp: Implementation of the functions that carry out the
//                            tensor product multiplications to compute the
//                            contribution from a K face to the residual.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"

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
void TensorProductResKFace(const int M,
                           const int N,
                           const int K,
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
  switch( K )
  {
    case  2: TensorProductResKFaceK2(M,N,ATx,ATy,AFacez,B,C); break;
    case  3: TensorProductResKFaceK3(M,N,ATx,ATy,AFacez,B,C); break;
    case  4: TensorProductResKFaceK4(M,N,ATx,ATy,AFacez,B,C); break;
    case  5: TensorProductResKFaceK5(M,N,ATx,ATy,AFacez,B,C); break;
    case  6: TensorProductResKFaceK6(M,N,ATx,ATy,AFacez,B,C); break;
    case  7: TensorProductResKFaceK7(M,N,ATx,ATy,AFacez,B,C); break;
    case  8: TensorProductResKFaceK8(M,N,ATx,ATy,AFacez,B,C); break;
    case  9: TensorProductResKFaceK9(M,N,ATx,ATy,AFacez,B,C); break;
    case 10: TensorProductResKFaceK10(M,N,ATx,ATy,AFacez,B,C); break;
    default:
      TerminateAll("TensorProductResKFace", __FILE__, __LINE__,
                   "Value of K not considered");
  }
}
