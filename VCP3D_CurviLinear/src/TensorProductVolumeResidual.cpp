//------------------------------------------------------------------------------
// TensorProductVolumeResidual.cpp: Implementation of the functions that carry
//                                  out the tensor product multiplications to
//                                  compute the volume residuals.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"

//------------------------------------------------------------------------------

// Function, which performs the tensor product multiplication to accumulate the
// volume residual. The C tensor is a rank 4 tensor of dimension C[K][K][K][N]
// stored in column major order, B is a rank 4 tensor of dimension B[M][M][M][N]
// stored in column major order, and ATx, ATy and ATz are rank 2 tensors (i.e.
// matrices) of dimension AT[K][M] in column major order where the first
// dimension is padded to a multiple of the required vector length.
// The tensor product carried out is
// C += ATz*(ATy*(ATx*B))
void TensorProductVolumeResidual(const int M,
                                 const int N,
                                 const int K,
                                 su2double *ATx,
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
    case  2: TensorProductVolumeResidualK2(M,N,ATx,ATy,ATz,B,C); break;
    case  3: TensorProductVolumeResidualK3(M,N,ATx,ATy,ATz,B,C); break;
    case  4: TensorProductVolumeResidualK4(M,N,ATx,ATy,ATz,B,C); break;
    case  5: TensorProductVolumeResidualK5(M,N,ATx,ATy,ATz,B,C); break;
    case  6: TensorProductVolumeResidualK6(M,N,ATx,ATy,ATz,B,C); break;
    case  7: TensorProductVolumeResidualK7(M,N,ATx,ATy,ATz,B,C); break;
    case  8: TensorProductVolumeResidualK8(M,N,ATx,ATy,ATz,B,C); break;
    case  9: TensorProductVolumeResidualK9(M,N,ATx,ATy,ATz,B,C); break;
    case 10: TensorProductVolumeResidualK10(M,N,ATx,ATy,ATz,B,C); break;
    default:
      TerminateAll("TensorProductVolumeResidual", __FILE__, __LINE__,
                   "Value of K not considered");
  }
}
