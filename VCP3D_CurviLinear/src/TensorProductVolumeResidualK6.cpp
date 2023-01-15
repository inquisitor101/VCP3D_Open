//------------------------------------------------------------------------------
// TensorProductVolumeResidualK6.cpp: Implementation of the functions that carries
//                                    out the tensor product multiplications to
//                                    compute the volume residuals.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"
#include "TensorProductVolumeResidual.hpp"

//------------------------------------------------------------------------------

// Function, which performs the tensor product multiplication to accumulate the
// volume residual. The C tensor is a rank 4 tensor of dimension C[K][K][K][N]
// stored in column major order, B is a rank 4 tensor of dimension B[M][M][M][N]
// stored in column major order, and ATx, ATy and ATz are rank 2 tensors (i.e.
// matrices) of dimension AT[K][M] in column major order where the first
// dimension is padded to a multiple of the required vector length.
// The tensor product carried out is
// C += ATz*(ATy*(ATx*B))
void TensorProductVolumeResidualK6(const int M,
                                   const int N,
                                   su2double *ATx,
                                   su2double *ATy,
                                   su2double *ATz,
                                   su2double **B,
                                   su2double **C)
{
  switch( M )
  {
    case  8: TPAVR<6, 8>(N,ATx,ATy,ATz,B,C); break;
    default: TerminateTensorProduct("TensorProductVolumeResidual", 6, M);
  }
}
