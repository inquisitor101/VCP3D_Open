//------------------------------------------------------------------------------
// TensorProductSolAndGradVolume.cpp: Implementation of the functions that carry
//                                    out the tensor product to compute the
//                                    solution and possibly gradients in the
//                                    integration points of the volume.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"

//------------------------------------------------------------------------------

// Function, which performs the tensor product multiplications to compute both
// the solution and its gradients in the volume. The C tensors are rank 4 tensors
// of dimension C[M][M][M][N] stored in column major order, B is a rank 4 tensor
// of dimension B[K][K][K][N] stored in column major order, and A and ADer are
// rank 2 tensors (i.e. matrices) of dimension A[M][K] in column major order
// where the first dimension is padded to a multiple of the required vector
// length. The tensor products carried out are
// C     = A*(A*(A*B)),
// CDerX = ADer*(A*(A*B))
// CDerY = A*(ADer*(A*B))
// CDerZ = A*(A*(ADer*B))
//
// By simultaneously computing the solution and the gradients, about 1/4 of the
// operations are saved.
void TensorProductSolAndGradVolume(const int M,
                                   const int N,
                                   const int K,
                                   su2double *A,
                                   su2double *ADer,
                                   su2double **B,
                                   su2double **C,
                                   su2double **CDerX,
                                   su2double **CDerY,
                                   su2double **CDerZ)
{
  // In order to obtain the best performance the sizes of M and K must be known at
  // compile time. Therefore determine the size of these variables and call the
  // corresponding function that carries out the actual multiplication.
  switch( K )
  {
    case  2: TensorProductSolAndGradVolumeK2(M,N,A,ADer,B,C,CDerX,CDerY,CDerZ); break;
    case  3: TensorProductSolAndGradVolumeK3(M,N,A,ADer,B,C,CDerX,CDerY,CDerZ); break;
    case  4: TensorProductSolAndGradVolumeK4(M,N,A,ADer,B,C,CDerX,CDerY,CDerZ); break;
    case  5: TensorProductSolAndGradVolumeK5(M,N,A,ADer,B,C,CDerX,CDerY,CDerZ); break;
    case  6: TensorProductSolAndGradVolumeK6(M,N,A,ADer,B,C,CDerX,CDerY,CDerZ); break;
    case  7: TensorProductSolAndGradVolumeK7(M,N,A,ADer,B,C,CDerX,CDerY,CDerZ); break;
    case  8: TensorProductSolAndGradVolumeK8(M,N,A,ADer,B,C,CDerX,CDerY,CDerZ); break;
    case  9: TensorProductSolAndGradVolumeK9(M,N,A,ADer,B,C,CDerX,CDerY,CDerZ); break;
    case 10: TensorProductSolAndGradVolumeK10(M,N,A,ADer,B,C,CDerX,CDerY,CDerZ); break;
    default:
      TerminateAll("TensorProductSolAndGradVolume", __FILE__, __LINE__,
                   "Value of K not considered");
  }
}
