//------------------------------------------------------------------------------
// TensorProductSolAndGradIFaceK5.cpp: Implementation of the function that carries
//                                     out the tensor product multiplications to
//                                     compute the solution and gradients in the
//                                     integration points of an I-face.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"
#include "TensorProductSolAndGradIFace.hpp"

//------------------------------------------------------------------------------

// Function, which performs the tensor product multiplications to compute both
// the solution and its gradients on a face in i-direction. The C tensors are
// rank 3 tensors of dimension C[M][M][N] stored in column major order, B is a
// rank 4 tensor of dimension B[K][K][K][N] stored in column major order, 
// A and ADer are rank 2 tensors (i.e. matrices) of dimension A[M][K] in column
// major order and AFace and ADerFace are rank 1 tensor of dimension AFace[K].
// The first dimension of the A tensors is padded to a multiple of the required
// vector length. The tensor products carried out are
// C     = A*(A*(AFace*B)),
// CDerX = A*(A*(ADerFace*B))
// CDerY = A*(ADer*(AFace*B))
// CDerZ = ADer*(A*(AFace*B))
void TensorProductSolAndGradIFaceK5(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ)
{
  // In order to obtain the best performance the sizes of M and K must be known at
  // compile time. Therefore determine the size of these variables and call the
  // corresponding template instantiation of the function that carries out the
  // actual multiplication.
  switch( M )
  {
    case  7: TPSAGI<5, 7>(N,A,ADer,AFace,ADerFace,B,C,CDerX,CDerY,CDerZ); break;
    default: TerminateTensorProduct("TensorProductSolAndGradIFace", 5, M);
  }
}
