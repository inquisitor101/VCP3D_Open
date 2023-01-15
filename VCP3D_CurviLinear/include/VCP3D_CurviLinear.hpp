//------------------------------------------------------------------------------
// VCP3D_CurviLinear.hpp: Main include file of the program VCP3D_CurviLinear.
//                        It is set up such that the source files only need to
//                        include this file. This may not be the most efficient
//                        way of getting things compiled, but it is certainly
//                        the easiest way.
//------------------------------------------------------------------------------

#pragma once

#define _USE_MATH_DEFINES

// Include the mostly needed C++ include files.
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cctype>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <algorithm>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <random>
#include <cassert>

// Include file for the signal handling.
#ifndef USE_NO_SIGNALS
#include <signal.h>
#endif

// Include file for runtime NaN catching.
#ifdef ENABLE_NAN_CHECK
#include <fenv.h>
#endif

// MPI include file, if needed.
#ifdef HAVE_MPI
#include <mpi.h>
#endif

// OpenMP include file, if needed.
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

// MKL include file, if needed.
#ifdef HAVE_MKL
#include "mkl.h"
#endif

//------------------------------------------------------------------------------
//            Parameters, which define the floating point accuracy.
//------------------------------------------------------------------------------

// Precision used in the code.
#ifdef USE_SINGLE_PRECISION  // Use single precision.
typedef float  su2double;

#ifdef HAVE_MPI
const MPI_Datatype MPI_SU2DOUBLE = MPI_FLOAT;
#endif

#else  // Use double precision.
typedef double su2double;

#ifdef HAVE_MPI
const MPI_Datatype MPI_SU2DOUBLE = MPI_DOUBLE;
#endif

#endif

//------------------------------------------------------------------------------
//                     Global constants.
//------------------------------------------------------------------------------

// Floating point values and functions.
// Make a distinction between the precisions used.
#ifdef USE_SINGLE_PRECISION

const su2double zero      =  0.0f;
const su2double one       =  1.0f;
const su2double two       =  2.0f;
const su2double three     =  3.0f;
const su2double four      =  4.0f;
const su2double five      =  5.0f;
const su2double six       =  6.0f;
const su2double seven     =  7.0f;
const su2double eight     =  8.0f;
const su2double nine      =  9.0f;
const su2double eleven    = 11.0f;
const su2double twentyOne = 21.0f;

const su2double half = 0.5f;

const su2double epsThreshold         = 1.e-5f;
const su2double convergenceThreshold = 1.e-6f;

const su2double epsSmall = 1.e-30f;
const su2double valLarge = 1.e+30f;

#define SQRT  sqrtf
#define POW   powf
#define SIN   sinf
#define COS   cosf
#define TAN   tanf
#define SINH  sinhf
#define COSH  coshf
#define TANH  tanhf
#define FABS  fabsf
#define LOG   logf
#define EXP   expf

#else   // Use double precision.

const su2double zero      =  0.0;
const su2double one       =  1.0;
const su2double two       =  2.0;
const su2double three     =  3.0;
const su2double four      =  4.0;
const su2double five      =  5.0;
const su2double six       =  6.0;
const su2double seven     =  7.0;
const su2double eight     =  8.0;
const su2double nine      =  9.0;
const su2double eleven    = 11.0;
const su2double twentyOne = 21.0;

const su2double half = 0.5;

const su2double epsThreshold         = 1.e-10;
const su2double convergenceThreshold = 1.e-12;

const su2double epsSmall = 1.e-30;
const su2double valLarge = 1.e+30;

#define SQRT  sqrt
#define POW   pow
#define SIN   sin
#define COS   cos
#define TAN   tan
#define SINH  sinh
#define COSH  cosh
#define TANH  tanh
#define FABS  fabs
#define LOG   log
#define EXP   exp

#endif

// Byte alignment in the memory allocation.
const size_t byteAlignment = 32;

// The vector length used in one space dimension must be a multiple of vecLen1D.
const size_t vecLen1D = 8;

// The vector length used in two space dimensions must be a multiple of vecLen2D.
const size_t vecLen2D = 8;

// The vector length used in three space dimensions must be a multiple of vecLen3D.
const size_t vecLen3D = 8;

// Number of conservative variables.
const int nVar = 5;

// Specific heat ratio, Prandtl number, Gas constant and specific
// heat at constant volume.
const su2double GamConstant  = (su2double) 1.4;
const su2double Prandtl_Lam  = (su2double) 0.72;
const su2double Prandtl_Turb = (su2double) 0.90;
const su2double RGas         = (su2double) 287.058;
const su2double Cv           = RGas/(GamConstant-one);
const su2double Cp           = GamConstant*Cv;
const su2double Karman       = (su2double) 0.41;

// Constant factor present in the heat flux vector and the ratio of
// the second viscosity and the viscosity itself.
const su2double factHeatFlux_Lam  =  GamConstant/Prandtl_Lam;
const su2double factHeatFlux_Turb =  GamConstant/Prandtl_Turb;
const su2double lambdaOverMu      = -two/three;

// Length of the strings to store the variables names in the restart files.
const int CGNS_STRING_SIZE = 33;

// Magic number in the restart file when using the SU2 format.
const int SU2_MagicNumber = 535532;

// VTK element type definitions.
const int LINE          =  3;
const int TRIANGLE      =  5;
const int QUADRILATERAL =  9;
const int TETRAHEDRON   = 10;
const int HEXAHEDRON    = 12;
const int PRISM         = 13;
const int PYRAMID       = 14;

//------------------------------------------------------------------------------
//                 Enumerated types.
//------------------------------------------------------------------------------

enum ENUM_DOF_LOCATION {
  NO_LOCATION = 0,
  EQUIDISTANT = 1,
  LGL_POINTS  = 2
};

enum ENUM_FEM_VARIABLES {
  NO_VARIABLES           = 0,
  CONSERVATIVE_VARIABLES = 1,
  ENTROPY_VARIABLES      = 2
};

enum ENUM_RIEMANN {
  NO_RIEMANN = 0,
  ROE        = 1,
  ISMAIL_ROE = 2
};

enum ENUM_BOUNDARY_CONDITIONS {
  INTERNAL_1TO1             =  0,
  PERIODIC_1TO1_TRANS       =  1,
  PERIODIC_1TO1_ROT         =  2,
  BC_FARFIELD               =  3,
  BC_ISOTHERMAL_WALL        =  4,
  BC_HEATFLUX_WALL          =  5,
  BC_INVISCID_WALL          =  6,
  BC_SYMMETRY               =  7,
  BC_INFLOW_SUBSONIC        =  8,
  BC_INFLOW_SUPERSONIC      =  9,
  BC_OUTFLOW_SUBSONIC       = 10,
  BC_OUTFLOW_SUPERSONIC     = 11,
	BC_INFLOW_CHARACTERISTIC  = 12,
	BC_OUTFLOW_CHARACTERISTIC = 13,
	BC_UNDEFINED_TYPE         = 99
};

enum ENUM_WALL_MODEL {
  NO_WALL_MODEL          = 0,
  EQUILIBRIUM_WALL_MODEL = 1,
  LOGARITHMIC_WALL_MODEL = 2
};

enum ENUM_SGS_MODEL {
  NO_SGS_MODEL = 0,
  WALE_MODEL   = 1,
  VREMAN_MODEL = 2
};

enum ENUM_FLUCTUATIONS {
  NO_FLUCTUATIONS = 0,
  RANDOM          = 1,
  SINE_AND_COSINE = 2,
  SYNTHETIC       = 3,
	PRESSURE_PULSE  = 4
};

enum ENUM_ELEM_TYPE {
  INTERNAL_ELEMENT  = 0,
  HALO_ELEMENT_KMIN = 1,
  HALO_ELEMENT_JMIN = 2,
  HALO_ELEMENT_IMIN = 3,
  HALO_ELEMENT_IMAX = 4,
  HALO_ELEMENT_JMAX = 5,
  HALO_ELEMENT_KMAX = 6
};

enum ENUM_HEX_FACE {
  KMIN = 0,
  JMIN = 1,
  IMIN = 2,
  IMAX = 3,
  JMAX = 4,
  KMAX = 5
};

enum ENUM_SIGNAL {
  NoSignal        = 0,
  SignalWrite     = 1,
  SignalWriteQuit = 2
};

//------------------------------------------------------------------------------
//                     Global variables.
//------------------------------------------------------------------------------

// Number of MPI ranks and the current MPI rank.
extern int rank, nRanks;

// The current value of the signal for signal handling.
extern ENUM_SIGNAL CurrentSignal;

// Dimensional values of the density and pressure.
extern su2double rhoDim;
extern su2double pDim;

// Reference values to make the equations dimensionless.
extern su2double pRef;
extern su2double uRef;
extern su2double rhoRef;

// The non-dimensional value of the viscosity.
// Note that the reference length is missing in this non-dimensionalization,
// because it is taken to be 1.0.
extern su2double mu;

//------------------------------------------------------------------------------
//                 Prototypes of the stand alone functions.
//------------------------------------------------------------------------------

void *AllocateMemory(const size_t sizeAlloc);

void CreateLowerCase(std::string &lineBuf);

void FreeMemory(void **memPointer);

double GetWallTime(void);

su2double GradNormJacobi(int       n,
                         int       alpha,
                         int       beta,
                         su2double x);

void GradVandermonde1D(const int                    nDOFs,
                       const std::vector<su2double> &r,
                             std::vector<su2double> &VDr);

void InverseMatrix(int                    nn,
                   std::vector<su2double> &A);

void LocationDOFs1D(const int               nDOFs1D,
                    const ENUM_DOF_LOCATION DOFLocation,
                    std::vector<su2double>  &rDOFs1D);

su2double NormJacobi(int       n,
                     int       alpha,
                     int       beta,
                     su2double x);

void RemoveBlanks(std::string &lineBuf);

void ReplaceTabsAndReturns(std::string &lineBuf);

void StartUpCode(int  *argc,
                 char ***argv);

void SwapBytes(void   *buffer,
               size_t nBytes,
               int    nItems);

void Vandermonde1D(const int                    nDOFs,
                   const std::vector<su2double> &r,
                         std::vector<su2double> &V);

void Terminate(const char        *functionName,
               const char        *fileName,
               const int         lineNumber,
               const std::string &errorMessageIn);

void TerminateAll(const char        *functionName,
                  const char        *fileName,
                  const int         lineNumber,
                  const std::string &errorMessageIn);

void TerminateTensorProduct(const char *functionName,
                            const int  K,
                            const int  M);

#ifdef HAVE_MPI
MPI_Datatype DetermineMPI_Type(const float dummy);

MPI_Datatype DetermineMPI_Type(const double dummy);
#endif

#ifndef USE_NO_SIGNALS
void ConnectSignals(void);

void SignalHandler(int sig);
#endif

//------------------------------------------------------------------------------
//                  Prototypes for TensorProductResIFace
//------------------------------------------------------------------------------

void TensorProductResIFace(const int M,
                           const int N,
                           const int K,
                           su2double *AFacex,
                           su2double *ATy,
                           su2double *ATz,
                           su2double **B,
                           su2double **C);

void TensorProductResIFaceK2(const int M,
                             const int N,
                             su2double *AFacex,
                             su2double *ATy,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResIFaceK3(const int M,
                             const int N,
                             su2double *AFacex,
                             su2double *ATy,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResIFaceK4(const int M,
                             const int N,
                             su2double *AFacex,
                             su2double *ATy,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResIFaceK5(const int M,
                             const int N,
                             su2double *AFacex,
                             su2double *ATy,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResIFaceK6(const int M,
                             const int N,
                             su2double *AFacex,
                             su2double *ATy,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResIFaceK7(const int M,
                             const int N,
                             su2double *AFacex,
                             su2double *ATy,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResIFaceK8(const int M,
                             const int N,
                             su2double *AFacex,
                             su2double *ATy,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResIFaceK9(const int M,
                             const int N,
                             su2double *AFacex,
                             su2double *ATy,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResIFaceK10(const int M,
                              const int N,
                              su2double *AFacex,
                              su2double *ATy,
                              su2double *ATz,
                              su2double **B,
                              su2double **C);

//------------------------------------------------------------------------------
//                  Prototypes for TensorProductResJFace
//------------------------------------------------------------------------------

void TensorProductResJFace(const int M,
                           const int N,
                           const int K,
                           su2double *ATx,
                           su2double *AFacey,
                           su2double *ATz,
                           su2double **B,
                           su2double **C);

void TensorProductResJFaceK2(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *AFacey,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResJFaceK3(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *AFacey,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResJFaceK4(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *AFacey,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResJFaceK5(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *AFacey,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResJFaceK6(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *AFacey,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResJFaceK7(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *AFacey,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResJFaceK8(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *AFacey,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResJFaceK9(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *AFacey,
                             su2double *ATz,
                             su2double **B,
                             su2double **C);

void TensorProductResJFaceK10(const int M,
                              const int N,
                              su2double *ATx,
                              su2double *AFacey,
                              su2double *ATz,
                              su2double **B,
                              su2double **C);

//------------------------------------------------------------------------------
//                  Prototypes for TensorProductResKFace
//------------------------------------------------------------------------------

void TensorProductResKFace(const int M,
                           const int N,
                           const int K,
                           su2double *ATx,
                           su2double *ATy,
                           su2double *AFacez,
                           su2double **B,
                           su2double **C);

void TensorProductResKFaceK2(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *ATy,
                             su2double *AFacez,
                             su2double **B,
                             su2double **C);

void TensorProductResKFaceK3(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *ATy,
                             su2double *AFacez,
                             su2double **B,
                             su2double **C);

void TensorProductResKFaceK4(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *ATy,
                             su2double *AFacez,
                             su2double **B,
                             su2double **C);

void TensorProductResKFaceK5(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *ATy,
                             su2double *AFacez,
                             su2double **B,
                             su2double **C);

void TensorProductResKFaceK6(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *ATy,
                             su2double *AFacez,
                             su2double **B,
                             su2double **C);

void TensorProductResKFaceK7(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *ATy,
                             su2double *AFacez,
                             su2double **B,
                             su2double **C);

void TensorProductResKFaceK8(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *ATy,
                             su2double *AFacez,
                             su2double **B,
                             su2double **C);

void TensorProductResKFaceK9(const int M,
                             const int N,
                             su2double *ATx,
                             su2double *ATy,
                             su2double *AFacez,
                             su2double **B,
                             su2double **C);

void TensorProductResKFaceK10(const int M,
                              const int N,
                              su2double *ATx,
                              su2double *ATy,
                              su2double *AFacez,
                              su2double **B,
                              su2double **C);

//------------------------------------------------------------------------------
//                  Prototypes for TensorProductSolAndGradIFace.
//------------------------------------------------------------------------------

void TensorProductSolAndGradIFace(const int M,
                                  const int N,
                                  const int K,
                                  su2double *A,
                                  su2double *ADer,
                                  su2double *AFace,
                                  su2double *ADerFace,
                                  su2double **B,
                                  su2double **C,
                                  su2double **CDerX,
                                  su2double **CDerY,
                                  su2double **CDerZ);

void TensorProductSolAndGradIFaceK2(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradIFaceK3(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradIFaceK4(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

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
                                    su2double **CDerZ);

void TensorProductSolAndGradIFaceK6(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradIFaceK7(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradIFaceK8(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradIFaceK9(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradIFaceK10(const int M,
                                     const int N,
                                     su2double *A,
                                     su2double *ADer,
                                     su2double *AFace,
                                     su2double *ADerFace,
                                     su2double **B,
                                     su2double **C,
                                     su2double **CDerX,
                                     su2double **CDerY,
                                     su2double **CDerZ);

//------------------------------------------------------------------------------
//                  Prototypes for TensorProductSolAndGradJFace.
//------------------------------------------------------------------------------

void TensorProductSolAndGradJFace(const int M,
                                  const int N,
                                  const int K,
                                  su2double *A,
                                  su2double *ADer,
                                  su2double *AFace,
                                  su2double *ADerFace,
                                  su2double **B,
                                  su2double **C,
                                  su2double **CDerX,
                                  su2double **CDerY,
                                  su2double **CDerZ);

void TensorProductSolAndGradJFaceK2(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradJFaceK3(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradJFaceK4(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradJFaceK5(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradJFaceK6(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradJFaceK7(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradJFaceK8(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradJFaceK9(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradJFaceK10(const int M,
                                     const int N,
                                     su2double *A,
                                     su2double *ADer,
                                     su2double *AFace,
                                     su2double *ADerFace,
                                     su2double **B,
                                     su2double **C,
                                     su2double **CDerX,
                                     su2double **CDerY,
                                     su2double **CDerZ);

//------------------------------------------------------------------------------
//                  Prototypes for TensorProductSolAndGradKFace.
//------------------------------------------------------------------------------

void TensorProductSolAndGradKFace(const int M,
                                  const int N,
                                  const int K,
                                  su2double *A,
                                  su2double *ADer,
                                  su2double *AFace,
                                  su2double *ADerFace,
                                  su2double **B,
                                  su2double **C,
                                  su2double **CDerX,
                                  su2double **CDerY,
                                  su2double **CDerZ);

void TensorProductSolAndGradKFaceK2(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradKFaceK3(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradKFaceK4(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradKFaceK5(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradKFaceK6(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradKFaceK7(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradKFaceK8(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradKFaceK9(const int M,
                                    const int N,
                                    su2double *A,
                                    su2double *ADer,
                                    su2double *AFace,
                                    su2double *ADerFace,
                                    su2double **B,
                                    su2double **C,
                                    su2double **CDerX,
                                    su2double **CDerY,
                                    su2double **CDerZ);

void TensorProductSolAndGradKFaceK10(const int M,
                                     const int N,
                                     su2double *A,
                                     su2double *ADer,
                                     su2double *AFace,
                                     su2double *ADerFace,
                                     su2double **B,
                                     su2double **C,
                                     su2double **CDerX,
                                     su2double **CDerY,
                                     su2double **CDerZ);

//------------------------------------------------------------------------------
//                  Prototypes for TensorProductSolAndGradVolume
//------------------------------------------------------------------------------

void TensorProductSolAndGradVolume(const int M,
                                   const int N,
                                   const int K,
                                   su2double *A,
                                   su2double *ADer,
                                   su2double **B,
                                   su2double **C,
                                   su2double **CDerX,
                                   su2double **CDerY,
                                   su2double **CDerZ);

void TensorProductSolAndGradVolumeK2(const int M,
                                     const int N,
                                     su2double *A,
                                     su2double *ADer,
                                     su2double **B,
                                     su2double **C,
                                     su2double **CDerX,
                                     su2double **CDerY,
                                     su2double **CDerZ);

void TensorProductSolAndGradVolumeK3(const int M,
                                     const int N,
                                     su2double *A,
                                     su2double *ADer,
                                     su2double **B,
                                     su2double **C,
                                     su2double **CDerX,
                                     su2double **CDerY,
                                     su2double **CDerZ);

void TensorProductSolAndGradVolumeK4(const int M,
                                     const int N,
                                     su2double *A,
                                     su2double *ADer,
                                     su2double **B,
                                     su2double **C,
                                     su2double **CDerX,
                                     su2double **CDerY,
                                     su2double **CDerZ);

void TensorProductSolAndGradVolumeK5(const int M,
                                     const int N,
                                     su2double *A,
                                     su2double *ADer,
                                     su2double **B,
                                     su2double **C,
                                     su2double **CDerX,
                                     su2double **CDerY,
                                     su2double **CDerZ);

void TensorProductSolAndGradVolumeK6(const int M,
                                     const int N,
                                     su2double *A,
                                     su2double *ADer,
                                     su2double **B,
                                     su2double **C,
                                     su2double **CDerX,
                                     su2double **CDerY,
                                     su2double **CDerZ);

void TensorProductSolAndGradVolumeK7(const int M,
                                     const int N,
                                     su2double *A,
                                     su2double *ADer,
                                     su2double **B,
                                     su2double **C,
                                     su2double **CDerX,
                                     su2double **CDerY,
                                     su2double **CDerZ);

void TensorProductSolAndGradVolumeK8(const int M,
                                     const int N,
                                     su2double *A,
                                     su2double *ADer,
                                     su2double **B,
                                     su2double **C,
                                     su2double **CDerX,
                                     su2double **CDerY,
                                     su2double **CDerZ);

void TensorProductSolAndGradVolumeK9(const int M,
                                     const int N,
                                     su2double *A,
                                     su2double *ADer,
                                     su2double **B,
                                     su2double **C,
                                     su2double **CDerX,
                                     su2double **CDerY,
                                     su2double **CDerZ);

void TensorProductSolAndGradVolumeK10(const int M,
                                      const int N,
                                      su2double *A,
                                      su2double *ADer,
                                      su2double **B,
                                      su2double **C,
                                      su2double **CDerX,
                                      su2double **CDerY,
                                      su2double **CDerZ);

//------------------------------------------------------------------------------
//                  Prototypes for TensorProductVolumeResidual.
//------------------------------------------------------------------------------

void TensorProductVolumeResidual(const int M,
                                 const int N,
                                 const int K,
                                 su2double *ATx,
                                 su2double *ATy,
                                 su2double *ATz,
                                 su2double **B,
                                 su2double **C);

void TensorProductVolumeResidualK2(const int M,
                                   const int N,
                                   su2double *ATx,
                                   su2double *ATy,
                                   su2double *ATz,
                                   su2double **B,
                                   su2double **C);

void TensorProductVolumeResidualK3(const int M,
                                   const int N,
                                   su2double *ATx,
                                   su2double *ATy,
                                   su2double *ATz,
                                   su2double **B,
                                   su2double **C);

void TensorProductVolumeResidualK4(const int M,
                                   const int N,
                                   su2double *ATx,
                                   su2double *ATy,
                                   su2double *ATz,
                                   su2double **B,
                                   su2double **C);

void TensorProductVolumeResidualK5(const int M,
                                   const int N,
                                   su2double *ATx,
                                   su2double *ATy,
                                   su2double *ATz,
                                   su2double **B,
                                   su2double **C);

void TensorProductVolumeResidualK6(const int M,
                                   const int N,
                                   su2double *ATx,
                                   su2double *ATy,
                                   su2double *ATz,
                                   su2double **B,
                                   su2double **C);

void TensorProductVolumeResidualK7(const int M,
                                   const int N,
                                   su2double *ATx,
                                   su2double *ATy,
                                   su2double *ATz,
                                   su2double **B,
                                   su2double **C);

void TensorProductVolumeResidualK8(const int M,
                                   const int N,
                                   su2double *ATx,
                                   su2double *ATy,
                                   su2double *ATz,
                                   su2double **B,
                                   su2double **C);

void TensorProductVolumeResidualK9(const int M,
                                   const int N,
                                   su2double *ATx,
                                   su2double *ATy,
                                   su2double *ATz,
                                   su2double **B,
                                   su2double **C);

void TensorProductVolumeResidualK10(const int M,
                                    const int N,
                                    su2double *ATx,
                                    su2double *ATy,
                                    su2double *ATz,
                                    su2double **B,
                                    su2double **C);

//------------------------------------------------------------------------------
//                            Class definitions.
//------------------------------------------------------------------------------

#include "adt_structure.hpp"
#include "InputParamClass.hpp"
#include "SolverClass.hpp"
#include "GaussJacobiQuadrature.hpp"
#include "SGSModel.hpp"
#include "WallModel.hpp"

//------------------------------------------------------------------------------
// Prototype of the functions that require some of the above classes.
//------------------------------------------------------------------------------

void CG_ConservativeVarMassMatrix(const StandardElementClass *standardHex,
                                  su2double                  **sol,
                                  const su2double            *Jac,
                                  const su2double            volElem,
                                  su2double                  **res,
                                  su2double                  **workArray);

void CG_EntropyVarMassMatrix(const StandardElementClass *standardHex,
                             su2double                  **sol,
                             const su2double            *Jac,
                             const su2double            volElem,
                             su2double                  **dUdV,
                             su2double                  **res,
                             su2double                  **workArray);

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
                        su2double                                                *forceCoef);

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
                        su2double             &EddyVisMax);

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
                        su2double            **pInt);

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
                            su2double             **fluxSymZ);

void PenAndSymFluxesSymmVar(const InputParamClass *inputParam,
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
                            su2double             **fluxSymZ);

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
                       su2double             *forceCoef);

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
                                su2double                                                *forceCoef);
