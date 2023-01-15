//------------------------------------------------------------------------------
// StandAloneFunctions.cpp: Implementation of the stand alone functions.
//                          These are relatively small functions and therefore
//                          they are put in one file.
//------------------------------------------------------------------------------

// General include file.
#include "VCP3D_CurviLinear.hpp"

// Function prototypes for the Lapack routines used.
// Not needed when the MKL is used.
#ifndef HAVE_MKL
extern "C"
{
#ifdef USE_SINGLE_PRECISION  // Use single precision

  // LU decomoposition of a general matrix
  void sgetrf_(int* M, int *N, float* A, int* lda, int* IPIV, int* INFO);

  // Generate inverse of a matrix given its LU decomposition
  void sgetri_(int* N, float* A, int* lda, int* IPIV, float* WORK,
               int* lwork, int* INFO);

#else // Use double precision

  // LU decomoposition of a general matrix
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);

  // Generate inverse of a matrix given its LU decomposition
  void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK,
               int* lwork, int* INFO);
#endif
}
#endif

//------------------------------------------------------------------------------
//               Implementation of stand alone functions.
//------------------------------------------------------------------------------

// Function, which carries out the aligned memory allocation.
void *AllocateMemory(const size_t sizeAlloc)
{
#ifdef HAVE_MKL
  return mkl_malloc(sizeAlloc, byteAlignment);
#else
  return aligned_alloc(byteAlignment, sizeAlloc);
#endif
}

//------------------------------------------------------------------------------

// Function, which converts a given string to lower case.
void CreateLowerCase(std::string &lineBuf)
{
  // Determine the length of the string and convert its elements to
  // lower case.
  std::string::size_type strLength = lineBuf.length();
  for(std::string::size_type i=0; i<strLength; ++i)
    lineBuf[i] = (char) tolower(lineBuf[i]);
}

//------------------------------------------------------------------------------

// Function, which frees the, aligned, allocated memory.
void FreeMemory(void **memPointer)
{
  if( *memPointer ) {
#ifdef HAVE_MKL
    mkl_free(*memPointer);
#else
    free(*memPointer);
#endif
    *memPointer = NULL;
  }
}

//------------------------------------------------------------------------------

// Function, which returns the wall clock time of the application.
double GetWallTime(void)
{
#ifdef HAVE_OPENMP
  return omp_get_wtime();
#elif HAVE_MPI
  return MPI_Wtime();
#else
  return double(clock())/double(CLOCKS_PER_SEC);
#endif
}

//------------------------------------------------------------------------------

// Function, which computes the value of the derivative to x of the
// normalized Jacobi polynomial with coefficients alpha and beta of
// order n for the given value of x.
su2double GradNormJacobi(int       n,
                         int       alpha,
                         int       beta,
                         su2double x)
{
  // Make a distinction for n == 0 and n > 0. For n == 0 the derivative is
  // zero, because the polynomial itself is constant.
  su2double grad;
  if(n == 0) grad = 0.0;
  else
  {
    const su2double tmp = n*(n+alpha+beta+one);
    grad = SQRT(tmp)*NormJacobi(n-1, alpha+1, beta+1, x);
  }

  // Return the gradient.
  return grad;
}

//------------------------------------------------------------------------------

// Function, which computes the gradient of the Vandermonde matrix
// for a standard edge.
void GradVandermonde1D(const int                    nDOFs,
                       const std::vector<su2double> &r,
                             std::vector<su2double> &VDr)
{
  // Determine the number or rows of the gradient of the Vandermonde matrix
  // and check if the dimension of VDr is correct.
  const int nRows = r.size();

  if(VDr.size() != (unsigned int) (nRows*nDOFs))
    Terminate("GradVandermonde1D", __FILE__, __LINE__,
              "Wrong size of the VDr matrix");

  // Compute the gradient of the Vandermonde matrix.
  int ii = 0;
  for(int i=0; i<nDOFs; ++i)
  {
    for(int k=0; k<nRows; ++k, ++ii)
      VDr[ii] = GradNormJacobi(i,0,0,r[k]);
  }
}

//------------------------------------------------------------------------------

// Function, which computes the inverse of the given matrix.
// Double precision version.
void InverseMatrix(int                    nn,
                   std::vector<su2double> &A)
{
  // Allocate the memory for the work arrays.
  std::vector<int>       ipivVec(nn);
  std::vector<su2double> workVec(nn);

  // Determine the pointers for the data in the vectors.
  int *ipiv = ipivVec.data();

  su2double *work = workVec.data();
  su2double *AA   = A.data();

  // Call the appropriate Lapack functions to compute the inverse.
  int ierr;

#ifdef USE_SINGLE_PRECISION
  sgetrf_(&nn, &nn, AA, &nn, ipiv, &ierr);
  if(ierr != 0)
    Terminate("InverseMatrix", __FILE__, __LINE__, "Matrix is singular");

  sgetri_(&nn, AA, &nn, ipiv, work, &nn, &ierr);
  if(ierr != 0)
    Terminate("InverseMatrix", __FILE__, __LINE__, "Matrix inversion failed");
#else
  dgetrf_(&nn, &nn, AA, &nn, ipiv, &ierr);
  if(ierr != 0)
    Terminate("InverseMatrix", __FILE__, __LINE__, "Matrix is singular");

  dgetri_(&nn, AA, &nn, ipiv, work, &nn, &ierr);
  if(ierr != 0)
    Terminate("InverseMatrix", __FILE__, __LINE__, "Matrix inversion failed");
#endif
}

//------------------------------------------------------------------------------

// Function, which determines the 1D location of the DOFs in the standard edge.
void LocationDOFs1D(const int               nDOFs1D,
                    const ENUM_DOF_LOCATION DOFLocation,
                    std::vector<su2double>  &rDOFs1D)
{
  // Allocate the memory for rDOFs1D.
  rDOFs1D.resize(nDOFs1D);

  // Make a distinction between the possibilities for the DOF locations.
  switch( DOFLocation )
  {
    case EQUIDISTANT:
    {
      // Equidistant spacing is used.
      const su2double dh = two/(nDOFs1D-1);

      for(int i=0; i<nDOFs1D; ++i)
        rDOFs1D[i] = -one + i*dh;

      break;
    }

    //--------------------------------------------------------------------------

    case LGL_POINTS:
    {
      // The Gauss Lobatto points are used as parametric coordinates of the
      // DOFs. There are algorithms to compute those locations, but below these
      // points are set explicitly, depending on the number of DOFs. Up till
      // 10 DOFs are considered, which should be enough.
      switch( nDOFs1D )
      {
        case 1:
        {
          rDOFs1D[0] = zero;
          break;
        }

        case 2:
        {
          rDOFs1D[0] = -one; rDOFs1D[1] = one;
          break;
        }

        case 3:
        {
          rDOFs1D[0] = -one; rDOFs1D[1] = zero; rDOFs1D[2] = one;
          break;
        }

        case 4:
        {
          rDOFs1D[0] = -one;           rDOFs1D[1] = -SQRT(one/five);
          rDOFs1D[2] = SQRT(one/five); rDOFs1D[3] = one;
          break;
        }

        case 5:
        {
          rDOFs1D[0] = -one;              rDOFs1D[1] = -SQRT(three/seven); rDOFs1D[2] = zero;
          rDOFs1D[3] = SQRT(three/seven); rDOFs1D[4] = one;
          break;
        }

        case 6:
        {
          const su2double t1 = two*SQRT(seven)/twentyOne;
          const su2double t2 = SQRT(one/three + t1);
          const su2double t3 = SQRT(one/three - t1);
          rDOFs1D[0] = -one; rDOFs1D[1] = -t2; rDOFs1D[2] = -t3;
          rDOFs1D[3] =  t3;  rDOFs1D[4] =  t2; rDOFs1D[5] =  one;
          break;
        }

        case 7:
        {
          const su2double t1 = two*SQRT(five/three)/eleven;
          const su2double t2 = SQRT(five/eleven + t1);
          const su2double t3 = SQRT(five/eleven - t1);
          rDOFs1D[0] = -one; rDOFs1D[1] = -t2; rDOFs1D[2] = -t3;  rDOFs1D[3] = zero;
          rDOFs1D[4] =  t3;  rDOFs1D[5] =  t2; rDOFs1D[6] =  one;
          break;
        }

        case 8:
        {
          rDOFs1D[0] = -one;
          rDOFs1D[1] = (su2double) -0.8717401485096066153375;
          rDOFs1D[2] = (su2double) -0.5917001814331423021445;
          rDOFs1D[3] = (su2double) -0.2092992179024788687687;
          rDOFs1D[4] = (su2double)  0.2092992179024788687687;
          rDOFs1D[5] = (su2double)  0.5917001814331423021445;
          rDOFs1D[6] = (su2double)  0.8717401485096066153375;
          rDOFs1D[7] =  one;
          break;
        }

        case 9:
        {
          rDOFs1D[0] = -one;
          rDOFs1D[1] = (su2double) -0.8997579954114601573124;
          rDOFs1D[2] = (su2double) -0.6771862795107377534459;
          rDOFs1D[3] = (su2double) -0.3631174638261781587108;
          rDOFs1D[4] =  zero;
          rDOFs1D[5] = (su2double)  0.3631174638261781587108;
          rDOFs1D[6] = (su2double)  0.6771862795107377534459;
          rDOFs1D[7] = (su2double)  0.8997579954114601573124;
          rDOFs1D[8] =  one;
          break;
        }

        case 10:
        {
          rDOFs1D[0] = -one;
          rDOFs1D[1] = -0.9195339081664588138289;
          rDOFs1D[2] = -0.7387738651055050750031;
          rDOFs1D[3] = -0.4779249498104444956612;
          rDOFs1D[4] = -0.1652789576663870246262;
          rDOFs1D[5] =  0.1652789576663870246262;
          rDOFs1D[6] =  0.4779249498104444956612;
          rDOFs1D[7] =  0.7387738651055050750031;
          rDOFs1D[8] =  0.9195339081664588138289;
          rDOFs1D[9] =  one;
          break;
        }

        default:
        {
          TerminateAll("LocationDOFs1D", __FILE__, __LINE__,
                       "Number of DOFs not accounted for in LGL points");
        }
      }

      break;
    }

    //--------------------------------------------------------------------------

    default:
      TerminateAll("LocationDOFs1D", __FILE__, __LINE__,
                   "This is to avoid a compiler warning and should not happen");
  }
}

//------------------------------------------------------------------------------

// Function, which computes the value of the normalized Jacobi polynomial
// with coefficients alpha and beta of order n for the given value of x.
su2double NormJacobi(int      n,
                     int      alpha,
                     int      beta,
                     su2double x)
{
  // Some abbreviations.
  const su2double ap1   = (su2double) (alpha + 1);
  const su2double bp1   = (su2double) (beta  + 1);
  const su2double apb   = (su2double) (alpha + beta);
  const su2double apbp1 = (su2double) (apb + 1);
  const su2double apbp2 = (su2double) (apb + 2);
  const su2double apbp3 = (su2double) (apb + 3);
  const su2double b2ma2 = (su2double) (beta*beta - alpha*alpha);

  // Determine the terms, which involves the gamma function. As the
  // arguments are integers, this term can be computed easily, because
  // Gamma(n+1) = n!.
  su2double Gamap1 = one, Gambp1 = one, Gamapbp2 = one;
  for(int i=2; i<=alpha; ++i)          Gamap1   *= i;
  for(int i=2; i<=beta; ++i)           Gambp1   *= i;
  for(int i=2; i<=(alpha+beta+1); ++i) Gamapbp2 *= i;

  // Initialize the normalized polynomials.
  su2double Pnm1 = SQRT(POW(half,apbp1)*Gamapbp2/(Gamap1*Gambp1));
  su2double Pn   = half*Pnm1*(apbp2*x + alpha - beta)*SQRT(apbp3/(ap1*bp1));

  // Take care of the special situation of n == 0.
  if(n == 0) Pn = Pnm1;
  else
  {
    // The value of the normalized Jacobi polynomial must be obtained
    // via recursion.
    for(int i=2; i<=n; ++i)
    {
      // Compute the coefficients a for i and i-1 and the coefficient bi.
      int j = i-1;
      su2double tmp  = 2*j + apb;
      su2double aim1 = two*SQRT(j*(j+apb)*(j+alpha)*(j+beta)/((tmp-one)*(tmp+one)))/tmp;

      su2double bi = b2ma2/(tmp*(tmp+two));

      tmp = 2*i + apb;
      su2double ai = two*SQRT(i*(i+apb)*(i+alpha)*(i+beta)/((tmp-one)*(tmp+one)))/tmp;

      // Compute the new value of Pn and make sure to store Pnm1 correctly.
      tmp  = Pnm1;
      Pnm1 = Pn;

      Pn = ((x-bi)*Pn - aim1*tmp)/ai;
    }
  }

  // Return Pn.
  return Pn;
}

//------------------------------------------------------------------------------

// Function, which removes leading and trailing blanks.
void RemoveBlanks(std::string &lineBuf)
{
  // Find the first non space character.
  int strLength = lineBuf.length();
  int posLeading = 0;
  while(posLeading < strLength && isspace(lineBuf[posLeading])) ++posLeading;

  // Find the last non space character.
  int posTrailing = strLength - 1;
  while(posTrailing >= 0 && isspace(lineBuf[posTrailing])) --posTrailing;

  // Determine the situation.
  if(posLeading == strLength || posTrailing < 0)
  {
    // No non-blanks in the string. Set lineBuf to an empty string.
    lineBuf = "";
  }
  else
  {
    // Non-blanks are present. Remove the blanks. First the trailing ones,
    // because otherwise the positions are screwed up.
    int eraseLenBack = strLength - posTrailing - 1;
    if( eraseLenBack ) lineBuf.erase(posTrailing+1, eraseLenBack);
    lineBuf.erase(0, posLeading);
  }
}

//------------------------------------------------------------------------------

// Function, which replaces the tabs and return characters.
void ReplaceTabsAndReturns(std::string &lineBuf)
{
  // Replace the tabs.
  for(;;)
  {
    std::string::size_type pos = lineBuf.find("\t");
    if(pos == std::string::npos) break;
    lineBuf.replace(pos, 1, " ");
  }

  // Replace the returns.
  for(;;)
  {
    std::string::size_type pos = lineBuf.find("\n");
    if(pos == std::string::npos) break;
    lineBuf.replace(pos, 1, " ");
  }
}

//------------------------------------------------------------------------------

void StartUpCode(int  *argc,
                 char ***argv)
{
  // Initialize the variables for the current MPI rank and the number of MPI ranks.
  rank   = 0;
  nRanks = 1;

#ifdef HAVE_MPI

#ifdef HAVE_OPENMP
  // As this program also uses OpenMP, we would like to have support for
  // MPI_THREAD_MULTIPLE.
  int providedThreadSupport;
  MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &providedThreadSupport);

  // Check if the provided thread support is ok.
  if(providedThreadSupport != MPI_THREAD_MULTIPLE)
  {
    Terminate("StartUpCode", __FILE__, __LINE__,
              "MPI implementation does not support MPI_THREAD_MULTIPLE, needed for hybrid MPI-OpenMP");
  }

  // Determine the number of MPI ranks and my rank.
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

  // Test for rank 0.
  if(rank == 0)
  {
#pragma omp parallel
#pragma omp single nowait
    {
    std::cout << "#" << std::endl;
    std::cout << "# This is a hybrid MPI-OpenMP parallel executable." << std::endl;
    std::cout << "# It is running on " << nRanks << " MPI ranks, "
              <<   "each using " << omp_get_num_threads() << " OpenMP threads." << std::endl;
    std::cout << "#" << std::endl;
    }
  }

#else   // HAVE_OPENMP
  // Only MPI is used for the parallelization. So MPI_Init is fine.
  MPI_Init(argc, argv);

  // Determine the number of MPI ranks and my rank.
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nRanks);

  // Write a message if this is rank 0.
  if(rank == 0)
  {
    std::cout << "#" << std::endl;
    std::cout << "# This is a pure MPI parallel executable running on " << nRanks << " ranks." << std::endl;
    std::cout << "#" << std::endl;
  }

#endif  // HAVE_OPENMP

#else   // HAVE_MPI

#ifdef HAVE_OPENMP

  // Only OpenMP is used for the parallelization. Write a message.
#pragma omp parallel
#pragma omp single nowait
  {
  std::cout << "#" << std::endl;
  std::cout << "# This is a pure OpenMP parallel executable "
            <<   "using " << omp_get_num_threads() << " OpenMP threads." << std::endl;
  std::cout << "#" << std::endl;
  }

#else   // HAVE_OPENMP

  // No parallelization is used. Write a message.
  std::cout << "#" << std::endl;
  std::cout << "# This is a sequential executable." << std::endl;
  std::cout << "#" << std::endl;

#endif  // HAVE_OPENMP

#endif  // HAVE_MPI

  // Write a message about the signals, if supported.
#ifdef USE_NO_SIGNALS
  if(rank == 0)
  {
    std::cout << "# No support for signals." << std::endl;
    std::cout << "#" << std::endl;
  }
#else
  if(rank == 0)
  {
    std::cout << "# The following user signals are supported." << std::endl;
    std::cout << "# kill -USR1 <pid>: Solution file is written after "
              <<   "the current time step." << std::endl;
    std::cout << "# kill -USR2 <pid>: Solution file is written after "
              <<   "the current time step" << std::endl
              << "#                   and the program will stop." << std::endl;
    std::cout << "#" << std::endl;
  }
#endif

  // Enable the runtime NaN catching, if desired.
#ifdef ENABLE_NAN_CHECK
#ifndef __APPLE__
    feenableexcept(FE_INVALID | FE_OVERFLOW);
#else
    feclearexcept(FE_INVALID | FE_OVERFLOW);
#endif
#endif

  // Check if the program is called correctly.
  if((rank == 0) && (*argc != 2))
    Terminate("main", __FILE__, __LINE__,
              "Parameter file not provided as comment line argument");

  // Check if the number of MPI ranks is a power of 2.
  int mRanks = nRanks;
  for(;;)
  {
    if(mRanks == 1) break;
    const int mRanks2 = mRanks/2;
    const int resDiv  = mRanks - 2*mRanks2;
    mRanks = mRanks2;

    if((rank == 0) && (resDiv != 0))
      Terminate("main", __FILE__, __LINE__,
                "Number of MPI ranks must be a power of 2");
  }
}

//------------------------------------------------------------------------------

// Function to swap the bytes for certain primitive types.
void SwapBytes(void   *buffer,
               size_t nBytes,
               int    nItems)
{
  // Store half the number of bytes in kk and cast the buffer
  // to a character buffer.
  char *buf = (char *) buffer;
  const size_t kk = nBytes/2;

  // Loop over the number of items in the buffer.
  for(int j=0; j<nItems; ++j)
  {
    // Initialize ii and jj, which are used to store the
    // indices of the bytes to be swapped.
    size_t ii = j*nBytes;
    size_t jj = ii + nBytes - 1;

    // Swap the bytes.
    for(size_t i=0; i<kk; ++i, ++ii, --jj)
    {
      const char tmp = buf[jj];
      buf[jj] = buf[ii];
      buf[ii] = tmp;
    }
  }
}

//------------------------------------------------------------------------------

// Function, which computes the Vandermonde matrix for a standard edge.
void Vandermonde1D(const int                    nDOFs,
                   const std::vector<su2double> &r,
                         std::vector<su2double> &V)
{
  // Determine the number or rows of the Vandermonde matrix and check
  // if the dimension of V is correct.
  const int nRows = r.size();

  if(V.size() != (unsigned int) (nRows*nDOFs))
    Terminate("Vandermonde1D", __FILE__, __LINE__,
              "Wrong size of the V matrix");

  // Compute the Vandermonde matrix.
  int ii = 0;
  for(int i=0; i<nDOFs; ++i)
  {
    for(int k=0; k<nRows; ++k, ++ii)
      V[ii] = NormJacobi(i,0,0,r[k]);
  }
}

//------------------------------------------------------------------------------

// Error function, which prints an error message and exits.
void Terminate(const char        *functionName,
               const char        *fileName,
               const int         lineNumber,
               const std::string &errorMessageIn)
{
  // Insert a comment character and a space on places where a newline
  // character occurs in the string, such that the error message
  // looks nicer.
  std::string  errorMessage  = errorMessageIn;
  std::string  tmpString     = errorMessage;
  std::string::size_type off = 1;

  for(;;)
  {
    std::string::size_type loc = tmpString.find("\n");
    if(loc == std::string::npos) break;
    errorMessage.insert(loc+off, "# ");
    off += loc+3;
    tmpString.erase(0,loc+1);
  }

  // Header of the error message.
  std::cout << "#" << std::endl
            << "#=========================== !!! Error !!! "
            << "============================" << std::endl
            << "#" << std::endl;

  // Write the function name, file name and line number.
  std::cout << "#* Run-time error in " << functionName << std::endl;
  std::cout << "#* File: " << fileName
            << ", Line: "  << lineNumber << std::endl
            << "#" << std::endl;

  // Write the error message and the terminating line.
  std::cout << "# " << errorMessage << std::endl
            << "#" << std::endl;
  std::cout << "#=========================================="
            << "============================" << std::endl
            << "#" << std::endl << std::flush;

  // And exit.
#ifdef HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif
  exit(1);
}

//------------------------------------------------------------------------------

// Error function, which prints an error message and exits. The error is
// handled by the master MPI rank to avoid multiple messages.
void TerminateAll(const char        *functionName,
                  const char        *fileName,
                  const int         lineNumber,
                  const std::string &errorMessageIn)
{
  // The master rank handles the error.
  if(rank == 0)
  {
#ifdef HAVE_OPENMP
#pragma omp single
#endif
    Terminate(functionName, fileName, lineNumber, errorMessageIn);
  }

  // For a MPI executable the remaining ranks wait to get killed.
#ifdef HAVE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

//------------------------------------------------------------------------------

// Function, which serves as an interface from the tensor product functions
// to the error handler. In this way a lot of very similar code is avoided.
void TerminateTensorProduct(const char *functionName,
                            const int  K,
                            const int  M)
{
#ifdef HAVE_OPENMP
#pragma omp single
#endif
  // Rank 0 prints a nice message.
  if(rank == 0)
  {
    std::cout << "#" << std::endl;
    std::cout << "#----------------------------------------------------------" << std::endl;
    std::cout << "# In function " << functionName << "." << std::endl;
    std::cout << "# Template not instantiated for K = " << K << " and M = " << M
              << "." << std::endl;
    std::cout << "#----------------------------------------------------------" << std::endl
              << std::flush;
  }

  // Call TerminateAll to stop the computation.
#ifdef HAVE_OPENMP
  if(omp_get_thread_num() == 0)
#endif
    TerminateAll("TerminateTensorProduct", __FILE__, __LINE__,
                 "Tensor product size not supported yet.");
}

//------------------------------------------------------------------------------

#ifdef HAVE_MPI
// Functions, which return the corresponding MPI datatype for the given argument.
MPI_Datatype DetermineMPI_Type(const float dummy){return MPI_FLOAT;}

MPI_Datatype DetermineMPI_Type(const double dummy){return MPI_DOUBLE;}
#endif

//------------------------------------------------------------------------------

#ifndef USE_NO_SIGNALS
// Function, which connects the user defined signals.
void ConnectSignals(void)
{
  // Connect the kill -USR1 and kill -USR2 to the signal handler function.
  signal(SIGUSR1, SignalHandler);
  signal(SIGUSR2, SignalHandler);
}

//------------------------------------------------------------------------------

// Function, which handles the given kill signal.
void SignalHandler(int sig)
{
  // The signal handling only needs to be carried out by one thread.
#ifdef HAVE_OPENMP
#pragma omp single nowait
#endif
  {
  // The user signals must be reconnected again, because the
  // connection is lost after a signal has been given.
  ConnectSignals();

  // Write a message that a signal was received.
  std::cout << "#" << std::endl
            << "#=========================================="
            << "============================" << std::endl
            << "#" << std::endl;

  if(nRanks > 1) std::cout << "# MPI rank " << rank << " received a";
  else           std::cout << "# Received a";

  switch( sig )
  {
    case SIGUSR1:
      std::cout << " USR1 signal." << std::endl;
      break;

    case SIGUSR2:
      std::cout << " USR2 signal." << std::endl;
      break;

    default:
      std::cout << "n unknown signal." << std::endl;
      break;
  }

  // Check if a signal was set previously.
  if(CurrentSignal != NoSignal)
  {
    std::cout << "# Signal was already set and the corresponding action "
              << "has not been carried out yet." << std::endl;
    std::cout << "# Signal will be overwritten to NoSignal. This means "
              << "that nothing will happen." << std::endl;
    CurrentSignal = NoSignal;
  }
  else
  {
    // Signal was not set. So set it now.
    switch( sig )
    {
      case SIGUSR1:
        std::cout << "# Solution file will be written after the current time step."
                  << std::endl;
        CurrentSignal = SignalWrite;
        break;

      case SIGUSR2:
        std::cout << "# Solution file will be written after the current time step "
                  <<   "and the program will stop" << std::endl;
        CurrentSignal = SignalWriteQuit;
        break;

      default:
        std::cout << "# Signal will be ignored." << std::endl;
        break;
    }
  }
    
  // End of the message.
  std::cout << "#" << std::endl
            << "#=========================================="
            << "============================" << std::endl
            << "#" << std::endl << std::flush;
  }
}

#endif
