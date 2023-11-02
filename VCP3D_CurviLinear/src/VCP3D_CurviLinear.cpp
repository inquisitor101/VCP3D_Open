//------------------------------------------------------------------------------
// Main program to solve the 3D NavierStokes equations using a high order
// Discontinuous Galerkin discretization for a 3D geometry for which a
// single block grid structure can be used.
//------------------------------------------------------------------------------

#include "VCP3D_CurviLinear.hpp"

// Definition of the global variables, see VCP3D_CurviLinear.hpp for more details.
int rank, nRanks;

ENUM_SIGNAL CurrentSignal;

su2double rhoDim;
su2double pDim;
su2double pRef;
su2double uRef;
su2double rhoRef;

su2double mu;

// Start of the main program.
int main(int argc, char **argv)
{
  // Start up the code.
  StartUpCode(&argc, &argv);

  // Create an object for the input parameters and read them.
  InputParamClass inputParam;
  inputParam.ReadInputParameters(argv[1]);

	// Read the sponge layer information, if specified.
	if( inputParam.mSpongeLayer || inputParam.mCharacteristicMatchingLayer ) inputParam.ReadSpongeLayerData(argv[1]);

  // Read the subface information and the prescribed boundary data.
  inputParam.ReadSubfaceInformation(argv[1]);
  inputParam.ReadBoundaryData(argv[1]);

  // Compute the reference conditions.
  inputParam.ComputeReferenceConditions();

  // Create an object of the class SolverClass.
  SolverClass solver(&inputParam);

  // Read the grid, compute the metric terms and determine
  // the bounding box for the fluctuations.
  solver.ReadGrid();
  solver.ComputeMetricTerms();
  solver.BoundingBoxFluctuations();

  // Write the SU2 grid file, if needed.
  if( inputParam.mWriteSU2GridFile ) solver.WriteSU2GridFile();

  // Initializes the solution. Overwrite when a restart is performed.
  solver.InitSol();
  if( inputParam.mRestart ) solver.ReadRestartSolution();

	// If a sponge layer is specified, initialize the layers and then 
	// proceed to read the averaged solution used as the damping state.
	if( inputParam.mSpongeLayer ) solver.InitSpongeLayer();	

  // Determine the prescribed boundary conditions in the
  // integration points of boundary faces.
  solver.PrescribedDataIntegrationPoints();

  // Add fluctuations to the initial solution, if desired.
  solver.AddFluctuations();

	// If NSCBC is specified, initialize the boundaries.
	if( inputParam.mNSCBC_Specified ) solver.InitNSCBC();

  // Connect the user signals, if supported.
  CurrentSignal = NoSignal;
#ifndef USE_NO_SIGNALS
  ConnectSignals();
#endif

  // Write the initial solution.
  solver.WriteRestartSolution(0);
  solver.WriteParaviewSolution(0);

  // Write the header for the time step history.
  inputParam.WriteHeaderTimeStepHistory();

  // Compute the reference time for the timings.
  const double time0 = GetWallTime();

	// Compute the average data required in the NSCBC, in case imposed.
	if( inputParam.mNSCBC_Specified ) solver.AverageDataNSCBC();

  // Define the array containing the monitoring data and compute the
  // residual of the current solution.
  su2double monitoringData[8];
  solver.Residual(monitoringData, true);

  // Write the time step data.
  solver.WriteTimeStepData(0, time0, monitoringData);

  // Initialize the global signal to no signal. Note that the
  // actual signal decision making is based on GlobalSignal
  // and not on CurrentSignal.
  ENUM_SIGNAL GlobalSignal = NoSignal;

  // Loop over the number of time synchronization steps to be taken.
  bool solSaved = true;
  for(int timeStep=1; timeStep<=inputParam.mNTimeSynchrMax; ++timeStep)
  {
    // Set solSaved to false.
    solSaved = false;

    // Rewrite the header for the time step history, if needed.
    if(timeStep%(50*inputParam.mNWriteTimeHistory) == 0)
      inputParam.WriteHeaderTimeStepHistory();

    // Carry out one time synchronization step.
    solver.ExecuteTimeSynchronizationStep(monitoringData);

    // Write the time history data, if needed. Also the all reduce for
    // the signal takes place here.
    if(timeStep%inputParam.mNWriteTimeHistory == 0)
    {
      solver.WriteTimeStepData(timeStep, time0, monitoringData);

      // Check for signal support.
#ifndef USE_NO_SIGNALS

      // Determine the global signal. For MPI an Allreduce must be
      // used, while in sequential mode it is simply a copy.
#ifdef HAVE_MPI
      short localSignal = (short) CurrentSignal, globalSignal;
      MPI_Allreduce(&localSignal, &globalSignal, 1, MPI_SHORT,
                    MPI_MAX, MPI_COMM_WORLD);
      GlobalSignal = (ENUM_SIGNAL) globalSignal;
#else
      GlobalSignal = CurrentSignal;
#endif
      // Set CurrentSignal back to no signal.
      CurrentSignal = NoSignal;
#endif
    }

    // Check if a solution file must be written.
    if((timeStep%inputParam.mNSave == 0) || (GlobalSignal == SignalWrite) ||
       (GlobalSignal == SignalWriteQuit))
    {
      // Write the solution files and set solSaved to true.
      solver.WriteRestartSolution(timeStep);
      solver.WriteParaviewSolution(timeStep);
      solver.WriteSpaceTimeAverageData();
      solSaved = true;

      // Break the for loop over the number of time steps
      // if the corresponding kill signal has been received.
      if(GlobalSignal == SignalWriteQuit) break;

      // Reset the global signal to no signal.
      GlobalSignal = NoSignal;
    }
  }

  // When this point is reached and the solution is not saved, that means
  // that the maximum number of time synchronization steps has been reached.
  if( !solSaved )
  {
    solver.WriteRestartSolution(inputParam.mNTimeSynchrMax);
    solver.WriteParaviewSolution(inputParam.mNTimeSynchrMax);
    solver.WriteSpaceTimeAverageData();
  }

  // Quit from MPI, if needed.
#ifdef HAVE_MPI
  solver.FreeCommRequests();
  MPI_Finalize();
#endif

  // Return 0 to indicate that everything went fine.
  return 0;
}
