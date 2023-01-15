//------------------------------------------------------------------------------
// ElementClass.hpp: Definition of the class, which contains the variables and
//                   functions of a hexahedral element used in the
//                   Discontinuous Galerkin method. This include file should not
//                   be included by any source code file directly. These files
//                   should just include VCP3D_CurviLinear.hpp.
//------------------------------------------------------------------------------

#pragma once

// Include the class to store the interpolation data for the exchange points
// for the wall modeling.
#include "ExchangeDataWallModelClass.hpp"

// Definition of the class ElementClass.
class ElementClass
{
public:
  //--------------------------------------------------
  // Constructor and destructor.
  //--------------------------------------------------

  // Overloaded constructor. Initialize the member variables.
  ElementClass(const InputParamClass      *inputParam,
               const StandardElementClass *standardHex,
               const int                  *globInd,
               const ENUM_ELEM_TYPE       elemType);

  // Destructor. Release the memory.
  ~ElementClass();

  //--------------------------------------
  // Public member variables.
  //--------------------------------------

  // What type of element.
  ENUM_ELEM_TYPE mElemType;

  // The global indices of the element.
  int mGlobalInd[3];

  // The local indices of the element.
  int mLocalInd[3];

  // Vector of arrays with the coordinates of the nodal grid DOFs.
  std::vector<su2double *> mCoorNodalGridDOFs;

  // Vector of arrays with the coordinates of the nodal solution DOFs.
  std::vector<su2double *> mCoorNodalSolDOFs;

  // Vector of arrays with the metric terms in the surface integration points
  // of the iMin boundary.
  std::vector<su2double *> mSurfMetricIntIMin;

  // Vector of arrays with the metric terms in the surface integration points
  // of the iMax boundary.
  std::vector<su2double *> mSurfMetricIntIMax;

  // Vector of arrays with the metric terms in the surface integration points
  // of the jMin boundary.
  std::vector<su2double *> mSurfMetricIntJMin;

  // Vector of arrays with the metric terms in the surface integration points
  // of the iMax boundary.
  std::vector<su2double *> mSurfMetricIntJMax;

  // Vector of arrays with the metric terms in the surface integration points
  // of the kMin boundary.
  std::vector<su2double *> mSurfMetricIntKMin;

  // Vector of arrays with the metric terms in the surface integration points
  // of the kMax boundary.
  std::vector<su2double *> mSurfMetricIntKMax;

  // Vector of arrays with the metric terms in the volume integration points.
  std::vector<su2double *> mVolMetricInt;

  // Vector of arrays with the metric terms in the volume solution DOFs.
  std::vector<su2double *> mVolMetricSolDOFs;

  // Vector of arrays of the solution for the modal form.
  std::vector<su2double *> mSol;

  // Vector of arrays of the average of the primitive variables in the DOFs.
  std::vector<su2double *> mAvePrim;

  // Vector of arrays of the average of the velocity products in the DOFs.
  std::vector<su2double *> mAveVelProd;

  // Array with the average of the eddy viscosity in the DOFs.
  su2double *mAveEddyVis;
  
  // Vector of arrays of the old solution for the modal form, needed for Runge Kutta.
  std::vector<su2double *> mSolOld;

  // Vector of arrays of the residuals for the modal form.
  std::vector<su2double *> mRes;

  // Vector of arrays of the residuals for the neighboring element on the iMin face.
  std::vector<su2double *> mResIMin;

  // Vector of arrays of the residuals for the neighboring element on the jMin face.
  std::vector<su2double *> mResJMin;

  // Vector of arrays of the residuals for the neighboring element on the kMin face.
  std::vector<su2double *> mResKMin;

  // Vector of arrays, which contains the transformation matrices dUdV in the integration
  // points. Only when entropy variables are used.
  std::vector<su2double *> mDUdVInt;

  // LES length scale of the element.
  su2double mLenScaleLES;

  // Volume length scale of the element, i.e. the cubic root of the volume.
  su2double mLenScaleVolume;

  // Length scale in i-direction of the element.
  su2double mLenScaleIDir;

  // Length scale in j-direction of the element.
  su2double mLenScaleJDir;

  // Length scale in k-direction of the element.
  su2double mLenScaleKDir;

  // Averaged unit normal in i-direction.
  su2double mAveNormIDir[3];

  // Averaged unit normal in j-direction.
  su2double mAveNormJDir[3];

  // Averaged unit normal in k-direction.
  su2double mAveNormKDir[3];

  // Volume of the element.
  su2double mVolElem;

  // Possible periodic translation to the neighbor on the iMin boundary.
  su2double mTransIMin[3];

  // Possible periodic translation to the neighbor on the iMax boundary.
  su2double mTransIMax[3];

  // Possible periodic translation to the neighbor on the jMin boundary.
  su2double mTransJMin[3];

  // Possible periodic translation to the neighbor on the jMax boundary.
  su2double mTransJMax[3];

  // Possible periodic translation to the neighbor on the kMin boundary.
  su2double mTransKMin[3];

  // Possible periodic translation to the neighbor on the kMax boundary.
  su2double mTransKMax[3];

  // The possible corresponding boundary subface for the iMin boundary of the element.
  // It is NULL if it is an internal boundary.
  SubfaceBaseClass *mBCIMin;

  // The possible corresponding boundary subface for the iMax boundary of the element.
  // It is NULL if it is an internal boundary.
  SubfaceBaseClass *mBCIMax;

  // The possible corresponding boundary subface for the jMin boundary of the element.
  // It is NULL if it is an internal boundary.
  SubfaceBaseClass *mBCJMin;

  // The possible corresponding boundary subface for the jMax boundary of the element.
  // It is NULL if it is an internal boundary.
  SubfaceBaseClass *mBCJMax;

  // The possible corresponding boundary subface for the kMin boundary of the element.
  // It is NULL if it is an internal boundary.
  SubfaceBaseClass *mBCKMin;

  // The possible corresponding boundary subface for the kMax boundary of the element.
  // It is NULL if it is an internal boundary.
  SubfaceBaseClass *mBCKMax;

  // The interpolation data for the exchange points on the iMin boundary if
  // a wall model is used on this boundary.
  ExchangeDataWallModelClass mExchangeDataIMin;

  // The interpolation data for the exchange points on the iMax boundary if
  // a wall model is used on this boundary.
  ExchangeDataWallModelClass mExchangeDataIMax;

  // The interpolation data for the exchange points on the jMin boundary if
  // a wall model is used on this boundary.
  ExchangeDataWallModelClass mExchangeDataJMin;

  // The interpolation data for the exchange points on the jMax boundary if
  // a wall model is used on this boundary.
  ExchangeDataWallModelClass mExchangeDataJMax;

  // The interpolation data for the exchange points on the kMin boundary if
  // a wall model is used on this boundary.
  ExchangeDataWallModelClass mExchangeDataKMin;

  // The interpolation data for the exchange points on the kMax boundary if
  // a wall model is used on this boundary.
  ExchangeDataWallModelClass mExchangeDataKMax;

  // Prescribed boundary data for the integration points on the iMin boundary
  // if this boundary is a physical boundary.
  std::vector<su2double *> mPrescribedDataIMin;

  // Prescribed boundary data for the integration points on the iMax boundary
  // if this boundary is a physical boundary.
  std::vector<su2double *> mPrescribedDataIMax;

  // Prescribed boundary data for the integration points on the jMin boundary
  // if this boundary is a physical boundary.
  std::vector<su2double *> mPrescribedDataJMin;

  // Prescribed boundary data for the integration points on the jMax boundary
  // if this boundary is a physical boundary.
  std::vector<su2double *> mPrescribedDataJMax;

  // Prescribed boundary data for the integration points on the kMin boundary
  // if this boundary is a physical boundary.
  std::vector<su2double *> mPrescribedDataKMin;

  // Prescribed boundary data for the integration points on the kMax boundary
  // if this boundary is a physical boundary.
  std::vector<su2double *> mPrescribedDataKMax;

  //--------------------------------------
  // Public member functions.
  //--------------------------------------

	// Function, which computes the surface area of each (physical) boundary
	void ComputeSurfaceAreaBoundary(const InputParamClass      *inputParam,
																	const StandardElementClass *standardHex,
																	su2double                  *surfArea);

	// Function, which sets the tuning parameters required for a NSCBC. 
	void TuningParamNSCBC(const InputParamClass 		 *inputParam,
												const StandardElementClass *standardHex);

	// Function, which computes the required averaged Mach on each NSCBC boundary.
	void AverageMachNSCBC(const InputParamClass       *inputParam,
											  const StandardElementClass	*standardHex,
											  su2double 								 **workArray,
											  su2double                   *LocalMachWeighted);

	// Function, which sets the averaged Mach number over the NSCBC boundary.
	void SetAverageBoundaryMachNSCBC(const su2double *value);

  // Function, which adds fluctuations to the dimensional primitive variables.
  void AddFluctuationsPrimVar(const InputParamClass      *inputParam,
                              const StandardElementClass *standardHex,
                              const unsigned long        timeStep,
                              su2double                  **primVar);

  // Function, which allocates the residual and halo solution arrays.
  void AllocateResidualAndHaloSolutionArrays(const StandardElementClass *standardHex);

  // Function, which computes the coordinates of the nodal solution DOFs.
  void ComputeCoorNodalSolDOFs(const int nDOFs1DGrid,
                               const int nDOFs1DSol,
                               su2double *lagrangeDOFs1D,
                               su2double *derLagrangeDOFs1D);

  // Function, which computes the different length scales of the element.
  void ComputeLengthScales(const InputParamClass      *inputParam,
                           const StandardElementClass *standardHex);

  // Function, which computes the metric terms.
  void ComputeMetricTerms(const int nDOFs1DGrid,
                          const int nInt1D,
                          su2double *lagrangeInt1D,
                          su2double *derLagrangeInt1D,
                          su2double *lagrangeMinFace1D,
                          su2double *lagrangeMaxFace1D,
                          su2double *derLagrangeMinFace1D,
                          su2double *derLagrangeMaxFace1D);

  // Function, which computes the dimensional primitive variables in
  // the nodal DOFs.
  void ComputePrimitiveVariablesNodalDOFs(const InputParamClass      *inputParam,
                                          const StandardElementClass *standardHex,
                                          su2double                  **primVar);

  // Function, which computes the allowable time step for the element.
  su2double ComputeTimeStep(const InputParamClass      *inputParam,
                            const StandardElementClass *standardHex,
                            const su2double            *f1,
                            const su2double            *f2,
                            su2double                  **solNodal,
                            su2double                  **gradSolNodal);

  // Function, which computes the velocity gradients in the nodal DOFs.
  void ComputeVelocityGradientsNodalDOFs(const InputParamClass      *inputParam,
                                         const StandardElementClass *standardHex,
                                         su2double                  **primVarNodal,
                                         su2double                  **gradSolNodal);

  // Function, which computes the velocity gradients from the primitive variables.
  void ComputeVelocityGradients(const StandardElementClass *standardHex,
                                su2double                  **primVar,
                                su2double                  **gradVel);

  // Function, which computes the working variables from the primitive variables.
  void ComputeWorkingVariables(const InputParamClass      *inputParam,
                               const StandardElementClass *standardHex,
                               su2double                  **primVar);

  // Function, which initializes the solution in the DOFs of the element.
  void InitSol(const InputParamClass      *inputParam,
               const StandardElementClass *standardHex,
               const su2double            *primVar);

  // Function, which determines the interpolation weights for the exchange location
  // of the wall model of the wall boundary faces, if needed.
  void InterpolationWeightsExchangeLocation(
                    const InputParamClass                                    *inputParam,
                    const StandardElementClass                               *standardHex,
                    std::vector<std::vector< std::vector<ElementClass *> > > &elements,
                    const int                                                nElemPerRankI,
                    const int                                                nElemPerRankJ,
                    const int                                                nElemPerRankK,
                    CADTElemClass                                            *localVolumeADT,
                    const std::vector<su2double>                             &rSample1D);

  // Function, which computes the residuals of the i face(s) of this element.
  void IFaceResiduals(const InputParamClass                                    *inputParam,
                      const StandardElementClass                               *standardHex,
                      std::vector<std::vector< std::vector<ElementClass *> > > &elements,
                      su2double                                                **workArray,
                      const bool                                               ComputeMonitoringData,
                      su2double                                                &EddyVisMax,
                      su2double                                                *forceCoef);

  // Function, which computes the residuals of the j face(s) of this element.
  void JFaceResiduals(const InputParamClass                                    *inputParam,
                      const StandardElementClass                               *standardHex,
                      std::vector<std::vector< std::vector<ElementClass *> > > &elements,
                      su2double                                                **workArray,
                      const bool                                               ComputeMonitoringData,
                      su2double                                                &EddyVisMax,
                      su2double                                                *forceCoef);

  // Function, which computes the residuals of the k face(s) of this element.
  void KFaceResiduals(const InputParamClass                                    *inputParam,
                      const StandardElementClass                               *standardHex,
                      std::vector<std::vector< std::vector<ElementClass *> > > &elements,
                      su2double                                                **workArray,
                      const bool                                               ComputeMonitoringData,
                      su2double                                                &EddyVisMax,
                      su2double                                                *forceCoef);

  // Function, which determines whether or not an LGL distribution of the
  // nodal grid DOFs is used.
  bool LGLDistribution(const InputParamClass      *inputParam,
                       const StandardElementClass *standardHex);

  // Function, which multiplies the residual with the inverse of
  // the mass matrix.
  void MultiplyResInverseMassMatrix(const InputParamClass      *inputParam,
                                    const StandardElementClass *standardHex,
                                    su2double                  **workArray);

  // Function, which determines the prescribed data in the integration
  // points of the boundary faces.
  void PrescribedDataIntegrationPoints(const InputParamClass      *inputParam,
                                       const StandardElementClass *standardHex,
                                       const su2double            *nonDimPrimVarFreeStream);

  // Function, which stores the local indices of the element.
  void StoreLocalIndices(const int i,
                         const int j,
                         const int k);

  // Function, which updates the data for the averages.
  void UpdateAverageData(const InputParamClass      *inputParam,
                         const StandardElementClass *standardHex,
                         const int                  nTimeSteps,
                         su2double                  **primVar,
                         su2double                  **gradVel);

  // Function, which computes the volume residuals of this element.
  void VolumeResidual(const InputParamClass      *inputParam,
                      const StandardElementClass *standardHex,
                      su2double                  **workArray,
                      const bool                 ComputeMonitoringData,
                      su2double                  &Mach2Max,
                      su2double                  &EddyVisMax);

private:
  //--------------------------------------
  // Private member variables.
  //--------------------------------------

  //------------------------------------------------
  // Private member functions.
  //------------------------------------------------

  // Function, which computes the dimensional primitive variables
  // for the given number of items.
  void ComputePrimitiveVariables(const InputParamClass *inputParam,
                                 const int             nItems,
                                 su2double             **primVar);

  // Function, which creates the final metric terms for a face in i-direction.
  void FinalMetricIFace(const int                nItems,
                        std::vector<su2double *> &surfMetricInt);

  // Function, which creates the final metric terms for a face in j-direction.
  void FinalMetricJFace(const int                nItems,
                        std::vector<su2double *> &surfMetricInt);

  // Function, which creates the final metric terms for a face in k-direction.
  void FinalMetricKFace(const int                nItems,
                        std::vector<su2double *> &surfMetricInt);

  // Function, which carries out a containment search in a high order element,
  // where an initial guess can be obtained from linear sub-elements.
  void HighOrderContainmentSearch(const su2double              *coorExchange,
                                  const unsigned short         subElem,
                                  const StandardElementClass   *standardHex,
                                  const su2double              *weightsInterpol,
                                  const std::vector<su2double> &VInv,
                                  const std::vector<su2double> &rSample1D,
                                  su2double                    *parCoor);

  // Function, which determines the interpolation weights for the exchange locations
  // corresponding to the integration points of the given face.
  void InterpolationWeightsExchangeLocationFace(
                    su2double                                                **faceMetric,
                    const su2double                                          factNorm,
                    const InputParamClass                                    *inputParam,
                    const StandardElementClass                               *standardHex,
                    std::vector<std::vector< std::vector<ElementClass *> > > &elements,
                    const int                                                nElemPerRankI,
                    const int                                                nElemPerRankJ,
                    const int                                                nElemPerRankK,
                    CADTElemClass                                            *localVolumeADT,
                    ExchangeDataWallModelClass                               &exchangeData,
                    const std::vector<su2double>                             &rSample1D);

  // Function, which determines the prescribed data in the integration
  // points of the given boundary face.
  void PrescribedDataIntegrationPointsFace(const std::vector<su2double *> &faceMetric,
                                           SubfaceBaseClass               *BC,
                                           const InputParamClass          *inputParam,
                                           const StandardElementClass     *standardHex,
                                           const su2double                *nonDimPrimVarFreeStream,
                                           std::vector<su2double *>       &prescribedData);

  //------------------------------------------------
  // Disabled constructors and operators.
  //------------------------------------------------

  // Default constructor, disabled.
  ElementClass();

  // Copy constructor, disabled.
  ElementClass(const ElementClass&);

  // Assignment operator, disabled.
  ElementClass& operator=(const ElementClass&);
};
