//------------------------------------------------------------------------------
// StandardElementClass.hpp: Definition of the class, which contains the
//                           variables and functions to store the data of a
//                           standard element. This include file should not be
//                           included by any source code file directly. These
//                           files should just include VCP3D_CurviLinear.hpp.
//------------------------------------------------------------------------------

#pragma once

// Definition of the class StandardElementClass.
class StandardElementClass
{
public:
  //------------------------------------------
  // Constructor and destructor.
  //------------------------------------------

  // Constructor.
  StandardElementClass();

  // Destructor.
  ~StandardElementClass();

  //--------------------------------------
  // Public member variables.
  //--------------------------------------

  // Polynomial degree of the element.
  int mNPoly;

  // Number of DOFs in 1 dimension.
  int mNDOFs1D;

  // Padded value for the number of DOFs in 1 dimension.
  // Padding is present to allow for vectorization.
  int mNDOFs1DPad;

  // Total number of DOFs of the element.
  int mNDOFs;

  // Padded value of mNDOFs. Padding is present to allow for vectorization.
  int mNDOFsPad;

  // Number of integration points in 1 dimension.
  int mNIntegration1D;

  // Padded value for the number of integration points in 1 dimension.
  // Padding is present to allow for vectorization.
  int mNIntegration1DPad;

  // Number of integration points in 2 dimensions.
  int mNIntegration2D;

  // Padded value for the number of integration points in 2 dimensions.
  // Padding is present to allow for vectorization.
  int mNIntegration2DPad;

  // Total number of points used in the numerical integration.
  int mNIntegration;

  // Padded value for mNIntegration. Padding is present to allow for vectorization.
  int mNIntegrationPad;

  // 1D parametric coordinates of the grid DOFs for this standard element.
  std::vector<su2double> mRDOFsGrid1D;

  // 1D parametric coordinates of the DOFs for this standard element.
  std::vector<su2double> mRDOFs1D;

  // 1D parametric coordinates of the integration points for this standard element.
  std::vector<su2double> mRIntegration1D;

  // 1D integration weights for the numerical integration.
  std::vector<su2double> mIntegrationWeights1D;

  // 1D integration weights for the numerical integration based on the DOFs.
  std::vector<su2double> mIntegrationWeights1DDOFs;

  // Integration weights for all integration points.
  su2double *mIntWeights;

  // Array, which contains the values of the 1D Vandermonde matrix in the 1D DOFs.
  // The size of the fastest running index is mNDOFs1DPad.
  su2double *mVandermonde1D;

  // Array, which contains the values of the inverse of the 1D Vandermonde matrix
  // in the 1D DOFs. The size of the fastest running index is mNDOFs1DPad.
  su2double *mVandermonde1DInverse;

  // Array, which contains the values of the transpose of the 1D Vandermonde
  // matrix in the 1D DOFs. The size of the fastest running index is mNDOFs1DPad.
  su2double *mVandermonde1DTranspose;

  // Array, which contains the values of the 1D Legendre basis functions in
  // the 1D DOFs. This array is actually a matrix, where the size of the fastest
  // running index is mNDOFs1DPad.
  su2double *mLegendreDOFs1D;

  // Array, which contains the values of the derivatives of the 1D Legendre basis
  // functions in the 1D DOFs. This array is actually a matrix, where the size of
  // the fastest running index is mNDOFs1DPad.
  su2double *mDerLegendreDOFs1D;

  // Array, which contains the values of the 1D Legendre basis functions in
  // the 1D integration points. This array is actually a matrix, where the
  // size of the fastest running index is mNIntegration1DPad.
  su2double *mLegendreInt1D;

  // Transpose of mLegendreInt1D. The fastest running index is mNDOFs1DPad.
  su2double *mLegendreInt1DTranspose;

  // Array, which contains the values of the derivatives of the 1D Legendre basis
  // functions in the 1D integration points. This array is actually a matrix,
  // where the size of the fastest running index is mNIntegration1DPad.
  su2double *mDerLegendreInt1D;

  // Transpose of mDerLegendreInt1D. The fastest running index is mNDOFs1DPad.
  su2double *mDerLegendreInt1DTranspose;

  // Integration weights for all integration points of a face.
  su2double *mIntWeightsFace;

  // Array, which contains the values of the 1D Legendre basis functions at the
  // min face, where the min face is located at r = -1.
  su2double *mLegendreMinFace1D;

  // Array, which contains the values of the 1D Legendre basis functions at the
  // max face, where the max face is located at r = 1.
  su2double *mLegendreMaxFace1D;

  // Array, which contains the values of the derivatives of the 1D Legendre basis
  // functions at the min face, where the min face is located at r = -1.
  su2double *mDerLegendreMinFace1D;

  // Array, which contains the values of the derivatives of the 1D Legendre basis
  // functions at the max face, where the max face is located at r = 1.
  su2double *mDerLegendreMaxFace1D;

  // Array, which contains the values of the 3D basis functions in the center of
  // the element. Needed to compute the value of the solution in the center of
  // the element, which is used in the preconditioning step of the CG algorithm.
  su2double *mBasisCenter;

	// Array, which contains the values of the 1D Lagrange basis functions at 
	// the integration points in 1D.
	su2double *mLagrangeInt1D;

	// Transpose of mLagrangeInt1D. The fastest running index is mNDOFs1DPad.
  su2double *mLagrangeInt1DTranspose;

	// Array, which contains the values of the derivatives of the 1D Lagrange basis 
	// functions at the integration points in 1D.
	su2double *mDerLagrangeInt1D;

	// Transpose of mDerLagrangeInt1D. The fastest running index is mNDOFs1DPad.
  su2double *mDerLagrangeInt1DTranspose;
	
	// Array, which contains the values of the derivatives of the 1D Lagrange basis
	// functions at the DOFs in 1D.
	su2double *mDerLagrangeDOFs1D;

	// Transpose of mDerLagrangeDOFs1D. The fastest running index is mNDOFs1DPad.
  su2double *mDerLagrangeDOFs1DTranspose;

	// Array, which contains the least-squares polynomial-correction matrix needed
	// in case an NSCBC is specified.
	su2double *mPolynomialCorrectionMatrix;

	// Array, which contains the values of the derivatives of the 1D Lagrange basis
	// functions at the integration points in 1D with the integration points as the
	// basis locations. 
	su2double *mDerLagrangeInt1D_BaseInt;

	// Transpose of mDerLagrangeInt1D_BaseInt. The fastest running index is mNIntegration1DPad.
  su2double *mDerLagrangeInt1D_BaseIntTranspose;


  //--------------------------------------
  // Public member functions.
  //--------------------------------------

  // Function, which determines the data for this standard element.
  void DataStandardElement(const InputParamClass *inputParam);

  // Function, which determines the 1D Lagrangian basis functions for the
  // given 1D location of the DOFs and points.
  void LagrangianBasisFunctions(const std::vector<su2double> &rDOFs1D,
                                const std::vector<su2double> &rPoints1D,
                                su2double                    *lagrangePoints1D,
                                su2double                    *derLagrangePoints1D,
                                su2double                    *derLagrangeMinFace1D,
                                su2double                    *derLagrangeMaxFace1D);

private:
  //--------------------------------------
  // Private member variables.
  //--------------------------------------

  //------------------------------------------------
  // Private member functions.
  //------------------------------------------------

  // Function, which computes the values of the basis functions in the
  // center of the element.
  void BasisFunctionsCenter(void);

  // Function, which determines the standard element data needed for the
  // computation of the face integrals.
  void DataStandardFace(void);

  // Function, which determines the 1D location of the integration points.
  void LocationIntegration1D(const int nPolyIntRule);

  // Function, which determines the values of the 1D Legendre basis functions
  // and its derivatives in the 1D DOFs.
  void ValuesLegendreBasisDOFs1D(void);

  // Function, which determines the values of the 1D Legendre basis functions
  // and its derivatives in the 1D integration points.
  void ValuesLegendreBasisInt1D(void);

  // Function, which determines the values of the, possibly padded, 1D
  // Vandermonde matrix in the DOFs.
  void VandermondeDOFs1D(void);

	// Function, which determines the values of the 1D Lagrange basis functions
	// and its derivatives in the 1D interpolation points and 1D DOFs.
	void ValuesLagrangeBasis1D(void);

  // Function, which determine the weights of the 1D integration rule
  // based on the DOFs.
  void Weights1DIntegrationRuleDOFs(void);

	// Function, which determines the values of the polynomial-correction
	// operators, in case an NSCBC is specified.
	void AssemblePolynomialCorrectionOperators(void);
};
