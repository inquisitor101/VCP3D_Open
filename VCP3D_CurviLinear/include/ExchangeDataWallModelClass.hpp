//------------------------------------------------------------------------------
// ExchangeDataWallModelClass.hpp: Definition of the class, which contains the
//                                 interpolation data for the corresponding
//                                 exchange points of the integration points on
//                                 a boundary face when a wall model is used.
//                                 This include file should not be included by
//                                 any source code file directly. These files
//                                 should just include VCP3D_CurviLinear.hpp.
//------------------------------------------------------------------------------

#pragma once

// Definition of the class ExchangeDataWallModelClass.
class ExchangeDataWallModelClass
{
public:
  //--------------------------------------------------
  // Constructor and destructor.
  //--------------------------------------------------

  // Constructor, nothing to be done.
  ExchangeDataWallModelClass();

  // Destructor, release the allocated memory.
  ~ExchangeDataWallModelClass();

  //--------------------------------------
  // Public member variables.
  //--------------------------------------

  // Vector, which contains the number of exchange points for donor element
  // in cumulative storage format.
  std::vector<int> mNExchangePointsPerDonorElement;

  // Vector, which contains the indices of the exchange points per donor
  // element. If multiple donor elements are present, this sequence may
  // not be contiguous.
  std::vector<int> mExchangePointsPerDonorElement;

  // Local indices (i,j,k) of the donor elements.
  std::vector<int> mIndicesDonorElements;

  // Interpolation coefficients of the exchange points. The first index
  // corresponds to the solution DOFs of the element, the second index
  // corresponds to the integration points. In this way the inner loop
  // to compute the interpolated data can be vectorized.
  std::vector<su2double *> mWeightsExchangePoints;
};
