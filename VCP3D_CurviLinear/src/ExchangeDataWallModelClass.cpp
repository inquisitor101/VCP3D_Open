//------------------------------------------------------------------------------
// File, which contains the implementation of the member functions of the
// class ExchangeDataWallModelClass.
//------------------------------------------------------------------------------

// Include files needed.
#include "VCP3D_CurviLinear.hpp"

//-----------------------------------------------------------------------------

// Constructor, nothing to be done.
ExchangeDataWallModelClass::ExchangeDataWallModelClass(){}

// Destructor, release the allocated memory.
ExchangeDataWallModelClass::~ExchangeDataWallModelClass()
{
   for(unsigned long i=0; i<mWeightsExchangePoints.size(); ++i)
     FreeMemory((void **) &mWeightsExchangePoints[i]);
}
