/**
 * @file specialcoef.cpp
 * @author  C.Lehrenfeld
 * @date August 2015
 * @brief 
 */
#include "specialcoef.hpp"
  
namespace bcip {
  using namespace ngfem;


  double MySpecialCoeff::Evaluate (const BaseMappedIntegrationPoint & ip) const
  {
    int onbound = 0;
    for (int d = 0; d < 3; ++d)
      if (ip.IP()(d) < 1e-12) onbound++;
    if (fabs(1-ip.IP()(0)-ip.IP()(1)-ip.IP()(2)) < 3e-12) onbound++;

    if (onbound>1)
      return 0.0;
    else
      return 1.0;
  }

  double ElementNrCoeff::Evaluate (const BaseMappedIntegrationPoint & mip) const
  {
    return mip.GetTransformation().GetElementNr();
  }

  double ElementIdxCoeff::Evaluate (const BaseMappedIntegrationPoint & mip) const
  {
    return mip.GetTransformation().GetElementIndex();
  }

  double ElementHCoeff::Evaluate (const BaseMappedIntegrationPoint & mip) const
  {
    return pow(mip.GetMeasure(),1.0/mip.Dim());
  }
  
}//end of namespace bcip
