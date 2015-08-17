/**
 * @file specialcoef.hpp
 * @author  C.Lehrenfeld
 * @date August 2015
 * @brief
 */

#pragma once 

#include <fem.hpp>

namespace bcip {
  using namespace ngfem;

  class MySpecialCoeff : public CoefficientFunction
  {
    ///
  public:
    ///
    MySpecialCoeff (){;};
    ///
    virtual ~MySpecialCoeff (){;};
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;
  };

}//end of namespace bcip



