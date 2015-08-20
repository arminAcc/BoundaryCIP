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

  class ElementNrCoeff : public CoefficientFunction
  {
    ///
  public:
    ///
    ElementNrCoeff (){;};
    ///
    virtual ~ElementNrCoeff (){;};
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;
  };

  class ElementIdxCoeff : public CoefficientFunction
  {
    ///
  public:
    ///
    ElementIdxCoeff (){;};
    ///
    virtual ~ElementIdxCoeff (){;};
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;
  };

  class ElementHCoeff : public CoefficientFunction
  {
    ///
  public:
    ///
    ElementHCoeff (){;};
    ///
    virtual ~ElementHCoeff (){;};
    ///
    virtual double Evaluate (const BaseMappedIntegrationPoint & ip) const;
  };
  
  
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



