/**
 * @file boundlaplace.cpp
 * @author  C.Lehrenfeld
 * @date August 2015
 * @brief 
 */

#include <fem.hpp>
  
namespace bcip {
  using namespace ngfem;

//DiffOp und Integrator für Normalkomponente
  class DiffOp_GradTang : public DiffOp<DiffOp_GradTang>
  {
  public:
    enum { DIM = 2 };
    enum { DIM_SPACE = 2 };
    enum { DIM_ELEMENT = 1 }; //DIM_ELEMENT < DIM_SPACE loest Abfrage nach boundary element aus
    enum { DIM_DMAT = 1 };
    enum { DIFFORDER = 1 };

    template <typename FEL, typename SIP, typename MAT>
    static void GenerateMatrix (const FEL & bfel, const SIP & sip,
                                MAT & mat, LocalHeap & lh)
    {
      const ScalarFiniteElement<DIM_ELEMENT> & fel_u = 
        dynamic_cast<const ScalarFiniteElement<DIM_ELEMENT>&> (bfel);
    
      FlatMatrixFixWidth<DIM_ELEMENT> dshape_tmp = fel_u.GetDShape (sip.IP(), lh);
      FlatMatrixFixWidth<DIM_SPACE> dshape(fel_u.GetNDof(), lh); // = dshape_tmp * trafo;
      for (int i = 0; i < fel_u.GetNDof(); ++i)
        dshape.Row(i) = Trans(sip.GetJacobianInverse()) * dshape_tmp.Row(i);
      mat = dshape;
    }
  };  


  class BoundaryLaplaceIntegrator 
    : public T_BDBIntegrator<DiffOp_GradTang, DiagDMat<2>, FiniteElement >   //DimSpace-1

  {
  public:
    BoundaryLaplaceIntegrator (shared_ptr<CoefficientFunction> coeff)
      : T_BDBIntegrator<DiffOp_GradTang, DiagDMat<2>, FiniteElement > (DiagDMat<2> (coeff))
    { ; }

    BoundaryLaplaceIntegrator (const Array<shared_ptr<CoefficientFunction>>& coeff)
      : T_BDBIntegrator<DiffOp_GradTang, DiagDMat<2>, FiniteElement > (DiagDMat<2> (coeff[0]))
    { ; }
 
 
    // static shared_ptr<Integrator> Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    // {
    //   return make_shared<BoundaryLaplaceIntegrator> (coeffs[0]);
    // }
 
    virtual bool BoundaryForm () const { return true; }
    virtual string Name () const { return "BoundLaplace"; }
  };
  
  static RegisterBilinearFormIntegrator<BoundaryLaplaceIntegrator> initmassasdfsdf ("BoundLaplace", 2, 1);

}

