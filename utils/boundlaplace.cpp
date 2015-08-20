/**
 * @file boundlaplace.cpp
 * @author  C.Lehrenfeld
 * @date August 2015
 * @brief 
 */

#include <fem.hpp>
#include "surffacetfes.hpp"

namespace bcip {
  using namespace ngfem;

//DiffOp und Integrator für Normalkomponente
  class DiffOp_GradTang : public DiffOp<DiffOp_GradTang>
  {
  public:
    enum { DIM = 1 };
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

      // cout << " sip.GetJacobianInverse() = " << sip.GetJacobianInverse() << endl;
      
      // for (int i = 0; i < fel_u.GetNDof(); ++i)
      //   dshape.Row(i) = Trans(sip.GetJacobianInverse()) * dshape_tmp.Row(i);

      double len = sip.GetJacobiDet();
      
      for (int i = 0; i < fel_u.GetNDof(); ++i)
        mat(0,i) = 1.0/len * dshape_tmp(i,0);

      // mat = dshape;
      // cout << " sip.GetPoint() = " << sip.GetPoint() << endl;
      // cout << " dshape_tmp = " << dshape_tmp << endl;
      // cout << " dshape = " << dshape << endl;
      // cout << " mat = " << mat << endl;
      // getchar();
    }
  };  


  class BoundaryLaplaceIntegrator 
    : public T_BDBIntegrator<DiffOp_GradTang, DiagDMat<1>, FiniteElement >   //DimSpace-1

  {
  public:
    BoundaryLaplaceIntegrator (shared_ptr<CoefficientFunction> coeff)
      : T_BDBIntegrator<DiffOp_GradTang, DiagDMat<1>, FiniteElement > (DiagDMat<1> (coeff))
    { ; }

    BoundaryLaplaceIntegrator (const Array<shared_ptr<CoefficientFunction>>& coeff)
      : T_BDBIntegrator<DiffOp_GradTang, DiagDMat<1>, FiniteElement > (DiagDMat<1> (coeff[0]))
    { ; }
 
 
    // static shared_ptr<Integrator> Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    // {
    //   return make_shared<BoundaryLaplaceIntegrator> (coeffs[0]);
    // }
 
    virtual bool BoundaryForm () const { return true; }
    virtual string Name () const { return "BoundLaplace"; }
  };
  
  static RegisterBilinearFormIntegrator<BoundaryLaplaceIntegrator> initmassasdfsdf ("BoundLaplace", 2, 1);









  class BCIPIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coef;
  public:
    BCIPIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs);
    virtual ~BCIPIntegrator(){ ; };
    virtual string Name () const { return "BCIPIntegrator"; }
    virtual int DimElement () const { return 1; }
    virtual int DimSpace () const { return 2; }
    virtual bool BoundaryForm () const { return true; }
    virtual bool IsSymmetric () const { return true; }
    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans,
                                    FlatMatrix<double> elmat,
                                    LocalHeap & lh) const;
  };


  BCIPIntegrator :: BCIPIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) : coef(coeffs[0])
  { ; }

  void BCIPIntegrator :: CalcElementMatrix (const FiniteElement & fel,
                                            const ElementTransformation & eltrans,
                                            FlatMatrix<double> elmat,
                                            LocalHeap & lh) const
  {
    static Timer time_fct ("BCIPIntegrator::CalcElementMatrix");
    RegionTimer reg (time_fct);
    
    elmat = 0.0;
    
    const CompoundFiniteElement & compfe = dynamic_cast<const CompoundFiniteElement &>(fel);
    
    const ScalarFiniteElement<1> & scafe
      = dynamic_cast<const ScalarFiniteElement<1> &>(compfe[2]);

    const MyDummyFE<ET_SEGM> * hybridfe
      = dynamic_cast<const MyDummyFE<ET_SEGM> *>(&(compfe[3]));

    if (hybridfe)
    {

      IntegrationPoint ip[] = {IntegrationPoint(0.0),IntegrationPoint(1.0)};

      FlatVector<> bmat(scafe.GetNDof()+2,lh);
      FlatMatrix<> innerelmat(scafe.GetNDof()+2,lh);
      innerelmat = 0.0;
      
      for (int j = 0; j < 2; ++j)
      {
        bmat = 0.0;
        MappedIntegrationPoint<1,2> mip(ip[j],eltrans);

        FlatMatrixFixWidth<1> dshape_tmp = scafe.GetDShape (mip.IP(), lh);
        double len = mip.GetJacobiDet();

        if (j==0)
          for (int i = 0; i < scafe.GetNDof(); ++i)
            bmat(i) = -1.0/len * dshape_tmp(i,0);
        else
          for (int i = 0; i < scafe.GetNDof(); ++i)
            bmat(i) = 1.0/len * dshape_tmp(i,0);
        
        if (hybridfe->pos[j])
          // if (j==0)
          bmat(scafe.GetNDof()+j) = 1.0; ///len;
          // else
          //   bmat(scafe.GetNDof()+j) = -1.0/len;
        else
          // if (j==0)
          bmat(scafe.GetNDof()+j) = -1.0; ///len;
          // else
          //   bmat(scafe.GetNDof()+j) = 1.0/len;
        innerelmat += coef->EvaluateConst() * pow(len,4.0) * bmat * Trans(bmat);

        // cout << " mip.GetPoint() = " << mip.GetPoint() << endl;
        // if (hybridfe->pos[0])
        //   std::cout << " left is pos " << std::endl;
        // else
        //   std::cout << " left is neg " << std::endl;
        // if (hybridfe->pos[1])
        //   std::cout << " right is pos " << std::endl;
        // else
        //   std::cout << " right is neg " << std::endl;
          
        // cout << " bmat = " << bmat << endl;
        
        // cout << " bmat * Trans(bmat) = \n" << bmat * Trans(bmat) << endl;
        // cout << " innerelmat = \n" << innerelmat << endl;
      }

      IntRange scarange(compfe.GetRange(2).First(),compfe.GetRange(2).Next()+2) ;
      elmat.Rows(scarange).Cols(scarange) = innerelmat;

      // cout << " elmat = " << elmat << endl;
      // getchar();

      // FlatMatrixFixWidth<DIM_ELEMENT> dshape_tmp = fel_u.GetDShape (sip.IP(), lh);
      // FlatMatrixFixWidth<DIM_SPACE> dshape(fel_u.GetNDof(), lh); // = dshape_tmp * trafo;

      // double len = sip.GetJacobiDet();
      
      // for (int i = 0; i < fel_u.GetNDof(); ++i)
      //   mat(0,i) = 1.0/len * dshape_tmp(i,0);


    }
    
    // if (!ElementInRelevantBand(coef_lset_p1, eltrans, lower_lset_bound, upper_lset_bound))
    //   return;
      
    // FlatVector<> shape (scafe.GetNDof(),lh);
      
    // IntegrationRule ir = SelectIntegrationRule (eltrans.GetElementType(), 2*scafe.Order());
    // for (int l = 0 ; l < ir.GetNIP(); l++)
    // {
    //   MappedIntegrationPoint<D,D> mip(ir[l], eltrans);
    //   const double coef_val = coef->Evaluate(mip);
    //   scafe.CalcShape(ir[l],shape);
    //   elmat += coef_val * mip.GetWeight() * shape * Trans(shape);
    // }      
  }

  static RegisterBilinearFormIntegrator<BCIPIntegrator > initbcip ("bcip", 2, 1);
  
}

