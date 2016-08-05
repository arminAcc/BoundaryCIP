/*********************************************************************/
/* File:   cipintegrator.cpp                                         */
/* Author: Christoph Lehrenfeld                                      */
/* Date:   August  2015                                              */
/*********************************************************************/

/*
  

 */


#include <fem.hpp>


namespace bcip {
  using namespace ngfem;

  void MySetRefBaryCenter(ELEMENT_TYPE eltype, FlatVector<> & point)
  {
    switch (eltype)
    {
    case ET_POINT: point = 0.0; break;
    case ET_SEGM: point = 0.5; break;
    case ET_TRIG: point = 1.0/3.0; break;
    case ET_QUAD: point = 0.5; break;
    case ET_TET: point = 1.0/4.0; break;
    case ET_PYRAMID: point(0) = 1.0/2.0; point(1) = 1.0/2.0; point(2) = 1.0/4.0; break;
    case ET_PRISM: point(0) = 1.0/3.0; point(1) = 1.0/3.0; point(2) = 1.0/2.0; break;
    case ET_HEX: point = 1.0/2.0; break;
    };
  };


  template <int D, int difforder=1>
  class dudnJumpIntegrator : public FacetBilinearFormIntegrator
  {
  protected:
    shared_ptr<CoefficientFunction> coef;
  public:
    dudnJumpIntegrator (const Array<shared_ptr<CoefficientFunction>> & coeffs) 
      : FacetBilinearFormIntegrator(),coef(coeffs[0])
    { 
    }

    virtual ~dudnJumpIntegrator () { ; }
    
    virtual bool BoundaryForm () const 
    { return 0; }

    virtual bool IsSymmetric () const 
    { return true; }
    
    static Integrator * Create (const Array<shared_ptr<CoefficientFunction>> & coeffs)
    {
      return new dudnJumpIntegrator (coeffs);
    }
    
    virtual void CalcElementMatrix (const FiniteElement & fel,
                                    const ElementTransformation & eltrans, 
                                    FlatMatrix<double> elmat,
                                    LocalHeap & lh) const
    {
      throw Exception("dudnJumpIntegrator::CalcElementMatrix - not implemented!");
    }
    
    virtual void CalcFacetMatrix (const FiniteElement & volumefel1, int LocalFacetNr1,
                                  const ElementTransformation & eltrans1, FlatArray<int> & ElVertices1,
                                  const FiniteElement & volumefel2, int LocalFacetNr2,
                                  const ElementTransformation & eltrans2, FlatArray<int> & ElVertices2,
                                  FlatMatrix<double> & elmat,
                                  LocalHeap & lh) const
    {
      static int timer = NgProfiler::CreateTimer ("dudnJumpIntegratorIntegrator");
      NgProfiler::RegionTimer reg (timer);

      if (LocalFacetNr2==-1) throw Exception("dudnJumpIntegrator: LocalFacetNr2==-1");

      const ScalarFiniteElement<D> * fel1_l2 = 
        dynamic_cast<const ScalarFiniteElement<D>*> (&volumefel1);
      ELEMENT_TYPE eltype1 = volumefel1.ElementType();
      int nd1 = fel1_l2->GetNDof();

      const ScalarFiniteElement<D> * fel2_l2 = NULL;
      ELEMENT_TYPE eltype2 = eltype1;

      fel2_l2 = dynamic_cast<const ScalarFiniteElement<D>*> (&volumefel2);
      eltype2 = volumefel2.ElementType();
      int nd2 = fel2_l2->GetNDof();
      double maxorder = max(fel1_l2->Order(),fel2_l2->Order());

      FlatVector<double> eltype1_ref_barcenter(3,lh);
      MySetRefBaryCenter(eltype1,eltype1_ref_barcenter);
      IntegrationPoint refbarycenterv(eltype1_ref_barcenter,0.0);
      MappedIntegrationPoint<D,D> barycenterv (refbarycenterv, eltrans1);
      Vec<D> pointv = barycenterv.GetPoint();
      FlatVector<double> eltype2_ref_barcenter(3,lh);
      MySetRefBaryCenter(eltype2,eltype2_ref_barcenter);
      IntegrationPoint refbarycenterv2(eltype2_ref_barcenter,0.0);
      MappedIntegrationPoint<D,D> barycenterv2 (refbarycenterv2, eltrans2);
      Vec<D> pointv2 = barycenterv2.GetPoint();
      Vec<D> pointsdiff = pointv - pointv2;
      
      // *testout << " pointv  = " << pointv << endl;
      // *testout << " pointv2 = " << pointv2 << endl;


      elmat = 0.0;

      FlatVector<> mat1_shape(nd1, lh);
      FlatVector<> mat1_dudn(nd1, lh);
      FlatVector<> mat2_shape(nd2, lh);
      FlatVector<> mat2_dudn(nd2, lh);
      
      FlatMatrixFixHeight<1> bmat(nd1+nd2, lh);
      FlatMatrixFixHeight<1> dbmat(nd1+nd2, lh);
      Mat<1> dmat;

      FlatMatrixFixWidth<D> dshape(nd1, lh);
      FlatMatrixFixWidth<D> fac_dshape(nd1, lh);
      
      Facet2ElementTrafo transform1(eltype1,ElVertices1); 
      Facet2ElementTrafo transform2(eltype2,ElVertices2); 

      const NORMAL * normals1 = ElementTopology::GetNormals(eltype1);
      const NORMAL * normals2 = ElementTopology::GetNormals(eltype2);

      HeapReset hr(lh);
      ELEMENT_TYPE etfacet = ElementTopology::GetFacetType (eltype1, LocalFacetNr1);

      Vec<D> normal_ref1, normal_ref2;
      for (int i=0; i<D; i++){
        normal_ref1(i) = normals1[LocalFacetNr1][i];
        normal_ref2(i) = normals2[LocalFacetNr2][i];
      }
      const IntegrationRule & ir_facet =
        SelectIntegrationRule (etfacet, 2*maxorder);
	
      if (maxorder==0) maxorder=1;
   
      bmat = 0.0;
      for (int l = 0; l < ir_facet.GetNIP(); l++)
      {
        IntegrationPoint ip1 = transform1(LocalFacetNr1, ir_facet[l]);
	  
        MappedIntegrationPoint<D,D> sip1 (ip1, eltrans1);
        // double lam = coef->Evaluate(sip1);

        // Mat<D> jac1 = sip1.GetJacobian();
        Mat<D> inv_jac1 = sip1.GetJacobianInverse();
        double det1 = sip1.GetJacobiDet();

        Vec<D> normal1 = det1 * Trans (inv_jac1) * normal_ref1;       
        double len1 = L2Norm (normal1);
        normal1 /= len1;

        fel1_l2->CalcShape(sip1.IP(), mat1_shape);
        Vec<D> invjac_normal1 = inv_jac1 * normal1;



        if (difforder == 1)
        {
          mat1_dudn = fel1_l2->GetDShape (sip1.IP(), lh) * invjac_normal1;
        }
        else
        {
          // DO num diff in normal direction
          double eps = 1e-7;

          IntegrationPoint ipl(sip1.IP());
          IntegrationPoint ipc(sip1.IP());
          IntegrationPoint ipr(sip1.IP());

          for (int d = 0; d < D; ++d)
            ipl(d) -= eps * invjac_normal1(d);
          for (int d = 0; d < D; ++d)
            ipr(d) += eps * invjac_normal1(d);

          // double len = L2Norm(invjac_normal2);

          FlatVector<> shapeleft = fel1_l2->GetShape (ipl, lh);
          FlatVector<> shapecenter = fel1_l2->GetShape (ipc, lh);
          FlatVector<> shaperight = fel1_l2->GetShape (ipr, lh);

          mat1_dudn = shaperight - 2 * shapecenter + shapeleft;
          mat1_dudn *= 1.0/(eps*eps);
        }

        
        IntegrationPoint ip2 = (LocalFacetNr2!=-1) ? transform2(LocalFacetNr2, ir_facet[l]) : ip1;
        MappedIntegrationPoint<D,D> sip2 (ip2, eltrans2);
        // double lam2 = coef_lam->Evaluate(sip2);
        // Mat<D> jac2 = sip2.GetJacobian();
        Mat<D> inv_jac2 = sip2.GetJacobianInverse();
        double det2 = sip2.GetJacobiDet();
	  
        Vec<D> normal2 = det2 * Trans (inv_jac2) * normal_ref2;       
        double len2 = L2Norm (normal2); 
        if(abs(len1-len2)>1e-6){
          std::cout << "len :\t" << len1 << "\t=?=\t" << len2 << std::endl;
          throw Exception ("dudnJumpIntegrator: len1!=len2");
        }
        normal2 /= len2;
        Vec<D> invjac_normal2;;
        fel2_l2->CalcShape(sip2.IP(), mat2_shape);
        invjac_normal2 = inv_jac2 * normal2;

        if (difforder == 1)
        {
          mat2_dudn = fel2_l2->GetDShape (sip2.IP(), lh) * invjac_normal2;
        }
        else
        {
          // DO num diff in normal direction
          double eps = 1e-7;

          IntegrationPoint ipl(sip2.IP());
          IntegrationPoint ipc(sip2.IP());
          IntegrationPoint ipr(sip2.IP());

          for (int d = 0; d < D; ++d)
            ipl(d) -= eps * invjac_normal2(d);
          for (int d = 0; d < D; ++d)
            ipr(d) += eps * invjac_normal2(d);

          // double len = L2Norm(invjac_normal2);

          FlatVector<> shapeleft = fel2_l2->GetShape (ipl, lh);
          FlatVector<> shapecenter = fel2_l2->GetShape (ipc, lh);
          FlatVector<> shaperight = fel2_l2->GetShape (ipr, lh);

          mat2_dudn = shaperight - 2 * shapecenter + shapeleft;
          mat2_dudn *= 1.0/(eps*eps);
        }

        
        bmat.Row(0).Range (0   , nd1)   = mat1_dudn;	    
        bmat.Row(0).Range (nd1   , nd1+nd2)   = mat2_dudn;

        dmat(0,0) = 1.0;

        const double orthdist = abs(InnerProduct(pointsdiff,normal1));

        *testout << " orthdist h = " << orthdist << endl;
        // *testout << "       len1 = " << len1 << endl;
        double scal = pow(orthdist, 2*difforder + 1);
        dmat *= coef->Evaluate(sip1) * len1 * scal * ir_facet[l].Weight();
        dbmat = dmat * bmat;
        elmat += Trans (bmat) * dbmat;
        // elmat = 1.0;
      }
      if (LocalFacetNr2==-1) elmat=0.0;
    }
  };


  static RegisterBilinearFormIntegrator<dudnJumpIntegrator<2,1> > init_gp_2d_1cip ("cip", 2, 1);
  static RegisterBilinearFormIntegrator<dudnJumpIntegrator<2,2> > init_gp_2d_1cip2 ("cip2ndorder", 2, 1);
}
