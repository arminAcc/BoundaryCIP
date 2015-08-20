#include <python_ngstd.hpp>

#include <comp.hpp>
#include "../utils/specialcoef.hpp"
#include "../utils/surffacetfes.hpp"
using namespace ngcomp;
using namespace bcip;

void ExportNgsBCIP() 
{
  std::string nested_name = "bcip";
  if( bp::scope() )
    nested_name = bp::extract<std::string>(bp::scope().attr("__name__") + ".bcip");
  
  bp::object module(bp::handle<>(bp::borrowed(PyImport_AddModule(nested_name.c_str()))));

  cout << "exporting bcip as " << nested_name << endl;
  bp::object parent = bp::scope() ? bp::scope() : bp::import("__main__");
  parent.attr("bcip") = module ;

  bp::scope local_scope(module);

  bp::def("GetSElEdge", FunctionPointer( [] (shared_ptr<MeshAccess> ma, int selnr)
                                              {
                                                Array<int> fnums(1);
                                                ma->GetSElFacets(selnr,fnums);
                                                return fnums[0];
                                              } ),
          (bp::arg("ma")=NULL,bp::arg("selnr")=0))
    ;


  bp::def("SetEdgeOrder", FunctionPointer( [] (shared_ptr<FESpace> fes, int edgenr, int order)
                                              {
                                                auto fesh1 = dynamic_pointer_cast<H1HighOrderFESpace>(fes);
                                                fesh1->SetEdgeOrder(edgenr, order);
                                              } ),
          (bp::arg("space")=NULL,bp::arg("edgenr")=0,bp::arg("order")=1))
    ;



  bp::def("UpdateDofTables", FunctionPointer( [] (shared_ptr<FESpace> fes)
                                              {
                                                auto fesh1 = dynamic_pointer_cast<H1HighOrderFESpace>(fes);
                                                fesh1->UpdateDofTables();
                                              } ),
          (bp::arg("space")=NULL))
    ;


  bp::class_<ElementHCoeff, shared_ptr<ElementHCoeff>, bp::bases<CoefficientFunction>, boost::noncopyable> ("h", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([]()
                           {
                             return make_shared<ElementHCoeff> ();
                           })));
    
  bp::implicitly_convertible 
    <shared_ptr<ElementHCoeff>, shared_ptr<CoefficientFunction> >();
  
  bp::class_<ElementIdxCoeff, shared_ptr<ElementIdxCoeff>, bp::bases<CoefficientFunction>, boost::noncopyable> ("idx", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([]()
                           {
                             return make_shared<ElementIdxCoeff> ();
                           })));

  bp::implicitly_convertible 
    <shared_ptr<ElementIdxCoeff>, shared_ptr<CoefficientFunction> >();
  
  bp::class_<ElementNrCoeff, shared_ptr<ElementNrCoeff>, bp::bases<CoefficientFunction>, boost::noncopyable> ("nr", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([]()
                           {
                             return make_shared<ElementNrCoeff> ();
                           })));

  bp::implicitly_convertible 
    <shared_ptr<ElementNrCoeff>, shared_ptr<CoefficientFunction> >(); 


  bp::class_<SurfaceFacetFESpace, shared_ptr<SurfaceFacetFESpace>, bp::bases<FESpace>, boost::noncopyable> ("SurfaceFacetFESpace", bp::no_init)
    .def("__init__", bp::make_constructor 
         (FunctionPointer ([](shared_ptr<MeshAccess> ma)
                           {
                             Flags flags;
                             auto fes = make_shared<SurfaceFacetFESpace> (ma,flags);

                             LocalHeap lh (1000000, "FESpace::Update-heap");
                             fes->Update(lh);
                             fes->FinalizeUpdate(lh);
                             
                             return fes;
                           }),
          bp::default_call_policies(),     // need it to use named arguments
          (bp::arg("mesh")=NULL))
      );

  bp::implicitly_convertible 
    <shared_ptr<SurfaceFacetFESpace>, shared_ptr<FESpace> >(); 


  

}

BOOST_PYTHON_MODULE(libbcip) 
{
  cout << " exporting bcip library " << endl;
  ExportNgsBCIP();
}
