#include <python_ngstd.hpp>

#include <comp.hpp>
using namespace ngcomp;

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

}

BOOST_PYTHON_MODULE(libbcip) 
{
  cout << " exporting bcip library " << endl;
  ExportNgsBCIP();
}
