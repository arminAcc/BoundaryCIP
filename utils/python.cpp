#include <python_ngstd.hpp>
//using namespace ngcomp;

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

}

BOOST_PYTHON_MODULE(libbcip) 
{
  ExportNgsBCIP();
}
