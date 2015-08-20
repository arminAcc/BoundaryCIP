add_library(bcip ${NGS_LIB_TYPE}
  specialcoef.cpp
  boundlaplace.cpp
  python.cpp
  cipintegrator.cpp
  surffacetfes.cpp
)

set_target_properties(bcip PROPERTIES SUFFIX ".so")

if(NOT WIN32)
  target_link_libraries(bcip  ${PYTHON_LIBS} ${NETGEN_LIBS})
  install( TARGETS bcip ${ngs_install_dir})
endif(NOT WIN32)

install( FILES
  specialcoef.hpp 
  surffacetfes.hpp
  DESTINATION include
)
