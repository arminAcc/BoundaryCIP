add_library(bcip_utils ${NGS_LIB_TYPE}
  specialcoef.cpp
)

set_target_properties(bcip_utils PROPERTIES SUFFIX ".so")

if(NOT WIN32)
  # target_link_libraries(bcip_utils ${NETGEN_LIBS})
  install( TARGETS bcip_utils ${ngs_install_dir})
endif(NOT WIN32)

install( FILES
  specialcoef.hpp 
  DESTINATION include
)
