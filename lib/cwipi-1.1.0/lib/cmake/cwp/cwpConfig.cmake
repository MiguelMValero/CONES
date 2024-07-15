set(CWP_ENABLE_SHARED                   ON)
set(CWP_ENABLE_STATIC                   ON)
set(CWP_ENABLE_Fortran                  OFF)
set(cwp_VERSION                         1.1.0)
set(CWP_ENABLE_SHARED_Fortran_INTERFACE )
set(CWP_ENABLE_STATIC_Fortran_INTERFACE )

include("${CMAKE_CURRENT_LIST_DIR}/cwpTargets.cmake")
if(CWP_ENABLE_Fortran)
  include("${CMAKE_CURRENT_LIST_DIR}/../cwpf/pdmfTargets.cmake")
endif()

if(CWP_ENABLE_SHARED)
  include("${CMAKE_CURRENT_LIST_DIR}/../cwp_shared/cwp_sharedTargets.cmake")
endif()
if(CWP_ENABLE_SHARED_Fortran_INTERFACE)
  include("${CMAKE_CURRENT_LIST_DIR}/../cwpf_shared/cwpf_sharedTargets.cmake")
endif()


if(CWP_ENABLE_STATIC)
  include("${CMAKE_CURRENT_LIST_DIR}/../cwp_static/cwp_staticTargets.cmake")
endif()
if(CWP_ENABLE_STATIC_Fortran_INTERFACE)
  include("${CMAKE_CURRENT_LIST_DIR}/../cwpf_static/cwpf_staticTargets.cmake")
endif()

