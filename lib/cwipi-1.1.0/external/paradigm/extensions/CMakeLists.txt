if(PDM_ENABLE_EXTENSION_PDMA)
  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/paradigma/CMakeLists.txt)
    add_subdirectory(paradigma)  
    set (PDM_LIST_TEST ${PDM_LIST_TEST} PARENT_SCOPE )
    set (PDM_LIST_NRANK ${PDM_LIST_NRANK} PARENT_SCOPE )
  else()
    message(
      FATAL_ERROR
        "Sources of ParaDiGMA extension are not into the extensions/paradigma directory. Please, donwload them from the git command :
     git submodule update --init")
  endif()
endif()
