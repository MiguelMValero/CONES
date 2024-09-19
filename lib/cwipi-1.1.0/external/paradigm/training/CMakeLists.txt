include(check_python_module)

if (PDM_ENABLE_FORTRAN)
  check_python_module(sphinxfortran REQUIRED)
endif()

find_package(Jupyter REQUIRED)
find_package(Jupytext REQUIRED)

set(TRAINING_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/training)
set(TRAINING_BUILD  ${CMAKE_CURRENT_BINARY_DIR}/training)

file(GLOB_RECURSE SOURCES_JUPITEXT
  *.md
)

file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/C DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/Fortran DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/Python DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file (COPY ${CMAKE_CURRENT_SOURCE_DIR}/utils DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


if (PDM_ENABLE_SHARED)
    set(JUPYTEXT_DEP pdm_shared)
else ()
    set(JUPYTEXT_DEP pdm_static)
endif()
if(PDM_ENABLE_Fortran)
  if (PDM_ENABLE_SHARED)
      set(JUPYTEXT_DEP ${JUPYTEXT_DEP} pdmf_shared)
  else ()
      set(JUPYTEXT_DEP ${JUPYTEXT_DEP} pdmf_static)
  endif()
endif()

set(SOURCES_JUPITER "")
foreach(_md_file IN ITEMS ${SOURCES_JUPITEXT})
  file(RELATIVE_PATH _md_file_rel ${CMAKE_CURRENT_SOURCE_DIR} ${_md_file})

  get_filename_component(_md_file_dir ${_md_file_rel} DIRECTORY )
  get_filename_component(_md_file_basename ${_md_file_rel} NAME_WE )

  set(_ipynb_file "${CMAKE_CURRENT_BINARY_DIR}/${_md_file_dir}/${_md_file_basename}.ipynb")

  set(SOURCES_JUPITER ${SOURCES_JUPITER} "${_ipynb_file}")

  file(REMOVE ${CMAKE_CURRENT_BINARY_DIR}/${_md_file_dir}/${_md_file_basename}.md)

  add_custom_command(OUTPUT ${_ipynb_file}
                     COMMAND ${JUPYTEXT_EXECUTABLE} --to ipynb ${_md_file} -o ${_ipynb_file} -q
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                     DEPENDS 
                     ${JUPYTEXT_DEP}
                     MAIN_DEPENDENCY ${_md_file}
                     VERBATIM)
endforeach()

add_custom_target(training ALL DEPENDS ${SOURCES_JUPITER} COMMENT "Generating Jupyter notebook files")

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/C
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share/training/pdm)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Fortran
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share/training/pdm)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Python
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share/training/pdm)
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/utils
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share/training/pdm)


