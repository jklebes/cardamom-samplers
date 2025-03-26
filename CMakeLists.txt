cmake_minimum_required(VERSION 3.25)

project(cardamom-samplers VERSION 0.0.1
	DESCRIPTION "fortran implementation of differential evolution" 
  LANGUAGES Fortran)
enable_language(Fortran)


set(PROJECT_NAME cardamom-samplers)

add_subdirectory(src)

#install(TARGETS cardamom-samplers DESTINATION "lib")

if(ENABLE_TESTING)
  # add the testsuite
  enable_testing()
  add_subdirectory("test")
endif()

find_package(OpenMP)
if(OpenMP_Fortran_FOUND)
  target_link_libraries("${PROJECT_NAME}-lib" PUBLIC OpenMP::OpenMP_Fortran)
endif()
