cmake_minimum_required(VERSION 3.28)

project(fDE VERSION 0.0.1
	DESCRIPTION "fortran implementation of differential evolution" 
  LANGUAGES Fortran)
enable_language(Fortran)

set(PROJECT_NAME fDEMCz)

add_subdirectory(src)

#install(TARGETS fDE DESTINATION "lib")

if(ENABLE_TESTING)
  # add the testsuite
  enable_testing()
  add_subdirectory("test")
endif()
