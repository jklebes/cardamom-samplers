
add_library(cardamom-samplers-lib SHARED
              math.f90
              random_uniform.f90
              samplers_io/samplers_io.f90
              samplers_shared.f90
              DEMCz.f90
              MHMCMC/MHMCMC.f90
)

find_package(OpenMP)
if(OpenMP_Fortran_FOUND)
  target_link_libraries(cardamom-samplers-lib PUBLIC OpenMP::OpenMP_Fortran)
endif()
