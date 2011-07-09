# Benchmarks.
add_subdirectory(smooth-iso)
add_subdirectory(smooth-aniso-x)
add_subdirectory(smooth-aniso-y)
add_subdirectory(lshape)
add_subdirectory(layer-interior)
add_subdirectory(layer-boundary)
add_subdirectory(moving-front-time-adapt)
add_subdirectory(moving-front-space-adapt)
# TODO: This does not work in Windows.
if((NOT NMAKE) AND (NOT MSVC))
add_subdirectory(neutronics-heat-conduction)
add_subdirectory(neutronics-2-group-adapt) 
endif((NOT NMAKE) AND (NOT MSVC))

add_subdirectory(stabilized-advection-diffusion)
add_subdirectory(stabilized-advection-reaction)
add_subdirectory(nonsym-check)

#if(NOT WITH_TRILINOS)
  add_subdirectory(screen)
#endif(NOT WITH_TRILINOS)
