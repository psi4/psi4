# Note: order is important as some projects might be interdependent
foreach(dir
  pybind11
  ambit
  chemps2
  cppe
  dkh
  libefp
  erd
  gau2grid
  gdma
  libint
  mdi
  pcmsolver
  pylibefp
  qcelemental
  qcengine
  simint
  libxc
  )
  add_subdirectory(${dir})
endforeach()

# Do the xtensor stack separately
# Order is important here too!
add_subdirectory(xtl)
if(USE_XSIMD)
  add_subdirectory(xsimd)
else()
  add_library(xsimd_external INTERFACE)  # dummy
endif()
add_subdirectory(xtensor)
add_subdirectory(xtensor-blas)
add_subdirectory(xtensor-python)
