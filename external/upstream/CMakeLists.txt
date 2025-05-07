# Note: order is important as some projects might be interdependent
foreach(dir
  pybind11
  ambit
  bse
  chemps2
  cppe
  dkh
  ecpint
  libefp
  einsums
  gau2grid
  gdma
  libint2
  mdi
  pcmsolver
  pylibefp
  qcelemental
  qcengine
  optking
  simint
  libxc
  brianqc
  ddx
  integratorxx
  gauxc  # keep this _after_ integratorxx as it's checking <project>_DIR
  )
  add_subdirectory(${dir})
endforeach()

# for rare case where resolved ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR} isn't nested, make libdir. may yet be empty
message("pre ${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" _plain_libdir)
message("NATIVE ${_plain_libdir}")
file(TO_CMAKE_PATH "${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR}" _cmake_libdir)
message("CMAKE  ${_cmake_libdir}")

file(MAKE_DIRECTORY ${_plain_libdir})
