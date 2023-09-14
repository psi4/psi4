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
  erd
  gau2grid
  gdma
  libint
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
  )
  add_subdirectory(${dir})
endforeach()

# for rare case where resolved ${CMAKE_INSTALL_LIBDIR}${PYMOD_INSTALL_LIBDIR} isn't nested, make libdir. may yet be empty
file(TO_NATIVE_PATH "${STAGED_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" _plain_libdir)
file(MAKE_DIRECTORY ${_plain_libdir})
