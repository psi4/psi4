list(APPEND sources
  a_spinad.cc
  b_spinad.cc
  c_sort.cc
  cache.cc
  cctransort.cc
  d_sort.cc
  d_spinad.cc
  denom.cc
  e_sort.cc
  e_spinad.cc
  f_sort.cc
  fock.cc
  memcheck.cc
  pitzer2qt.cc
  scf_check.cc
  sort_tei_rhf.cc
  sort_tei_uhf.cc
  )

psi4_add_module(bin cctransort sources)

if(TARGET PCMSolver::pcm)
    target_link_libraries(cctransort
      PUBLIC
        PCMSolver::pcm
      )
endif()
