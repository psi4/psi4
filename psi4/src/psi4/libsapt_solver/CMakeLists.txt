list(APPEND sources
  amplitudes.cc
  disp20.cc
  disp21.cc
  disp22sdq.cc
  disp22t.cc
  disp2ccd.cc
  disp30.cc
  elst10.cc
  elst12.cc
  elst13.cc
  exch-disp20.cc
  exch-disp30.cc
  exch-ind-disp30.cc
  exch-ind20.cc
  exch-ind30.cc
  exch10.cc
  exch11.cc
  exch12.cc
  fdds_disp.cc
  ind-disp30.cc
  ind20.cc
  ind22.cc
  ind30.cc
  sapt.cc
  sapt0.cc
  sapt2.cc
  sapt2p.cc
  sapt2p3.cc
  usapt0.cc
  utils.cc
  )
psi4_add_module(lib sapt_solver sources)

if(TARGET gauxc::gauxc)
  target_include_directories(sapt_solver
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
      #/storage/home/hhive1/dpoole34/data/Codes/GauXC/include/ 
    )
  target_link_libraries(sapt_solver
    PUBLIC 
      gauxc::gauxc 
    )
endif()


