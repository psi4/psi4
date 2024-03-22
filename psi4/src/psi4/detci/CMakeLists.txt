list(APPEND sources
  b2brepl.cc
  calc_d.cc
  civect.cc
  ciwave.cc
  compute_cc.cc
  compute_mpn.cc
  detci.cc
  diag_h.cc
  form_ov.cc
  get_mo_info.cc
  h0block.cc
  import_vector.cc
  ints.cc
  misc.cc
  mitrush_iter.cc
  odometer.cc
  og_addr.cc
  olsengraph.cc
  olsenupdt.cc
  opdm.cc
  params.cc
  printing.cc
  s1v.cc
  s2v.cc
  s3_block_bz.cc
  s3v.cc
  sem.cc
  set_ciblks.cc
  sigma.cc
  slater.cc
  slater_matel.cc
  slaterd.cc
  stringlist.cc
  tpdm.cc
  vector.cc
  )

set_source_files_properties(
  s1v.cc
  s2v.cc
  s3v.cc
  vector.cc
  tpdm.cc
  PROPERTIES
    COMPILE_FLAGS -O3
  )
psi4_add_module(bin detci sources)

if(TARGET gauxc::gauxc)
  target_include_directories(detci
    PUBLIC
      $<TARGET_PROPERTY:gauxc::gauxc,INTERFACE_INCLUDE_DIRECTORIES>
    )
  target_link_libraries(detci
    PUBLIC 
      gauxc::gauxc 
    )
endif()

target_link_libraries(detci
  PRIVATE
    pybind11::module
  )
