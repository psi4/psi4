set(sources_list c_sort.cc d_sort.cc e_sort.cc fock.cc scf_check.cc a_spinad.cc cache.cc d_spinad.cc e_spinad.cc memcheck.cc sort_tei_rhf.cc b_spinad.cc cctransort.cc denom.cc f_sort.cc pitzer2qt.cc sort_tei_uhf.cc)
psi4_add_binary(cctransort sources_list dpd iwl trans options qt psio mints)
