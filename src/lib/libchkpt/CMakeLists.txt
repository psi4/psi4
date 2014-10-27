set(SRC am2canon_shell_order.cc atom_dummy.cc atom_pos.cc cartrep.cc
       ccvecs.cc cdsalc2cd.cc cdsalcpi.cc close.cc clsdpi.cc contr.cc
       contr_full.cc disp.cc disp_irrep.cc e_labeled.cc e_t.cc
       eccsd.cc ecorr.cc efzc.cc emp2.cc enuc.cc eref.cc escf.cc
       etot.cc evals.cc exist_add_prefix.cc exist.cc exps.cc felement.cc 
       fgeom.cc fock.cc
       fragment_coeff.cc frzcpi.cc frzvpi.cc geom.cc grad.cc ict.cc init.cc
       iopen.cc irr_labs.cc keyword.cc label.cc lagr.cc max_am.cc
       nallatom.cc nallatom_per_fragment.cc nao.cc natom.cc
       natom_per_fragment.cc ncalcs.cc nfragment.cc nfzc.cc nfzv.cc
       nirreps.cc nmo.cc nprim.cc nref_per_fragment.cc nshell.cc nso.cc
       nsymhf.cc num_unique_atom.cc num_unique_shell.cc openpi.cc orbspi.cc
       override_occ.cc phase_check.cc prefix.cc puream.cc ref.cc
       rot_symm_num.cc rotconst.cc rottype.cc rref.cc scf.cc shell_transm.cc
       shells_per_am.cc sloc.cc sloc_new.cc snuc.cc snumg.cc sopi.cc sprim.cc
       statespi.cc stype.cc sym_label.cc symoper.cc ua2a.cc us2s.cc usotao.cc
       usotbf.cc zmat.cc zvals.cc)
add_library(chkpt ${SRC})
if(BUILD_CUSTOM_BOOST)
   add_dependencies(chkpt custom_boost)
endif()   
