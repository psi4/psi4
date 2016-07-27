foreach(dir_name adc ccdensity ccenergy cceom cchbar cclambda ccresponse
                 ccsort cctransort cctriples dcft deriv_wrapper detci dfmp2
                 dfocc efp_interface findif fisapt fnocc mcscf mints_wrapper mrcc
                 occ optking psimrcc sapt scf scfgrad thermo transqt2
                 gdma_interface dmrg psi4_main
)
   add_subdirectory(${dir_name})
endforeach()
