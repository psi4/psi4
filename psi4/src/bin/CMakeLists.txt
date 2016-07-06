foreach(dir_name adc ccdensity ccenergy cceom cchbar cclambda ccresponse
                 ccsort cctransort cctriples dcft deriv_wrapper detci dfmp2
                 dfocc efp_interface findif fisapt fnocc mcscf mints_wrapper mrcc
                 occ optking psi4_main psimrcc sapt scf scfgrad thermo transqt2)
   add_subdirectory(${dir_name})
endforeach()

if(ENABLE_CHEMPS2)
    add_subdirectory(dmrg)
endif()


