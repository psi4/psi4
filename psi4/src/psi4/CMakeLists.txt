foreach (dir_name libpsi4util liboptions libparallel libciomr libpsio libqt
         libiwl libtrans libdpd libefp_solver libmints libdisp
         libfunctional libdiis libscf_solver libthce libfock lib3index
         libcubeprop libmoinfo libplugin libsapt_solver libpsipcm libfilesystem)
    add_subdirectory(${dir_name})
endforeach ()

foreach (dir_name adc ccdensity ccenergy cceom cchbar cclambda ccresponse
         cctransort cctriples dcft detci dfmp2
         dfocc efp_interface findif fisapt fnocc mcscf mrcc
         occ optking psimrcc sapt scfgrad thermo
         gdma_interface dmrg
        )
    add_subdirectory(${dir_name})
endforeach ()
