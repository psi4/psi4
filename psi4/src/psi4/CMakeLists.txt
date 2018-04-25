foreach (dir_name libpsi4util liboptions libciomr libpsio libqt
         libiwl libtrans libdpd libefp_solver libmints libdisp
         libfunctional libdiis libscf_solver libfock lib3index
         libcubeprop libmoinfo libplugin libsapt_solver libpsipcm libfilesystem)
    add_subdirectory(${dir_name})
endforeach ()

foreach (dir_name adc ccdensity ccenergy cceom cchbar cclambda ccresponse
         cctransort cctriples dcft detci dfmp2 dfep2
         dfocc efp_interface findif fisapt fnocc mcscf mrcc
         occ optking psimrcc sapt scfgrad
         gdma_interface dmrg ambit_interface
        )
    add_subdirectory(${dir_name})
endforeach ()
