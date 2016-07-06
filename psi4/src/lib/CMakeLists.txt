foreach(dir_name libpsi4util liboptions libparallel libciomr libpsio libqt
                 libiwl libtrans libdpd libefp_solver libmints libdisp
                 libfunctional libdiis libscf_solver libthce libfock lib3index
                 libcubeprop libmoinfo libplugin libsapt_solver)
   add_subdirectory(${dir_name})
endforeach()
if(ENABLE_PCMSOLVER)
   add_subdirectory(libpsipcm)
endif()
