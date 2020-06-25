foreach(dir pasture
            adcc
            cct3
            psi4fockci
            gpu_dfcc
            snsmp2
            v2rdm_casscf
)

   add_subdirectory(${dir})
endforeach()

