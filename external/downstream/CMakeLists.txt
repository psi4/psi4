foreach(dir pasture
            gpu_dfcc
            snsmp2
            v2rdm_casscf
            adcc
            psi4fockci
)

   add_subdirectory(${dir})
endforeach()

