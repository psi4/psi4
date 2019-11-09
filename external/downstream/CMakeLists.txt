foreach(dir pasture
            gpu_dfcc
            snsmp2
            v2rdm_casscf
            adcc
)

   add_subdirectory(${dir})
endforeach()

