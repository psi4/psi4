foreach(dir pasture
            gpu_dfcc
            qcelemental
            snsmp2
            v2rdm_casscf
)

   add_subdirectory(${dir})
endforeach()

