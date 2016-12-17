foreach(dir ambit
            chemps2
            dkh
            libefp
            erd
            gdma
            libint
            pcmsolver
            pybind11
)
   add_subdirectory(${dir})
endforeach()

