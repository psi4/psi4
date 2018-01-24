foreach(dir pybind11
            ambit
            chemps2
            dkh
            libefp
            erd
            gau2grid
            gdma
            libint
            pcmsolver
            simint
            libxc
)
   add_subdirectory(${dir})
endforeach()

