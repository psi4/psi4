foreach(dir pybind11
            ambit
            chemps2
            dkh
            libecpint
            libefp
            erd
            gau2grid
            gdma
            libint
            pcmsolver
            pylibefp
            qcelemental
            simint
            libxc
)
   add_subdirectory(${dir})
endforeach()

