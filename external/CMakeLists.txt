foreach(dir pybind11 libint ambit libefp pcmsolver chemps2 
            liberd dkh gdma
)
   add_subdirectory(${dir})
endforeach()

