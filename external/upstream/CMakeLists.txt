foreach(dir pybind11 libint ambit libefp pcmsolver chemps2 
            erd dkh gdma
)
   add_subdirectory(${dir})
endforeach()

