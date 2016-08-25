foreach(dir pybind11 libint ambit libefp PCMSolver chemps2 
            liberd dkh gdma
)
   add_subdirectory(${dir})
endforeach()

