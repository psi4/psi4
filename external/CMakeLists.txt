foreach(dir boost libint libderiv ambit libefp PCMSolver chemps2 
            liberd dkh gdma
)
   add_subdirectory(${dir})
endforeach()

