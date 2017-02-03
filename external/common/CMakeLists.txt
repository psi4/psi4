foreach(dir lapack
)
   add_subdirectory(${dir})
endforeach()

if(ENABLE_Ambit OR ENABLE_CheMPS2)
    add_subdirectory(hdf5)
endif()
