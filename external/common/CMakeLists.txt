add_subdirectory(lapack)
if(ENABLE_Ambit OR ENABLE_CheMPS2)
    add_subdirectory(hdf5)
endif()
