foreach(dir lapack
            hdf5
)
    add_subdirectory(${dir})
endforeach()
