set(sources_list mxmb.cc 
                 eigout.cc 
                 init_array.cc 
                 init_matrix.cc 
                 block_matrix.cc 
                 tqli.cc 
                 flin.cc 
                 sq_rsp.cc 
                 add_arr.cc 
                 int_array.cc 
                 ffile.cc 
                 zero.cc 
                 print_mat.cc 
                 tri_to_sq.cc 
                 rsp.cc 
                 add_mat.cc 
                 mmult.cc 
                 dot.cc 
                 tred2.cc 
                 tstart.cc 
                 eivout.cc 
                 eigsort.cc 
                 ludcmp.cc 
                 print_array.cc 
                 long_int_array.cc 
                 sq_to_tri.cc 
                 lubksb.cc
)

add_library(ciomr ${sources_list})
target_link_libraries(ciomr parallel)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        DESTINATION psi4/src/lib/libciomr
        FILES_MATCHING PATTERN "*.h")
