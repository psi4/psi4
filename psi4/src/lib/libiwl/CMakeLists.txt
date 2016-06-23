set(sources_list buf_rd_arr2.cc 
                 buf_wrt_mp2.cc 
                 buf_wrt_mp2r12a.cc 
                 buf_rd_all_mp2r12a.cc 
                 rdtwo.cc 
                 buf_fetch.cc 
                 rdone.cc 
                 buf_wrt_val.cc 
                 buf_wrt_arr2.cc 
                 buf_rd_arr.cc 
                 buf_wrt_val_SI.cc 
                 buf_wrt_arr_SI_nocut.cc 
                 buf_rd_all_act.cc 
                 buf_close.cc 
                 buf_toend.cc 
                 buf_wrt_mat.cc 
                 buf_wrt_all.cc 
                 buf_wrt_arr.cc 
                 buf_wrt.cc 
                 buf_flush.cc 
                 buf_wrt_arr_SI.cc 
                 buf_rd.cc 
                 sortbuf.cc 
                 wrtone.cc 
                 buf_init.cc 
                 wrttwo.cc 
                 buf_put.cc 
                 buf_rd_all.cc 
)
add_library(iwl ${sources_list})
target_include_directories(iwl PUBLIC options)
target_link_libraries(iwl options)

if(APPLE)
   set_property(GLOBAL APPEND PROPERTY LIBLIST -Wl,-force_load;iwl)
else()
   set_property(GLOBAL APPEND PROPERTY LIBLIST -Wl,--whole-archive;iwl;-Wl,--no-whole-archive)
endif()
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} 
        DESTINATION psi4/src/lib/libiwl
        FILES_MATCHING PATTERN "*.h")
