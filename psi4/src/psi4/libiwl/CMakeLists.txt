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
                 #buf_toend.cc 
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
psi4_add_module(lib iwl sources_list)
target_include_directories(iwl PUBLIC options)
#if(APPLE)
#   target_link_libraries(iwl PRIVATE options "-Wl,-force_load")
#else()
#   target_link_libraries(iwl PRIVATE options "-Wl,--whole-archive")
#endif()
