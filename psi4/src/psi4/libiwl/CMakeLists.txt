set(sources_list buf_close.cc
                 buf_flush.cc
                 buf_put.cc
                 buf_wrt_mat.cc
                 wrtone.cc
                 buf_fetch.cc
                 buf_init.cc
                 buf_wrt.cc
                 buf_wrt_val.cc
                 rdone.cc
)
psi4_add_module(lib iwl sources_list)
target_include_directories(iwl PUBLIC options)
#if(APPLE)
#   target_link_libraries(iwl PRIVATE options "-Wl,-force_load")
#else()
#   target_link_libraries(iwl PRIVATE options "-Wl,--whole-archive")
#endif()
