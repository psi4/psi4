set(sources_list PBE_Cfunctional.cc
                 PZ81_Cfunctional.cc
                 PW91_Cfunctional.cc
                 wpbec_functional.cc
                 VWN3_Cfunctional.cc
                 superfunctional.cc
                 LYP_Cfunctional.cc
                 FT97B_Xfunctional.cc
                 wpbex_functional.cc
                 FT97_Cfunctional.cc
                 VWN5_Cfunctional.cc
                 PW92_Cfunctional.cc
                 P86_Cfunctional.cc
                 LibXCfunctional.cc
                 cfunctional.cc
                 factory.cc
                 utility.cc
                 functional.cc
                 xfunctional.cc
)
psi4_add_module(lib functional sources_list disp)
target_link_libraries(functional PUBLIC Libxc::xc)
#target_link_libraries(functional PUBLIC ${LIBXC_LIBRARIES})
#target_include_directories(functional PUBLIC ${LIBXC_INCLUDE_DIRS})
