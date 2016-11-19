set(sources_list superfunctional.cc
                 LibXCfunctional.cc
                 factory.cc
                 utility.cc
                 functional.cc
)
psi4_add_module(lib functional sources_list disp)
target_link_libraries(functional PUBLIC Libxc::xc)
#target_link_libraries(functional PUBLIC ${LIBXC_LIBRARIES})
#target_include_directories(functional PUBLIC ${LIBXC_INCLUDE_DIRS})
