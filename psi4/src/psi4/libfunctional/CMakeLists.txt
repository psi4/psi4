set(sources_list superfunctional.cc
                 LibXCfunctional.cc
                 factory.cc
                 functional.cc
)
psi4_add_module(lib functional sources_list disp)
target_link_libraries(functional PUBLIC Libxc::xc)
