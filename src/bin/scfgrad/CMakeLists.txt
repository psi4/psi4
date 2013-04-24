set(SRC jk_grad.cc response.cc scf_grad.cc wrapper.cc)
add_library(scfgrad ${SRC})
add_dependencies(scfgrad mints)
