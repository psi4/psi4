set(SRC corr_grad.cc mp2.cc wrapper.cc)
add_library(dfmp2 ${SRC})
add_dependencies(dfmp2 mints)
