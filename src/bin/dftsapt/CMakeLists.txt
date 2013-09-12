set(SRC asapt.cc atomic.cc  dftsapt.cc infsapt.cc vis.cc     wrapper.cc)
add_library(dftsapt ${SRC})
add_dependencies(dftsapt mints)
