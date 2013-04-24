set(SRC dftsapt.cc infsapt.cc wrapper.cc)
add_library(dftsapt ${SRC})
add_dependencies(dftsapt mints)
