set(sources_list wrapper.cc )
psi4_add_binary(sapt sources_list)
target_link_libraries(sapt mints sapt_solver)
