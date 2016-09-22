set(sources_list fd_misc.cc 
                 fd_geoms_freq_1.cc 
                 fd_geoms_1_0.cc 
                 fd_1_0.cc 
                 fd_geoms_freq_0.cc 
                 fd_freq_0.cc 
                 fd_freq_1.cc 
)
psi4_add_module(bin findif sources_list mints)
