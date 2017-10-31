set(sources_list fd_misc.cc 
                 fd_geoms_freq_1.cc 
                 fd_geoms_1_0.cc 
                 fd_1_0.cc 
                 fd_geoms_freq_0.cc 
                 fd_freq_0.cc 
                 fd_freq_1.cc 
)
# plugin dependency here is fake, just so it gets linked
psi4_add_module(binlib findif sources_list mints plugin)
