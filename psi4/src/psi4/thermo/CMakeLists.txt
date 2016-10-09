set(sources_list thermo.cc )
# plugin dependency here is fake, just so it gets linked
psi4_add_module(bin thermo sources_list mints plugin)
