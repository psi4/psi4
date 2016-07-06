set(sources_list model_space_build.cc 
                 slater_determinant.cc 
                 moinfo_model_space.cc 
                 moinfo_mappings.cc 
                 moinfo_slaterdeterminant.cc 
                 moinfo.cc 
                 model_space.cc 
                 moinfo_scf.cc 
                 moinfo_base.cc 
)
psi4_add_library(moinfo sources_list mints)
