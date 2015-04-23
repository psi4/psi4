set_property(GLOBAL PROPERTY INTERFACES)
###### START internally stored, psi4-dependency projects ########

# <<<  libefp  >>>

add_subdirectory(libefp)

###### END internally stored, psi4-dependency projects ##########

###### START externally stored, no-dependency projects ##########

# <<<  dummyplugin  >>>

if(${USEEXT_DUMMY_PLUGIN})
    add_subdirectory(dummyplugin)
endif()

# <<<  CheMPS2  >>>
# handled by ConfigChemps2

###### END externally stored, no-dependency projects ############
