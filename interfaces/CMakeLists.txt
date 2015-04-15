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

if(${ENABLE_CHEMPS2})
    add_subdirectory(chemps2)
endif()

###### END externally stored, no-dependency projects ############
