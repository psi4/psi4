###### START internally stored, psi4-dependency projects ########

# <<<  libefp  >>>

add_subdirectory(libefp)

###### END internally stored, psi4-dependency projects ##########

###### START externally stored, no-dependency projects ##########

# <<<  dummyplugin  >>>

if(${USEEXT_DUMMY_PLUGIN})
    add_subdirectory(dummyplugin)
endif()

###### END externally stored, no-dependency projects ############

# send new dependency info up the directory tree
set(PSILIB ${PSILIB} PARENT_SCOPE)

