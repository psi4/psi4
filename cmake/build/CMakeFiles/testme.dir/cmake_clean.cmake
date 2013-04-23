FILE(REMOVE_RECURSE
  "CMakeFiles/testme.dir/demo.cpp.o"
  "testme.pdb"
  "testme"
)

# Per-language clean rules from dependency scanning.
FOREACH(lang CXX)
  INCLUDE(CMakeFiles/testme.dir/cmake_clean_${lang}.cmake OPTIONAL)
ENDFOREACH(lang)
