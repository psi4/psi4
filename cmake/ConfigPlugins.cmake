# Include this module **after** the ConfigCompilerFlags.cmake module
# compiler flags are set there, according to build type.
# If extra flags were provided they are included in CMAKE_CXX_FLAGS
# while CMAKE_CXX_FLAGS_DEBUG, CMAKE_CXX_FLAGS_RELEASE and CMAKE_CXX_FLAGS_PROFILE
# are set to the defaults.
# If custom flags were provided, they are set into the CMAKE_CXX_FLAGS variable
# while the build type-dependent flags are empty.
if(CMAKE_BUILD_TYPE STREQUAL "debug")
   set(CXX_FLAGS_PLUGIN "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
elseif(CMAKE_BUILD_TYPE STREQUAL "release")
   set(CXX_FLAGS_PLUGIN "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE STREQUAL "profile")
   # CMAKE_CXX_FLAGS_PROFILE contains CMAKE_CXX_FLAGS_RELEASE
   set(CXX_FLAGS_PLUGIN "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_PROFILE}")
else()
   set(CXX_FLAGS_PLUGIN "${CMAKE_CXX_FLAGS}")
endif()
