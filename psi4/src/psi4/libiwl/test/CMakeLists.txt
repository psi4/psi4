add_executable(iwl_roundtrip_test.x iwl_roundtrip_test.cc stubs.cc)
target_link_libraries(iwl_roundtrip_test.x PRIVATE iwl psio psi4util options ciomr)
target_include_directories(iwl_roundtrip_test.x PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_BINARY_DIR}/src
)
add_test(NAME iwl_roundtrip
         COMMAND iwl_roundtrip_test.x
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_tests_properties(iwl_roundtrip PROPERTIES LABELS "libiwl;quicktests")
