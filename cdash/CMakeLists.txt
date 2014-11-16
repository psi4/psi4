#
# Generate the valgrind suppressions file.
#

file(GLOB suppression_files ${CMAKE_SOURCE_DIR}/cdash/*-suppressions.txt)

file(REMOVE ${CMAKE_BINARY_DIR}/valgrind-suppressions.txt)

foreach(fname ${suppression_files})
    file(READ ${fname} vgs)
    file(APPEND ${CMAKE_BINARY_DIR}/valgrind-suppressions.txt ${vgs})
endforeach()

# vim:et:sw=4:ts=4:
