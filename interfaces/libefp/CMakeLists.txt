
# Directory path spacer between this CMakeLists.txt file and the code that it's building 
set(SD "libefp/src")

set(headers_list "")
# List of headers
set(tmp balance.h bvec.h clapack.h compat.h efp.h elec.h int.h int_shift.h log.h math_util.h private.h stream.h swf.h terms.h util.h)

foreach(fl ${tmp})
    list(APPEND headers_list ${SD}/${fl})
endforeach(fl)

list(SORT headers_list)

set(sources_list "")
# List of sources
set(tmp ai_disp.c balance.c bvec.c clapack.c disp.c efp.c elec.c elec_terms.c int.c log.c parse.c pol.c pol_direct.c stream.c swf.c util.c xr.c)

foreach(fl ${tmp})
    list(APPEND sources_list ${SD}/${fl})
endforeach(fl)

# Write list of files to be passed to cloc for counting lines of code.
# Only files that are actually compiled are counted.
set(to_count "${sources_list}" "${headers_list}")
write_to_cloc_list("${to_count}")

# Build static library
add_library(interface_libefp STATIC ${sources_list})

# Specify dependencies for the library (if any)
set(libs_to_merge interface_libefp ${libs_to_merge} PARENT_SCOPE)
install(TARGETS interface_libefp ARCHIVE DESTINATION lib)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/interface_libefp)

# Make sure psi4 links in
set(PSILIB ${PSILIB} ${CMAKE_BINARY_DIR}/lib/libinterface_libefp.a PARENT_SCOPE)

