
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

# Build static library
add_library(interface_libefp STATIC ${sources_list})
set_property(TARGET interface_libefp PROPERTY COMPILE_FLAGS "-O1")

# Specify dependencies for the library (if any)
# Make sure psi4 links in
set_property(GLOBAL APPEND PROPERTY INTERFACES interface_libefp)

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/interface_libefp)

