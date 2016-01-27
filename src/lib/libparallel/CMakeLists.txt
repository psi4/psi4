set(headers_list "")
# List of headers
list(APPEND headers_list PsiFileImpl.h TableSpecs.h StreamBase.h ParallelScanner.h local.h elem.h openmp.h parallel.h PsiOutStream.h ParallelPrinter.h BasesBase.h PsiInStream.h mad.h serialize.h threaded_storage.h mpi_wrapper.h libparallel.h )

# If you want to remove some headers specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM headers_list "")
else()
   list(REMOVE_ITEM headers_list "")
endif()
if(NOT MPI_FOUND)
   list(REMOVE_ITEM headers_list "mpi_wrapper.h")
endif()
# Sort alphabetically
list(SORT headers_list)

set(sources_list "")
# List of sources
list(APPEND sources_list ParallelScanner.cc p_fprintf.cc PsiOutStream.cc BasesBase.cc process.cc ParallelPrinter.cc parallel.cc mpi_wrapper.cc PsiInStream.cc )

# If you want to remove some sources specify them explictly here
if(DEVELOPMENT_CODE)
   list(REMOVE_ITEM sources_list "")
else()
   list(REMOVE_ITEM sources_list "")
endif()
if(NOT MPI_FOUND)
   list(REMOVE_ITEM sources_list "mpi_wrapper.cc")
endif()


# Build static library
add_library(parallel STATIC ${sources_list})
# Specify dependencies for the library (if any)
add_dependencies(parallel mints)
set_property(GLOBAL APPEND PROPERTY LIBLIST parallel)
if(BUILD_CUSTOM_BOOST)
   add_dependencies(parallel custom_boost)
endif()

# Sets install directory for all the headers in the list
install_list_FILES("${headers_list}" include/psi4/libparallel)
