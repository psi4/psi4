#At this point these variables should have been set:
#
#HAVE_GLOBAL_ARRAYS True if the user specified a global array package
#                    If true will assume that we also know JK_GA_INC
#                    and JK_GA_LIB, which respectively are the includes
#                    and library directories for GA
#
#JKROOT The path to the libJKFactory folder (absolute preferred)
#
#JK_BLAS_INC Where the header files for blas are (I'm assuming a lot
#             here, in particular that lapack, scalapack, and blas all
#             have header files given by the same include
#JK_BLAS_LIB The equivalent of JK_BLAS_INC, except for the libraries


#Build Global Arrays (if needed)
include(cmake/BuildGA.txt)
set(JK_GA_LIBS ${JK_GA_LIB}/libga.a
               ${JK_GA_LIB}/libarmci.a
)
#I'm going to assume the user doesn't have GTFock
#so let's build that now
add_library(GTFockInt STATIC MinimalInterface.h
                          MinimalInterface.cpp )
set_property(GLOBAL APPEND PROPERTY LIBLIST GTFockInt)
add_subdirectory(gtfock)
add_dependencies(GTFockInt GTFock mints)

set(JK_LIBS ${CMAKE_BINARY_DIR}/libGTFock.a)
set(JK_LIBS ${JK_LIBS} ${JK_GA_LIBS} )
add_library(GTFockCInt STATIC CifiedFxns.h
                             CifiedFxns.cc)
add_dependencies(GTFockCInt qt mints)
set_property(GLOBAL APPEND PROPERTY LIBLIST GTFockCInt)
