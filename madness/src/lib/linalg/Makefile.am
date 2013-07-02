include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = test  test_solvers
lib_LIBRARIES = libMADlinalg.a

thisincludedir = $(includedir)/linalg
thisinclude_HEADERS = tensor_lapack.h clapack.h  lapack_functions.h \
                         solvers.cc solvers.h gmres.h

# cblas.h and cblas.cc temporarily removed

libMADlinalg_a_SOURCES = lapack.cc \
                         \
                         tensor_lapack.h clapack.h  lapack_functions.h \
                         solvers.cc solvers.h

test_SOURCES = test.cc
test_LDADD = libMADlinalg.a $(LIBTENSOR) $(LIBMISC) $(LIBWORLD)

test_solvers_SOURCES = test_solvers.cc
test_solvers_LDADD =  libMADlinalg.a $(LIBTENSOR) $(LIBMISC) $(LIBWORLD)

libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done
