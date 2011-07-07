SUBDIRS = nanophoto

include $(top_srcdir)/config/MakeGlobal.am

bin_PROGRAMS = depth embedded_dirichlet

depth_SOURCES = depth.cc

embedded_dirichlet_SOURCES = embedded_dirichlet.cc test_problems.h

LDADD = $(MRALIBS)

