libcomex_la_SOURCES += src-mpi/comex.c
libcomex_la_SOURCES += src-mpi/comex_impl.h
libcomex_la_SOURCES += src-mpi/groups.c
libcomex_la_SOURCES += src-mpi/groups.h

AM_CPPFLAGS += -I$(top_srcdir)/src-mpi
