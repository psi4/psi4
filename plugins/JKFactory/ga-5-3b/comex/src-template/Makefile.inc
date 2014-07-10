# NOTE: Replace 'template' with appropriate name of new src-xxx directory.
libcomex_la_SOURCES += src-template/comex.c
libcomex_la_SOURCES += src-template/comex_impl.h
libcomex_la_SOURCES += src-template/groups.c
libcomex_la_SOURCES += src-template/groups.h

AM_CPPFLAGS += -I$(top_srcdir)/src-template
