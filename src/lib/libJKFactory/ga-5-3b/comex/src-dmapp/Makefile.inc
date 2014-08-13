libcomex_la_SOURCES += src-dmapp/comex.c
libcomex_la_SOURCES += src-dmapp/comex_impl.h
#libcomex_la_SOURCES += src-dmapp/clusterinfo.c
#libcomex_la_SOURCES += src-dmapp/clusterinfo.h
libcomex_la_SOURCES += src-dmapp/groups.c
libcomex_la_SOURCES += src-dmapp/groups.h
libcomex_la_SOURCES += src-dmapp/reg_cache.c
libcomex_la_SOURCES += src-dmapp/reg_cache.h

AM_CPPFLAGS += -I$(top_srcdir)/src-dmapp
