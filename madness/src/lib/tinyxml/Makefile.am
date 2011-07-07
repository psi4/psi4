include $(top_srcdir)/config/MakeGlobal.am

lib_LIBRARIES = libMADtinyxml.a

thisincludedir = $(includedir)/tinyxml
thisinclude_HEADERS = tinystr.h  tinyxml.h

libMADtinyxml_a_SOURCES = tinystr.cc  tinyxml.cc  tinyxmlerror.cc  tinyxmlparser.cc \
                          \
                          tinystr.h  tinyxml.h

AM_CPPFLAGS += -I.

libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done
