include $(top_srcdir)/config/MakeGlobal.am

lib_LIBRARIES = libMADmisc.a

bin_PROGRAMS = interp3

thisincludedir = $(includedir)/misc/
thisinclude_HEADERS = misc.h ran.h phandler.h cfft.h


libMADmisc_a_SOURCES = checksum_file.cc position_stream.cc gprofexit.cc ran.cc cfft.cc\
                       misc.h ran.h cfft.h phandler.h

interp3_SOURCES = interp3.cc

libraries:	$(lib_LIBRARIES)
	for dir in $(SUBDIRS) ; do $(MAKE) -C $$dir $@ ; done
