lib_LTLIBRARIES += libarmci.la

libarmci_la_SOURCES =
libarmci_la_SOURCES += src-armci/armci.c
libarmci_la_SOURCES += src-armci/groups.c
libarmci_la_SOURCES += src-armci/message.c
libarmci_la_SOURCES += src-armci/iterator.c
libarmci_la_SOURCES += src-armci/iterator.h

libarmci_la_LIBADD = libcomexi.la

include_HEADERS += src-armci/armci.h
include_HEADERS += src-armci/message.h

noinst_LTLIBRARIES += libarmcii.la

libarmcii_la_SOURCES = $(libarmci_la_SOURCES)
