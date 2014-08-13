libarmci_la_SOURCES += src-gemini/acc.h
libarmci_la_SOURCES += src-gemini/aggregate.c
libarmci_la_SOURCES += src-gemini/armci.c
libarmci_la_SOURCES += src-gemini/armci-onesided.c
libarmci_la_SOURCES += src-gemini/armci-onesided.h
libarmci_la_SOURCES += src-gemini/armcip.h
libarmci_la_SOURCES += src-gemini/atomics-i386.h
libarmci_la_SOURCES += src-gemini/buffers.c
libarmci_la_SOURCES += src-gemini/caccumulate.c
libarmci_la_SOURCES += src-gemini/ccopy.c
libarmci_la_SOURCES += src-gemini/clusterinfo.c
libarmci_la_SOURCES += src-gemini/copy.h
libarmci_la_SOURCES += src-gemini/ds-shared.c
libarmci_la_SOURCES += src-gemini/fence.c
libarmci_la_SOURCES += src-gemini/groups.c
libarmci_la_SOURCES += src-gemini/kr_malloc.c
libarmci_la_SOURCES += src-gemini/kr_malloc.h
libarmci_la_SOURCES += src-gemini/locks.c
libarmci_la_SOURCES += src-gemini/locks.h
libarmci_la_SOURCES += src-gemini/memlock.c
libarmci_la_SOURCES += src-gemini/memlock.h
libarmci_la_SOURCES += src-gemini/memory.c
libarmci_la_SOURCES += src-gemini/message.c
libarmci_la_SOURCES += src-gemini/mutex.c
libarmci_la_SOURCES += src-gemini/pack.c
libarmci_la_SOURCES += src-gemini/pendbufs.h
libarmci_la_SOURCES += src-gemini/request.c
libarmci_la_SOURCES += src-gemini/request.h
libarmci_la_SOURCES += src-gemini/rmw.c
libarmci_la_SOURCES += src-gemini/rtinfo.c
libarmci_la_SOURCES += src-gemini/semaphores.c
libarmci_la_SOURCES += src-gemini/semaphores.h
libarmci_la_SOURCES += src-gemini/shmalloc.h
libarmci_la_SOURCES += src-gemini/shmem.c
libarmci_la_SOURCES += src-gemini/shmem.h
libarmci_la_SOURCES += src-gemini/shmlimit.c
libarmci_la_SOURCES += src-gemini/shmlimit.h
libarmci_la_SOURCES += src-gemini/signaltrap.c
libarmci_la_SOURCES += src-gemini/signaltrap.h
libarmci_la_SOURCES += src-gemini/sockets.h
libarmci_la_SOURCES += src-gemini/spawn.c
libarmci_la_SOURCES += src-gemini/spinlock.h
libarmci_la_SOURCES += src-gemini/strided.c
libarmci_la_SOURCES += src-gemini/utils.h
libarmci_la_SOURCES += src-gemini/vector.c
if ARMCI_ENABLE_GPC_CALLS
libarmci_la_SOURCES += src-gemini/gpc.c
endif
if THREAD_SAFE
libarmci_la_SOURCES += src-gemini/threads.c
libarmci_la_SOURCES += src-gemini/utils.c
endif

include_HEADERS += src-gemini/armci.h
include_HEADERS += src-gemini/gpc.h
include_HEADERS += src-gemini/message.h

AM_CPPFLAGS += -I$(top_srcdir)/src-gemini
AM_CPPFLAGS += -I$(top_srcdir)/src/include
LDADD += -lnumatoolkit
LDADD += -lonesided
