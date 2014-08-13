##############################################################################
# Makefile.am for ARMCI.
#
# Rationale:
# This Makefile.am follows many of the suggestions outlined in the paper
# "Recursive Make Considered Harmful".  We do not use Automake's
# 'include' feature (instead preferring a single, large Makefile.am).
#
# Additional targets:
# Besides the traditional make targets supplied by Automake, we have added the
# "examples" and "checkprogs" targets to build example programs and test
# programs, respectively.
#
# Notes:
# In general, each subdirectory has a corresponding section of this
# Makefile.am with the notable exception being the many examples getting
# rolled up into the examples section.
#

# The usual aclocal nonsense to get include paths right.
ACLOCAL_AMFLAGS = -I m4 -I ../m4

# All public headers, installed programs, test programs, and example programs
# are listed in these variables.  Appended to throughout.  These are the
# automake variables used.
include_HEADERS =
bin_PROGRAMS =
check_PROGRAMS =
check_LTLIBRARIES =
lib_LTLIBRARIES =
EXTRA_DIST = README
BUILT_SOURCES =
MOSTLYCLEANFILES =
CLEANFILES =
DISTCLEANFILES =
MAINTAINERCLEANFILES =
AM_FFLAGS =
AM_CFLAGS =
AM_CXXFLAGS =
AM_CPPFLAGS =
AM_LDFLAGS =
LDADD =

if ARMCI_NETWORK_ARMCI
else
lib_LTLIBRARIES += libarmci.la
endif
libarmci_la_SOURCES =
nodist_libarmci_la_SOURCES =
libarmci_la_LIBADD =

AM_FFLAGS += $(ARMCI_FOPT)
AM_FFLAGS += $(GA_F_WARN)
AM_FFLAGS += $(FFLAG_INT)
AM_CFLAGS += $(ARMCI_COPT)
AM_CFLAGS += $(GA_C_WARN)
AM_CXXFLAGS += $(ARMCI_CXXOPT)
AM_CXXFLAGS += $(GA_CXX_WARN)
AM_CPPFLAGS += $(GA_MP_CPPFLAGS)
AM_CPPFLAGS += $(ARMCI_NETWORK_CPPFLAGS)
if MSG_COMMS_MPI
else
AM_CPPFLAGS += -I$(top_srcdir)/tcgmsg
endif
LDFLAGS_BASE =
LDFLAGS_BASE += $(GA_MP_LDFLAGS)
LDFLAGS_BASE += $(ARMCI_NETWORK_LDFLAGS)
AM_LDFLAGS += $(LDFLAGS_BASE)

if ARMCI_NETWORK_ARMCI
LDADD += $(ARMCI_NETWORK_LIBS)
LDADD += $(GA_MP_LIBS)
else
LDADD += libarmci.la
LDADD += $(LTLIBOBJS)
LDADD += $(GA_MP_LIBS)
LDADD += $(ARMCI_NETWORK_LIBS)
endif

if ARMCI_NETWORK_ARMCI
else
libarmci_la_LIBADD += $(GA_MP_LIBS)
libarmci_la_LIBADD += $(ARMCI_NETWORK_LIBS)
endif

##############################################################################
# compiler and linker flags
#
# Important for external tools wanting to know how to link to ARMCI.
SED_NORMALIZE_WHITESPACE = $(SED) 's/ [ ]*/ /g;s/" /"/g;s/ "/"/g'
.PHONY: flags
flags:
	@echo ' =========================================================================== '
	@echo ' Suggested compiler/linker options are as follows.'
	@echo ' ARMCI libraries are installed in $(libdir)' 
	@echo ' ARMCI headers are installed in $(includedir)' 
	@echo ' '
	@echo ' CPPFLAGS="$(GA_MP_CPPFLAGS) $(ARMCI_NETWORK_CPPFLAGS) -I$(includedir)"' | $(SED_NORMALIZE_WHITESPACE)
	@echo ' '
	@echo ' LDFLAGS="$(GA_MP_LDFLAGS) $(ARMCI_NETWORK_LDFLAGS) -L$(libdir)"' | $(SED_NORMALIZE_WHITESPACE)
	@echo ' '
	@echo ' For C/C++ Programs: '
	@echo ' LIBS="-larmci $(GA_MP_LIBS) $(ARMCI_NETWORK_LIBS) $(MAYBE_FLIBS)"' | $(SED_NORMALIZE_WHITESPACE)
	@echo ' =========================================================================== '

##############################################################################
# src-portals
#
if ARMCI_NETWORK_PORTALS
include $(top_srcdir)/src-portals/Makefile.inc

##############################################################################
# src-gemini
#
else # !ARMCI_NETWORK_PORTALS
if ARMCI_NETWORK_GEMINI
include $(top_srcdir)/src-gemini/Makefile.inc

##############################################################################
# src
#
else # !ARMCI_NETWORK_GEMINI
libarmci_la_SOURCES += src/collectives/message.c
libarmci_la_SOURCES += src/common/aggregate.c
libarmci_la_SOURCES += src/common/armci.c
libarmci_la_SOURCES += src/common/ccopy.c
libarmci_la_SOURCES += src/common/clusterinfo.c
libarmci_la_SOURCES += src/common/pack.c
libarmci_la_SOURCES += src/common/iterator.c
libarmci_la_SOURCES += src/include/iterator.h
libarmci_la_SOURCES += src/ft/armci_chkpt.h
libarmci_la_SOURCES += src/ft/armci_storage.h
libarmci_la_SOURCES += src/include/acc.h
libarmci_la_SOURCES += src/include/armcip.h
libarmci_la_SOURCES += src/include/asm-ppc.h
libarmci_la_SOURCES += src/include/atomic_ops_ia64.h
libarmci_la_SOURCES += src/include/atomics-i386.h
libarmci_la_SOURCES += src/include/copy.h
libarmci_la_SOURCES += src/include/fujitsu-vpp.h
libarmci_la_SOURCES += src/include/kr_malloc.h
libarmci_la_SOURCES += src/include/locks.h
libarmci_la_SOURCES += src/include/memlock.h
libarmci_la_SOURCES += src/include/request.h
libarmci_la_SOURCES += src/include/semaphores.h
libarmci_la_SOURCES += src/include/shmalloc.h
libarmci_la_SOURCES += src/include/shmem.h
libarmci_la_SOURCES += src/include/shmlimit.h
libarmci_la_SOURCES += src/include/signaltrap.h
libarmci_la_SOURCES += src/include/spinlock.h
libarmci_la_SOURCES += src/include/utils.h
libarmci_la_SOURCES += src/locks/locks.c
libarmci_la_SOURCES += src/locks/memlock.c
libarmci_la_SOURCES += src/locks/mutex.c
libarmci_la_SOURCES += src/locks/semaphores.c
libarmci_la_SOURCES += src/memory/kr_malloc.c
libarmci_la_SOURCES += src/memory/memory.c
libarmci_la_SOURCES += src/progress/fence.c
libarmci_la_SOURCES += src/progress/wait.c
libarmci_la_SOURCES += src/xfer/caccumulate.c
libarmci_la_SOURCES += src/xfer/rmw.c
libarmci_la_SOURCES += src/xfer/strided.c
libarmci_la_SOURCES += src/xfer/vector.c
if MSG_COMMS_MPI
libarmci_la_SOURCES += src/common/groups.c
endif
if ARMCI_NETWORK_BGML
AM_CPPFLAGS += -I$(top_srcdir)/src/devices/bgml
libarmci_la_SOURCES += src/devices/bgml/bgmlacc.h
libarmci_la_SOURCES += src/devices/bgml/bgml.c
libarmci_la_SOURCES += src/devices/bgml/bgmldefs.h
libarmci_la_SOURCES += src/devices/bgml/bgml.h
endif
if ARMCI_NETWORK_CRAY_SHMEM
AM_CPPFLAGS += -I$(top_srcdir)/src/devices/cray-shmem
libarmci_la_SOURCES += src/memory/shmalloc.c
endif
if ARMCI_NETWORK_DCMF
AM_CPPFLAGS += -I$(top_srcdir)/src/devices/dcmf
libarmci_la_SOURCES += src/devices/dcmf/armcix.h
if DCMF_VER_2
libarmci_la_SOURCES += src/devices/dcmf/dcmf2/armcix_acc.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf2/armcix_fence.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf2/armcix_get.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf2/armcix_impl.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf2/armcix_lock.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf2/armcix_put.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf2/armcix_rmw.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf2/armcix_wait.c
endif
if DCMF_VER_0_2
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.2.0/armcix_acc.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.2.0/armcix_fence.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.2.0/armcix_get.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.2.0/armcix_impl.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.2.0/armcix_lock.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.2.0/armcix_put.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.2.0/armcix_rmw.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.2.0/armcix_wait.c
endif
if DCMF_VER_0_3
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.3.0/armcix_acc.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.3.0/armcix_fence.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.3.0/armcix_get.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.3.0/armcix_impl.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.3.0/armcix_lock.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.3.0/armcix_put.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.3.0/armcix_rmw.c
libarmci_la_SOURCES += src/devices/dcmf/dcmf-0.3.0/armcix_wait.c
endif
endif # ARMCI_NETWORK_DCMF
if ARMCI_NETWORK_LAPI
AM_CPPFLAGS += -I$(top_srcdir)/src/devices/lapi
libarmci_la_SOURCES += src/common/async.c
libarmci_la_SOURCES += src/common/request.c
libarmci_la_SOURCES += src/devices/lapi/lapi.c
libarmci_la_SOURCES += src/devices/lapi/lapidefs.h
libarmci_la_SOURCES += src/memory/buffers.c
endif
if ARMCI_NETWORK_MPI_MT
AM_CPPFLAGS += -I$(top_srcdir)/src/devices/mpi-mt
libarmci_la_SOURCES += src/common/ds-shared.c
libarmci_la_SOURCES += src/common/request.c
libarmci_la_SOURCES += src/common/spawn.c
libarmci_la_SOURCES += src/devices/mpi-mt/mpi2_client.c
libarmci_la_SOURCES += src/devices/mpi-mt/mpi2.h
libarmci_la_SOURCES += src/devices/mpi-mt/mpi2_server.c
libarmci_la_SOURCES += src/memory/buffers.c
endif
if ARMCI_NETWORK_MPI_SPAWN
AM_CPPFLAGS += -I$(top_srcdir)/src/devices/mpi-spawn
libarmci_la_SOURCES += src/common/ds-shared.c
libarmci_la_SOURCES += src/common/request.c
libarmci_la_SOURCES += src/devices/mpi-spawn/mpi2_client.c
libarmci_la_SOURCES += src/devices/mpi-spawn/mpi2.h
libarmci_la_SOURCES += src/devices/mpi-spawn/mpi2_server.c
libarmci_la_SOURCES += src/memory/buffers.c
endif
if ARMCI_NETWORK_OPENIB
AM_CPPFLAGS += -I$(top_srcdir)/src/devices/openib
libarmci_la_SOURCES += src/common/async.c
libarmci_la_SOURCES += src/common/ds-shared.c
libarmci_la_SOURCES += src/common/regions.c
libarmci_la_SOURCES += src/common/request.c
libarmci_la_SOURCES += src/common/spawn.c
libarmci_la_SOURCES += src/devices/openib/armci-vapi.h
libarmci_la_SOURCES += src/devices/openib/cbuf.c
libarmci_la_SOURCES += src/devices/openib/cbuf.h
libarmci_la_SOURCES += src/devices/openib/openib.c
libarmci_la_SOURCES += src/devices/openib/pendbufs.c
libarmci_la_SOURCES += src/devices/openib/pendbufs.h
libarmci_la_SOURCES += src/devices/openib/rtinfo.c
libarmci_la_SOURCES += src/memory/buffers.c
endif
if ARMCI_NETWORK_SOCKETS
AM_CPPFLAGS += -I$(top_srcdir)/src/devices/sockets
libarmci_la_SOURCES += src/common/async.c
libarmci_la_SOURCES += src/common/ds-shared.c
libarmci_la_SOURCES += src/common/request.c
libarmci_la_SOURCES += src/common/spawn.c
libarmci_la_SOURCES += src/devices/sockets/dataserv.c
libarmci_la_SOURCES += src/devices/sockets/sockets.c
libarmci_la_SOURCES += src/devices/sockets/sockets.h
libarmci_la_SOURCES += src/memory/buffers.c
endif
if ENABLE_CHECKPOINT
libarmci_la_SOURCES += src/ft/armci_chkpt.c
libarmci_la_SOURCES += src/ft/armci_storage.c
endif
if ARMCI_ENABLE_GPC_CALLS
libarmci_la_SOURCES += src/common/gpc.c
endif
if MACX
libarmci_la_SOURCES += src/common/signaltrap.c
libarmci_la_SOURCES += src/memory/shmalloc.c
libarmci_la_SOURCES += src/memory/winshmem.c
AM_CPPFLAGS += -DSHMEM -DMMAP
else
if CYGWIN
libarmci_la_SOURCES += src/common/signaltrap.c
libarmci_la_SOURCES += src/memory/shmalloc.c
libarmci_la_SOURCES += src/memory/winshmem.c
AM_CPPFLAGS += -DSHMEM -DMMAP
else
if SYSV
libarmci_la_SOURCES += src/common/signaltrap.c
libarmci_la_SOURCES += src/memory/shmem.c
libarmci_la_SOURCES += src/memory/shmlimit.c
endif
endif
endif
if NB_NONCONT
libarmci_la_SOURCES += src/common/noncont.c
endif
if THREAD_SAFE
libarmci_la_SOURCES += src/common/utils.c
libarmci_la_SOURCES += src/util/threads.c
endif
if NEC
libarmci_la_SOURCES += src/locks/tas-sx.s

src/locks/tas-sx.o: src/locks/tas-sx.s
	$(AM_V_CCAS) \
	$(CCAS) $(AM_CCASFLAGS) $(CCASFLAGS) \
	-c -o src/locks/tas-sx.o $(top_srcdir)/src/locks/tas-sx.s

src/locks/tas-sx.lo: src/locks/tas-sx.s
	$(AM_V_CCAS) \
	$(LIBTOOL) $(AM_V_lt) --tag=CC $(AM_LIBTOOLFLAGS) \
	$(LIBTOOLFLAGS) --mode=compile $(CCAS) $(AM_CCASFLAGS) \
	$(CCASFLAGS) -c -o src/locks/tas-sx.lo $(top_srcdir)/src/locks/tas-sx.s
endif

include_HEADERS += src/include/armci.h
include_HEADERS += src/include/gpc.h
include_HEADERS += src/include/message.h

AM_CPPFLAGS += -I$(top_srcdir)/src/include

endif # ARMCI_NETWORK_GEMINI
endif # ARMCI_NETWORK_PORTALS

##############################################################################
# profiling
#
if ENABLE_PROFILING
if HAVE_SYS_WEAK_ALIAS_PRAGMA
lib_LTLIBRARIES += libarmci_prof.la
libarmci_prof_la_SOURCES = 
libarmci_prof_la_SOURCES += tools/armci_prof.c
libarmci_la_SOURCES += src/common/capi.c
else # HAVE_SYS_WEAK_ALIAS_PRAGMA
libarmci_la_SOURCES += tools/armci_prof.c
endif # HAVE_SYS_WEAK_ALIAS_PRAGMA
else # ENABLE_PROFILING
libarmci_la_SOURCES += src/common/capi.c
endif

include_HEADERS += src/include/parmci.h

##############################################################################
# tcgmsg
#
EXTRA_DIST += tcgmsg/README

##############################################################################
# tcgmsg/ipcv4.0
#
if MSG_COMMS_TCGMSG4
AM_CPPFLAGS += -I$(top_srcdir)/tcgmsg/ipcv4.0

libarmci_la_SOURCES += tcgmsg/ipcv4.0/brdcst.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/checkbyte.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/cluster.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/cluster.h
libarmci_la_SOURCES += tcgmsg/ipcv4.0/defglobals.h
libarmci_la_SOURCES += tcgmsg/ipcv4.0/drand48.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/error.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/evlog.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/evlog.h
libarmci_la_SOURCES += tcgmsg/ipcv4.0/evon.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/globalop.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/llog.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/mdtob.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/mdtoi.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/mitob.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/mitod.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/mtime.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/niceftn.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/nnodes.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/nodeid.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/nxtval.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/pbegin.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/pbeginf.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/pfilecopy.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/sema.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/sema.h
libarmci_la_SOURCES += tcgmsg/ipcv4.0/setdbg.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/shmem.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/signals.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/signals.h
libarmci_la_SOURCES += tcgmsg/ipcv4.0/snd.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/sndrcv.h
libarmci_la_SOURCES += tcgmsg/ipcv4.0/sndrcvP.h
libarmci_la_SOURCES += tcgmsg/ipcv4.0/sockets.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/srftoc.h
libarmci_la_SOURCES += tcgmsg/ipcv4.0/srmover.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/stats.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/synch.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/usleep.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/waitall.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/waitcom.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/xdrstuff.c
libarmci_la_SOURCES += tcgmsg/ipcv4.0/xdrstuff.h

bin_PROGRAMS += parallel

parallel_SOURCES = tcgmsg/ipcv4.0/parallel.c

endif # MSG_COMMS_TCGMSG4

##############################################################################
# tcgmsg/ipcv5.0
#
if MSG_COMMS_TCGMSG5
AM_CPPFLAGS += -I$(top_srcdir)/tcgmsg/ipcv5.0

libarmci_la_SOURCES += tcgmsg/ipcv5.0/drand48.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/evlog.h
libarmci_la_SOURCES += tcgmsg/ipcv5.0/globalop.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/llog.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/mdtob.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/mdtoi.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/misc.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/mitob.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/mitod.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/mtime.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/niceftn.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/nnodes.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/nodeid.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/pbegin.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/pfilecopy.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/queues.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/signals.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/sndrcv.h
libarmci_la_SOURCES += tcgmsg/ipcv5.0/snd_rcv_probe.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/srftoc.h
libarmci_la_SOURCES += tcgmsg/ipcv5.0/tcgmsgP.h
libarmci_la_SOURCES += tcgmsg/ipcv5.0/usleep.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/waitall.c
if ARMCI_NETWORK_LAPI
libarmci_la_SOURCES += tcgmsg/ipcv5.0/async_send_lapi.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/lapi_putget.c
else
libarmci_la_SOURCES += tcgmsg/ipcv5.0/async_send.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/busy.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/error.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/nxtval.shm.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/pbeginf.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/shmem.c
libarmci_la_SOURCES += tcgmsg/ipcv5.0/synch.c
endif # ARMCI_NETWORK_LAPI

endif # MSG_COMMS_TCGMSG5

##############################################################################
# testing
#
check_PROGRAMS += testing/msgcheck
check_PROGRAMS += testing/perf
check_PROGRAMS += testing/perf2
check_PROGRAMS += testing/perf_aggr
check_PROGRAMS += testing/perf_nb
check_PROGRAMS += testing/shmclean
check_PROGRAMS += testing/shmtest
check_PROGRAMS += testing/simplelock
check_PROGRAMS += testing/test
check_PROGRAMS += testing/test2
if HAVE_ARMCI_NOTIFY
check_PROGRAMS += testing/testnotify
endif
if ENABLE_CHECKPOINT
check_PROGRAMS += testing/fttest
endif
if ARMCI_ENABLE_GPC_CALLS
check_PROGRAMS += testing/gpctest
endif
if MSG_COMMS_MPI
check_PROGRAMS += testing/simple
check_PROGRAMS += testing/test_groups
if SYSV
if ARMCI_NETWORK_ARMCI
else
check_PROGRAMS += testing/ipctest
endif
endif
endif
if THREAD_SAFE
check_PROGRAMS += testing/test_mt
endif
if ARMCI_NETWORK_PORTALS
check_PROGRAMS += testing/clone
check_PROGRAMS += testing/fork
check_PROGRAMS += testing/origptl
check_PROGRAMS += testing/ptltest
check_PROGRAMS += testing/te
endif
if HAVE_ARMCI_STRIDE_INFO_INIT
check_PROGRAMS += testing/testitr
endif

atsrc =
if HAVE_ARMCI_MSG_INIT
else
atsrc += testing/rpl_armci_msg_init.c
endif
if HAVE_ARMCI_MSG_FINALIZE
else
atsrc += testing/rpl_armci_msg_finalize.c
endif

ARMCI_SERIAL_TESTS =
ARMCI_SERIAL_TESTS_XFAIL =
ARMCI_PARALLEL_TESTS =
ARMCI_PARALLEL_TESTS_XFAIL =
ARMCI_TESTS = $(ARMCI_SERIAL_TESTS) $(ARMCI_PARALLEL_TESTS)
ARMCI_TESTS_XFAIL = $(ARMCI_SERIAL_TESTS_XFAIL) $(ARMCI_PARALLEL_TESTS_XFAIL)

ARMCI_PARALLEL_TESTS += testing/msgcheck$(EXEEXT)
ARMCI_PARALLEL_TESTS += testing/perf$(EXEEXT)
#ARMCI_PARALLEL_TESTS += testing/perf2(EXEEXT) # only a benchmark
ARMCI_PARALLEL_TESTS += testing/perf_aggr$(EXEEXT)
#ARMCI_PARALLEL_TESTS += testing/perf_nb$(EXEEXT) # needs bug fixes
#ARMCI_SERIAL_TESTS += testing/shmclean$(EXEEXT) # not sure we need this
#ARMCI_SERIAL_TESTS += testing/shmtest$(EXEEXT) # not sure we need this
ARMCI_PARALLEL_TESTS += testing/simplelock$(EXEEXT)
ARMCI_PARALLEL_TESTS += testing/test$(EXEEXT)
ARMCI_PARALLEL_TESTS += testing/test2$(EXEEXT)
if HAVE_ARMCI_NOTIFY
ARMCI_PARALLEL_TESTS += testing/testnotify$(EXEEXT)
endif
if ENABLE_CHECKPOINT
ARMCI_PARALLEL_TESTS += testing/fttest$(EXEEXT)
endif
if ARMCI_ENABLE_GPC_CALLS
ARMCI_PARALLEL_TESTS += testing/gpctest$(EXEEXT)
endif
if MSG_COMMS_MPI
ARMCI_PARALLEL_TESTS += testing/simple$(EXEEXT)
ARMCI_PARALLEL_TESTS += testing/test_groups$(EXEEXT)
if SYSV
if ARMCI_NETWORK_ARMCI
else
#ARMCI_PARALLEL_TESTS += testing/ipctest$(EXEEXT) # needs bug fixes
endif
endif
endif
if THREAD_SAFE
ARMCI_PARALLEL_TESTS += testing/test_mt$(EXEEXT)
endif
if ARMCI_NETWORK_PORTALS
#ARMCI_PARALLEL_TESTS += testing/clone$(EXEEXT)
#ARMCI_PARALLEL_TESTS += testing/fork$(EXEEXT)
#ARMCI_PARALLEL_TESTS += testing/origptl$(EXEEXT)
#ARMCI_PARALLEL_TESTS += testing/ptltest$(EXEEXT)
#ARMCI_PARALLEL_TESTS += testing/te$(EXEEXT)
endif
if HAVE_ARMCI_STRIDE_INFO_INIT
ARMCI_SERIAL_TESTS += testing/testitr$(EXEEXT)
endif

testing_clone_SOURCES       = testing/clone.c $(atsrc)
testing_fork_SOURCES        = testing/fork.c $(atsrc)
testing_fttest_SOURCES      = testing/fttest.c $(atsrc)
testing_gpctest_SOURCES     = testing/gpctest.c $(atsrc)
testing_ipctest_SOURCES     = testing/ipctest.c $(atsrc)
testing_msgcheck_SOURCES    = testing/msgcheck.c $(atsrc)
testing_origptl_SOURCES     = testing/origptl.c $(atsrc)
testing_perf_aggr_SOURCES   = testing/perf_aggr.c $(atsrc)
testing_perf_nb_SOURCES     = testing/perf_nb.c $(atsrc)
testing_perf_SOURCES        = testing/perf.c $(atsrc)
testing_perf2_SOURCES       = testing/perf2.c $(atsrc)
testing_ptltest_SOURCES     = testing/ptltest.c $(atsrc)
testing_shmclean_SOURCES    = testing/shmclean.c $(atsrc)
testing_shmtest_SOURCES     = testing/shmtest.c $(atsrc)
testing_simple_SOURCES      = testing/simple.c $(atsrc)
testing_simplelock_SOURCES  = testing/simplelock.c $(atsrc)
testing_te_SOURCES          = testing/te.c $(atsrc)
testing_test2_SOURCES       = testing/test2.c $(atsrc)
testing_test_groups_SOURCES = testing/test_groups.c $(atsrc)
testing_testitr_SOURCES     = testing/testitr.c $(atsrc)
testing_test_mt_SOURCES     = testing/test_mt.c $(atsrc)
testing_testnotify_SOURCES  = testing/testnotify.c $(atsrc)
testing_test_SOURCES        = testing/test.c $(atsrc)

testing_testitr_CPPFLAGS    = $(AM_CPPFLAGS) -I$(top_srcdir)/src/include

##############################################################################
# examples
#
check_PROGRAMS += examples/benchmarks/cg/armci_sharedmemory/cg
check_PROGRAMS += examples/benchmarks/lu/armci_blocking/lu-b-bc
check_PROGRAMS += examples/benchmarks/lu/armci_blocking/lu-block
check_PROGRAMS += examples/benchmarks/lu/armci_blocking/lu
check_PROGRAMS += examples/benchmarks/lu/armci_nonblocking/lu_nb_get
check_PROGRAMS += examples/benchmarks/lu/armci_nonblocking/lu_nb_put
check_PROGRAMS += examples/benchmarks/lu/pthreads/lu-thread
check_PROGRAMS += examples/features/aggregation/simple/simple
check_PROGRAMS += examples/features/aggregation/sparse_matvecmul/sparse_matvecmul
check_PROGRAMS += examples/features/non-blocking/simple/simple
check_PROGRAMS += examples/features/symmetric_memory/simple/simple
check_PROGRAMS += examples/simpleApps/transp1D-c
if HAVE_ARMCI_NOTIFY
check_PROGRAMS += examples/features/concurrency/simple/comdegree
check_PROGRAMS += examples/features/notification/simple/testnotify
endif
if ARMCI_ENABLE_GPC_CALLS
check_PROGRAMS += examples/features/gpc/hashtable/test_hashtable
if HAVE_BLAS
if HAVE_LAPACK
check_PROGRAMS += examples/features/proc_affinity/computation_impact/computation_impact
endif # HAVE_LAPACK
endif # HAVE_BLAS
check_PROGRAMS += examples/features/proc_affinity/simple/simple
endif # ARMCI_ENABLE_GPC_CALLS
if MSG_COMMS_MPI
check_PROGRAMS += examples/features/concurrency/multidma/multidma
check_PROGRAMS += examples/features/non-blocking/overlap/overlap
endif
if RANDOM_ACCESS
check_PROGRAMS += examples/benchmarks/RandomAccess/simple/simple
check_PROGRAMS += examples/benchmarks/RandomAccess/smp_bucket/smp_bucket
endif
if THREAD_SAFE
check_PROGRAMS += examples/benchmarks/lu/armci_multithreaded/lu-block-th
endif

#ARMCI_PARALLEL_TESTS += examples/benchmarks/cg/armci_sharedmemory/cg$(EXEEXT) # needs input file
ARMCI_PARALLEL_TESTS += examples/benchmarks/lu/armci_blocking/lu-b-bc$(EXEEXT)
ARMCI_PARALLEL_TESTS += examples/benchmarks/lu/armci_blocking/lu-block$(EXEEXT)
ARMCI_PARALLEL_TESTS += examples/benchmarks/lu/armci_blocking/lu$(EXEEXT)
ARMCI_PARALLEL_TESTS += examples/benchmarks/lu/armci_nonblocking/lu_nb_get$(EXEEXT)
ARMCI_PARALLEL_TESTS += examples/benchmarks/lu/armci_nonblocking/lu_nb_put$(EXEEXT)
ARMCI_PARALLEL_TESTS += examples/benchmarks/lu/pthreads/lu-thread$(EXEEXT)
ARMCI_PARALLEL_TESTS += examples/features/aggregation/simple/simple$(EXEEXT)
#ARMCI_PARALLEL_TESTS += examples/features/aggregation/sparse_matvecmul/sparse_matvecmul$(EXEEXT) # needs input file
ARMCI_PARALLEL_TESTS += examples/features/non-blocking/simple/simple$(EXEEXT)
if HAVE_ARMCI_NOTIFY
ARMCI_PARALLEL_TESTS += examples/features/concurrency/simple/comdegree$(EXEEXT)
ARMCI_PARALLEL_TESTS += examples/features/notification/simple/testnotify$(EXEEXT)
endif
ARMCI_PARALLEL_TESTS += examples/features/symmetric_memory/simple/simple$(EXEEXT)
ARMCI_PARALLEL_TESTS += examples/simpleApps/transp1D-c$(EXEEXT)
if ARMCI_ENABLE_GPC_CALLS
ARMCI_PARALLEL_TESTS += examples/features/gpc/hashtable/test_hashtable$(EXEEXT)
if HAVE_BLAS
if HAVE_LAPACK
ARMCI_PARALLEL_TESTS += examples/features/proc_affinity/computation_impact/computation_impact$(EXEEXT)
endif # HAVE_LAPACK
endif # HAVE_BLAS
ARMCI_PARALLEL_TESTS += examples/features/proc_affinity/simple/simple$(EXEEXT)
endif # ARMCI_ENABLE_GPC_CALLS
if MSG_COMMS_MPI
#ARMCI_PARALLEL_TESTS += examples/features/concurrency/multidma/multidma$(EXEEXT) # needs bug fixes
#ARMCI_PARALLEL_TESTS += examples/features/non-blocking/overlap/overlap$(EXEEXT) # needs bug fixes
endif
if RANDOM_ACCESS
ARMCI_PARALLEL_TESTS += examples/benchmarks/RandomAccess/simple/simple$(EXEEXT)
ARMCI_PARALLEL_TESTS += examples/benchmarks/RandomAccess/smp_bucket/smp_bucket$(EXEEXT)
endif
if THREAD_SAFE
ARMCI_PARALLEL_TESTS += examples/benchmarks/lu/armci_multithreaded/lu-block-th$(EXEEXT)
endif

examples_benchmarks_cg_armci_sharedmemory_cg_SOURCES = \
examples/benchmarks/cg/armci_sharedmemory/cg.c         \
examples/benchmarks/cg/armci_sharedmemory/compute.c    \
examples/benchmarks/cg/armci_sharedmemory/read_input.c \
examples/benchmarks/cg/armci_sharedmemory/timing.c $(atsrc)

examples_benchmarks_lu_armci_blocking_lu_SOURCES = \
examples/benchmarks/lu/armci_blocking/lu.c         \
examples/benchmarks/lu/armci_blocking/timing.c $(atsrc)

examples_benchmarks_lu_armci_blocking_lu_block_SOURCES = \
examples/benchmarks/lu/armci_blocking/lu-block.c         \
examples/benchmarks/lu/armci_blocking/timing.c $(atsrc)

examples_benchmarks_lu_armci_blocking_lu_b_bc_SOURCES =  \
examples/benchmarks/lu/armci_blocking/lu-b-bc.c          \
examples/benchmarks/lu/armci_blocking/timing.c $(atsrc)

examples_benchmarks_lu_armci_multithreaded_lu_block_th_SOURCES = \
examples/benchmarks/lu/armci_multithreaded/lu-block-th.c         \
examples/benchmarks/lu/armci_multithreaded/timing.c $(atsrc)

examples_benchmarks_lu_armci_nonblocking_lu_nb_put_SOURCES = \
examples/benchmarks/lu/armci_nonblocking/lu_nb_put.c         \
examples/benchmarks/lu/armci_nonblocking/timing.c $(atsrc)

examples_benchmarks_lu_armci_nonblocking_lu_nb_get_SOURCES = \
examples/benchmarks/lu/armci_nonblocking/lu_nb_get.c         \
examples/benchmarks/lu/armci_nonblocking/timing.c $(atsrc)

examples_benchmarks_lu_pthreads_lu_thread_SOURCES = \
examples/benchmarks/lu/pthreads/barrier.c           \
examples/benchmarks/lu/pthreads/barrier.h           \
examples/benchmarks/lu/pthreads/errors.h            \
examples/benchmarks/lu/pthreads/lu-thread.c         \
examples/benchmarks/lu/pthreads/timing.c $(atsrc)
examples_benchmarks_lu_pthreads_lu_thread_LDADD = -lpthread

examples_benchmarks_RandomAccess_simple_simple_SOURCES = \
examples/benchmarks/RandomAccess/simple/simple.c         \
examples/benchmarks/RandomAccess/timing.c $(atsrc)

examples_benchmarks_RandomAccess_smp_bucket_smp_bucket_SOURCES = \
examples/benchmarks/RandomAccess/smp_bucket/smp_bucket.c         \
examples/benchmarks/RandomAccess/timing.c $(atsrc)

examples_features_aggregation_simple_simple_SOURCES = \
examples/features/aggregation/simple/simple.c $(atsrc)

examples_features_aggregation_sparse_matvecmul_sparse_matvecmul_SOURCES = \
examples/features/aggregation/sparse_matvecmul/sparse_matvecmul.c $(atsrc)

examples_features_concurrency_multidma_multidma_SOURCES = \
examples/features/concurrency/multidma/multidma.c $(atsrc)

examples_features_concurrency_simple_comdegree_SOURCES = \
examples/features/concurrency/simple/comdegree.c $(atsrc)

examples_features_gpc_hashtable_test_hashtable_SOURCES = \
examples/features/gpc/hashtable/DistHashmap.cc           \
examples/features/gpc/hashtable/DistHashmap.h            \
examples/features/gpc/hashtable/GPCHashmap.cc            \
examples/features/gpc/hashtable/GPCHashmap.h             \
examples/features/gpc/hashtable/GPCHashmapHandler.cc     \
examples/features/gpc/hashtable/Hash_common.h            \
examples/features/gpc/hashtable/HashFunctions.cc         \
examples/features/gpc/hashtable/hash_map.h               \
examples/features/gpc/hashtable/HashUtil.cc              \
examples/features/gpc/hashtable/test_hashtable.cc        \
examples/features/gpc/hashtable/Util.h $(atsrc)
examples_features_gpc_hashtable_test_hashtable_LDFLAGS = $(LDFLAGS_BASE)
examples_features_gpc_hashtable_test_hashtable_LDADD = $(LDADD) $(MAYBE_FLIBS)

examples_features_non_blocking_overlap_overlap_SOURCES = \
examples/features/non-blocking/overlap/overlap.c $(atsrc)

examples_features_non_blocking_simple_simple_SOURCES = \
examples/features/non-blocking/simple/simple.c $(atsrc)

examples_features_notification_simple_testnotify_SOURCES = \
examples/features/notification/simple/testnotify.c $(atsrc)

examples_features_proc_affinity_computation_impact_computation_impact_SOURCES = \
examples/features/proc_affinity/computation_impact/computation_impact.c $(atsrc)
examples_features_proc_affinity_computation_impact_computation_impact_LDFLAGS = \
$(AM_LDFLAGS) \
$(LAPACK_LDFLAGS) \
$(BLAS_LDFLAGS)
examples_features_proc_affinity_computation_impact_computation_impact_LDADD = \
$(LDADD) \
$(LAPACK_LIBS) \
$(BLAS_LIBS)

examples_features_proc_affinity_simple_simple_SOURCES = \
examples/features/proc_affinity/simple/simple.c $(atsrc)

examples_features_symmetric_memory_simple_simple_SOURCES = \
examples/features/symmetric_memory/simple/simple.c $(atsrc)

examples_simpleApps_transp1D_c_SOURCES = \
examples/simpleApps/transp1D-c.c $(atsrc)

EXTRA_DIST += examples/benchmarks/lu/README
EXTRA_DIST += examples/features/aggregation/README
EXTRA_DIST += examples/features/concurrency/README
EXTRA_DIST += examples/features/gpc/hashtable/README
EXTRA_DIST += examples/features/non-blocking/README
EXTRA_DIST += examples/features/notification/README
EXTRA_DIST += examples/features/proc_affinity/README
EXTRA_DIST += examples/README

##############################################################################
# the end
#
.PHONY: checkprogs
checkprogs: $(check_PROGRAMS)

# support verbose/silent make rules for additional programs
# sed
SED_V = $(SED__v_$(V))
SED__v_ = $(SED__v_$(AM_DEFAULT_VERBOSITY))
SED__v_0 = @echo "  SED   " $@;
# ARMCI_TAS_AS
ARMCI_TAS_AS_V = $(ARMCI_TAS_AS__v_$(V))
ARMCI_TAS_AS__v_ = $(ARMCI_TAS_AS__v_$(AM_DEFAULT_VERBOSITY))
ARMCI_TAS_AS__v_0 = @echo "  GCC   " $@;
# ARMCI_X86COPY_AS
ARMCI_X86COPY_AS_V = $(ARMCI_X86COPY_AS__v_$(V))
ARMCI_X86COPY_AS__v_ = $(ARMCI_X86COPY_AS__v_$(AM_DEFAULT_VERBOSITY))
ARMCI_X86COPY_AS__v_0 = @echo "  GCC   " $@;

.PHONY: pfiles clean-pfiles
if MSG_COMMS_TCGMSG4
pfiles: $(check_PROGRAMS)
	for p in $(check_PROGRAMS); do \
	    echo "`whoami` `hostname` $(NPROCS) `pwd`/$$p /tmp" > $$p.p; \
	done
check: pfiles
clean-pfiles:
	find . -name '*.p' -exec rm {} \; ;
clean-local: clean-pfiles
endif # MSG_COMMS_TCGMSG4

##############################################################################
# test suite
#
# Some tests were commented out either because they required an input file or
# they were failing for unknown reasons and we didn't want to further debug.
#
SERIAL_TESTS =
SERIAL_TESTS += $(GAF2C_SERIAL_TESTS)
SERIAL_TESTS += $(TCGMSG_SERIAL_TESTS)
SERIAL_TESTS += $(ARMCI_SERIAL_TESTS)

SERIAL_TESTS_XFAIL =
SERIAL_TESTS_XFAIL += $(GAF2C_SERIAL_TESTS_XFAIL)
SERIAL_TESTS_XFAIL += $(TCGMSG_SERIAL_TESTS_XFAIL)
SERIAL_TESTS_XFAIL += $(ARMCI_SERIAL_TESTS_XFAIL)

PARALLEL_TESTS =
PARALLEL_TESTS += $(GAF2C_PARALLEL_TESTS)
PARALLEL_TESTS += $(TCGMSG_PARALLEL_TESTS)
PARALLEL_TESTS += $(ARMCI_PARALLEL_TESTS)

PARALLEL_TESTS_XFAIL =
PARALLEL_TESTS_XFAIL += $(GAF2C_PARALLEL_TESTS_XFAIL)
PARALLEL_TESTS_XFAIL += $(TCGMSG_PARALLEL_TESTS_XFAIL)
PARALLEL_TESTS_XFAIL += $(ARMCI_PARALLEL_TESTS_XFAIL)

TESTS =
TESTS += $(SERIAL_TESTS)
TESTS += $(PARALLEL_TESTS)

XFAIL_TESTS =
XFAIL_TESTS += $(SERIAL_TESTS_XFAIL)
XFAIL_TESTS += $(PARALLEL_TESTS_XFAIL)

if MSG_COMMS_TCGMSG4
LOG_COMPILER = \
maybe_parexec=`if echo "$(SERIAL_TESTS)" | $(GREP) "$$p" > /dev/null; then echo ""; else echo "$(TCGEXEC)"; fi`; eval $$maybe_parexec
else
if MSG_COMMS_TCGMSG5
LOG_COMPILER = tst="$$dir$$f -np $(NPROCS)"; eval
else
if CROSS_COMPILING
maybe_mpiexec=`if echo "$(SERIAL_TESTS)" | $(GREP) "$$p" > /dev/null; then echo "$(MPIEXEC)" | $(SED) 's/%NP%/1/'; else echo "$(MPIEXEC)" | $(SED) 's/%NP%/$(NPROCS)/'; fi`; eval $$maybe_mpiexec
else
LOG_COMPILER = \
maybe_mpiexec=`if echo "$(SERIAL_TESTS)" | $(GREP) "$$p" > /dev/null; then echo ""; else if echo "$(MPIEXEC)" | $(GREP) "%NP%" > /dev/null; then echo "$(MPIEXEC)" | $(SED) 's/%NP%/$(NPROCS)/'; else echo "$(MPIEXEC)"; fi; fi`; eval $$maybe_mpiexec
endif # CROSS_COMPILING
endif # MSG_COMMS_TCGMSG5
endif # MSG_COMMS_TCGMSG4

.PHONY: check-gaf2c
check-gaf2c: $(all-am) $(check_LTLIBRARIES) $(GAF2C_TESTS) $(MAYBE_PFILES)
	$(MAKE) $(AM_MAKEFLAGS) TESTS="$(GAF2C_TESTS)" XFAIL_TESTS="$(GAF2C_TESTS_XFAIL)" check-TESTS

.PHONY: check-tcgmsg
check-tcgmsg: $(all-am) $(check_LTLIBRARIES) $(TCGMSG_TESTS) $(MAYBE_PFILES)
	$(MAKE) $(AM_MAKEFLAGS) TESTS="$(TCGMSG_TESTS)" XFAIL_TESTS="$(TCGMSG_TESTS_XFAIL)" check-TESTS

.PHONY: check-armci
check-armci: $(all-am) $(check_LTLIBRARIES) $(ARMCI_TESTS) $(MAYBE_PFILES)
	$(MAKE) $(AM_MAKEFLAGS) TESTS="$(ARMCI_TESTS)" XFAIL_TESTS="$(ARMCI_TESTS_XFAIL)" check-TESTS
