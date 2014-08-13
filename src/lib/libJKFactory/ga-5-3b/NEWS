5.0   November 2010
    * Now built using GNU autotools (autoconf,automake,libtool)
    * Restricted arrays (see user manual)
    * ARMCI runtime enhancements
        + On-demand connection management
        + Improved scalability for fence
    * New Python interface
    * Task Scheduling Library (tascel)

5.0b  July 2010
    * Now built using GNU autotools (autoconf,automake,libtool)

4.3   May 2010
    * Optimized portals port to scale upto 200K procs
    * Optimized OpenIB port
    * BlueGene/P
    * Support for Sparse Data Operations
        + (See GA user manual - Chapter 11 for more details)

4.2   July 2009
    * Support for several new platforms
    * Optimized portals port for Cray XT5
    * BlueGene/P
    * Optimized OpenIB port
    * Support for Sparse Data Operations
        + (See GA user manual - Chapter 11 for more details)

4.1   May 2008
    * Support for several new platforms
        + Cray XT4
        + BlueGene/L, BlueGene/P
        + OpenIB network
    * Optimized one-sided non-blocking operations
    * New networks. i.e. ARMCI_NETWORK
        + OPENIB
        + PORTALS
        + MPI-SPAWN (one-sided communication thru' MPI2 Dynamic
          Process management and Send/Recv)

4.0   April 2006
    * Support for multilevel parallelism: processor group awareness
    * GA_Dgemm matrix multiplication based on SRUMMA
    * Support for large arrays (a terabyte Global Array now possible)
    * Optimized one-sided non-blocking operations
    * Supports various platforms (Crays, IBM SPs, SGI Altix, ...) and
      interconnects (Myrinet, Quadrics, Infiniband, ...)

3.4   2004 (beta)
    * Initial support for multilevel parallelism
    * processor group awareness
    * nonblocking one-sided interfaces
    * optimized port to Cray X1, SGI Altix

3.3   July, 2003
    * nonblocking one-sided interfaces
    * optimized matrix multiplication (ga_dgemm)
    * mirrored arrays
    * unoptimized port to Cray X1
    * support for infiniband network (Mellanox VAPI)
    * GA CCA component
    
3.2   October, 2002
    * support for ghost cells
    * port to Linux64/Itanium, Hitachi, Cray SV1, MAC X, and Cygwin 
    * support two underscores in fortran names on Linux/Cygwin
      (enabled by setting F2C_TWO_UNDERSCORES)
    * C++ bindings to Global Arrays
    * element-wise and matrix operations 
    * synchronization control in collective operations
    * additional datatypes: float & long in C, real in Fortran
    * cluster information calls
    * n-dim DRA

3.1   November 6th, 2000
    * ports to 64-bit IBM, HP, Sun, Fujitsu, Linux/Alpha platforms
    * new operations:
        + periodic version of put,get,accumulate
        + nga_scatter_acc
        + nga_select_elem
        + ga_nblock
        + nga_select_elem
    * restructured makefiles to simplify changes and enable
      export of make definitions to application makefiles
    * new platforms supported by ARMCI (and thus GA)
        + Myrinet on Linux (Intel, Sparc) and Solaris/Sparc
        + Giganet on Linux/Intel
        + Compaq/Quadrics clusters 
        + clusters of Windows NT under TCP/IP

3.0   November 15th, 1999
    * This is the most significant update of the GA package to date: 
      all source code for the core functionality of the library 
      has been rewritten, ~90% of the library code is new/changed; 
    * added support for n-dimensional arrays (n<=7);
    * added several n-dimensional test programs;
    * updated utility tools: ga_print, ga_summarize, which
      support n-dimensional operations;
    * improved interoperability with MPI on clusters of
      workstations
    * ports to Cray J90 and Windows NT (single machine)
    * ga_create uses default distribution for block size <1
      rather than <2 (the old scheme can be restored see
      OLD_DEFAULT_BLK in src/GNUmakefile)
    * added a new, easier-to-use C interface (see ntestc.c) 
    * added User's Manual

2.4
    * port to IBM SP running LAPI active messages
      (significant performance improvement).
      Version >=2.3 of the IBM PSSP parallel environment  and
      threaded version of MPI are required to run with LAPI.
    * port to Fujitsu VX-4 and VPP-700.
    * minor platform specific updates. 

2.3
    * added lock operations: ga_lock, ga_unlock,
      ga_create_mutex, ga_destroy_mutex
    * ga_acc works for integer data type
    * port to the CRAY T3E (with streams disabled)
    * environment variable TARGET_CPU can be use to guide
      compilation for particular CPU (R10000 and R8000 are
      currently recognized on the SGI systems)
    * a new TARGET, SGI_N32 used for the 64-bit SGI systems

2.2
    * added support for double complex datatype;
    * port to HP and Convex SPP;
    * latency optimizations for IBM SP running AIX 4.X;
    * relaxed restriction on data types in ga_dscal,
      ga_dscal_patch, ga_dadd, ga_dadd_patch,
      ga_ifill_patch, ga_dfill_patch - these functions has
      been renamed to: ga_scale, ga_scale_patch, ga_add,
      ga_add_patch, ga_fill_patch;
    * new functions: ga_fence and ga_init_fence, ga_zdot,
      and ga_zdot_patch, ga_print_stats;
    * new interface functions to ScaLAPACK: ga_cholesky,
      ga_llt_solve, ga_spd_invert, ga_solve
    * several configuration parameters are exposed to the
      user in global/src/config.h to allow customization of
      the package.

2.1
    * selection of message-passing (TCGMSG or MPI) library
    * new operations: ga_mpi_communicator and ga_proc_topology;
    * performance optimizations;
    * alternative C interface and C test programs;
    * port to Sun Solaris and PC running Linux.

2.0
    * support for networks of multiprocessors (SMP);
    * (tightly) integrated shared and distributed memory versions;
    * ability to restrict memory usage in global arrays;
    * ga_scatter/gather supports integers;
    * a program to test performance of GA primitives;
    * different numbering of global array handles (< 0);
    * and other internal improvements, like randomized order
      of remote operations, prototypes of C-callable ops, etc.

1.3.1
    * port to CRAY-T3D;
    * new operation ga_summarize;
    * improved ga_dgemm.

1.3
    * introduced patch versions of ga_ddot, ga_dscal, ga_dadd;
    * also ga_matmul_patch, ga_ifill_patch, ga_dfill_patch, and ga_duplicate. 
    * ga_copy_patch now has a transpose argument.

1.2
    * (loosely) integrated shared and distributed memory versions
