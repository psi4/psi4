.. comment Note: This document contains light reStructuredText mark-up. 
   (Ignore the symbols .. :: ``.) It can be read here as plain-text or viewed in html at 
   http://sirius.chem.vt.edu/psi4manual/latest/installfile.html .


Installation Instructions for PSI4
==================================

* :ref:`I.   Compilation Prerequisites                                     <sec:install_I>`
* :ref:`II.  Brief Summary of Configuration, Compilation, and Installation <sec:install_II>`
* :ref:`III. Detailed Installation Instructions                            <sec:install_III>`
* :ref:`IV.  Recommendations for BLAS and LAPACK libraries                 <sec:install_IV>`
* :ref:`V.   Miscellaneous architecture-specific notes                     <sec:install_V>`
* :ref:`VI.  Common Problems with PSI Compilation                          <sec:install_VI>`


.. _`sec:install_I`:

I. Compilation Prerequisites
----------------------------

* C++ Compiler

* F77 Compiler (the F95 compiler, gfortran, with gcc-4.X will work)

  .. note:: The F77 compiler is only used to determine the symbol-naming
     convention of and some system routines for the BLAS and LAPACK libraries
     on a few architectures.  It is optional in a few cases (e.g. Mac OS X
     systems).

* Optimized BLAS library (preferably NOT one supplied by a standard
  Linux distribution; see recommendations at :ref:`Section IV <sec:install_IV>` below)

* Optimized LAPACK library (preferably NOT one supplied by a standard
  Linux distribution; see recommendations at :ref:`Section IV <sec:install_IV>` below)

* POSIX threads (Pthreads) library (optional)

* Python interpreter (2.6 or higher; psi4 is Python3 compatible)

  * psi4 is Python3 compatible intermittently every few months when tested.

* Python developer libraries corresponding to your interpreter.

  .. note:: To check to see if you have the Python developer libraries
     installed look for the Python config program. If your Python interpreter
     is named ``python`` look for the config program ``python-config``,
     likewise if your interpreter is ``python2.6``, the config program is
     ``python2.6-config``. If you cannot find the config program the
     developer libraries will not be detected and the PSI4 configure script
     will fail. The library is called ``python-devel`` for Fedora and
     ``python-dev`` for Ubuntu.

* GNU utilities: (see http://www.gnu.org)

  * make
  * autoconf (version 2.52 or higher)

    .. note:: autoconf is only needed for special
       architectures or to compile the development branch.

  * aclocal
  * fileutils (esp. install)

* For documentation: (available pre-built off http://www.psicode.org)

  * latex
  * sphinx (version 1.1 or higher)
  * dvipng (for latex math in html)

For Ubuntu users, you will need the following packages installed:
gfortran [for linking to BLAS/LAPACK], g++, autoconf, python-dev 


.. _`sec:install_II`:

II. Brief Summary of Configuration, Compilation, and Installation
-----------------------------------------------------------------

This section outlines the main steps of configuring, compiling, and
installing PSI.  More detail is given below in :ref:`Section III <sec:install_III>`.

A. Autoconf

   1. For unusual architectures (or for developers working in the development
      branch), one needs to first run autoconf to generate
      the file "configure" in the top-level psi4 directory.  For most Linux
      and Mac compilations, this should not be necessary because the configure
      file provided with PSI4 should be sufficient.  To replace the general
      configure file with one specific to your architecture, in the top-level
      psi4 directory, run autoconf::

         >>> autoconf

   2. Distributed-parallel compilation.

      Not recommended at this time except for developers.  Shared-memory
      parallelization is already enabled by default in the standard
      compilation.
      
      Distributed-parallel versions of PSI4 require madness. If you select
      mpicxx as the compiler, the distributed-parallel version (including
      madness) will compile. For distributed-parallel compilation, you must
      run the following command in the madness directory, otherwise the PSI4
      configure script will fail (autoreconf is provided by package autoconf,
      but it calls another program provided by libtool, so that package must
      also be installed)::

         >>> cd madness
         >>> autoreconf
         >>> cd ..

B. Configuration and Compilation

   Make an object directory in which you can compile the code ::

      >>> mkdir obj

   Next you need to configure the code. Find a configuration
   :ref:`option line <sec:install_III_1_configurelines>` or combination of
   configuration options at :ref:`Section III(1)A <sec:install_III_1_A>`.

   * Either, use the line directly::

        >>> cd obj
        >>> ../configure [your compilation configuration options here]

   * or, save your configuration options for a future compilation.
     In the top-level psi4 directory, create a file like "do-configure" with 
     the configure command and options on one line. ::

        >>> vi do-configure
        ../configure [your compilation configuration options here]
        >>> chmod u+x do-configure
        >>> cd obj
        >>> ../do-configure

   Compile the code, run the tests, and (if tests pass) install it. ::

      >>> make
      >>> make tests
      >>> make install

That's it!  The details about final user configuration are given below in 
:ref:`Section III(7) <sec:install_III_7>`.  If something goes wrong, 
check :ref:`Section VI <sec:install_VI>` about common compilation problems.


.. _`sec:install_III`:

III. Detailed Installation Instructions
---------------------------------------

This section provides a more detailed explanation of the procedure for
compiling and installing the PSI4 package.

* Step 1: Configuration

  A. General Information about Configuration

     First, we recommend that you choose for the top-level psi4 source
     directory something other than ``/usr/local/psi``; ``$HOME/psi4`` or
     ``/usr/local/src/psi4`` are convenient choices.  Next, in the top-level psi4
     source directory you've chosen, first run autoconf to generate the configure
     script from configure.ac.  It is best to keep the source code separate
     from the compilation area, so you must first choose a subdirectory for
     compilation of the codes.  A simple option is ``psi4/objdir``, which should
     work for most environments.  However, if you need executables for several
     architectures, you should choose more meaningful subdirectory names.

     .. note:: The compilation directory will be referred to as $objdir for the
        remainder of these instructions.

     In $objdir, run the configure script found in the PSI4 top-level source
     directory.  This script will scan your system to locate certain libraries,
     header files, etc. needed for complete compilation.  The script accepts a
     number of options, all of which are listed above.  The most important of
     these is the ``--prefix`` option, which selects the installation directory for
     the executables, the libraries, header files, basis set data, and other
     administrative files.  The default ``--prefix`` is ``/usr/local/psi``.

     .. note:: The configure script's ``--prefix`` directory will be referred to as
        $prefix for the remainder of these instructions.

     .. _`sec:install_III_1_A`:

     Besides ``--prefix``, PSI often needs a few additional options for the
     configure script.  To make it easy to recompile later (especially if
     you're a developer), it can be convenient (but not necessary) to to put
     the configure options in a small executable file, so you can re-do the
     configuration later very easily. Let us assume that we will be putting
     the configure options in a file named do-configure, in the top-level
     psi4 directory (we'll keep it up there instead of down in the compilation
     directory $objdir, so that if we delete the compilation directory later,
     we'll still have the do-configure file). All configure options must be
     on one line in the do-configure script.
     
     .. note:: The configure options below are for the most common architectures and
        compilers. The developers would appreciate it if you would share any special
        configuration options that might be needed for less commonly encountered
        situations. 
     
     For g++, if you have BLAS and LAPACK in standard locations (like ``/usr/lib64``),
     configuration is very easy. Pick one of the following scenarios, and place the
     text given in the psi4/do-configure file (all on one long line). Replace the
     text after prefix with whatever directory you want to use for your
     installation directory.

     .. _`sec:install_III_1_configurelines`:

     * Intel compiler with MKL math library [highly recommended; if you don't use
       this, then at least make sure you have a threaded BLAS (see BLAS
       recommendations at :ref:`Section IV <sec:install_IV>` below)] ::

          ../configure --prefix=/usr/local/psi4 --with-blas='-mkl' --with-cc=icc --with-cxx=icpc --with-fc=ifort  --with-opt='-O2 -static -no-prec-div' --with-incdirs=-mkl

       .. note:: It's ``-mkl``, not ``-lmkl``.

       .. warning:: A few users have reported errors with MKL 10.  Use at
          least version 11.

       .. warning:: There seems to be a problem with icpc 12.0.2 and possibly earlier
          12.0 versions, giving an error like::

             error: identifier "__is_trivial" is undefined.

          Use at least version 12.0.4.

     * Gnu compiler with ACML math library (better than MKL for AMD processors) ::

          ../configure --prefix=/usr/local/psi4 --with-cc=gcc --with-cxx=g++ --with-fc=gfortran --with-opt=-O2 --with-blas="-L/opt/acml5.2.0/gfortran64_mp/lib -lacml_mp" --with-lapack="-L/opt/acml5.2.0/gfortran64_mp/lib -lacml_mp"

     * g++, optimized ::

         ../configure --prefix=/usr/local/psi4
      
     * g++, for debugging ::

         ../configure --prefix=/usr/local/psi4 --without-opt --with-debug


     * Compiling for Mac

       PSI4 has been compiled on OS X 10.7 (Lion) and 10.8 (Mountain Lion). 
       To get the compilers needed, it's easiest to install Xcode.
       However, Xcode does not provide a Fortran compiler. Although
       Fortran compilers are not needed to compile Psi, a broken one can
       prevent Psi from configuring properly. Do not download the latest
       version of GFortran from the HPC website; this is unlikely to be
       compatible with your version of GCC. Instead, you should run ``gcc
       -v`` to find out what version of GCC you're using, and then
       download the corresponding GFortran from
       http://r.research.att.com/tools/.  If you configure Psi on a Mac
       without any Fortran compiler it will set itself up correctly, so
       this is only necessary if you want a Fortran compiler for other
       purposes. You can configure Psi by adding something like ::

          ../configure --with-plugins

       to the do-configure script. If you want to use the new LLVM compilers that
       ship with Xcode 4 (they compile quicker than GCC), use ::
       
          ../configure --with-plugins --with-cxx=llvm-g++

       .. warning:: If you still happen to encounter an error like::

             checking Fortran symbols... giving up
             configure: error: could not determine fortran symbol names

          adding the following tag to your configure may help ::

             --with-f77symbol=lcu

       .. warning:: An error like the one below has been seen
          when inadvertantly linking to 32-bit libraries ::

             Undefined symbols:
             "_omp_get_num_threads", referenced from:
                 __ZN3psi5dfmp26UDFMP28form_AiaEv.omp_fn.4 in libPSI_dfmp2.a(mp2.o)
                 ...

  B. List of Specific Configuration Options

     The example configuration options in the previous subsection are usually
     sufficient.  However, if not, you may need to make use of one or more
     of the following options to the configure script:

     * ``--prefix=directory`` --- Use this option if you wish to install the
       PSI4 package somewhere other than the default directory, ``/usr/local/psi``.
  
     * ``--with-cxx=compiler`` --- Use this option to specify a C++ compiler.
       One should use compilers that generate reentrant code, if possible.
       The default search order for compilers is: xlC_r (AIX only), g++, c++,
       icpc, cxx.  
  
     * ``--with-fc=compiler`` --- Use this option to specify a Fortran-77 compiler,
       which is used to determine linking coventions for BLAS and LAPACK libraries
       and to provide system routines for those libraries.  Note that no fortran
       compiler is necessary on Mac OS X systems (see below).  The default search
       order for compilers is: xlf_r (AIX only), gfortran, g77, ifort, f77, f2c.
  
     * ``--with-f77symbol=value`` --- This option allows manual assignment of the 
       FORTRAN77 symbol convention, which is necessary for C programs to link
       Fortran-interface libraries such as BLAS and LAPACK. This option should
       only be used by experts and even then should almost never be necessary. 
       Allowed values are:

       * lc  : lower-case
       * lcu : lower-case with underscore (default)
       * uc  : upper-case
       * ucu : upper-case with underscore
  
     * ``--with-ld=linker`` --- Use this option to specify a linker.  The
       default is 'ld'.
  
     * ``--with-ar=archiver`` --- Use this option to specify an archiver.  The
       default is to look for 'ar' automatically.
  
     * ``--with-ar-flags=flags`` --- Use this option to specify additional archiver 
       flags. The default is 'r'.
  
     * ``--with-incdirs=directories`` --- Use this option to specify extra
       directories where to look for header files. Directories should be specified
       prepended by ``-I``, i.e. ``-Idir1 -Idir2``, etc. If several directories are 
       specified, enclose the list with single right-quotes, e.g., ::

          --with-incdirs='-I/usr/local/include -I/home/psi4/include'
  
     * ``--with-libs=libraries`` --- Use this option to specify extra
       libraries which should be used during linking. Libraries should be 
       specified by their full names or in the usual ``-l`` notation, e.g. 
       ``-lm /usr/lib/libm.a``.  If several libraries are specified, enclose 
       the list with single right-quotes, e.g., ::

          --with-libs='-libm -lgcc_s'
  
     * ``--with-libdirs=directories`` --- Use this option to specify extra
       directories where to look for libraries. Directories should be specified
       prepended by ``-L``, e.g., ``-Ldir1 -Ldir2``. If several directories are 
       specified, enclose the list with single right-quotes, e.g., ::

          --with-libdirs='-L/usr/local/lib -I/home/psi4/lib'
  
     * ``--with-blas=library`` --- Use this option to specify a BLAS library.
       (Many BLAS libraries can be detected automatically.)
       If your BLAS library has multiple components, enclose the file list
       with single right-quotes, e.g., ::

          --with-blas='-lf77blas -latlas'
  
     * ``--with-lapack=library`` --- Use this option to specify a LAPACK library.
       (Many LAPACK libraries can be detected automatically.)
       If your LAPACK library has multiple components, enclose the file list
       with single right-quotes, e.g., ::

          --with-lapack='-llapack -lcblas -latlas'
  
     * ``--with-max-am-eri=integer`` --- Specifies the maximum angular momentum
       level for the primitive Gaussian basis functions when computing
       electron repulsion integrals.  This is set to h-type functions (AM=5)
       by default.
  
     * ``--with-max-am-deriv1=integer`` --- Specifies the maximum angular
       momentum level for first derivatives of the primitive Gaussian
       basis functions.  This is set to g-type functions (AM=4) by default.
  
     * ``--with-max-am-deriv2=integer`` --- Specifies the maximum angular
       momentum level for second derivatives of the primitive Gaussian
       basis functions.  This is set to f-type functions (AM=3) by default.
  
     * ``--with-debug=yes/no`` --- Turns on debugging flags (-g) if yes.  This is
       set to no by default.
  
     * ``--with-opt=yes/no`` --- Turns off compiler optimizations (-OX) if no.
       This is set to yes by default.
  
     * ``--with-strict=yes`` --- Turns on strict compiler warnings.

  C. Python interpreter

     Usually Python will be detected automatically.  If this fails, or if
     you have multiple versions installed and want to specify a particular
     one, set the PYTHON environmental variable to the full path name
     of the Python interpreter you want to use.  This defaults to the
     ``python`` in your path. For example, if you want to use
     ``python2.6`` located in /usr/bin set the environmental variable to be::

        PYTHON=/usr/bin/python2.6

     .. note:: If the variable PYTHON is set, the config program must be 
        present with a similar name. For instance, in the above example 
        the following must exist::

           /usr/bin/python2.6-config

     You either set the environmental variable before you call configure, or
     tell configure about it::

        ../configure PYTHON=/usr/bin/python2.6

  D. Boost Libraries

     PSI4 can use a user-provided boost C++ library, or, alternatively,
     build the boost version 1.53.0 that comes bundled with the distribution.
     By default, PSI4 will look in your include/library paths for
     a compatible and complete boost installation (boost 1.46 or newer). A
     boost installation in a nonstandard location can be specified by the
     ``--with-boost=PATH`` and ``--with-boost-libdir=PATH`` configure flags. If a
     default or user-specified boost installation is found to be incomplete,
     incompatible, or nonexistent, boost 1.53.0 will be unpacked automatically
     and built as part of the PSI4 build process.

     Required Compiled Boost Modules (all Boost 1.46.0 or later): 

     * Filesystem
     * Python
     * Regex
     * Serialization
     * System
     * Thread

     Relevant Configure Options:

     * ``--with-boost[=value]`` --- Use Boost library from a standard location
       if yes (default), from the specified location if <path>, or disable
       it if no.

     * ``--with-boost-libdir=directory`` ---
       Force given directory for boost libraries. Note that this will override
       library path detection, so use this parameter only if default library
       detection fails and you know exactly where your boost libraries are
       located. 
 
     * ``--with-boost-filesystem[=special-lib]`` ---
       Use the Filesystem library from boost. It is possible to specify a 
       certain library for the linker e.g., ::

          --with-boost-filesystem=boost_filesystem-gcc-mt

     * ``--with-boost-python`` --- Specify the boost python library or suffix to use.

     * ``--with-boost-regex[=special-lib]`` ---
       Use the Regex library from boost. It is possible to specify a certain
       library for the linker e.g., ::

          --with-boost-regex=boost_regex-gcc-mt-d-1_33_1

     * ``--with-boost-serialization[=special-lib]`` ---
       Use the Serialization library from boost. It is possible to specify a
       certain library for the linker e.g., ::

          --with-boost-serialization=boost_serialization-gcc-mt-d-1_33_1

     * ``--with-boost-system[=special-lib]`` ---
       Use the System library from boost. It is possible to specify a certain
       library for the linker e.g., ::

          --with-boost-system=boost_system-gcc-mt

     * ``--with-boost-thread[=special-lib]`` ---
       Use the Thread library from boost. It is possible to specify a certain
       library for the linker e.g., ::

          --with-boost-thread=boost_thread-gcc-mt


* Step 2: Compilation

  Running ``make`` (which must be GNU's 'make' utility) in $objdir will compile
  the PSI4 libraries and executable modules.

* Step 3: Testing

  To execute automatically the ever-growing number of test cases after
  compilation, simply execute ``make tests`` in the $objdir directory.
  This will run each (relatively small) test case and report the results.
  Failure of any of the test cases should be reported to the developers.
  By default, any such failure will stop the testing process.  If you desire
  to run the entire testing suit without interruption, execute ``make tests
  TESTFLAGS='-u -q'``. Note that you must do a ``make testsclean`` in $objdir
  to run the test suite again.

* Step 4: Installation

  Once testing is complete, installation into $prefix is accomplished by
  running ``make install`` in $objdir. Executable modules are installed in
  $prefix/bin, include files in $prefix/include, libraries in $prefix/lib, and 
  basis set data and various control structures in $prefix/share.

* Step 5: Building Documentation

  This is not recommended because all of the documentation should be
  available at http://sirius.chem.vt.edu/psi4manual/latest/index.html
  (link "docs" off http://www.psicode.org), and it is automatically updated.  However,
  if your system has the appropriate utilities (notably the sphinx package
  and LaTeX), you may build the package documentation from the top-level
  $objdir by running ``make doc``.  The resulting files will appear in the
  $prefix/doc area.

* Step 6: Cleaning

  All object files and libraries can be removed to save disk space by running
  ``make clean`` in $objdir.


.. _`sec:install_III_7`:

* Step 7: User Configuration

  After the PSI4 package has been successfully installed, the user will need
  to add the installation directory into his/her path.  If the package has
  been installed in the default location ``/usr/local/psi``, then in C shell,
  the user should add something like the following to their ``.cshrc`` file::

     setenv PSI /usr/local/psi
     set path = ($path $PSI/bin)

  Next, the user needs to tell the PSI4 I/O manager how to handle scratch files.
  Identify the path to a fast scratch disk for which the user has write access.  
  If the local ``/tmp`` volume is large enough, it might be used.
  However, a dedicated scratch volume (using RAID0 striping for speed) is
  recommended.

  .. warning:: Scratch should NOT be a NFS-mounted volume, as writes to a
     remote disk over the network can be very slow and can tie up the network
     and negatively impact other users.

  Specify scratch location by editing the ``.cshrc`` file to set the scratch 
  environment variable :envvar:`PSI_SCRATCH`. If the selected location is 
  ``/scratch/user``, add something like the following::

     setenv PSI_SCRATCH /scratch/user

  In a bash shell, the corresponding commands to be added to ``.bashrc`` is
  the following::

     export PSI=/usr/local/psi
     PATH=$PSI/bin:$PATH ; export PATH
     export PSI_SCRATCH=/scratch/user

  More advanced control of scratch files and is handled through a
  ``.psi4rc`` file, which is discussed at section :ref:`sec:psirc`.

  .. note:: For developers: during compilation and testing, PSI4 finds its basis sets,
     grids, etc., in ``psi4/lib``.  After installation, PSI4 will look in 
     $prefix/share/psi.  If you want to specify a non-standard location for this
     information, you can do this by setting the environmental variable
     $PSI4DATADIR to the directory containg the basis, grids, etc.,
     subdirectories.


.. _`sec:install_IV`:

IV. Recommendations for BLAS and LAPACK libraries
-------------------------------------------------

Much of the speed and efficiency of the PSI4 programs depends on the
corresponding speed and efficiency of the available BLAS and LAPACK libraries
(especially the former).  In addition, the most common compilation problems
involve these libraries.  Users may therefore wish to consider the following
BLAS and LAPACK recommendations when building PSI4:

(1) It is NOT wise to use the stock BLAS library provided with many
    Linux distributions like RedHat. This library is usually just the
    netlib distribution and is completely unoptimized. PSI4's
    performance will suffer if you choose this route. 

    The choice of LAPACK is less critical, and so the unoptimized
    netlib distribution is acceptable.  If you do choose to use the
    RedHat/Fedora stock BLAS and LAPACK, make sure that the blas-devel
    and lapack-devel packages are installed.

(2) Perhaps the best choice, if you have it available, is
    Intel's MKL library, which includes BLAS and LAPACK (note: use
    version 11 or later, we had reports of occasional errors using version 
    10).  MKL is efficient and works well in threaded mode.

    Otherwise, the simplest choice is to use ATLAS
    (http://math-atlas.sourceforge.net/), which is readily available
    on all Linux distributions. Another alternative is OpenBLAS
    (https://github.com/xianyi/OpenBLAS, formerly GotoBLAS). These
    work well on nearly every achitecture to which the PSI4 developers
    have access, though we have identified at least one case in which
    the Goto libraries yielded faulty DGEMM calls.  On Mac OS X
    systems, the vecLib package that comes with Xcode works well.

    If you prefer to use the ACML
    (http://developer.amd.com/tools/cpu-development/amd-core-math-library-acml/)
    we highly recommend using the latest version. Older versions
    of ACML have been known to cause problems.

.. _`sec:install_IV_3`:

(3) PSI4 does not require a Fortran compiler, unless the resident BLAS
    and LAPACK libraries require Fortran-based system libraries.  If you see
    compiler complaints about missing symbols like "do_fio" or "e_wsfe", then
    your libraries were most likely compiled with g77 or gfortran, which
    require ``-lg2c`` to resolve the Fortran I/O calls.  Use of the same gcc
    package for PSI4 should normally resolve this problem.

(4) The PSI4 configure script can often identify and use several
    different BLAS and LAPACK libraries, but its ability to do this
    automatically depends on a number of factors, including correspondence
    between the compiler used for PSI4 and the compiler used to build
    BLAS/LAPACK, placement of the libraries in commonly searched directories,
    etc. PSI4's configure script will find your BLAS and LAPACK if any of the
    the following are installed in standard locations (e.g. ``/usr/local/lib``):

    (a) ATLAS: ``libf77blas.a`` and ``libatlas.a``, plus netlib's ``liblapack.a``
    (b) MKL 8: ``libmkl.so`` and ``libmkl_lapack64.a`` (with the corresponding
        Intel compilers)
    (c) Goto: ``libgoto.a`` and netlib's ``liblapack.a``
    (d) Cray SCSL (e.g. on SGI Altix): ``libscs.so`` (NB: No Fortran compiler
        is necessary in this case, so ``--with-fc=no`` should work.)
    (e) ESSL (e.g. on AIX systems): ``libessl.a``


(5) If configure cannot identify your BLAS and LAPACK libraries
    automatically, you can specify them on the command-line using the
    ``--with-blas`` and ``--with-lapack`` arguments described above.  Here are a few
    examples that work on the PSI4 developers' systems:

    (a) Linux with ATLAS::

        --with-blas='-lf77blas -latlas' --with-lapack='-llapack -lcblas'

    (b) Mac OS X with vecLib::

        --with-blas='-altivec -framework vecLib' --with-lapack=' '

    (c) Linux with MKL 8.1 and icc/icpc/ifort 9.1::

        --with-libdirs=-L/usr/local/opt/intel/mkl/8.0.2/lib/32 --with-blas=-lmkl --with-lapack=-lmkl_lapack32

    (d) Linux on ia32 with MKL 10.1 and icc/icpc 11.0::

        --with-blas='-Wl,--start-group -L/usr/local/opt/intel/mkl/10.1.0.015/lib/32 -l mkl -Wl,--end-group -lguide -lpthread'

* Compilation notes for ATLAS

  These shortcut notes might be helpful if you are using Linux.  However,
  we recommend reading and following the full ATLAS installation notes.

  You'll need a Fortran compiler installed.   

  Unpack the source code, then make a compilation directory (could
  be an obj subdirectory in the source directory, or elsewhere).

  Turn off CPU throttling so the auto-tuning capabilities have a chance
  to work.  On Linux, this can be tune using ::

     /usr/bin/cpufreq-selector -g performance

  cd into the compilation directory and run the source
  directory configure script there, with any necessary flags, e.g., ::
    
     /usr/local/src/atlas/configure --prefix=/usr/local/atlas

  where prefix gives the installation directory.
  It should automatically detect if you're on an x86_64

  Then make and check using ::

     make; make check; make ptcheck

  And install ::

     make install
   
* Compilation notes for netlib's LAPACK

  These shortcut notes might be helpful if you are using Linux.  However,
  we recommend reading and following the full LAPACK installation notes.

  You'll need a Fortran compiler installed.

  If you decide to compile LAPACK from source, it may be obtained from 
  http://www.netlib.org/lapack/.  Unpack the source code, and in the
  top-level source directory, you need to create a make.inc file with
  the appropriate options for your machine.  For Linux/gfortran,
  simply ::
 
     cp make.inc.example make.inc

  Next, edit BLASLIB in make.inc to point to your BLAS library
  (full pathnames are recommended)::

     BLASLIB = /home/david/software/atlas3.9.25/lib/libf77blas.a /home/david/software/atlas3.9.25/lib/libatlas.a

  Edit Makefile as necessary (probably not needed). ::

     make

  Copy the resulting file [lapack_($ARCH).a] where you want it
  (a standard location like /usr/local/lib is easier for PSI to find).
  It is probably helpful to rename the file liblapack.a.
     

.. _`sec:install_V`:

V. Miscellaneous Architecture-Specific Notes
--------------------------------------------

* Linux on x86 and x86_64

  (1) Intel compilers: We had trouble with icpc 12.0.x.  Use 12.1 or
      later.

.. _`sec:install_VI`:

VI. Common Problems with PSI Compilation
----------------------------------------

* No rule to make target foo.h, needed by bar.d. Stop.

  This commonly happens after pulling updates from the repository. It happens
  when a library header file is removed or renamed by the update, but there are
  still old dependency files in the object directory, which think that they
  still need to know about that header. There's a simple remedy, just run ::

     >>> make DODEPEND=no dclean

  in the object directory.

* Make gets stuck in an infinite loop

  This means that the makefiles have not been properly updated. Running ::

     >>> autoconf

  in the top-level Psi directory, followed by ::

     >>> ./config.status --recheck
     >>> ./config.status

  in the object directory should fix it. This procedure will need to be run
  whenever an update changes the directory structure. 

* Incompatible g++/icpc

  The Intel compilers require an installed set of C++ headers. Unfortunately,
  the GNU compilers tend to be more cutting-edge than the Intel compilers,
  meaning that Intel is always playing catch-up to new features in g++. This
  means the two are often incompatible, leading to trouble if one wants to use
  icpc to compile PSI4 (or anything else...). Your best bet in general is to not
  upgrade Linux too fast, and always keep the very latest Intel compilers
  around.

* Missing symbols like "do_fio" or "e_wsfe"

  See :ref:`Section IV(3) <sec:install_IV_3>` above.


