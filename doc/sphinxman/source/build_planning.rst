

.. _`faq:cmakeviasetup`:

=====================================================
Planning: how to configure |PSIfour| and invoke CMake
=====================================================

|PSIfour| is built through CMake. Prior to 1.1, |PSIfour| had a Python
script ``setup`` as a frontend to CMake, but this is no more, and
``cmake`` is now invoked directly. An abbreviated build guide can be found
within the source itself at
https://github.com/psi4/psi4/blob/master/CMakeLists.txt#L13-L106.

CMake does a good job scanning your computer to locate libraries, header
files, and executables needed for compilation. So it's very possible that
from :samp:`{$top-level-psi4-dir}` source directory, you can call :samp:`cmake -H.
-B{objdir}` without any further arguments, and it will invoke ``cmake``,
detect some appropriate defaults, configure the build, create a build
directory :samp:`{objdir}`, and complete, leaving you to only invoke
``make`` in the new build directory.

Should this happy scenario not come to pass, or if the default build
options are not to your taste, use the links within [core
dependencies](2_Planning#coredepend) and [add-on
dependencies](2_Planning#addondepend) to plan a set of arguments to
``setup`` tailored to your computer. The following topics may also be
helpful.

* [How to see what build configuration options are available](2_Planning#setuphelp)
* [How to compile elsewhere than ``$top-level-psi4-dir/objdir``](2_Planning#setupobjdir)
* :ref:`bfaq:setupprefix`
* [How to compile for debugging](2_Planning#setuptype)
* [How to configure code to use high angular momentum basis sets](2_Planning#setupmaxameri)
* [How to set CMake and Preprocessor options through the ``setup`` script](2_Planning#setupd)
* [How to set up a profiling build](2_Planning#profiling)


.. _`faq:coredepend`:

What are the tools and dependencies strictly required for building |PSIfour|
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The core Psi4 build requires the software below. Note that CMake and Python can be satisfied through the [Conda Psi4dependencies Package](zArchive_Conda#getting-and-using-the-psi4dependencies-package). The links below give examples of how to configure that software for Psi4 and any notes and warnings pertaining to it.

* [C++ Compiler](9_CXX)

* [F77 Compiler](9_Fortran)

  * only used to determine the symbol-naming convention for BLAS and LAPACK libraries
  * optional for Mac OS X

* Optimized [BLAS and LAPACK libraries](9_BlasLapack) (preferably NOT one supplied by a standard
  Linux distribution)

* Selected [Boost libraries](9_Boost) (1.55 or higher; optional in that bundled with Psi4 distribution)

* [Python interpreter and corresponding developer libraries](9_Python) (2.7)

* [CMake](9_CMake) [(3.0 or higher)](http://www.cmake.org/download/)

* [NumPy](9_NumPy) (needed at runtime, not buildtime)

* System utilities

  * GNU make
  * GNU install
  * POSIX threads (Pthreads) library



.. _`faq:addondepend`:

What are the add-on capabilities for |PSIfour| and what are their dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each of the items below is an independent additional capability that can
be built with Psi4. Sub-items below are the respective additional
dependencies of the add-on. Select which, if any, you want, and examine
the links for appropriate enable/disable arguments to ``setup``. Note that
the default for each may be enabled, disabled, or enabled if dependencies
already satisfied, so you may need to explicitly disable add-ons in to
tailor your build.

* Psi4 Testing

  * [CTest](9_CTest)
  * Perl (for some coupled-cluster tests) <http://perl.org>

* Psi4 Documentation (available pre-built at http://www.psicode.org/psi4manual/master/index.html)

  * LaTeX <http://latex-project.org>
  * Sphinx (1.1, 1.2, or 1.3) <http://sphinx-doc.org>
  * dvipng (for LaTeX math in html) <http://savannah.nongnu.org/projects/dvipng>
  * Perl (for some auto-documentation scripts) <http://perl.org>

* [ERD](9_ERD) (in place of LibInt) [_what is this?_](http://psicode.org/psi4manual/master/erd.html)

  * [Fortran Compiler](9_Fortran)

* [PCMSolver](9_PCMSolver) [_what is this?_](http://psicode.org/psi4manual/master/pcmsolver.html)

  * [Fortran Compiler](9_Fortran)
  * [Zlib Library](9_ZLib)

* [CheMPS2](9_CheMPS2) [_what is this?_](http://psicode.org/psi4manual/master/chemps2.html)

  * [GSL Library](9_GSL)
  * [HDF5](9_HDF5)


.. _`faq:setupmaxameri`:

How to configure code to use high angular momentum basis sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :ref:`libint <addon:libint>` integral code handles arbitrary order
angular momentum, but compiling that is prohibitive. The default of ``5``
is generally good. ``6`` has met all of a research group's needs for
years. ``4`` is handy for quickly testing other parts of the build.

* Build with Higher Angular Momentum

  .. code-block:: bash

   >>> cmake -DMAX_AM_ERI=6

* Relevant CMake Options:

  .. code-block:: bash

   MAX_AM_ERI=N        # The maximum angular momentum level (1=p, 2=d, 3=f,
                       # etc.) for the libint integrals and derivative
                       # integrals. A value of N implies a maximum first
                       # derivative of N-1, and maximum second derivative of
                       # N-2, so for an atom such as Neon, the default 5 gets
                       # you conventional cc-pV5Z for energies, cc-pVQZ for
                       # gradients, cc-pVTZ for frequencies and density-fitted
                       # cc-pVQZ for energies, cc-pVTZ for gradients, cc-pVDZ
                       # for frequencies. [default: 5]



.. .. code-block:: cmake
.. 
..     cmake_minimum_required(VERSION 2.8.12)
..     project(example)
.. 
..     add_subdirectory(pybind11)
..     pybind11_add_module(example example.cpp)




.. _`faq:setuphelp`:
###<a name="setuphelp"></a> How to see what build configuration options are available

Run ``./setup --help`` to get the summary below (ca. June 2015) and more.

    >>> cd $top-level-psi4-dir
    >>> ./setup --help
    usage: setup [-h] [--cc STRING] [--cxx STRING] [--fc STRING]
             [--max-am-eri MAX_ANGULAR_MOMENTUM]
             [--type [{release,debug,profile}]] [--prefix PATH] [--show]
             [--cmake STRING] [--boost-incdir PATH] [--boost-libdir PATH]
             [--python PYTHON] [--mpi] [--sgi-mpt] [--omp]
             [--mkl [{sequential,parallel,cluster}]]
             [--blas [{auto,builtin,none,/full/path/lib.a}]]
             [--lapack [{auto,builtin,none,/full/path/lib.a}]]
             [--extra-math-flags STRING] [--accelerate] [--cray] [--csr]
             [--scalapack] [--scalasca] [--cxx11 [{on,off}]]
             [--plugins [{on,off}]] [--suffix STRING] [--check] [--memcheck]
             [--coverage] [--static] [--unit-tests] [--vectorization]
             [-D STRING] [--host STRING] [--generator STRING] [--timings]
             [--asan | --msan | --tsan | --ubsan] [--erd {on,off}]
             [--jkfactory {on,off}] [--gpu-dfcc {on,off}]
             [--dummy-plugin {on,off}] [--pcmsolver {on,off}]
             [--chemps2 {on,off}] [--chemps2-dir PATH] [--zlib-dir PATH]
             [--gsl-dir PATH] [--hdf5-dir PATH] [--extra-cc-flags STRING]
             [--extra-cxx-flags STRING] [--extra-fc-flags STRING]
             [--custom-cc-flags STRING] [--custom-cxx-flags STRING]
             [--custom-fc-flags STRING]
             [OBJDIR]
             ...


---
.. _`faq:setupd`:
###<a name="setupd"></a> How to set CMake and Preprocessor options through the ``setup`` script

CMake can always be invoked directly to build Psi4 [](see active cmake). But more often you have a working ``setup`` configuration and just need to convey a couple CMake or Preprocessor variables.

* ###### Build with Hint Variable to CMake

    ```
    setup -DGSL_ROOT_DIR=$CONDA/envs/boostenv
    ```

* ###### Relevant ``setup`` Options:

    ```
    -D STRING             forward directly to cmake (example: -D ENABLE_THIS=1
                          -D ENABLE_THAT=1); you can also forward CPP definitions
                          all the way to the program (example: -D CPP="-DDEBUG"); 
                          also handle multi-word arguments
                          (example: -D MORELIBS="-L/path/to/lib /path/to/lib2")
                          (default: [])
    ```

* ###### Relevant ``cmake`` Options:

    ```
    -DSTRING              -express to cmake
    ```


.. _`bfaq:setupprefix`:

How to install elsewhere than :samp:`/usr/local/psi4`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The installation directory is the filesystem location for the executable
script, the Python module, basis set data, and other administrative files.
Unless using the conda package, which is relocatable, the installation
directory must be specified with CMake variable ``CMAKE_INSTALL_PREFIX``
before compiling.

* Build with Specific Install Directory

  .. code-block:: bash

   cmake -DCMAKE_INSTALL_PREFIX=/nfs/common/software/psi4

* Relevant CMake Options:

  .. code-block:: bash

   CMAKE_INSTALL_PREFIX=PATH  # Location to which Psi4 and internally built
                              # add-ons are installed (default: /usr/local/psi4)


---
.. ###<a name="profiling"></a> 
.. 
.. How to set up a profiling build
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. 
.. Specifying build type ``setup --type profile`` prepares a release build
.. type with the addition of extra flags for linking against the profiler
.. `gprof`.



.. _`faq:setuptype`:

How to compile for debugging
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Flags to turn optimizations off and debugging on can be set across the
project and plugins with CMake variable ``CMAKE_BUILD_TYPE`` before
compiling. Note that these flags *will not* propagate to any add-ons that
are detected pre-built rather than built.

* Build without optimization

  .. code-block:: bash

    cmake -DCMAKE_BUILD_TYPE=debug

                                  set the CMake build type [default: release]

* Relevant CMake Options:

  .. code-block:: bash

    CMAKE_BUILD_TYPE=[debug|release]  # Build type (Release or Debug)" [default: release]



.. _`faq:setupobjdir`:

How to compile elsewhere than ``$top-level-psi4-dir/objdir``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

[How to choose the compilation directory, ``$objdir``](2_Planning#chooseobjdir)

* Build in Specific Directory

  .. code-block:: bash

   cd $top-level-psi4-dir
   cmake -H. -Bobj-gcc
   cd obj-gcc


.. _`faq:erroreriam`:

How to fix error ``RuntimeError: value for ERI``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You will need to rebuild the code. Start with a fresh $objdir and adjust
``N`` in ``setup --max-am-eri`` according to
[here](2_Planning#how-to-configure-code-to-use-high-angular-momentum-basis-sets)


.. _`faq:chooseobjdir`:

How to choose the compilation directory, ``$objdir``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* there is no default
* common choices are ``objdir`` or ``build`` under :samp:`{$top-level-psi4-dir}`
  * ``cd $top-level-psi4-dir && cmake -H. -Bobjdir``
  * ``cd $top-level-psi4-dir && cmake -H. -Bbuild``
* in-source builds (``*.cc`` and ``*.o`` in same directory) are disallowed
* builds *outside* :samp:`{$top-level-psi4-dir}` are permitted


---
###<a name="doconfigure"></a>How to save configuration settings for a future compilation
.. _`faq:doconfigure`:

Create a file like ``do-configure`` with the ``cmake`` command and options
*on one line*. ::

 >>> cd $top-level-psi4-dir
 >>> cat do-configure
     cmake -H. -Bobjdir5 \
         -DCMAKE_INSTALL_PATH="/Users/me/psi4" \
         -DCMAKE_PREFIX_PATH="/Users/me/externals/install-libint" \
         -DMAX_AM_ERI=6 \
         -DENABLE_gdma=ON \
         -DBUILD_SHARED_LIBS=ON
 >>> chmod u+x do-configure
 >>> ./do-configure

---

.. _`faq:dirlayoutinstall`:

What is the directory layout of the staged or installed |PSIfour|
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After compilation (:samp:`cd {objdir} && make`), a directory structure like the
below will exist at :samp:`{objdir}/stage/{prefix}`. This may be tested and used
just like a full installation, see #runningfromcompdir.

After installation (:samp:`cd {objdir} && make && make install`), a directory
structure like the below will exist at :samp:`/{prefix}`. This is a full
installation, see #runningfrominstall to run.

* :envvar:`PATH` pointing to ``bin``
* :envvar:`PYTHONPATH` pointing to ``lib``
* :envvar:`PSIDATADIR` pointing to ``share/psi4``

.. code-block:: bash

 /
 bin/
 bin/psi4                                    (py script masquerading as c++ exe)
 external/                                         (external projects installed)
 external/bin
 external/include
 external/lib
 external/share
 include/psi4/                                   (header files for #include-ing)
 include/psi4/psi4-dec.h                                   (primary psi4 header)
 include/psi4/masses.h                                (project-wide psi4 header)
 include/psi4/libmints/                                   (psi4 library headers)
 include/psi4/libfock/                                                       (")
 share/                                             (read-only arch-indep files)
 share/cmake/psi4/                       (files for detecting installed targets)
 share/cmake/psi4/psi4Config.cmake                     (psi4 build/install info)
 share/cmake/psi4/psi4ConfigVersion.cmake                    (psi4 version info)
 share/doc/psi4/html/                                (sphinx html documentation)
 share/psi4/                                         (text files needed by psi4)
 share/psi4/basis                                                   (basis sets)
 share/psi4/plugins                                      (plugin template files)
 share/psi4/fsapt                                                (fsapt scripts)
 # ordinary
 lib/psi4/                                                        (object files)
 lib/psi4/driver/                                          (py-side, uncompiled)
 lib/psi4/header.py                                         (prints file header)
 lib/psi4/__init__.py                       (module marker/loader for psi4.core)
 lib/psi4/core.so                       (c-side, compiled and bound by pybind11)
 # conda
 lib/pythonX.X/site-packages/psi4/


.. ^^^^^^^^^^^^^^^^^^^^

###<a name="runfromobjdir"></a> 
###<a name="runfromprefix"></a> 

How to run |PSIfour| as executable after compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Substituting the installation directory :samp:`{prefix}` (or
:samp:`{objdir}/stage/{prefix}` to run from a compilation directory
(staged installation)) and a suitable scratch directory, issue the
following commands directly in your terminal or place them into your "rc"
file and open a new terminal.

.. code-block:: tcsh

    # csh, tcsh: add to shell or ~/.tcshrc file
    setenv PATH {prefix}/bin:$PATH
    setenv PSI_SCRATCH /path/to/existing/writable/local-not-network/directory/for/scratch/files

.. code-block:: bash

    # sh, bash: add to shell or ~/.bashrc (Linux) or ~/.bash_profile (Mac) file
    export PATH={prefix}/bin:$PATH
    export PSI_SCRATCH=/path/to/existing/writable/local-not-network/directory/for/scratch/files

Run |PSIfour|. ::

    >>> cat sample.in
    molecule {
    He
    }
    energy('hf/cc-pvdz')
    >>> psi4 sample.in

todo how to check if current py is compatible with compilation

How to run |PSIfour| as Python module after compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Substituting the installation directory :samp:`{prefix}` (or
:samp:`{objdir}/stage/{prefix}` to run from a compilation directory
(staged installation)) and a suitable scratch directory, issue the
following commands directly in your terminal or place them into your "rc"
file and open a new terminal.

.. code-block:: tcsh

    # csh, tcsh: add to shell or ~/.tcshrc file
    setenv PYTHONPATH {prefix}/lib:$PYTHONPATH
    setenv PSI_SCRATCH /path/to/existing/writable/local-not-network/directory/for/scratch/files

.. code-block:: bash

    # sh, bash: add to shell or ~/.bashrc (Linux) or ~/.bash_profile (Mac) file
    export PYTHONPATH={prefix}/bin:$PYTHONPATH
    export PSI_SCRATCH=/path/to/existing/writable/local-not-network/directory/for/scratch/files

Run |PSIfour|. ::

    >>> cat sample.py
    mol = psi4.geometry("""
    He
    """)
    psi4.energy('hf/cc-pvdz')
    >>> python sample.py

---

How to run |PSIfour| as executable from conda installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Substituting the installation directory :samp:`{prefix}` (or
:samp:`{objdir}/stage/{prefix}` to run from a compilation directory
(staged installation)) and a suitable scratch directory, issue the
following commands directly in your terminal or place them into your "rc"
file and open a new terminal.

.. code-block:: tcsh

    # csh, tcsh: add to shell or ~/.tcshrc file
    unsetenv PSIDATADIR
    setenv PATH {prefix}/bin:$PATH
    setenv PSI_SCRATCH /path/to/existing/writable/local-not-network/directory/for/scratch/files

.. code-block:: bash

    # sh, bash: add to shell or ~/.bashrc (Linux) or ~/.bash_profile (Mac) file
    unset PSIDATADIR
    export PATH={prefix}/bin:$PATH
    export PSI_SCRATCH=/path/to/existing/writable/local-not-network/directory/for/scratch/files

* |PSIfour| installed into the main Anaconda or Miniconda or Psi4conda

If you installed the Psi4conda distribution or installed the |PSIfour|
conda package into the main environment of an Anaconda or Miniconda
distribution and added that to your :envvar:`PATH`, as prompted, then
``which psi4`` likely yields :samp:`$HOME/{ana|mini|psi4}conda/bin/psi4`
and the ``PATH`` setting lines above are redundant. (prefix = home/ana|

If you installed into a conda environment :samp:`{p4env}` of
:samp:`{condadist}` and performed :samp:`source activate {p4env}`, then
the ``which psi4`` likely yields :samp:`{condadist}/envs/{p4env}/bin/psi4`
and the ``PATH`` setting lines above are redundant.

Run |PSIfour|. ::

    >>> cat sample.in
    molecule {
    He
    }
    energy('hf/cc-pvdz')
    >>> psi4 sample.in

How to run |PSIfour| as Python module from conda installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^




.. _`faq:psidatadir`:

How to set :envvar:`PSIDATADIR` and why
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``$PSIDATADIR`` is an environment variable containing the location of the
text resource 
parts of the |PSIfour| codebase (*e.g.*, basis sets,
databases, EFP fragments). When Psi4 is _installed_, the location of these
components is known relative to the executable, so the location is set
internally. When Psi4 is run from the _compilation directory_, the
relative location is not known, so the value must be explicitly set:

  * in the shell
    * ``csh``/``tcsh``: ``setenv PSIDATADIR $top-level-psi4-dir/share``
    * ``sh``/``bash``: ``export PSIDATADIR=$top-level-psi4-dir/share``
  * or at runtime: ``psi4 -p $top-level-psi4-dir/share``







