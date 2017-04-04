.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2017 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This program is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU General Public License as published by
.. # the Free Software Foundation; either version 2 of the License, or
.. # (at your option) any later version.
.. #
.. # This program is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU General Public License for more details.
.. #
.. # You should have received a copy of the GNU General Public License along
.. # with this program; if not, write to the Free Software Foundation, Inc.,
.. # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.. #
.. # @END LICENSE
.. #


.. index:: prerequisites, compiling, installing
.. _`sec:installFile`:

====================================
Compiling and Installing from Source
====================================

.. _`faq:cmakeviasetup`:

Planning: how to configure Psi4 and invoke CMake
------------------------------------------------

|PSIfour| is built through CMake. Prior to 1.1, |PSIfour| had a Python
script ``setup`` as a frontend to CMake, but this is no more, and
``cmake`` is now invoked directly. An abbreviated build guide can be found
`within the source itself
<https://github.com/psi4/psi4/blob/master/CMakeLists.txt#L13-L122>`_.

CMake does a good job scanning your computer to locate libraries, header
files, and executables needed for compilation. So it's very possible that
from :samp:`{top-level-psi4-dir}` source directory, you can call :samp:`cmake -H.
-B{objdir}` without any further arguments, and it will invoke ``cmake``,
detect some appropriate defaults, configure the build, create a build
directory :samp:`{objdir}`, and complete, leaving you to only invoke
``make`` in the new build directory.

Should this happy scenario not come to pass, or if the default build
options are not to your taste, use the links within :ref:`core
dependencies <faq:coredepend>` and :ref:`add-on dependencies
<faq:addondepend>` to plan a set of arguments to ``cmake`` tailored to
your computer. Proceed to :ref:`quick build <faq:buildquick>` or
:ref:`detailed build <faq:builddetailed>`.

.. The following topics may also be helpful.

.. * :ref:`faq:setuphelp`
.. * :ref:`faq:chooseobjdir`
.. * :ref:`faq:setupprefix`
.. * :ref:`faq:setuptype`
.. * :ref:`faq:setupmaxameri`


.. _`faq:buildquick`:

How to build and install Psi4, the compact version
--------------------------------------------------

This section outlines the main steps of configuring, compiling, and
installing |PSIfour|. More detail is given :ref:`here
<faq:builddetailed>`. ::

    >>> cd {top-level-psi4-dir}
    >>> cmake -H. -Bobjdir [your configuration options]
    >>> cd objdir
    >>> make -j`getconf _NPROCESSORS_ONLN`
    >>> make install


.. _`faq:builddetailed`:

How to build, test, and install Psi4, in detail
-----------------------------------------------

**1. Plan Directories**

   Get ahold of the |PSIfour| codebase, and navigate to the top level source
   directory, hereafter :samp:`{top-level-psi4-dir}`.

   * :ref:`faq:obtainpsi4`

   ::

    >>> cd {top-level-psi4-dir}

   Choose a compilation directory, hereafter :samp:`{objdir}`

   * :ref:`faq:chooseobjdir`

   Choose an installation directory, hereafter :samp:`{prefix}`

   * :ref:`faq:setupprefix`

**2. Plan Configuration**

   Examine the strict and optional software requirements to make sure the
   target computer has all the necessary dependencies installed.

   * :ref:`faq:coredepend`
   * :ref:`faq:addondepend`

   Prepare any necessary or desired configuration options for ``cmake``,
   hereafter ``[your configuration options]``

   * :ref:`faq:setuphelp`
   * :ref:`faq:cmakeviasetup`

**3. Configure**

   Run CMake with planned options and directories, as below. It reports on
   software found or unfound as it scans the computer, then (upon success)
   creates :samp:`{objdir}` ready for compilation.

   ::

    >>> cmake -H. -B{objdir} -DCMAKE_INSTALL_PREFIX={prefix} [your configuration options]

**4. Compile**

   Compile the code (optional ``-j`` triggers parallel compilation).

   ::

    >>> cd {objdir}
    >>> make -j`getconf _NPROCESSORS_ONLN`

**5. Test**

   Optionally, use CTest (thorough) or pytest (cursory) to test the build.

   * :ref:`faq:minutetests`
   * :ref:`faq:subsettests`
   * :ref:`faq:testsoutput`
   * :ref:`faq:pytest`

   ::

   >>> ctest -j`getconf _NPROCESSORS_ONLN`

   >>> make pytest

**6. Install**

   If tests pass, install the code.

   ::

   >>> make install

**7. Configure Runtime**

   To run Psi4 after installation, you need to configure a few variables:

   * :ref:`faq:runordinaryexe`
   * :ref:`faq:runordinarymodule`


.. _`faq:coredepend`:

What are the tools and dependencies strictly required for building Psi4
-----------------------------------------------------------------------

The core |PSIfour| build requires the software below. Note that Python,
CMake, NumPy, and libint (and even C++ compilers on Linux) can be
satisfied through conda. The links below give examples of how to configure
that software for |PSIfour| and any notes and warnings pertaining to it.

* :ref:`C++ and C Compilers <cmake:cxx>` (C++11 compliant)

* :ref:`Optimized BLAS and LAPACK libraries <cmake:lapack>` (preferably NOT one supplied by a standard
  Linux distribution)

* :ref:`Python interpreter and headers <cmake:python>` (2.7, 3.5, or 3.6) https://www.python.org/

* CMake (3.3+) http://www.cmake.org/download/

* NumPy (needed at runtime, not buildtime) http://www.numpy.org/

* System utilities: GNU make, GNU install, POSIX threads (Pthreads) library

The following are also required for |PSIfour|, but if not detected, the
build system will automatically download and build.

* :ref:`libint <cmake:libint>` |w---w| :ref:`[what is this?] <sec:libint>` `[min version] <https://github.com/psi4/psi4/blob/master/external/upstream/libint/CMakeLists.txt#L1>`_

* pybind11 |w---w| `[what is this?] <https://pybind11.readthedocs.io/en/master/>`_ `[min version] <https://github.com/psi4/psi4/blob/master/external/upstream/pybind11/CMakeLists.txt#L1>`_


.. _`faq:addondepend`:

What are the add-on capabilities for Psi4 and what are their dependencies
-------------------------------------------------------------------------

Each of the items below is an independent additional capability that can
be built with |PSIfour|. Sub-items below are the respective additional
dependencies of the add-on. Select which, if any, you want, and examine
the links for appropriate enabling arguments to ``cmake``. Note that many
are available pre-built from conda.

* |PSIfour| Testing

  * CTest http://www.cmake.org/download/
  * Perl (for some coupled-cluster CTest tests) http://perl.org
  * pytest (for installed testing) http://doc.pytest.org/en/latest/

* |PSIfour| Documentation (available pre-built at http://www.psicode.org/psi4manual/master/index.html)

  * Sphinx (1.4+) http://sphinx-doc.org
  * Perl (for some auto-documentation scripts) http://perl.org
  * nbsphinx (for converting Jupyter notebooks) http://nbsphinx.readthedocs.io/en/jupyter-theme/

* :ref:`CheMPS2 <cmake:chemps2>` |w---w| :ref:`[what is this?] <sec:chemps2>` `[min version] <https://github.com/psi4/psi4/blob/master/external/upstream/chemps2/CMakeLists.txt#L2>`_

  * HDF5 https://support.hdfgroup.org/HDF5/
  * zlib http://www.zlib.net/

* :ref:`libefp <cmake:libefp>` |w---w| :ref:`[what is this?] <sec:libefp>` `[min version] <https://github.com/psi4/psi4/blob/master/external/upstream/libefp/CMakeLists.txt#L1>`_

* :ref:`erd <cmake:erd>` |w---w| :ref:`[what is this?] <sec:erd>` `[min version] <https://github.com/psi4/psi4/blob/master/external/upstream/erd/CMakeLists.txt#L2>`_

  * :ref:`Fortran Compiler <cmake:fortran>`

* :ref:`dkh <cmake:dkh>` |w---w| :ref:`[what is this?] <sec:dkh>` `[min version] <https://github.com/psi4/psi4/blob/master/external/upstream/dkh/CMakeLists.txt#L2>`_

  * :ref:`Fortran Compiler <cmake:fortran>`

* :ref:`gdma <cmake:gdma>` |w---w| :ref:`[what is this?] <sec:gdma>` `[min version] <https://github.com/psi4/psi4/blob/master/external/upstream/gdma/CMakeLists.txt#L2>`_

  * :ref:`Fortran Compiler <cmake:fortran>`

* :ref:`PCMSolver <cmake:pcmsolver>` |w---w| :ref:`[what is this?] <sec:pcmsolver>`

  * :ref:`Fortran Compiler <cmake:fortran>`
  * zlib http://www.zlib.net/

* :ref:`simint <cmake:simint>` |w---w| :ref:`[what is this?] <sec:simint>` `[min version] <https://github.com/psi4/psi4/blob/master/external/upstream/simint/CMakeLists.txt#L2>`_

Additionally, there are runtime-only capabilities:

* cfour |w---w| :ref:`[what is this?] <sec:cfour>`

* dftd3 |w---w| :ref:`[what is this?] <sec:dftd3>`

* gcp |w---w| :ref:`[what is this?] <sec:gcp>`

* mrcc |w---w| :ref:`[what is this?] <sec:mrcc>`


.. _`faq:setupmaxameri`:

How to configure code to use high angular momentum basis sets
-------------------------------------------------------------

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

Note that since |PSIfour| 1.1, it is possible to build libint
independently (or install just the libint conda package), then have
any/all |PSIfour| builds detect that installation at compile-time.

* :ref:`cmake:libint`


.. _`faq:setuphelp`:

How to see what build configuration options are available
---------------------------------------------------------

CMake doesn't provide a summary for this (unless you want to try the CMake
GUI, which the developers have never looked at). However, the top half of
the main CMakeLists.txt is a passable summary:

.. literalinclude:: @SFNX_INCLUDE@CMakeLists.txt
   :lines: 13-122

Note that external projects will have their own sets of build
configuration options. Only the most-common user knobs of those are
mentioned above.

.. .. _`faq:setupd`:
.. ###<a name="setupd"></a> How to set CMake and Preprocessor options through the ``setup`` script
.. 
.. CMake can always be invoked directly to build Psi4 [](see active cmake). But more often you have a working ``setup`` configuration and just need to convey a couple CMake or Preprocessor variables.
.. 
.. * ###### Build with Hint Variable to CMake
.. 
..     ```
..     setup -DGSL_ROOT_DIR=$CONDA/envs/boostenv
..     ```
.. 
.. * ###### Relevant ``setup`` Options:
.. 
..     ```
..     -D STRING             forward directly to cmake (example: -D ENABLE_THIS=1
..                           -D ENABLE_THAT=1); you can also forward CPP definitions
..                           all the way to the program (example: -D CPP="-DDEBUG");
..                           also handle multi-word arguments
..                           (example: -D MORELIBS="-L/path/to/lib /path/to/lib2")
..                           (default: [])
..     ```
.. 
.. * ###### Relevant ``cmake`` Options:
.. 
..     ```
..     -DSTRING              -express to cmake
..     ```


.. _`faq:setupprefix`:

How to install elsewhere than :samp:`/usr/local/psi4`
-----------------------------------------------------

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
----------------------------

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



.. .. _`faq:setupobjdir`:
.. 
.. How to compile elsewhere than ``{top-level-psi4-dir}/objdir``
.. -------------------------------------------------------------
.. 
.. [How to choose the compilation directory, ``$objdir``](2_Planning#chooseobjdir)
.. 
.. * Build in Specific Directory
.. 
..   .. code-block:: bash
.. 
..    cd $top-level-psi4-dir
..    cmake -H. -Bobj-gcc
..    cd obj-gcc


.. _`faq:erroreriam`:

How to fix error "``RuntimeError: value for ERI``"
--------------------------------------------------

You will need to rebuild libint. Reissue ``cmake`` or edit
``CMakeCache.txt`` with larger ``MAX_AM_ERI`` and rebuilt.

* :ref:`faq:setupmaxameri`

.. _`faq:chooseobjdir`:

How to choose the compilation directory, ``{objdir}``
-----------------------------------------------------

* there is no default
* common choices are ``objdir`` or ``build`` under :samp:`{top-level-psi4-dir}`

  * ``cd {top-level-psi4-dir} && cmake -H. -Bobjdir``
  * ``cd {top-level-psi4-dir} && cmake -H. -Bbuild``

* in-source builds (``*.cc`` and ``*.o`` in same directory) are disallowed
* builds *outside* :samp:`{top-level-psi4-dir}` are permitted


.. _`faq:doconfigure`:

How to save configuration settings for a future compilation
-----------------------------------------------------------

Create a file like ``do-configure`` with the ``cmake`` command and options
*on one line*. ::

 >>> cd {top-level-psi4-dir}
 >>> cat do-configure
     cmake -H. -B{objdir} \
         -DCMAKE_INSTALL_PATH="/Users/me/psi4" \
         -DCMAKE_PREFIX_PATH="/Users/me/externals/install-libint" \
         -DMAX_AM_ERI=6 \
         -DENABLE_gdma=ON \
         -DBUILD_SHARED_LIBS=ON
 >>> chmod u+x do-configure
 >>> ./do-configure


.. _`faq:dirlayoutinstall`:

What is the directory layout of the installed or staged Psi4
------------------------------------------------------------

After compilation (:samp:`cd {objdir} && make`), a directory structure like the
below will exist at :samp:`{objdir}/stage/{prefix}`. This may be tested and used
just like a full installation.

After installation (:samp:`cd {objdir} && make && make install`), a directory
structure like the below will exist at :samp:`/{prefix}`. This is a full
installation.

.. code-block:: bash

 /
 bin/                                                   (executables for psi4 + any external proj)
 bin/psi4                                             (psi4 executable, actually just a py script)
 include/                                         (installed headers for psi4 + any external proj)
 include/psi4/                                                     (header files for #include-ing)
 include/psi4/psi4-dec.h                                                     (primary psi4 header)
 include/psi4/masses.h                                                (a project-wide psi4 header)
 include/psi4/libmints/                                                     (psi4 library headers)
 include/psi4/libfock/                                                                     (ditto)
 share/                                  (read-only arch-indep files for psi4 + any external proj)
 share/cmake/psi4/                                         (files for detecting installed targets)
 share/cmake/psi4/psi4Config.cmake                                       (psi4 build/install info)
 share/cmake/psi4/psi4ConfigVersion.cmake                                      (psi4 version info)
 share/doc/psi4/html/                                                  (sphinx html documentation)
 share/psi4/                                                           (text files needed by psi4)
 share/psi4/basis                                                                     (basis sets)
 share/psi4/plugins                                                        (plugin template files)
 share/psi4/fsapt                                                                  (fsapt scripts)
 share/psi4/samples/                                                          (sample input files)
 lib/                               (shared libraries and py modules for psi4 + any external proj)
 # ordinary
 lib/psi4/                                                                          (object files)
 lib/psi4/driver/                                                            (py-side, uncompiled)
 lib/psi4/header.py                                                           (prints file header)
 lib/psi4/__init__.py                                         (module marker/loader for psi4.core)
 lib/psi4/core.so                                         (c-side, compiled and bound by pybind11)
 # conda
 lib/pythonX.X/site-packages/psi4/

The following environment variables point to certain places in the above
directory structure. None to few need to be set; see for details:
:ref:`running compiled executable <faq:runordinaryexe>`,
:ref:`running compiled Python module <faq:runordinarymodule>`,
:ref:`running conda binary <faq:runfrombinary>`.

* :envvar:`PATH` pointing to ``bin``
* :envvar:`PYTHONPATH` pointing to ``lib`` (ordinary) or ``lib/pythonX.X/site-packages`` (conda)
* :envvar:`PSIDATADIR` pointing to ``share/psi4``


.. _`faq:runordinaryexe`:

How to run Psi4 as executable after compilation
-----------------------------------------------

Substituting the full installation directory :samp:`{prefix}` and a
suitable scratch directory, issue the following commands directly in your
terminal or place them into your "rc" file and open a new terminal. (To
use a staged installation directory, substitute
:samp:`{objdir}/stage/{prefix}` for :samp:`{prefix}`.)

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
    compare_values(-2.85518839, get_variable('current energy'), 5, 'SCF E')
    >>> psi4 sample.in
    SCF E.............................................................PASSED


todo how to check if current py is compatible with compilation


.. _`faq:runordinarymodule`:

How to run Psi4 as Python module after compilation
--------------------------------------------------

Substituting the full installation directory :samp:`{prefix}` and a
suitable scratch directory, issue the following commands directly in your
terminal or place them into your "rc" file and open a new terminal. (To
use a staged installation directory, substitute
:samp:`{objdir}/stage/{prefix}` for :samp:`{prefix}`.)

.. code-block:: tcsh

    # csh, tcsh: add to shell or ~/.tcshrc file
    setenv PYTHONPATH {prefix}/lib:$PYTHONPATH
    setenv PSI_SCRATCH /path/to/existing/writable/local-not-network/directory/for/scratch/files

.. code-block:: bash

    # sh, bash: add to shell or ~/.bashrc (Linux) or ~/.bash_profile (Mac) file
    export PYTHONPATH={prefix}/lib:$PYTHONPATH
    export PSI_SCRATCH=/path/to/existing/writable/local-not-network/directory/for/scratch/files

.. note:: If you know the location of the |PSIfour| executable
   (``bin/psi4``) and want to know the corresponding location to add to
   :envvar:`PYTHONPATH`, execute ``psi4 --psiapi-path``. It will return
   bash commands to set :envvar:`PATH` (for correct python interpreter)
   and :envvar:`PYTHONPATH` (to find psi4 module) correctly, after which
   ``python -c "import psi4"`` will work.

Run |PSIfour|. ::

    >>> cat sample.py
    import psi4
    mol = psi4.geometry("""
    He
    """)
    psi4.energy('hf/cc-pvdz')
    psi4.compare_values(-2.85518839, psi4.core.get_variable('current energy'), 5, 'SCF E')
    >>> python sample.py
    SCF E.............................................................PASSED


.. _`faq:runfrombinary`:

How to run Psi4 as executable or Python module from conda installation
----------------------------------------------------------------------

The configuration commands below are generic versions of the ones printed
to your screen as advice upon installing |PSIfour| into a Anaconda,
Miniconda, or Psi4conda distribution, :samp:`{condadist} =
{ana|mini|psi4}conda`. If ``which conda python psi4`` points to your
:samp:`{condadist}` and ``echo $PSI_SCRATCH`` is set, skip ahead to the
"Run |PSIfour|\" commands below. Otherwise, issue the following
commands directly in your terminal or place them into your "rc" file and
open a new terminal.

If you installed the Psi4conda distribution or installed the |PSIfour|
conda package into the main environment of an Anaconda or Miniconda
distribution and added that to your :envvar:`PATH`, as prompted, then
``which psi4`` likely yields :samp:`{condadist}/bin/psi4` and the ``PATH``
setting lines below are redundant.

If you installed into a conda environment :samp:`{p4env}` and performed
:samp:`source activate {p4env}`, then ``which psi4`` likely yields
:samp:`{condadist}/envs/{p4env}/bin/psi4` and the ``PATH`` setting lines
below are redundant.

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

.. If you installed the Psi4conda distribution or installed the |PSIfour|
.. conda package into the main environment of an Anaconda or Miniconda
.. distribution and added that to your :envvar:`PATH`, as prompted, then
.. :samp:`{condadist} = $HOME/{ana|mini|psi4}conda` and ``which psi4`` likely
.. yields :samp:`{condadist}/bin/psi4` and the ``PATH`` setting lines above
.. are redundant.

Run |PSIfour| as executable. ::

    >>> cat sample.in
    molecule {
    He
    }
    energy('hf/cc-pvdz')
    compare_values(-2.85518839, get_variable('current energy'), 5, 'SCF E')
    >>> psi4 sample.in
    SCF E.............................................................PASSED

*or* Run |PSIfour| as Python module. ::

    >>> cat sample.py
    import psi4
    mol = psi4.geometry("""
    He
    """)
    psi4.energy('hf/cc-pvdz')
    psi4.compare_values(-2.85518839, psi4.core.get_variable('current energy'), 5, 'SCF E')
    >>> python sample.py
    SCF E.............................................................PASSED


.. _`faq:inplace`:

How to run Psi4 as executable after compilation using driver from source
------------------------------------------------------------------------

When developing python driver code, it can be annoying to keep `make`\
ing to test the code. |PSIfour| can be run "inplace" through the
following procedure. To be clear, this is running compiled C++ from the
build directory and python from the source directory. This is an expert
option for development, and not all functionality will be available. ::

    >>> cd {objdir}
    >>> ln -s {top-level-psi4-dir}/{objdir}/stage/{prefix}/lib/psi4/core.so ../psi4/core.so
    >>> python ../psi4/run_psi4.py --inplace input.dat


.. _`faq:psidatadir`:

How to set :envvar:`PSIDATADIR` and why
---------------------------------------

:envvar:`PSIDATADIR` is an environment variable containing the location of the
text resource parts of the |PSIfour| codebase (*e.g.*, basis sets,
databases, EFP fragments). It is for developer use only. In |PSIfour| 1.1
and beyond, the program *always* knows where its resources are, and the
only reason to set this variable is to point to another location.
Previously in |PSIfour| 1.0 and previous, only installed executables knew
the location, so it always needed to be explicitly set when run from the
compilation directory.

At runtime

.. code-block:: bash

   >>> psi4 -p {top-level-psi4-dir}/psi4/share/psi4

Or in the shell

.. code-block:: tcsh

    # csh, tcsh: add to shell or ~/.tcshrc file
    setenv PSIDATADIR {top-level-psi4-dir}/psi4/share/psi4

.. code-block:: bash

    # sh, bash: add to shell or ~/.bashrc (Linux) or ~/.bash_profile (Mac) file
    export PSIDATADIR={top-level-psi4-dir}/psi4/share/psi4




.. _`cmake:cxx`:

How to configure C++ and C compilers for building Psi4
------------------------------------------------------

**Role and Dependencies**

* Role |w---w| In |PSIfour|, a C++ compiler is vital for building the code.

* Downstream Dependencies |w---w| |PSIfour| |dr| C++ Compiler

**CMake Variables**

* :makevar:`CMAKE_CXX_COMPILER` |w---w| CMake variable to specify name or full path to C++ compiler.
* :makevar:`CMAKE_C_COMPILER` |w---w| CMake variable to specify name or full path to C compiler.
* :makevar:`CMAKE_CXX_FLAGS` |w---w| CMake variable to specify any additional custom compiler flags for C++ source.
* :makevar:`CMAKE_C_FLAGS` |w---w| CMake variable to specify any additional custom compiler flags for C source.

**Examples**

A. Build with detected compilers from :envvar:`PATH`

  .. code-block:: bash

    >>> cmake

B. Build with specific (Intel) compilers from :envvar:`PATH`

  .. code-block:: bash

    >>> cmake -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc

C. Build with specific (GNU) compilers *not* in :envvar:`PATH`

  .. code-block:: bash

    >>> cmake -DCMAKE_CXX_COMPILER=/path/to/gcc6.2/bin/g++ -DCMAKE_C_COMPILER=/path/to/gcc6.2/bin/gcc

D. Build with specific (Intel) compilers from :envvar:`PATH` based on GCC *not* in :envvar:`PATH`

  .. code-block:: bash

    >>> cmake -DCMAKE_C_COMPILER=icc \
              -DCMAKE_CXX_COMPILER=icpc \
              -DCMAKE_C_FLAGS="-gcc-name=${GCC5}/bin/gcc" \
              -DCMAKE_CXX_FLAGS="-gcc-name=${GCC5}/bin/gcc -gxx-name=${GCC5}/bin/g++"

E. Build with specific (Intel) compilers from :envvar:`PATH` based on GCC *not* in :envvar:`PATH` and also building Fortran Add-Ons

  .. code-block:: bash

    >>> cmake -DCMAKE_C_COMPILER=icc \
              -DCMAKE_CXX_COMPILER=icpc \
              -DCMAKE_Fortran_COMPILER=ifort \
              -DCMAKE_C_FLAGS="-gcc-name=${GCC5}/bin/gcc" \
              -DCMAKE_CXX_FLAGS="-gcc-name=${GCC5}/bin/gcc -gxx-name=${GCC5}/bin/g++" \
              -DCMAKE_Fortran_FLAGS="-gcc-name=${GCC5}/bin/gcc -gxx-name=${GCC5}/bin/g++"


.. _`faq:approvedcxx`:

What C and C++ compilers and versions are approved
--------------------------------------------------

On Linux, the following work nicely.

  * GNU: ``gcc``, ``g++``
  * Intel: ``icc``, ``icpc``
  * Clang: ``clang``, ``clang++``

On Mac, the following work nicely.

  * Apple Clang: ``clang``, ``clang++``
  * Intel: ``icc``, ``icpc``

|PSIfour| requires *full* C++11 compliance, meaning, most importantly, GCC
>= 4.9. This compliance is checked for at build-time with file
:source:`cmake/custom_cxxstandard.cmake`, so either consult that file or
try a test build to ensure your compiler is approved. Note that Intel
compilers also rely on GCC, so both ``icpc`` and ``gcc`` versions are checked.

* :ref:`faq:modgcc`


.. _`faq:macxcode`:

How to obtain C and C++ compilers for Mac without Fink, MacPorts, or Homebrew
-----------------------------------------------------------------------------

The easiest compiler to obtain is ``clang`` which is a drop-in replacement
for ``gcc`` and ``g++``. Just install `XCode
<https://itunes.apple.com/us/app/xcode/id497799835>`_.  Some old versions
of XCode can't handle some of the advanced C++ language features, but this
is a *software* not *hardware* limitation. Checks for version compliance
performed at build-time.


.. _`faq:modgcc`:

How to satisfy the GCC >= 4.9 requirement on Linux without updating the OS
--------------------------------------------------------------------------

.. code-block:: bash

   # See if GCC too old (in this case, yes)
   >>> gcc --version
   gcc (GCC) 4.8.5 20150623 (Red Hat 4.8.5-4)

Even if you're not using GCC as a compiler, your Intel compiler seeks
`gcc` to generate code compatible with your current GCC version. If your
GCC is too old (like above), you can update your system GCC through your
package manager *or* install an acceptable version elsewhere on your
system. The latter route, tested on Linux with Intel compilers, is below.

.. code-block:: bash

   # Install GCC 5.2 into a non-primary conda environment
   >>> conda create -n gcc52 -c psi4 gcc-5-mp

   # To Build, either:

   #    (A) activate environment (prepends PATH)
   >>> source activate gcc52

   #    (B) tell CMake to tell the compiler which GCC to use
   >>> vi do-configure
      GCC5=/path/to/miniconda/envs/gcc52
      cmake -H. -Bbuild \
        -DCMAKE_C_COMPILER=icc \
        -DCMAKE_CXX_COMPILER=icpc \
        -DCMAKE_C_FLAGS="-gcc-name=${GCC5}/bin/gcc" \
        -DCMAKE_CXX_FLAGS="-gcc-name=${GCC5}/bin/gcc -gxx-name=${GCC5}/bin/g++" \
        ...
        # if Fortran active ...
        -DCMAKE_Fortran_COMPILER=icpc \
        -DCMAKE_Fortran_FLAGS="-gcc-name=${GCC5}/bin/gcc -gxx-name=${GCC5}/bin/g++" \

   # Configure and build

   # To Run:
   >>> export LD_LIBRARY_PATH=${GCC5}/lib:$LD_LIBRARY_PATH


.. _`faq:cray`:

How to configure a Psi4 build on Cray
-------------------------------------

Cray systems strongly prefer to build static libraries, but |PSIfour|
needs to be dynamic to function as a Python module. Courtesy of @misha
at the forum and various supercomputer guides, building |PSIfour| on
Cray requires setting environment variables :envvar:`CRAYPE_LINK_TYPE`
and :envvar:`CRAY_ADD_RPATH` before running `cmake`. ::

    CRAYPE_LINK_TYPE=dynamic CRAY_ADD_RPATH=yes cmake ...


.. _`cmake:fortran`:

How to configure Fortran compilers for building Psi4
----------------------------------------------------

**Role and Dependencies**

* Role |w---w| In |PSIfour|, a Fortran compiler in unneeded for core
  features but may be required for add-ons.

* Downstream Dependencies

  * |PSIfour| (\ |dr| optional) Fortran Compiler
  * erd, dkh, gdma, PCMSolver |dr| Fortran Compiler

**CMake Variables**

* :makevar:`CMAKE_Fortran_COMPILER` |w---w| CMake variable to specify name or full path to Fortran compiler.
* :makevar:`CMAKE_Fortran_FLAGS` |w---w| CMake variable to specify any additional custom compiler flags for Fortran source.

**Examples**

A. Build with detected compiler from :envvar:`PATH`

  .. code-block:: bash

    >>> cmake

B. Build with specific (Intel) compiler from :envvar:`PATH`

  .. code-block:: bash

    >>> cmake -DCMAKE_Fortran_COMPILER=ifort


.. _`faq:approvedfc`:

What Fortran compilers are approved
-----------------------------------

On Linux and Mac, the following work nicely.

  * GNU: ``gfortran``
  * Intel: ``ifort``

* Packages to install for specific OS or package managers:

  * Ubuntu ``gfortran``
  * conda ``gcc`` or ``gcc-5`` to get ``gfortran``


.. _`faq:macgfortran`:

How to obtain a Fortran compiler for Mac without Fink, MacPorts, or Homebrew
----------------------------------------------------------------------------

Xcode does not provide a Fortran compiler. A way to get one is to download
the ``gcc`` conda package. This provides ``gcc``, ``g++``, and
``gfortran`` compilers for Mac. The two former are 4.8.5 and so are too
old to compile |PSIfour|, but the Fortran compiler will work.

.. Xcode does not provide a Fortran compiler. Although a Fortran compiler is
.. not required for Psi4, a broken one can prevent correct configuration. Do
.. not download the latest version of GFortran from the HPC website; this is
.. unlikely to be compatible with your version of GCC. Instead, you should
.. run ``gcc -v`` to find out what version of GCC you're using, and then
.. download the corresponding GFortran from
.. <http://r.research.att.com/tools/>.  If you configure Psi on a Mac without
.. any Fortran compiler it will set itself up correctly, so this is only
.. necessary if you want a Fortran compiler for other purposes.


.. _`cmake:lapack`:

How to configure BLAS/LAPACK for building Psi4
----------------------------------------------

**Role and Dependencies**

* Role |w---w| In |PSIfour|, BLAS and LAPACK control much of the speed
  and efficiency of the code since computational chemistry is essentially
  linear algebra on molecular systems.

* Downstream Dependencies |w---w| |PSIfour| |dr| LAPACK Libraries

**CMake Variables**

* :makevar:`BLAS_TYPE` |w---w| CMake variable to specify which BLAS libraries to look for among ``MKL|OPENBLAS|ESSL|ATLAS|ACML|SYSTEM_NATIVE``.
* :makevar:`LAPACK_TYPE` |w---w| CMake variable to specify which LAPACK libraries to look for among ``MKL|OPENBLAS|ESSL|ATLAS|ACML|SYSTEM_NATIVE``.
* :envvar:`MKL_ROOT` |w---w| Environment variable set by Intel compilervars scripts. Sufficient to trigger math detection of MKL at this location.
* :envvar:`MATH_ROOT` |w---w| Environment variable to specify root directory in which BLAS/LAPACK libraries should be detected (*e.g.*, ``${MATH_ROOT}/lib64/libblas.so`` and ``${MATH_ROOT}/lib64/liblapack.so``).
* :makevar:`LAPACK_LIBRARIES` |w---w| CMake variable to specify BLAS/LAPACK libraries explicitly, bypassing math detection. Should be semicolon-separated list of full paths.
* :makevar:`LAPACK_INCLUDE_DIRS` |w---w| CMake variable to specify BLAS/LAPACK header location explicitly, bypassing math detection. Only needed for MKL.

**Examples**

A. Build with any LAPACK in standard location

  .. code-block:: bash

    >>> cmake

B. Build with native Accelerate LAPACK on Mac (MKL *not* also present)

  .. code-block:: bash

    >>> cmake

C. Build with native Accelerate LAPACK on Mac (MKL also present)

  .. code-block:: bash

    >>> cmake -DBLAS_TYPE=SYSTEM_NATIVE -DLAPACK_TYPE=SYSTEM_NATIVE

D. Build with Intel MKL

  .. code-block:: bash

    >>> source /path/to/intel/vers/linux/mkl/bin/mklvars.sh intel64  # adjust sh/csh and arch as needed
    >>> cmake

  .. code-block:: bash

    >>> MATH_ROOT=/path/to/intel/vers/linux/mkl/ cmake

E. Build with Intel MKL from conda

  .. code-block:: bash

    # won't work, as mkl.h header also needed

F. OpenBLAS

  .. code-block:: bash

    >>> MATH_ROOT=/path/to/openblas/0.2.13_seq/x86_64/gcc_5.2.0/lib cmake

G. Build with explicit MKL LAPACK

  .. code-block:: bash

    >>> cmake -DLAPACK_LIBRARIES="/path/to/lib/intel64/libmkl_lapack95_lp64.a;/path/to/lib/intel64/libmkl_rt.so" -DLAPACK_INCLUDE_DIRS="/path/to/mkl-h-include/"

H. Build with explicit non-MKL LAPACK

  .. code-block:: bash

    >>> cmake -DLAPACK_LIBRARIES="/path/to/lib/liblapack.so;/path/to/lib/libblas.a"

**Notes**

* Much of |PSIfours| speed and efficiency depends on the corresponding
  speed and efficiency of the linked BLAS and LAPACK libraries
  (especially the former). Consider the following recommendations:

  * It is NOT wise to use the stock BLAS library provided with many
    Linux distributions like RedHat, as it is usually just the completely
    unoptimized netlib distribution. The choice of LAPACK is less
    critical, and so the unoptimized netlib distribution is acceptable.

  * Perhaps the best choice, if available, is Intel's MKL library,
    which includes efficient threaded BLAS and LAPACK (as of |PSIfour|
    v1.1, earliest known working version is MKL 2013). On Mac, the
    native Accelerate libraries are also recommended.

  * For open-source LAPACK distributions, OpenBLAS (formerly GotoBLAS)
    is known to work, while ATLAS is known
    (https://github.com/psi4/psi4/issues/391) to have stability issues
    with the DFOCC module.

  * ACML libraries are known to work with |PSIfour| v1.1 at ACML 6.

* The BLAS/LAPACK detected for |PSIfour| are also linked into any
  Add-Ons (*e.g.*, libefp) that require them, rather than relying on
  those packages' native math detection.

* The separation between BLAS and LAPACK seen in detection printing
  and CMake variables is purely formal. In practice, they get run
  together and linked as ``${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}``.

* Sometimes the CMake's library search capabilites falter at SONAMEs
  (*e.g.*, ``libblas.so.3`` *vs.* ``libblas.so``), extensions (static
  *vs.* dynamic), or suffixes (*e.g.*, ``libacml_mp.so`` *vs.*
  ``libacml.so``). The developers would be interested in hearing
  of such problems to expand the math detection capabilities. The
  immediate solution, however, is to form symlinks between the
  library names that exist and the names expected. Consult file
  :source:`cmake/math/MathLibs.cmake` for the library patterns being
  sought.

* The BLAS/LAPACK interface is standardized, so only libraries, not
  headers, need to be detected. The exception is MKL, where the ``mkl.h``
  header defines additional functionality; it must be located to use
  BLAS threading.


.. _`cmake:python`:

How to configure Python for building Psi4
-----------------------------------------

**Role and Dependencies**

* Role |w---w| In |PSIfour|, Python allows the core compiled C++ code to
  be flexibly accessed for manipulation and extension in an interpreted
  language.

* Downstream Dependencies |w---w| |PSIfour| |dr| Python Interpreter

**CMake Variables**

* :makevar:`PYTHON_EXECUTABLE` |w---w| specify name or full path to Python interpreter.
* :makevar:`PYTHON_LIBRARY` |w---w| specify path to Python library.
* :makevar:`PYTHON_INCLUDE_DIR` |w---w| specify directory of Python headers. Contains ``Python.h``.

**Examples**

A. Build with detected Python from :envvar:`PATH`

  .. code-block:: bash

    >>> cmake

B. Build with specific Python

  .. code-block:: bash

    >>> cmake -DPYTHON_EXECUTABLE=/path/to/interp/python2.7

C. Build with full Python specification to root directory ``${PFXC}``

  .. code-block:: bash

    >>> cmake -DPYTHON_EXECUTABLE="${PFXC}/bin/python" \
              -DPYTHON_LIBRARY="${PFXC}/lib/libpython3.5m.so" \
              -DPYTHON_INCLUDE_DIR="${PFXC}/include/python3.5m"


.. _`faq:runtimepython`:

What Python is Psi4 running
---------------------------

The Python detected at build-time is embedded into the |PSIfour|
executable. That is, the top line of ``bin/psi4`` is something like
``#!/path/to/miniconda/envs/p4deps/bin/python3.5``, and that's the Python
through which |PSIfour| is running, *not* the Python of ``which python``.
To use a different Python with |PSIfour| in the short term, just
``path/to/desired/python psi4`` on the command line to override the
shebang line. To use a different Python with |PSIfour| in the long term,
edit the shebang line.

If you're using |PSIfour| as a Python module, then |PSIfour| *is* running
the Python of ``which python``.


.. _`faq:wrongpyfalse`:

How to fix "``undefined symbol: _Py_FalseStruct``"
--------------------------------------------------

You're probably loading a Py3-compiled Psi4 in Py2. Switch interpreters
and re-run. A python of proper Py2 or Py3-ness is baked into the |PSIfour|
"executable", so you'll see this error only for Psi4 as Python module.


.. _`faq:gdblldb`:

How to use ``gdb`` and ``lldb`` with Psi4
-----------------------------------------

Debugging |PSIfour| has gotten a little confusing now that it's running through Python. Here's the syntax ::

  >>> cd {objdir}
  >>> lldb -- python stage/{prefix}/bin/psi4 ../tests/tu1-h2o-energy/input.dat
  >>> (lldb) run

::

  >>> cd {objdir}
  >>> gdb --args python stage/{prefix}/bin/psi4 ../tests/tu1-h2o-energy/input.dat
  >>> (gdb) run


.. .. _`faq:valgrindpsi`:
.. 
.. How to use ``valgrind`` with Psi4
.. ---------------------------------
.. 
.. When you naively use Valgrind with Psi4, you're likely to get incomprehensible mess of garbage or it may just crash with a boost overflow error. This happens because the boost python layer looks really really bad as far as Valgrind is concerned, i.e. it looks like a ton of memory leaks.  It really isn't, so we want to ignore all such errors/warnings. Valgrind has a mechanism for this in the way of suppression files.  Calling Valgrind as:
.. 
.. ```bash
.. valgrind --suppressions=<file_name>
.. ```
.. 
.. will run valgrind with the suppression file located on disk at "file_name". Lucky for you, Psi4 comes with such a suppression file at [``$top-level-psi4-dir/lib/valgrind-python.supp``](../blob/master/lib/valgrind-python.supp). This should remove all the boost python errors.
..  
.. The other error, boost overflow error arises from ``src/lib/libmints/sieve.cc`` where the inverse
.. of the complementary error function is being called.  The internet seems to claim that this is a
.. bug the arises only in debugging mode and has something to do with the exponent boost chooses for
.. the default zero tolerance.  Anyways, commenting out lines 47 to 49, for valgrind purposes, should
.. allow you to run valgrind.  The consequence of commenting out these lines are you get no integeral
.. screening, so make sure you uncomment them when you actually run.


.. _`faq:cmakeverbose`:

How to see the actual compiling commands (or errors) with ``cmake``
-------------------------------------------------------------------

CMake by default hides a lot of useful debugging information to make the
compilation cleaner. Issue ``make VERBOSE=1`` to display the full
compilation commands and errors.


.. _`faq:vigitmerge`:

How to highlight git merge conflicts in ``vi``
----------------------------------------------

Edit your ``~/.vimrc`` file to include the lines below. Hitting the ``F7``
key will toggle highlighting of git's conflict markers.

.. code-block:: bash

   >>> cat ~/.vimrc
   set hlsearch
   map <F7> :/\(<<<<<<<\\|=======\\|>>>>>>>\)<CR>


.. _`faq:libmwcondapy`:

How to handle "runtime library may be hidden" when building with Anaconda Python
--------------------------------------------------------------------------------

When building against Ana/Miniconda python (e.g., ``cmake
-DPYTHON_EXECUTABLE=/path/to/conda/bin/python``), the warning below often
appears. It is harmless, proceed.

.. code-block:: bash

   CMake Warning at src/bin/psi4/CMakeLists.txt:58 (add_executable):
     Cannot generate a safe runtime search path for target psi4 because files in
     some directories may conflict with libraries in implicit directories:

       runtime library [libm.so.6] in /usr/lib64 may be hidden by files in:
         /theoryfs2/common/software/anaconda/lib

   Some of these libraries may not be found correctly.


.. _`faq:psi4scratch`:

How to set up the scratch directory
-----------------------------------

The scratch directory is where Psi4 stores potentially large files during
computation. It should thus be on a local, fast disk to minimize any
computational inefficiencies caused by I/O. The scratch directory is
commonly set up through the :envvar:`PSI_SCRATCH` environment variable:

.. code-block:: tcsh

    # csh, tcsh: add to shell or ~/.tcshrc file
    setenv PSI_SCRATCH /path/to/existing/writable/local-not-network/directory/for/scratch/files

.. code-block:: bash

    # sh, bash: add to shell or ~/.bashrc (Linux) or ~/.bash_profile (Mac) file
    export PSI_SCRATCH=/path/to/existing/writable/local-not-network/directory/for/scratch/files

See also the more general :ref:`scratch documentation <sec:Scratch>`.


.. _`faq:psi4fileretention`:

How do I retain specific Psi4 scratch files
-------------------------------------------

You can set up a specific path for |PSIfour| scratch file and keep them
for later use through the :ref:`psi4_io <sec:Scratch>` handler.


.. _`faq:psi4PBS`:

How to use Psi4 within a PBS queue
----------------------------------

You will usually need to set up a PBS job script that is setting all
necessary environment variables, making sure the scratch directories are
set up, and invokes the executable. An :ref:`example <sec:PBS>` PBS script
is provided in the manual, but make sure to also consult your own PBS
documentation for appropriate setup.


.. _`faq:recompile`:

How to update and rebuild Psi4
------------------------------

Obtain code updates as appropriate from LINKTOVARMODIES. Move into
:samp:`{objdir}` and reissue ``make``, whereupon CMake may reconfigure but
will only rebuild objects and libraries depending on changed files. It is
scarcely ever necessary for the user to reinvoke ``cmake`` to update
:samp:`{objdir}`.


.. _`faq:minutetests`:

How to run a minute's worth of tests
------------------------------------

When you want to do a very minimal test of the build and have
CTest installed, the following command can be useful. ::

    >>> ctest -L smoke -j`getconf _NPROCESSORS_ONLN`

If you have pytest installed, very similar coverage is obtained through::

    >>> make pytest

.. _`faq:subsettests`:

How to run a subset of tests
----------------------------

CTest allows flexibly partitioned running of the test suite. In
the examples below, *testname* are regex of :source:`test names <tests>`,
and *testlabel* are regex of labels (*e.g.*, ``cc``, ``mints``,
``libefp``).

* Run tests in parallel with ``-j`` flag. For maximum parallelism: :samp:`ctest -j\`getconf _NPROCESSORS_ONLN\`\ `
* Run full test suite: ``ctest``
* Run about a third of the tests in 10--20 minutes, the so-called *quicktests*: ``ctest -L quick``
* Run the minimal number of tests to ensure Psi4 and any add-ons in working order: ``ctest -L smoke``
* Run tests matching by name: ``ctest -R testname``
* Run tests excluding those by name: ``ctest -E testname``
* Run tests matching by label: ``ctest -L testlabel``
* Run tests excluding those by label: ``ctest -LE testlabel``


.. _`faq:testsoutput`:

How to see CTest testing errors
-------------------------------

::

 >>> ctest
 Test project /your/path/2/psi4/build/directory/tests
     Start 248: tu1-h2o-energy
 1/2 Test #248: tu1-h2o-energy ...................   Passed    1.73 sec
      Start  6: cc1
 2/2  Test  #6: cc1 ..............................***Failed    0.07 sec
 ...

When ``ctest`` reports that some (or all) tests have failed, look in your
build directory for file
:samp:`{objdir}/tests/Testing/Temporary/LastTest.log`. It may have a
``.tmp`` extension, depending on whether the last test was interrupted and
a few other factors. Either way, this file should contain CMake's testing
output, as well as everything that was printed to the screen.


.. _`faq:pytest`:

How to test a Psi4 installation
-------------------------------

``ctest`` requires a connection to source files and ``cmake``
machinery and so can only be performed from :samp:`{objdir}`
(staged installation). To test an installed |PSIfour| (full or staged
installation), a limited number of "smoke" tests are available to be
run via pytest.

  * From the executable::

    psi4 --test

  * From the library (|PSIfour| must be detectable as a Python
    module. See the setup note at :ref:`faq:runordinarymodule`
    if needed.)::

    python -c "import psi4; psi4.test()"


.. _`faq:writepsi4`:

How to refer to Psi4
--------------------

Ways to refer to |PSIfour| in text, in order of decreasing goodness:

  * as ``Psi4`` in Optima regular font with "si" in custom small caps
    according to :source:`media/README.md`.

    * html: ``<span style="font-family: Optima, sans-serif; color: #273896;">P<span style="font-size: 82%;">SI</span>4</span>``

  * as ``Psi4`` with "si" in generated small caps

    * html: ``<span style="font-variant: small-caps;">Psi4</span>``

  * as ``Psi4`` with "si" in lowercase

  * as ``psi4`` in code

  * **NOT** ``PSI4`` or ``PSI``

