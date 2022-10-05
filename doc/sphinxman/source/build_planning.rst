.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2022 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This file is part of Psi4.
.. #
.. # Psi4 is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU Lesser General Public License as published by
.. # the Free Software Foundation, version 3.
.. #
.. # Psi4 is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU Lesser General Public License for more details.
.. #
.. # You should have received a copy of the GNU Lesser General Public License along
.. # with Psi4; if not, write to the Free Software Foundation, Inc.,
.. # 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
.. #
.. # @END LICENSE
.. #


.. index:: prerequisites, compiling, installing
.. _`sec:installFile`:

====================================
Compiling and Installing from Source
====================================

This section fulfills the duty of every Linux project to have an INSTALL
file detailing how to build it from source. Few (mostly packagers and
sysadmins) are expected to read this file as it is a broad guide to the
many choices in building this project. For a specific, well-tested,
and performance-tuned build path, see :psicode:`installs/latest`
(select "source"). For pre-built binaries again well-tested,
performance-tuned, and available for all common operating systems,
see :psicode:`installs/latest` (select "conda").


.. _`faq:cmakeviasetup`:

Planning: how to configure Psi4 and invoke CMake
------------------------------------------------

|PSIfour| is built through CMake. An abbreviated build guide can be found
:source:`within the source itself <CMakeLists.txt#L22>` .

CMake does a good job scanning your computer to locate libraries, header
files, and executables needed for compilation. So it's very possible that
from :samp:`{top-level-psi4-dir}` source directory, you can call :samp:`cmake -S.
-B{objdir}` without any further arguments, and it will invoke ``cmake``,
detect some appropriate defaults, configure the build, create a build
directory :samp:`{objdir}`, and complete, leaving you to only invoke
``make`` or Ninja in the new build directory.

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
<faq:builddetailed>`.

.. code-block:: console

    >>> cd {top-level-psi4-dir}
    >>> cmake -S. -Bobjdir [your configuration options]
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

    >>> cmake -S. -B{objdir} -DCMAKE_INSTALL_PREFIX={prefix} [your configuration options]

**4. Compile**

   Compile the code (optional ``-j`` triggers parallel compilation).

.. code-block:: console

    >>> cd {objdir}
    >>> make -j`getconf _NPROCESSORS_ONLN`

**5. Test**

   Optionally, use CTest and Pytest to test the build.

   * :ref:`faq:minutetests`
   * :ref:`faq:subsettests`
   * :ref:`faq:testsoutput`
   * :ref:`faq:pytest`

.. code-block:: console

   >>> ctest -j`getconf _NPROCESSORS_ONLN`

   >>> make pytest

**6. Install**

   If tests pass, install the code.

   ::

   >>> make install

**7. Configure Runtime**

   To run |PSIfour| after installation, you need to configure a few variables:

   * :ref:`faq:runordinaryexe`
   * :ref:`faq:runordinarymodule`


.. _`faq:coredepend`:

What are the tools and dependencies strictly required for building Psi4
-----------------------------------------------------------------------

The core |PSIfour| build requires the software below. Note that
practically everything (including Python, CMake, NumPy, BLAS/LAPACK,
Libint, and even C++ compilers on Linux and Mac) can be
satisfied through conda. The links below give examples of how to configure
that software for |PSIfour| and any notes and warnings pertaining to it.

* :ref:`C++ and C Compilers <cmake:cxx>` (C++17 compliant)

* :ref:`Optimized BLAS and LAPACK libraries <cmake:lapack>` (preferably NOT one supplied by a standard
  Linux distribution)

* :ref:`Python interpreter and headers <cmake:python>` (3.8+) https://www.python.org/

* CMake (3.15+) https://cmake.org/download/

* NumPy (needed at runtime *and* buildtime) https://numpy.org/

* System utilities: GNU make and GNU install or Ninja, POSIX threads (Pthreads) library

The following are also required for |PSIfour|, but if not detected, the
build system will automatically download and build.

* :ref:`gau2grid <cmake:gau2grid>` |w---w| :ref:`[what is gau2grid?] <sec:gau2grid>` :source:`[gau2grid min version] <external/upstream/gau2grid/CMakeLists.txt#L1>`

* :ref:`Libint <cmake:libint>` |w---w| :ref:`[what is Libint?] <sec:libint>` :source:`[Libint min version] <external/upstream/libint/CMakeLists.txt#L1>` (Libint2 as of Nov 2020; added by v1.4)

  * Eigen https://eigen.tuxfamily.org/index.php?title=Main_Page

* :ref:`Libxc <cmake:libxc>` |w---w| :ref:`[what is Libxc?] <sec:libxc>` :source:`[Libxc min version] <external/upstream/libxc/CMakeLists.txt#L1>`
* pybind11 |w---w| `[what is Pybind11?] <https://pybind11.readthedocs.io/en/stable/>`_ :source:`[Pybind11 min version] <external/upstream/pybind11/CMakeLists.txt#L1>`
* QCElemental |w---w| `[what is QCElemental?] <https://qcelemental.readthedocs.io/en/latest/>`_

* QCEngine |w---w| `[what is QCEngine?] <https://qcengine.readthedocs.io/en/latest/>`_ (March 2019; added by v1.4)

Additionally, there are runtime-only dependencies:

* msgpack-python for MessagePack https://msgpack.org/ (transitive dependency of QCElemental; optional for QCElemental but required for |PSIfour|)

* NumPy https://numpy.org/

* networkx https://github.com/networkx/networkx (transitive dependency of QCElemental; optional for QCElemental but required for |PSIfour|)

* pint https://pint.readthedocs.io/en/latest/ (transitive dependency of QCElemental)

* pydantic https://pydantic-docs.helpmanual.io/# (transitive dependency of QCElemental)

* SciPy https://scipy.org/ (avoidable if necessary by following directions in runtime error message)

Dropped Dependencies:

* Boost (September 2016; dropped by v1.1)

* deepdiff https://github.com/seperman/deepdiff (May 2019; dropped by v1.4)

* MPFR https://www.mpfr.org/ (Mar 2022; dropped by v1.6) no longer needed to compile against Libint2

* Perl for testing (April 2022; dropped by v1.6)


.. _`faq:addondepend`:

What are the add-on capabilities for Psi4 and what are their dependencies
-------------------------------------------------------------------------

Each of the items below is an independent additional capability that can
be built with |PSIfour|. Sub-items below are the respective additional
dependencies of the add-on. Select which, if any, you want, and examine
the links for appropriate enabling arguments to ``cmake``. Note that many
are available pre-built from conda.

* |PSIfour| Testing

  * CTest https://cmake.org/download/
  * pytest (for installed testing) http://doc.pytest.org/en/latest/
  * pytest-xdist (for installed testing many tests at once) https://github.com/pytest-dev/pytest-xdist

* |PSIfour| Documentation (available pre-built at :psicode:`psi4manual/master/index.html`)

  * Sphinx (3.5+) http://sphinx-doc.org
  * Perl (for some auto-documentation scripts) https://www.perl.org/
  * nbsphinx (for converting Jupyter notebooks) http://nbsphinx.readthedocs.io/en/jupyter-theme/
  * sphinx-psi-theme https://github.com/psi4/sphinx-psi-theme
  * See `["message" lines] :source:`doc/sphinxman/CMakeLists.txt` for advice on obtaining docs dependencies
  * See :source:`.github/workflows/docs.yml` for full docs building procedure to follow

* Ambit |w---w| https://github.com/jturney/ambit

* :ref:`CheMPS2 <cmake:chemps2>` |w---w| :ref:`[what is CheMPS2?] <sec:chemps2>` :source:`[CheMPS2 min version] <external/upstream/chemps2/CMakeLists.txt#L2>`

  * HDF5 https://support.hdfgroup.org/HDF5/
  * zlib http://www.zlib.net/

.. * :ref:`erd <cmake:erd>` |w---w| :ref:`[what is this?] <sec:erd>` `[min version] <https://github.com/psi4/psi4/blob/master/external/upstream/erd/CMakeLists.txt#L2>`_

..  * :ref:`Fortran Compiler <cmake:fortran>`

* :ref:`ecpint <cmake:ecpint>` |w---w| :ref:`[what is ecpint?] <sec:ecpint>` :source:`[ecpint min version] <external/upstream/ecpint/CMakeLists.txt/#L2>`

* :ref:`dkh <cmake:dkh>` |w---w| :ref:`[what is dkh?] <sec:dkh>` :source:`[dkh min version] <external/upstream/dkh/CMakeLists.txt#L2>`

  * :ref:`Fortran Compiler <cmake:fortran>`

* :ref:`gdma <cmake:gdma>` |w---w| :ref:`[what is gdma?] <sec:gdma>` :source:`[gdma min version] <external/upstream/gdma/CMakeLists.txt#L2>`

  * :ref:`Fortran Compiler <cmake:fortran>`

* :ref:`PCMSolver <cmake:pcmsolver>` |w---w| :ref:`[what is PCMSolver?] <sec:pcmsolver>`

  * :ref:`Fortran Compiler <cmake:fortran>`
  * zlib http://www.zlib.net/

* :ref:`simint <cmake:simint>` |w---w| :ref:`[what is simint?] <sec:simint>` :source:`[simint min version] <external/upstream/simint/CMakeLists.txt#L2>`

Additionally, there are runtime-loaded capabilities:

* :ref:`PylibEFP & libefp <cmake:libefp>` |w---w| :ref:`[what is LibEFP?] <sec:libefp>` :source:`[LibEFP min version] <external/upstream/libefp/CMakeLists.txt#L1>`

* cfour |w---w| :ref:`[what is CFOUR?] <sec:cfour>`

* dftd3 |w---w| :ref:`[what is DFTD3?] <sec:dftd3>`

* gcp |w---w| :ref:`[what is gCP?] <sec:gcp>`

* mrcc |w---w| :ref:`[what is MPCC?] <sec:mrcc>`

* v2rdm_casscf |w---w| :ref:`[what is v2rdm_casscf?] <sec:v2rdm_casscf>`

* :ref:`sns-mp2 <sec:snsmp2>` |w---w| https://github.com/DEShawResearch/sns-mp2/commits/master

* resp |w---w| https://github.com/cdsgroup/resp

* gpu_dfcc |w---w| https://github.com/edeprince3/gpu_dfcc/commits/master

* OpenFermion-Psi4 |w---w| https://github.com/quantumlib/OpenFermion

* :ref:`cppe <sec:cppe>` |w---w| https://github.com/maxscheurer/cppe

* mp2d |w---w| https://github.com/Chandemonium/MP2D

* :ref:`FockCI <sec:fockci>` |w---w| https://github.com/shannonhouck/psi4fockci

* cct3 |w---w| https://github.com/piecuch-group/cct3

* :ref:`adcc <sec:adcc>` |w---w| https://github.com/adc-connect/adcc

* :ref:`BrianQC <sec:brianqc>` |w---w| https://www.brianqc.com/

* i-PI |w---w| http://ipi-code.org/

* psixas |w---w| https://github.com/Masterluke87/psixas

* MDI |w---w| https://github.com/MolSSI-MDI/MDI_Library


.. _`faq:condapsi4dev`:

How to use conda to get a Psi4 development environment
------------------------------------------------------

|PSIfour| has a lot of dependencies, some of which are complicated
to build. Fortunately, conda can supply them all. A conda
package ``psi4-dev`` is available to gather the dependencies
and build tools needed to compiler |PSIfour| from source. Go to
:psicode:`installs/latest`, select "source" and "nightly build" along
with Linux/Mac and a Python version. Follow its instructions to create
a conda environment that has a script ``psi4-path-advisor``. Running
``psi4-path-advisor --help`` shows more options. Running the script
without arguments gives a cmake command that supplies paths to dependency
conda packages to CMake's cache. You can add further cmake options to
this command when you execute it. Once the cmake command has built your
:samp:`{objdir}`, you can build |PSIfour|. Remember to activate this
conda environment during |PSIfour| running or development. There are a
few extra options to ``psi4-path-advisor`` to set up using conda's GNU
compilers (on Linux), conda's Clang compilers (on Mac; less reliably),
or your own Intel compilers backed by conda's GNU compilers (on Linux).

.. code-block:: console

   >>> git clone https://github.com/psi4/psi4.git && cd psi4
   >>> conda create -n p4dev psi4-dev python=3.8 -c psi4/label/dev
   >>> conda activate p4dev
   >>> `psi4-path-advisor --gcc`
   >>> cd objdir && make -j`getconf _NPROCESSORS_ONLN`


.. _`faq:setupmaxameri`:

How to configure code to use high angular momentum basis sets
-------------------------------------------------------------

The :ref:`Libint <sec:libint>` integral code handles
arbitrary order angular momentum (AM), but compiling that is
prohibitive. The build process for Libint2 takes longer than
for Libint1, so it's recommended to :ref:`use the conda packages
<faq:condamaxameri>`. The AM controlling keyword ``MAX_AM_ERI``
has definition changed a little from Libint1. It is used to ensure
that conventional energies can use at least ``MAX_AM_ERI`` zeta. See
:source:`/external/upstream/libint2/CMakeLists.txt` for a fuller
discussion of AM capabilities. Unless :makevar:`BUILD_Libint2_GENERATOR`
is set to ``ON``, an internal build of Libint2 will download a
pregenerated tarball of at least :makevar:`MAX_AM_ERI` AM. If you
want to build completely from generator source, turn on the keyword,
but parts can only run in serial and parts benefit from parallel, so
it's not recommended. All the AM generation settings would need to be
edited in the above ``CMakeLists.txt`` rather than passing them into
the |PSIfour| build directly.

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

Note that since |PSIfour| 1.1, it is possible to build Libint
independently (or install just the libint conda package), then have
any/all |PSIfour| builds detect that installation at compile-time.

To switch a finished build to a larger Libint, enter ``<objdir>``,
remove the result of ``make install`` (if ever invoked), remove
``external/upstream/libint/``. If earlier Libint was built
internally (as opposed to detecting an external installation),
remove ``stage/include/libint/``, ``stage/share/cmake/Libint/``,
``stage/lib/libint2.so.2``. Edit MAX_AM_ERI
value in ``CMakeCache.txt``. Reinvoke ``make`` and check the number
beside ``Found Libint`` in the CMake output.

Note that the same keyword MAX_AM_ERI controls Libint and simint. simint
can't compile above AM7, so you'll have to turn off simint to get an
internal compile of AM>=8 Libint. The keyword is only a lower bound for
detection, so you can set it at 7, have both Libint and simint enabled,
and point CMake toward a pre-built Libint AM8 just fine.

* :ref:`cmake:libint`


.. _`faq:condamaxameri`:

How to get high angular momentum integrals from conda
-----------------------------------------------------

Libint2 conda packages are available from conda for Linux, Mac, and Windows.

The Linux package has base ``MAX_AM_ERI=7`` and can run 7z energies
and gradients conventional, 6z energies and gradients density-fitted,
and qz Hessians both conventional and density fitted.

The Mac and Windows packages have base ``MAX_AM_ERI=5`` and can run
5z energies, qz gradients, and tz Hessians, both conventional and
density-fitted.

Details about angular momentum settings are available here
:source:`external/upstream/libint2/CMakeLists.txt`
.

.. Since February 2019, the |PSIfour| conda package on Linux has been the
.. large AM8. Likewise, this is the package you get upon ::
..
..     conda install libint -c psi4
..
.. There is no need for the extra ``am8`` argument previously documented here.
..
.. On Mac, the Libint conda package itself and the Libint compiled into
.. the |PSIfour| conda package remain AM6.
..
.. Beware this issue: https://github.com/psi4/psi4/issues/1533


.. _`faq:setuphelp`:

How to see what build configuration options are available
---------------------------------------------------------

CMake doesn't provide a summary for this (unless you want to try the CMake
GUI, which the developers have never looked at). However, the top half of
the main CMakeLists.txt is a passable summary:

.. literalinclude:: @SFNX_INCLUDE@CMakeLists.txt
   :lines: 14-142
   :language: none

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

.. note:: It's not guaranteed, but if, in a pinch, you need to install a
   built Psi4 to a location *not* configured by ``CMAKE_INSTALL_PREFIX``,
   recursively copy the folders under :samp:`{objdir}/stage/{prefix}` to
   the desired location, ``chown`` them if needed, edit the shebang in
   ``bin/psi4`` if needed, and recursively delete all the ".pyc" files. It
   may just run.

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
..    cmake -S. -Bobj-gcc
..    cd obj-gcc


.. .. _`faq:erroreriam`:
..
.. How to fix error "``RuntimeError: value for ERI``"
.. --------------------------------------------------
..
.. You will need to rebuild Libint. Reissue ``cmake`` or edit
.. ``CMakeCache.txt`` with larger ``MAX_AM_ERI`` and rebuild.
..
.. * :ref:`faq:setupmaxameri`
.. * :ref:`faq:condamaxameri`

.. _`faq:chooseobjdir`:

How to choose the compilation directory, ``{objdir}``
-----------------------------------------------------

* there is no default
* common choices are ``objdir`` or ``build`` under :samp:`{top-level-psi4-dir}`

  * ``cd {top-level-psi4-dir} && cmake -S. -Bobjdir``
  * ``cd {top-level-psi4-dir} && cmake -S. -Bbuild``

* in-source builds (``*.cc`` and ``*.o`` in same directory) are disallowed
* builds *outside* :samp:`{top-level-psi4-dir}` are permitted


.. _`faq:doconfigure`:

How to save configuration settings for a future compilation
-----------------------------------------------------------

Create a file like ``do-configure`` with the ``cmake`` command and options
*on one line*. ::

 >>> cd {top-level-psi4-dir}
 >>> cat do-configure
     cmake -S. -B{objdir} \
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
below will exist at :samp:`{objdir}/stage`. This may be tested and used
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
 share/cmake/psi4/psi4ConfigVersion.cmake                                (psi4 cmake version info)
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
 lib/psi4/metadata.py                                                          (psi4 version info)
 lib/psi4/__init__.py                                         (module marker/loader for psi4.core)
 lib/psi4/core.cpython-*.so                               (c-side, compiled and bound by pybind11)
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
:samp:`{objdir}/stage` for :samp:`{prefix}`.)

.. code-block:: tcsh

    # csh, tcsh: add to shell or ~/.tcshrc file
    setenv PATH {prefix}/bin:$PATH
    setenv PSI_SCRATCH /path/to/existing/writable/local-not-network/directory/for/scratch/files

.. code-block:: bash

    # sh, bash: add to shell or ~/.bashrc (Linux/Windows) or ~/.bash_profile (Mac) file
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


.. _`faq:modulenotfounderror`:

How to solve ``ModuleNotFoundError: No module named 'psi4'``
------------------------------------------------------------

You're using |PSIfour| in PsiAPI mode (``python input.py`` where
``input.py`` contains ``import psi4``) but the Python interpreter
can't find the |PSIfour| Python module. Follow the directions at
:ref:`faq:psi4psiapipath`.


.. _`faq:psi4psiapipath`:

How to configure paths for PsiAPI
---------------------------------

If you know the location of the |PSIfour| executable (``bin/psi4``)
for Psithon mode and want to know the corresponding location to add to
:envvar:`PYTHONPATH` for PsiAPI mode, execute ``psi4 --psiapi``. It
will return bash commands to set :envvar:`PATH` (for correct python
interpreter) and :envvar:`PYTHONPATH` (to find psi4 module) correctly,
after which ``import psi4`` will work.

.. code-block:: bash

    >>> psi4 --psiapi-path
    export PATH=/path/to/dir/of/python/interpreter/against/which/psi4/compiled:$PATH
    export PYTHONPATH=/path/to/dir/of/psi4/core-dot-so:$PYTHONPATH

    >>> export PATH=/path/to/dir/of/python/interpreter/against/which/psi4/compiled:$PATH
    >>> export PYTHONPATH=/path/to/dir/of/psi4/core-dot-so:$PYTHONPATH

    >>> python -c "import psi4"

    >>> python -c "import psi4;print(psi4.__file__, psi4.__version__)"


.. _`faq:runordinarymodule`:

How to run Psi4 as Python module after compilation
--------------------------------------------------

Substituting the full installation directory :samp:`{prefix}` and a
suitable scratch directory, issue the following commands directly in your
terminal or place them into your "rc" file and open a new terminal. (To
use a staged installation directory, substitute
:samp:`{objdir}/stage` for :samp:`{prefix}`.)

.. code-block:: tcsh

    # csh, tcsh: add to shell or ~/.tcshrc file
    setenv PYTHONPATH {prefix}/lib:$PYTHONPATH
    setenv PSI_SCRATCH /path/to/existing/writable/local-not-network/directory/for/scratch/files

.. code-block:: bash

    # sh, bash: add to shell or ~/.bashrc (Linux/Windows) or ~/.bash_profile (Mac) file
    export PYTHONPATH={prefix}/lib:$PYTHONPATH
    export PSI_SCRATCH=/path/to/existing/writable/local-not-network/directory/for/scratch/files

* :ref:`faq:psi4psiapipath`

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
{ana|mini|psi4}conda`. To see the message again after initial installation,
with the conda environment active, run ``.psi4-post-link.sh``.
If ``which conda python psi4`` points to your
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
:samp:`conda activate {p4env}`, then ``which psi4`` likely yields
:samp:`{condadist}/envs/{p4env}/bin/psi4` and the ``PATH`` setting lines
below are redundant.

.. code-block:: tcsh

    # csh, tcsh: add to shell or ~/.tcshrc file
    unsetenv PSIDATADIR
    setenv PATH {prefix}/bin:$PATH
    setenv PSI_SCRATCH /path/to/existing/writable/local-not-network/directory/for/scratch/files

.. code-block:: bash

    # sh, bash: add to shell or ~/.bashrc (Linux/Windows) or ~/.bash_profile (Mac) file
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
    >>> ln -s {top-level-psi4-dir}/{objdir}/stage/lib/psi4/core.cpython-{ext_will_vary}.so ../psi4/core.cpython-{ext_will_vary}.so
    >>> python ../psi4/run_psi4.py --inplace input.dat


.. _`faq:psidatadir`:

Why not to set :envvar:`PSIDATADIR`
-----------------------------------

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

    # sh, bash: add to shell or ~/.bashrc (Linux/Windows) or ~/.bash_profile (Mac) file
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

F. Build with specific (Intel) compilers from :envvar:`PATH` based on GCC
   with prefix and *not* in :envvar:`PATH`
   (``GCCPFX=/full/path/to/bin/prefix-`` compiler is ``$GCCPFX-gcc``)

  .. code-block:: bash

    >>> cmake -DCMAKE_C_COMPILER=icc \
              -DCMAKE_CXX_COMPILER=icpc \
              -DCMAKE_C_FLAGS="-gnu-prefix=${GCCPFX}" \
              -DCMAKE_CXX_FLAGS="-gnu-prefix=${GCCPFX}"

G. Build on Linux with specific (Intel) compilers from :envvar:`PATH`
   based on GCC from conda in **activated** environment
   (:envvar:`CONDA_PREFIX` and :envvar:`HOST` are defined upon
   activation)

  .. code-block:: bash

    >>> cmake -DCMAKE_C_COMPILER=icc \
              -DCMAKE_CXX_COMPILER=icpc \
              -DCMAKE_C_FLAGS="-gnu-prefix=${CONDA_PREFIX}/bin/${HOST} --sysroot=${CONDA_PREFIX}/${HOST}/sysroot" \
              -DCMAKE_CXX_FLAGS="-gnu-prefix=${CONDA_PREFIX}/bin/${HOST} --sysroot=${CONDA_PREFIX}/${HOST}/sysroot"

H. Build on Linux with specific (GCC) compilers from
   from conda in **activated** environment
   (:envvar:`CONDA_PREFIX` and :envvar:`HOST` are defined upon
   activation)

  .. code-block:: bash

    >>> cmake -DCMAKE_C_COMPILER=${GCC} \
              -DCMAKE_CXX_COMPILER=${GXX} \
              -DCMAKE_Fortran_COMPILER=${GFORTRAN}


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
compilers on Linux also rely on GCC, so both ``icpc`` and ``gcc`` versions are checked.
Intel OneAPI Classic compilers work fine. OneAPI beta compilers build but have
been only minimally tested.

* :ref:`faq:modgcc`


.. _`faq:macxcode`:

How to obtain C and C++ compilers for Mac without Fink, MacPorts, or Homebrew
-----------------------------------------------------------------------------

The easiest compiler to obtain is ``clang`` which is a drop-in replacement
for ``gcc`` and ``g++``. Just install `XCode
<https://apps.apple.com/us/app/xcode/id497799835>`_. Some old versions
of XCode can't handle some of the advanced C++ language features, but this
is a *software* not *hardware* limitation. Checks for version compliance
performed at build-time. Note that this "AppleClang" will not be compatible
with conda Mac packages using C++11, nor can it make use of OpenMP directives.

Another route to obtaining ``clang`` compilers without the above limitations
is through conda.

.. code-block:: bash

   # Install Clang 4.0.1 into a non-primary conda environment
   >>> conda create -n clang401 clangxx_osx-64 clang_osx-64 llvm-openmp intel-openmp

   # To Build, activate environment (prepends PATH and defines environment variables CLANG, CLANGXX, HOST, etc):
   >>> conda activate clang401
   >>> echo ${CLANGXX}
   /path/to/miniconda/envs/clang401/bin/x86_64-apple-darwin13.4.0-clang++
   >>> echo ${HOST}
   x86_64-apple-darwin13.4.0

   # build with Clang
   >>> cmake -S. -Bbuild \
        -DCMAKE_C_COMPILER=${CLANG} \
        -DCMAKE_CXX_COMPILER=${CLANGXX} \
        -DCMAKE_CXX_FLAGS="-stdlib=libc++" \
        -DOpenMP_CXX_FLAG="-fopenmp=libiomp5"

   # build with Intel
   >>> cmake -S. -Bbuild \
        -DCMAKE_C_COMPILER=icc \
        -DCMAKE_CXX_COMPILER=icpc \
        -DCMAKE_C_FLAGS="-clang-name=${CLANG}" \
        -DCMAKE_CXX_FLAGS="-clang-name=${CLANG} -clangxx-name=${CLANGXX} -stdlib=libc++ -I${CONDA_PREFIX}/include/c++/v1"

   # Configure and build

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

   # Install GCC 7.2 into a non-primary conda environment
   >>> conda create -n gcc72 gxx_linux-64 gcc_linux-64

   # To Build, either:

   # (A) activate environment (prepends PATH and defines environment variables CC, CXX, HOST, etc)
         >>> conda activate gcc72
         >>> echo ${CXX}
         /path/to/miniconda/envs/gcc72/bin/x86_64-conda_cos6-linux-gnu-g++
         >>> echo ${HOST}
         x86_64-conda_cos6-linux-gnu

         # build with GNU
         >>> cmake -S. -Bbuild \
              -DCMAKE_C_COMPILER=${CC} \
              -DCMAKE_CXX_COMPILER=${CXX} \

         # build with Intel
         >>> cmake -S. -Bbuild \
              -DCMAKE_C_COMPILER=icc \
              -DCMAKE_CXX_COMPILER=icpc \
              -DCMAKE_C_FLAGS="-gnu-prefix=${HOST}-" \
              -DCMAKE_CXX_FLAGS="-gnu-prefix=${HOST}-" \

   # (B) tell CMake to tell the compiler which GCC to use
         >>> GCC7=/path/to/miniconda/envs/gcc72
         >>> cmake -S. -Bbuild \
              -DCMAKE_C_COMPILER=icc \
              -DCMAKE_CXX_COMPILER=icpc \
              -DCMAKE_C_FLAGS="-gnu-prefix=${GCC7}/bin/x86_64-conda_cos6-linux-gnu-" \
              -DCMAKE_CXX_FLAGS="-gnu-prefix=${GCC7}/bin/x86_64-conda_cos6-linux-gnu-" \
              ...
              # if Fortran active ...
              -DCMAKE_Fortran_COMPILER=ifort \
              -DCMAKE_Fortran_FLAGS="-gnu-prefix=${GCC7}/bin/x86_64-conda_cos6-linux-gnu-" \

   # Configure and build

   # To Run:
   >>> export LD_LIBRARY_PATH=${GCC7}/lib:$LD_LIBRARY_PATH


.. _`faq:cray`:

How to configure a Psi4 build on Cray
-------------------------------------

Cray systems strongly prefer to build static libraries, but |PSIfour|
needs to be dynamic to function as a Python module. Courtesy of @misha
at the forum and various supercomputer guides, building |PSIfour| on
Cray requires setting environment variables before running `cmake`. ::

    CRAYPE_LINK_TYPE=dynamic CRAY_ADD_RPATH=yes cmake ...

.. envvar:: CRAYPE_LINK_TYPE

   Set to allow |PSIfour| to build as shared library on Cray supercomputers.

.. envvar:: CRAY_ADD_RPATH

   Set to true to allow |PSIfour| to build on Cray supercomputers.


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
  * conda ``gfortran_linux-64`` or ``gfortran_osx-64`` to get ``gfortran``


.. _`faq:macgfortran`:

How to obtain a Fortran compiler for Mac without Fink, MacPorts, or Homebrew
----------------------------------------------------------------------------

Xcode does not provide a Fortran compiler. A way to get one is to download
the ``gfortran_osx-64`` conda package. This provides
``gfortran`` compilers for Mac. The version is 4.8.5, which is quite old,
but the Fortran compiler will work.

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
* :makevar:`LAPACK_LIBRARIES` |w---w| CMake variable to specify BLAS/LAPACK libraries explicitly, bypassing math detection. Should be ";"-separated list of full paths.
* :makevar:`LAPACK_INCLUDE_DIRS` |w---w| CMake variable to specify BLAS/LAPACK header location explicitly, bypassing math detection. Only needed for MKL.
* :makevar:`OpenMP_LIBRARY_DIRS` |w---w| CMake variable to specify OpenMP library (iomp5/gomp/omp) directories explicitly. Should be ";"-separated list of full directory paths. Usually the solution to error ``Could NOT find MathOpenMP``.

.. envvar:: MKL_ROOT

  Environment variable set by Intel compilervars scripts. Sufficient to trigger math detection of MKL at this location.

.. envvar:: MATH_ROOT

   Environment variable to specify root directory in which BLAS/LAPACK libraries should be detected (*e.g.*, ``${MATH_ROOT}/lib64/libblas.so`` and ``${MATH_ROOT}/lib64/liblapack.so``).

**Examples**

A. Build with any LAPACK in standard location

  .. code-block:: bash

    >>> cmake

B. Build with native Accelerate LAPACK on Mac (MKL *not* also present).
   If NumPy *not* using native Accelerate LAPACK, then directing Psi4
   to use it is Bad Idea!

  .. code-block:: bash

    >>> cmake

C. Build with native Accelerate LAPACK on Mac (MKL also present)
   If NumPy *not* using native Accelerate LAPACK, then directing Psi4
   to use it is Bad Idea!

  .. code-block:: bash

    >>> cmake -DBLAS_TYPE=SYSTEM_NATIVE -DLAPACK_TYPE=SYSTEM_NATIVE

D. Build with Intel MKL

  .. code-block:: bash

    >>> source /path/to/intel/vers/linux/mkl/bin/mklvars.sh intel64  # adjust sh/csh and arch as needed
    >>> cmake

  .. code-block:: bash

    >>> MATH_ROOT=/path/to/intel/vers/linux/mkl/ cmake

E. Build with Intel MKL from conda (install ``mkl-devel`` package from defaults channel)

  .. code-block:: bash

    >>> cmake -DLAPACK_LIBRARIES="${CONDA_PREFIX}/lib/libmkl_rt.so" -DLAPACK_INCLUDE_DIRS="${CONDA_PREFIX}/include"

F. OpenBLAS - see note below.

  .. code-block:: bash

    >>> MATH_ROOT=/path/to/openblas/0.2.13_seq/x86_64/gcc_5.2.0/lib cmake

G. Build with explicit MKL LAPACK

  .. code-block:: bash

    >>> cmake -DLAPACK_LIBRARIES="/path/to/lib/intel64/libmkl_lapack95_lp64.a;/path/to/lib/intel64/libmkl_rt.so" -DLAPACK_INCLUDE_DIRS="/path/to/mkl-h-include/"

H. Build with explicit non-MKL LAPACK

  .. code-block:: bash

    >>> cmake -DLAPACK_LIBRARIES="/path/to/lib/liblapack.so;/path/to/lib/libblas.a"

I. Build with MKL and GCC (iomp5 needed instead of gomp for threading. use OpenMP_LIBRARY_DIRS to hint location.)

  .. code-block:: bash

    >>> cmake -DLAPACK_LIBRARIES=/opt/intel/mkl/lib/intel64/libmkl_rt.so -DLAPACK_INCLUDE_DIRS=/opt/intel/mkl/include -DOpenMP_LIBRARY_DIRS=/opt/intel/compiler/lib/intel64/

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
    v1.1, earliest known working version is MKL 2013; as of v1.4,
    should use at leat MKL 2019). MKL, which is
    freely available through conda, is the only threaded BLAS/LAPACK
    distribution fully supported by |PSIfour|.

  * On Mac, the native Accelerate libraries are very nice and would
    be recommended but for the potential conflict between |PSIfour|
    BLAS and NumPy BLAS. Unless you've a special NumPy, avoid!

  * The open-source LAPACK distributions OpenBLAS (formerly GotoBLAS)
    mostly works. Use it at your own risk and after testing your
    particular distribution, including tests run multithreaded,
    if you intend to run |PSIfour| so. Use at least 0.2.15, and
    pay attention to how it was compiled - unthreaded seems safe,
    openmp-threaded is mostly safe, default pthreaded is *not* safe. See
    https://github.com/psi4/psi4/issues/1009 for recent analysis.

  * Another open-source LAPACK distribution, ATLAS had
    stability issues with the DFOCC module at last testing,
    https://github.com/psi4/psi4/issues/391.

  * ACML libraries are known to work with |PSIfour| v1.1 at ACML 6.

* Because of how link loaders work, at runtime, the BLAS of |PSIfour|
  and the BLAS of NumPy are not independent. There can be unpredictable
  but reproducible numerical and thread-scaling errors if |PSIfour|
  and NumPy BLAS don't match down to the library name (that is,
  ``libmkl_rt``, ``libmkl_core.so``, ``libmkl_core.a`` are *not*
  interchangeable). See https://github.com/psi4/psi4/issues/1007,
  https://github.com/psi4/psi4/issues/748,
  https://github.com/psi4/psi4/issues/755 for gory discussions.
  Choose your NumPy and |PSIfour| compile conditions to use the same
  BLAS distribution.

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

Note that in v1.4, |PSIfour| switched from Numpy-written CMake Python
detection to new CMake-written CMake Python detection. Usually only
:makevar:`Python_EXECUTABLE` is now needed to hint Python's location. See
https://cmake.org/cmake/help/latest/module/FindPython.html for full
guidance.

* :makevar:`Python_EXECUTABLE` |w---w| specify name or full path to Python interpreter.
* :makevar:`Python_LIBRARY` |w---w| specify path to Python library.
* :makevar:`Python_INCLUDE_DIR` |w---w| specify directory of Python headers. Contains ``Python.h``.

**Examples**

A. Build with detected Python from :envvar:`PATH`

  .. code-block:: bash

    >>> cmake

B. Build with specific Python

  .. code-block:: bash

    >>> cmake -DPython_EXECUTABLE=/path/to/interp/python3.9

C. Build with full Python specification to root directory ``${PFXC}``

  .. code-block:: bash

    >>> cmake -DPython_EXECUTABLE="${PFXC}/bin/python" \
              -DPython_LIBRARY="${PFXC}/lib/libpython3.9.so" \
              -DPython_INCLUDE_DIR="${PFXC}/include/python3.9"


.. _`faq:runtimepython`:

What Python is Psi4 running
---------------------------

The Python detected at build-time is embedded into the |PSIfour|
executable. That is, the top line of ``bin/psi4`` is something like
``#!/path/to/miniconda/envs/p4deps/bin/python3.5``, and that's the Python
through which |PSIfour| is running, rather than the Python of ``which python``.
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
  >>> lldb -- python stage/bin/psi4 ../tests/tu1-h2o-energy/input.dat
  >>> (lldb) run

::

  >>> cd {objdir}
  >>> gdb --args python stage/bin/psi4 ../tests/tu1-h2o-energy/input.dat
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
.. will run valgrind with the suppression file located on disk at "file_name". Lucky for you, Psi4 comes with such a suppression file at [``$top-level-psi4-dir/psi4/share/psi4/scripts/valgrind-python.supp``](../blob/master/psi4/share/psi4/scripts/valgrind-python.supp). This should remove all the python errors.
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
-DPython_EXECUTABLE=/path/to/conda/bin/python``), the warning below often
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

    # sh, bash: add to shell or ~/.bashrc (Linux/Windows) or ~/.bash_profile (Mac) file
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
set up, and invoking the executable. An :ref:`example <sec:PBS>` PBS script
is provided in the manual, but make sure to also consult your own PBS
documentation for appropriate setup.


.. _`faq:recompile`:

How to update and rebuild Psi4
------------------------------

Obtain code updates as appropriate from :ref:`faq:binary`,
:ref:`faq:clonepsi4public`, or :ref:`faq:forkpsi4public`.  Move into
:samp:`{objdir}` and reissue ``make`` or ``cmake --build .``, whereupon CMake may reconfigure but
will only rebuild objects and libraries depending on changed files. It is
scarcely ever necessary for the user to reinvoke ``cmake`` to update
:samp:`{objdir}` based on changes to the |PSIfour| repository.
Upon switching conda environments, though, one should
make a new objdir since there are many full paths baked into
:samp:`{objdir}/CMakeCache.txt`.


.. _`faq:minutetests`:

How to run a minute's worth of tests
------------------------------------

When you want to do a very minimal test of the build and have
CTest installed, the following command can be useful.

.. code-block:: console

    >>> ctest -L smoke -j`getconf _NPROCESSORS_ONLN`

If you have pytest installed, very similar coverage is obtained through::

    >>> make pytest

.. _`faq:subsettests`:

How to run a subset of tests
----------------------------

CTest allows flexibly partitioned running of the test suite. In
the examples below, *testname* are regex of :source:`test names <tests>`,
and *testlabel* are regex of labels (*e.g.*, ``cc``, ``mints``,
``libefp`` defined :source:`[here, for example] <tests/ci-property/CMakeLists.txt#L3>` .

* Run tests in parallel with ``-j`` flag. For maximum parallelism (not Windows): :samp:`ctest -j\`getconf _NPROCESSORS_ONLN\`\ `
* Run full test suite: ``ctest``
* Run about a third of the tests in 5 minutes, the so-called *quicktests*: ``ctest -L quick``
* Run the same subset of tests that continuous integration checks (not the full test suite): ``ctest -L quick``
* Run the minimal number of tests to ensure Psi4 and any add-ons in working order: ``ctest -L smoke``
* Run tests matching by name: ``ctest -R testname``
* Run tests excluding those by name: ``ctest -E testname``
* Run tests matching by label: ``ctest -L testlabel``
* Run tests excluding those by label: ``ctest -LE testlabel``

Pytest has similar commands:

* Run tests in parallel with ``-n`` flag (if Python package ``pytest-xdist`` installed). For maximum parallelism: :samp:`pytest -n auto`
* Run full test suite: ``pytest``
* Run the quick tests: ``pytest -m quick``
* Run the minimal number of tests to ensure Psi4 and any add-ons in working order: ``pytest -m smoke``
* Run tests matching by name: ``pytest -k testname``
* Run tests excluding those by name: ``pytest -k "not testname"``
* Run tests matching by label: ``pytest -m testlabel``
* Run tests excluding those by label: ``pytest -m "not testlabel"``
* Run tests with complicated filtering by name: ``pytest -k "cc and not cc4 and not stdsuite"``
* Run tests and see test names: ``pytest -v``
* Run tests and see test names and sometimes more verbose errors: ``pytest -vv``

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
:samp:`{objdir}/Testing/Temporary/LastTest.log`. It may have a
``.tmp`` extension, depending on whether the last test was interrupted and
a few other factors. Either way, this file should contain CMake's testing
output, as well as everything that was printed to the screen.


.. _`faq:pytest`:

How to test a Psi4 installation
-------------------------------

``ctest`` requires a connection to source files and ``cmake``
machinery and so can only be performed from :samp:`{objdir}`
(staged installation). To test an installed |PSIfour| (full or staged
installation), a supplementary test suite is available to be
run via pytest. Running a limited number of "smoke" tests is
generally sufficient to verify an installation.

  * From the executable

    .. code-block:: bash

     # smoke tests only in serial
     psi4 --test
     # smoke tests only in parallel
     psi4 --test -n`getconf _NPROCESSORS_ONLN`
     # most tests
     psi4 --test full

  * From the library (|PSIfour| must be detectable as a Python
    module. See setup at :ref:`faq:psi4psiapipath`
    if needed.)

    .. code-block:: bash

     # smoke tests only
     python -c "import psi4; psi4.test('smoke')"
     # most tests
     python -c "import psi4; psi4.test()"

  * From pytest directly. If package ``pytest-xdist`` is installed, can run in parallel.

    .. code-block:: bash

     # smoke tests only
     pytest {prefix}/lib/{PYMOD_INSTALL_LIBDIR}/psi4/tests/ -m smoke -n auto
     # all tests
     pytest {prefix}/lib/{PYMOD_INSTALL_LIBDIR}/psi4/tests/ -n auto

Output looks something like the below. ``PASSED`` in green is good
(means test ran correctly); ``SKIPPED`` in yellow is good (means that
not all software required for test is available); ``XPASS`` or ``XFAIL``
in yellow is fine (unexpected pass or expected fail happens when we
include tests that need particular conditions (*e.g.*, multiple cores)
to run correctly); ``FAILED`` in red is bad. ::

    test_addons.py::test_gdma PASSED
    test_addons.py::test_mrcc SKIPPED
    test_addons.py::test_chemps2 PASSED
    test_addons.py::test_dftd3 PASSED
    test_addons.py::test_libefp PASSED
    test_addons.py::test_pcmsolver PASSED
    test_addons.py::test_erd PASSED
    test_addons.py::test_simint PASSED
    test_addons.py::test_json PASSED
    test_addons.py::test_cfour SKIPPED
    test_addons.py::test_v2rdm_casscf PASSED
    test_addons.py::test_grimme_3c PASSED
    test_addons.py::test_dkh PASSED
    test_psi4.py::test_psi4_basic PASSED
    test_psi4.py::test_psi4_cc PASSED
    test_psi4.py::test_psi4_cas PASSED
    test_psi4.py::test_psi4_dfmp2 PASSED
    test_psi4.py::test_psi4_sapt PASSED
    test_psi4.py::test_psi4_scfproperty PASSED


.. _`faq:writepsi4`:

How to refer to Psi4
--------------------

Ways to refer to |PSIfour| in text, in order of decreasing goodness:

  * as ``Psi4`` in Optima regular font with "si" in custom (82%) small caps
    according to :source:`media/README.md`.

    * html: ``<span style="font-family: Optima, sans-serif; color: #273896;">P<span style="font-size: 82%;">SI</span>4</span>``

  * as ``Psi4`` with "si" in generated small caps

    * html: ``<span style="font-variant: small-caps;">Psi4</span>``

  * as ``Psi4`` with "si" in lowercase

  * as ``psi4`` in code

  * **NOT** ``PSI4`` or ``PSI``


.. _`faq:psi4logos`:

How to get a Psi4 logo file
---------------------------

All image files are stored in https://github.com/psi4/psi4media


.. _`faq:localaddon`:

How to use a local Add-On repository in the Psi4 build
------------------------------------------------------

For each Add-On, |PSIfour| pulls source from a specific online Git
repository and a specific tag/branch/commit in it. This ensures success
of the |PSIfour| build, reproducibility of the runtime results, and
freedom for continued upstream development. Sometimes, you're the one
doing that development, and you need the CMake superbuild to pull source
from a local path rather than the approved codeset.

Find the ``CMakeLists.txt`` governing the target Add-On in
:source:`external` and make changes analogous to the below::

    #GIT_REPOSITORY https://github.com/jturney/ambit
    #GIT_TAG 1.0
    DOWNLOAD_COMMAND ""
    SOURCE_DIR "/path/to/ambit-directclone"

If you're changing the |PSIfour| repo codebase between compiles, there's
nothing more to do as CMake will handle the code rebuild deps for you.

If you're changing the local Add-On repo codebase between compiles,
CMake *does not* know when ``libaddon.[a|so|dylib]`` needs rebuilding. It
is recommended that the |PSIfour| build be initially configured with
``-DBUILD_SHARED_LIBS=ON`` (easier to notice changes). And to trigger
Add-On library rebuild, ``rm -rf {objdir}/external/upstream/addon/``
and ``rm -rf {objdir}/stage/share/cmake/AddOn``. This should
re-clone the Add-On, rebuild and install it, rebuild any parts of
|PSIfour| that interface to it, and relink the main ``core.so``.
If you're modifying the Add-On's file or directory structure, be
smart and ``rm`` all traces of it within ``{objdir}/stage/``,
especially any ``*.pyc`` files.

Alternatively and preferentially to the above, you can instead build and install the
Add-On library yourself, external to the |PSIfour| repository. This
is especially useful if you want to avoid full recompiles of the
Add-On at each change to the Add-On's source. Build the Add-On
library dynamically (``-DBUILD_SHARED_LIBS=ON``) and mind any
"Psi4 wants" in the Add-On's top-level CMakeLists.txt. Install the
Add-On and note the full path to ``AddOnConfig.cmake``. Pass
the path containing that file to |PSIfours| CMake as
``-DAddon_DIR=/path/to/config/usually/ending/in/share/cmake/AddON``
and build |PSIfour|. The main ``core.so`` should be dynamically linked
to your dev AddOn dynamic lib and update automatically when you rebuild
the AddOn lib. Naturally, you may need to delete ``core.so`` and remake
as needed.


.. _`faq:logging`:

How to use logging in Psi4
--------------------------

A log file is generated alongside the output file and has extension
``.log``. Selected key information is sent to it. The log file contents
are provisional and haphazard, so don't start relying upon them. It is
tentatively envisioned that log files will have much of the same results as output
files but more as structured data than narrative, nicely formatted data.
In some cases where external calls are involved, the most detailed
information may only be collectable in the log file. It is
worthwhile to check the log file if expected results are missing from
the usual output file. If the log file is empty, the logging level may
be set too low.

To change the logging level to DEBUG from the command-line::

  >>> psi4 --logging 10

To change the logging level to DEBUG in PsiAPI::

  >>> psi4.set_output_file(<filebase>, loglevel=10)


.. _`faq:findexemod`:

How to find the Psi4 module from the executable and vice versa
--------------------------------------------------------------

Because of the different conventions of C++/Linux and Python installation
layouts, the |PSIfour| executable and the |PSIfour| Python module aren't
at fixed locations relative to each other. Accessors have been added
so that the executable can show where the module is, and the module
can show where the executable is.

  .. code-block:: bash

    > /path/to/psi4/install/bin/psi4 --module
    /path/to/psi4/install/lib//
    > python -c "import psi4;print(psi4.executable)"
    /path/to/psi4/install/bin/psi4

