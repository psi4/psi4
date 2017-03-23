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

.. include:: autodoc_abbr_options_c.rst

.. _`sec:conda`:

Conda Binary Distribution
=========================

|PSIfour| is available as a pre-compiled binary for Linux and Mac architectures
through `Continuum Analytics
<https://store.continuum.io/cshop/anaconda/>`_, the company that produces
`Anaconda Python <http://docs.continuum.io/anaconda/index.html>`_ (a
full-fledged scientific python environment with package manager `conda
<http://conda.pydata.org/index.html>`_) and, more particularly, `Miniconda
<http://conda.pydata.org/miniconda.html>`_ (a lightweight python
distribution with same package manger `conda
<http://conda.pydata.org/index.html>`_). Some nice features for us:

* cross-platform (Linux only at present)

* no root, administrator, or sudo access required

* built with high-performance math libraries

* lightweight software stack (<100 MB w/o |PSIfour|; ~1 GB including |PSIfour|, numpy, and MKL)

* updated nightly so new features accessible

* standardizes python distribution so no need to find/install libpython packages

* add-ons (plugins, extra features requiring Fortran compiler, etc.) can be made available as conda packages

* develop |PSIfour| through plugins without a pre-existing development environment, see :ref:`sec:condaplugins`.

The |PSIfour| binary repository is at `Anaconda (formerly Binstar) <https://anaconda.org/psi4>`_.

For commands to get a default installation, go to :ref:`sec:psi4conda`
or the `psicode downloads page <http://psicode.org/downloads2.html>`_.
Users proficient with conda may prefer to consult :ref:`sec:condadetails`.
For more flexibility and a detailed explanation, go to
:ref:`sec:slowconda` and :ref:`sec:slowpsi4`.

.. _`sec:psi4conda`:

Psi4conda Installer
^^^^^^^^^^^^^^^^^^^

Sequence of commands to get you to a working |PSIfour| on Linux
or Mac. Installs Miniconda+Psi4 into ``$HOME/psi4conda`` and
the |PSIfour| executable into the main conda environment at
``$HOME/psi4conda/bin/psi4``.

.. code-block:: bash

    # Linux
    >>> curl -O "http://www.psicode.org/downloads/Psi4conda2-latest-Linux.sh" --keepalive-time 2
    >>> bash
    >>> bash Psi4conda-latest-Linux.sh -b -p $HOME/psi4conda  # agrees to license terms
    >>> echo "export PATH=$HOME/psi4conda/bin:\$PATH" >> ~/.bashrc
    # log out, log back in so conda and psi4 in path
    >>> psi4 "$(dirname $(which psi4))"/../share/psi4/samples/sapt1/test.in  # test installation. works b/c PSI_SCRATCH defaults to /tmp

.. code-block:: bash

    # Mac
    >>> curl -O "http://www.psicode.org/downloads/Psi4conda2-latest-MacOSX.sh" --keepalive-time 2
    >>> bash
    >>> bash Psi4conda-latest-MacOSX.sh -b -p $HOME/psi4conda  # agrees to license terms
    >>> echo "export PATH=$HOME/psi4conda/bin:\$PATH" >> ~/.bash_profile
    # log out, log back in so conda and psi4 in path
    >>> psi4 "$(dirname $(which psi4))"/../share/psi4/samples/sapt1/test.in  # test installation. works b/c PSI_SCRATCH defaults to /tmp

That last command tested that ``psi4`` is in your path, and it's finding
all the libraries it needs. Now you need only specify a scratch directory
(see :ref:`sec:Scratch`) by replacing the placeholder in the following:

.. code-block:: bash

    >>> echo "export PSI_SCRATCH=/path/to/existing/writable/local-not-network/directory/for/scratch/files" >> ~/.bashrc
    # log out, log back in so variable takes effect

All done!

.. note:: Above commands use bash for installation and set up your environment for bash at runtime. To use csh at runtime, follow the on-screen directions at the end of the installation or consult step 7 below.

.. _`sec:condadetails`:

Conda Proficients
^^^^^^^^^^^^^^^^^

The :ref:`sec:psi4conda` uses a `conda constructor
<https://github.com/conda/constructor>`_ to package up Miniconda,
the psi4 conda packages, the psi4 add-on conda packages, dependencies
thereof (possibly from particular channels), and the psi4 channel
as a default.  This is very convenient for novice users and robust
against differing channel settings in ``~/.condarc``. But proficient
conda users may prefer to treat ``psi4`` as a normal conda package and
not have another large Miniconda installation (including the hefty MKL)
lying around just for |PSIfour|. Installing just the ``psi4`` package
itself will get you |PSIfour|, whatever add-ons require linking in to
|PSIfour| (*e.g.*, CheMPS2 and PCMSolver), and the correct versions of
packages. However, just the ``psi4`` package won't get you add-ons that
don't need linking (*e.g.*, DFTD3 and v2rdm_casscf) or dependencies
from the "right" channels, which can be important for issues of fPIC
and libc++ vs. libstdc++. So ``conda create -c psi4 -n p4env psi4 dftd3
v2rdm_casscf`` *should* be equivalent to running the psi4conda installer,
but I wouldn't count on it. Instead, an `explicit environment spec
<http://conda.pydata.org/docs/using/envs.html#build-identical-conda-environments-with-urls>`_
will be available for download.

.. code-block:: bash

    # Linux
    >>> curl -o explicit-latest.sh "https://repo.continuum.io/miniconda/explicit2-latest-Linux-x86_64.txt"
    >>> conda create --name p4env --file explicitenv2-latest-Linux-x86_64.txt
    >>> source activate p4env

.. code-block:: bash

    # Mac
    >>> curl -o explicit-latest.sh "https://repo.continuum.io/miniconda/explicit2-latest-MacOSX-x86_64.txt"
    >>> conda create --name p4env --file explicitenv2-latest-MacOSX-x86_64.txt
    >>> source activate p4env

.. _`sec:quickconda`:

Quick Installation
^^^^^^^^^^^^^^^^^^

Sequence of commands to get you to a working |PSIfour| on Linux. Installs
Miniconda into ``$HOME/miniconda`` and the |PSIfour| executable into the
main conda environment at ``$HOME/miniconda/bin/psi4``.

.. code-block:: bash

    # Linux or Mac: select between next two lines
    >>> curl -o Miniconda-latest.sh "https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh"
    >>> curl -o Miniconda-latest.sh "https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh"

.. code-block:: bash

    >>> bash
    >>> bash Miniconda-latest.sh -b -p $HOME/miniconda  # agrees to conda's license terms
    >>> echo "export PATH=$HOME/miniconda/bin:\$PATH" >> ~/.bashrc
    # log out, log back in so conda in path
    >>> conda update --yes --all
    >>> conda config --add channels http://conda.anaconda.org/psi4
    >>> conda install --yes psi4
    >>> psi4 "$(dirname $(which psi4))"/../share/psi4/samples/sapt1/test.in  # test installation. works b/c PSI_SCRATCH defaults to /tmp

That last command tested that ``psi4`` is in your path, and it's finding
all the libraries it needs. Now you need only specify a scratch directory
(see :ref:`sec:Scratch`) by replacing the placeholder in the following:

.. code-block:: bash

    >>> echo "export PSI_SCRATCH=/path/to/existing/writable/local-not-network/directory/for/scratch/files" >> ~/.bashrc
    # log out, log back in so variable takes effect

All done!

.. note:: Above commands use bash for installation and set up your environment for bash at runtime. To use csh at runtime, follow the on-screen directions at the end of the installation or consult step 7 below.

.. _`sec:slowconda`:

Detailed Installation of Miniconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

0. Sanity check. If you already have Miniconda or Anaconda, skip to step 5. The whole installation takes ~5 min; reading this page takes far longer.

1. Get ``bzip2``. You'll need this slightly exotic command so run ``which`` to test for availability, and install from ``yum``, source, *etc.* if unavailable. You'll also need an internet connection for downloading; computers behind a firewall or with restricted login domains are eligible. So long as you can ssh *into* the computer to an account with write permissions and can connect to the internet *from* the computer, all is well.

.. code-block:: bash

    # check
    >>> which bzip2
    /usr/bin/bzip2
    >>> curl -O "http://psicode.org/psi4manual/master/introduction.html"
    >>> ls -1
    introduction.html

2. Get Miniconda installer script. Either issue the command below or download from http://conda.pydata.org/miniconda.html by clicking on the appropriate link for your OS. If you already have or would prefer to use Anaconda rather than Miniconda, that's fine. Locate or install Anaconda, check that ``conda`` is in your path, and skip to step 4.

.. code-block:: bash

    >>> curl -O "http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh"
    # check
    >>> ls -1
    Miniconda-latest-Linux-x86_64.sh

3. Install Miniconda. Execute the script and answer its questions, particularly your choice of installation location. You may need to replace the filename below with the correct filename for the OS/version of installer you downloaded. Execute with ``bash`` regardless of ``csh``/``bash`` shell. If you're a ``bash`` user, it's convenient to agree to its offer to prepend ``conda`` commands to your :envvar:`PATH` in ``~/.bashrc``. If you're a ``csh``/``tcsh`` user, it's convenient to do the same by hand to your ``~/.tcshrc``: ``setenv PATH /path/to/miniconda/bin:${PATH}``. Further directions assume that the ``conda`` command is in your path; you may have to log out and log back in for ``which conda`` to return correctly.

.. code-block:: bash

    >>> bash Miniconda-latest-Linux-x86_64.sh
    # check
    >>> which conda
    /path/to/miniconda/bin/conda

4. Update conda. This updates the package manager itself.

.. code-block:: bash

    >>> conda update conda

.. _`sec:slowpsi4`:

Detailed Installation of |PSIfour|
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

5. Subscribe to |PSIfour|. Subscribe to our channel at http://anaconda.org/psi4 that contains the |PSIfour| package and several dependency packages. Make sure this shows up in your ``~/.condarc`` file.

.. code-block:: bash

    >>> conda config --add channels http://conda.anaconda.org/psi4
    # check
    >>> cat ~/.condarc
    channels:
      - http://conda.anaconda.org/psi4
      - defaults

6. Install |PSIfour|. You can install into the main conda environment so that whenever commands ``conda`` or (Ana/Miniconda's) ``python`` are available, then ``psi4`` is available, too. 

.. code-block:: bash

    >>> conda install psi4
    # check
    >>> which psi4
    /path/to/miniconda/bin/psi4

Or, you can install into a `conda environment <http://conda.pydata.org/docs/faq.html#environments>`_ that places |PSIfour| and its dependencies (including python) into a sandbox unaffected by any other software installed in Ana/Miniconda. This is recommended for developers to avoid interference between multiple versions (including github/conda versions) or to test python versions, *etc.*. If your main conda is not python=2.7, then |PSIfour| *must* be installed into a conda environment. In practical terms, installing into a conda environment means you can turn |PSIfours| availability on/off by switching conda environments without turning on/off the whole Ana/Miniconda installation. Below, |PSIfour| is installed into an environment called ``p4env``. Then the environment is activated, removing the main Ana/Miniconda ``bin`` and adding ``envs/p4env/bin`` to :envvar:`PATH`. The activate command only works in ``bash``, so ``csh``/``tcsh`` will need corresponding adjustments.

.. code-block:: bash

    >>> conda create -n p4env psi4
    >>> source activate p4env
    # check
    >>> which psi4
    /path/to/miniconda/envs/p4env/bin/psi4

The output for either of the installation commands above looks like the following. It checks what packages are needed, gets your approval for downloading them, fetches and installs them, prints out some useful information, and runs a |PSIfour| test case to check that all's well.

.. code-block:: bash

    >>> conda install psi4
    Using Anaconda Cloud api site https://api.anaconda.org
    Fetching package metadata: ......
    Solving package specifications: .........
    
    Package plan for installation in environment /theoryfs2/ds/cdsgroup/miniconda/envs/tpsi4:
    
    The following packages will be downloaded:
    
        package                    |            build
        ---------------------------|-----------------
        psi4-0.4.322               |    py27_g84b3aa1        44.4 MB  http://conda.anaconda.org/psi4/linux-64/
    
    The following NEW packages will be INSTALLED:
    
        psi4: 0.4.322-py27_g84b3aa1 http://conda.anaconda.org/psi4/linux-64/
    
    Proceed ([y]/n)? y
    
    Fetching packages ...
    psi4-0.4.322-p 100% |####################################################################################| Time: 0:00:08   5.77 MB/s
    Extracting packages ...
    [      COMPLETE      ]|#######################################################################################################| 100%
    Linking packages ...
    
    
      Thank you for installing psi4. Additional resources:
        Website: www.psicode.org
        Inputs:  /theoryfs2/ds/cdsgroup/miniconda/envs/tpsi4/share/psi4/samples
        Manual:  http://psicode.org/psi4manual/master/index.html
        GitHub:  https://github.com/psi4/psi4/wiki
        Binary:  https://anaconda.org/psi4
        Youtube: https://www.youtube.com/user/psitutorials
    
      For csh/tcsh command-line use, add to shell or ~/.tcshrc file:
        unsetenv PSIDATADIR
        setenv PATH /theoryfs2/ds/cdsgroup/miniconda/envs/tpsi4/bin:$PATH
        setenv PSI_SCRATCH /path/to/existing/writable/local-not-network/disk/for/scratch/files
    
      For sh/bash command-line use, add to shell or ~/.bashrc file:
        unset PSIDATADIR
        export PATH=/theoryfs2/ds/cdsgroup/miniconda/envs/tpsi4/bin:$PATH
        export PSI_SCRATCH=/path/to/existing/writable/local-not-network/disk/for/scratch/files
    
      Report problems at http://forum.psicode.org/t/report-conda-update-psi4-oddities-here/32
    
    
        Nuclear Repulsion Energy..........................................PASSED
        SAPT0 Eelst.......................................................PASSED
        SAPT0 Eexch.......................................................PASSED
        SAPT0 Eind........................................................PASSED
        SAPT0 Edisp.......................................................PASSED
        SAPT0 Etotal......................................................PASSED
    
    [      COMPLETE      ]|#######################################################################################################| 100%

7. Configure environment. Preceding steps have placed ``conda`` and ``psi4`` in your :envvar:`PATH`, either permanently through rc-files or temporarily in this terminal session. You can keep or undo these changes. For general psi4 use, you must enable the ``psi4`` executable to be found through any of:

  #. prepending to :envvar:`PATH` in shell, ``~/.bashrc``, ``~/.tcshrc``, or PBS ``cmd`` file
  #. activating the conda environment (p4env above) in shell, ``~/.bashrc``, or PBS ``cmd`` file
  #. supplying full path to executable (shell or PBS ``cmd`` file)

Similarly, the scratch directory (see :ref:`sec:Scratch`) must be specified through:

  #. defining :envvar:`PSI_SCRATCH` in shell, ``~/.bashrc``, ``~/.tcshrc``, or PBS ``cmd`` file

Suitable values for these variables have been printed to screen during installation (see last codeblock in step 6).

Useful Commands
^^^^^^^^^^^^^^^

* Update to latest |PSIfour| version ::

    >>> conda update psi4

* Install into a conda environment "p4env" instead of "root". Second command only works on bash; for csh/tsch, ``setenv PATH /path/to/miniconda/envs/p4env/bin:$PATH`` instead. This creates a sandbox with just |PSIfour| and python (loaded as dependency). ::

    >>> conda create -y -n p4env psi4
    >>> source activate p4env

* Install a particular |PSIfour| version ::

    >>> conda install psi4=0.1.12

* Uninstall |PSIfour| from current environment ::

    >>> conda remove psi4

Troubleshooting
^^^^^^^^^^^^^^^

* If the target computer doesn't have libc >= 2.7 (released c.2007; for reference, 2.10 is newer than 2.7; unlike most libraries, libc generally not available in multiple versions on a computer), the |PSIfour| conda package won't work. ::

    # unsuitable computer
    >>> ldd --version
    ldd (GNU libc) 2.5
    # suitable computer
    >>> ldd --version
    ldd (GNU libc) 2.17

* It is of greatest importance that the |PSIfour| executable be linked against conda libpython.so *not* against any system libpython.so. This is arranged by setting ``RPATH`` to seek libraries relative to executable (thanks, conda binary relocation routine!). The conda |PSIfour| executable is not vulnerable to interference from your ``LD_LIBRARY_PATH`` settings. Below shows a well-linked executable.

    * no libraries "not found"
    * fundamental libraries like libc, ld-linux, pthreads found system libraries to link against
    * libpython linked against conda python *not* system python
    * libm is linked against conda *or* system
    * blas, c++, and gcc libraries are absent because statically linked ::


    >>> conda install conda-build  # needed for next command
    >>> conda inspect linkages psi4
    python-2.7.9-2:
      libpython2.7.so.1.0 (lib/libpython2.7.so.1.0)
    system-5.8-1:
      libm.so.6 (lib/libm.so.6)
    system:
      libc.so.6 (/lib64/libc.so.6)
      libdl.so.2 (/lib64/libdl.so.2)
      libpthread.so.0 (/lib64/libpthread.so.0)
      librt.so.1 (/lib64/librt.so.1)
      libutil.so.1 (/lib64/libutil.so.1)
      linux-vdso.so.1 ()
    not found:


.. comment find out about the current environment.
.. comment pythonhome should be empty
.. comment pythonpath should be empty or set to non-interfering packages (*e.g.*, qcdb)
.. comment ld_library_path shouldn't contain anything with a libpython
.. comment >>> conda info -a

