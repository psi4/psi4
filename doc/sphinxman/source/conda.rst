.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2023 The Psi4 Developers.
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

.. include:: autodoc_abbr_options_c.rst

.. _`sec:conda`:

Conda Binary Distribution
=========================

.. warning:: As of v1.2rc1, new (conda build 3; updated compilers)
   conda packages are available for Linux but not Mac).
   Psi4conda installers are not ready for either platform.

|PSIfour| is available as a pre-compiled binary for Mac and Linux (and
Windows, through the Ubuntu shell) and native Windows architectures
through `Anaconda (formerly Continuum Analytics
<https://www.anaconda.com/products/individual>`_, the company that produces
`Anaconda Python <http://docs.continuum.io/anaconda/index.html>`_ (a
full-fledged scientific python environment with package manager `conda
<https://conda.io/projects/conda/en/latest/user-guide/getting-started.html>`_) and, more particularly, `Miniconda
<https://docs.conda.io/en/latest/miniconda.html>`_ (a lightweight Python
distribution with same package manager `conda
<https://conda.io/projects/conda/en/latest/user-guide/getting-started.html>`_). Some nice features for us:

* cross-platform

* no root, administrator, or sudo access required

* built with high-performance math libraries

* lightweight software stack (<100 MB w/o |PSIfour|; ~1 GB including |PSIfour|, numpy, and MKL)

* updated nightly so new features accessible

* standardizes python distribution so no need to find/install libpython packages

* add-ons (plugins, extra features requiring Fortran compiler, etc.) can be made available as conda packages

* develop |PSIfour| through plugins without a pre-existing development environment, see :ref:`sec:condaplugins`.

The |PSIfour| binary repository is at `Anaconda (formerly Binstar) <https://anaconda.org/psi4>`_.

For commands to get a default installation, go to :ref:`sec:psi4conda`
or the :psicode:`psicode downloads page <installs/latest/>` .
Users proficient with conda may prefer to consult :ref:`sec:condadetails`.
For more flexibility and a detailed explanation, go to
:ref:`sec:slowconda` and :ref:`sec:slowpsi4`.


.. _`faq:psicodedownload`:

How to install a Psi4 binary with the Psi4conda installer, download site
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Download one of the nine installers
<http://psicode.org/downloads.html>`_ (Linux/Mac/Windows; Py38/39/310).
``bash`` it. Follow the prompts and *do* make the adjustments to
:envvar:`PATH` and :envvar:`PSI_SCRATCH` that it suggests at the end. Test
with ``psi4 --test`` (green and yellow good; red bad). Done. Explicit commands at :ref:`sec:psi4conda`.


.. _`sec:psi4conda`:

How to install a Psi4 binary with the Psi4conda installer, command-line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sequence of commands to get you to a working |PSIfour| on Linux,
Mac, Windows (Ubuntu shell for Windows accepts Linux installers), or native Windows.
Installs Miniconda+Psi4+AddOns into ``$HOME/psi4conda`` and
the |PSIfour| executable into the main conda environment at
``$HOME/psi4conda/bin/psi4``.

.. code-block:: bash

    # Linux or WSL (Windows Subsystem for Linux)
    # py38|py39|py310 for alternate python versions
    >>> curl "http://vergil.chemistry.gatech.edu/psicode-download/Psi4conda-1.4rc1-py38-Linux-x86_64.sh" -o Psi4conda-latest-py38-Linux-x86_64.sh --keepalive-time 2
    >>> bash Psi4conda-latest-py38-Linux-x86_64.sh -b -p $HOME/psi4conda  # agrees to license terms
    >>> (bash) echo $'. $HOME/psi4conda/etc/profile.d/conda.sh\nconda activate' >> ~/.bashrc
    >>> (tcsh) echo "source $HOME/psi4conda/etc/profile.d/conda.csh\nconda activate" >> ~/.tcshrc
    # log out, log back in so conda and psi4 in path
    >>> psi4 --test

.. code-block:: bash

    # Mac
    # py38|py39|py310 for alternate python versions
    >>> curl -O "http://vergil.chemistry.gatech.edu/download/Psi4conda-latest-py35-MacOSX-x86_64.sh" --keepalive-time 2
    >>> curl "http://vergil.chemistry.gatech.edu/psicode-download/Psi4conda-1.4rc1-py38-MacOSX-x86_64.sh" -o Psi4conda-latest-py38-MacOSX-x86_64.sh --keepalive-time 2
    >>> bash Psi4conda-latest-py38-MacOSX-x86_64.sh -b -p $HOME/psi4conda  # agrees to license terms
    >>> (bash) echo $'. $HOME/psi4conda/etc/profile.d/conda.sh\nconda activate' >> ~/.bash_profile
    >>> (tcsh) echo "source $HOME/psi4conda/etc/profile.d/conda.csh\nconda activate" >> ~/.tcshrc
    # log out, log back in so conda and psi4 in path
    >>> psi4 --test

.. code-block:: bash

    # Windows
    # py38 only python version
    # download via button at https://psicode.netlify.app/installs/latest with "Windows", "Installer", and "Stable Release" selected
    >>> # install via GUI by double-clicking downloaded `.exe` file analogous to https://conda.io/projects/conda/en/latest/user-guide/install/windows.html
    >>> # -OR- install via following line
    >>>  start /wait "" Psi4conda-1.4rc1-py38-Windows-x86_64.exe /InstallationType=JustMe /RegisterPython=0 /S /D=%UserProfile%\psi4conda
    >>>  psi4 --test

That last command tested that ``psi4`` is in your path, and it's finding
all the libraries it needs. It works because :envvar:`PSI_SCRATCH`
defaults to ``/tmp``. Now you need only specify a permanent scratch
directory (see :ref:`sec:Scratch`) by replacing the placeholder in the
following:

.. code-block:: bash

    >>> echo "export PSI_SCRATCH=/path/to/existing/writable/local-not-network/directory/for/scratch/files" >> ~/.bashrc
    # log out, log back in so variable takes effect

All done!

Configuration for this set-up is summarized at :ref:`faq:runfrombinary`.

.. note:: |PSIfour| installs a Python distribution alongside, so you should choose an installer based on the Python version you *want*, irrespective of any Python version you *have*.


.. _`faq:psi4pkg`:

How to install a Psi4 binary into an Ana/Miniconda distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Conda command to install the latest nightly build of |PSIfour| + compiled
add-ons + runtime add-ons into an existing Anaconda or Miniconda
distribution.

.. code-block:: bash

    # Linux or Mac or Windows
    # substitute x.x by 3.6|3.7|3.8|3.9 for alternate python versions
    # remove `-c psi4/label/dev` to get stable releases instead of nightly builds
    >>> conda create -n p4env python=x.x psi4 -c psi4/label/dev

Activate environment and make the adjustments to :envvar:`PATH` and
:envvar:`PSI_SCRATCH` that it suggests at the end. Test with ``psi4
--test``. Configuration for this set-up is summarized at
:ref:`faq:runfrombinary`.

**Details:**

* It is advised to place |PSIfour| into a conda
  environment where its libraries can't interfere with other programs
  rather than the main
  Anaconda or Miniconda environment. Hence the creation of the environment
  above, but the environment name (:samp:`{p4env}` above) can be
  substituted.

* The ``psi4-rt`` package can be added to the package list to get the
  QC runtime add-ons; could say any combination of ``v2rdm_casscf snsmp2
  resp`` etc. instead of ``psi4-rt``.
  As of |PSIfour| v1.7, the ``psi4-rt`` package is being slowly retired
  due to more optional dependencies being on conda-forge. Similar collections
  of dependencies can be obtained from environment spec files like :source:`devtools/conda-envs` .

* Grab a Miniconda through one of the below, selecting OS.

  >>> curl -O "https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Linux-x86_64.sh"
  >>> curl -O "https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-MacOSX-x86_64.sh"
  >>> curl -O "https://repo.anaconda.com/miniconda/Miniconda3-py39_4.9.2-Windows-x86_64.exe"


.. _`faq:updatepsi4`:

How to update a Psi4 binary
^^^^^^^^^^^^^^^^^^^^^^^^^^^

A. Conda command to update an existing |PSIfour| conda installation to the
newest stable release (roughly annually). It's often a better idea to create
a new environment rather than updating the old one.

.. code-block:: bash

    >>> # Linux/MacOS
    >>> conda update psi4 -c psi4
    >>> # Windows
    >>> conda update psi4 -c psi4 -c conda-forge

    # if psi4 channel in defaults (true for Psi4conda installers)
    >>> conda update psi4

B. Conda command to update an existing |PSIfour| conda installation to the
latest development head (roughly nightly).

.. code-block:: bash

    >>> # Linux/MacOS
    >>> conda update psi4 -c psi4/label/dev
    >>> # Windows
    >>> conda update psi4 -c psi4/label/dev -c conda-forge

C. Conda command to install a very specific package, including version,
build string, and subchannel. The final `-c psi4` represents any
additional channels or subchannels needed to locate all dependencies.

.. code-block:: bash

    >>> conda install psi4=1.2a1.dev249+623ad64=py36_sse41_0 -c psi4/label/subchannel -c psi4


.. _`faq:psi4deps`:

How to use conda to compile Psi4 faster and easier
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: console

    # Linux # c. v1.2rc1  ###or Mac or Windows
    # substitute x.x by 3.5|3.6|3.7 for alternate python versions
    >>> conda create -n p4dev python=x.x psi4-dev -c psi4/label/dev -c psi4
    >>> conda activate p4dev

    >>> cd {top-level-psi4-dir}
    >>> psi4-path-advisor --help
    usage: psi4-path-advisor [-h] [--psi4-compile] [--disable-addons]
                             [--disable-mkl] [--intel | --intel-multiarch | --gcc]
    
    Build and Run path advisor for Psi4
    
    optional arguments:
      -h, --help         show this help message and exit
      --psi4-compile     (Command Default) Generates a minimal CMake command for building Psi4 against
                             this psi4-dev conda metapackage.
                         >>> git clone https://github.com/psi4/psi4.git
                         >>> cd {top-level-psi4-dir}
                         >>> conda create -n p4dev python={3.6} psi4-dev [-c psi4/label/dev] -c psi4
                         >>> conda activate p4dev
                         >>> psi4-path-advisor
                         # execute or adapt `cmake` commands above; DepsCache handles python & addons;
                         #   DepsMKLCache handles math; further psi4-path-advisor options handle compilers.
                         >>> cd objdir && make -j`getconf _NPROCESSORS_ONLN`
                         >>> make install
      --disable-addons   Disengage building against the psi4-dev-provided _optional_ link-time Add-Ons like CheMPS2.
      --disable-mkl      Disengage building against the psi4-dev-provided MKL libraries (`libmkl_rt`).
      --intel            Engage self-provided icc/icpc/ifort compilers backed by conda's psi4-dev-provided gcc/g++.
      --intel-multiarch  Engage self-provided icc/icpc/ifort compilers backed by conda's psi4-dev-provided gcc/g++ PLUS compile for multiple architectures (useful for cluster deployments).
      --gcc              Engage conda's psi4-dev-provided gcc/g++/gfortran compilers.

    # execute or adapt `cmake` commands above; DepsCache handles python & addons;
    #   DepsMKLCache handles math; further psi4-path-advisor options handle compilers.
    >>> `psi4-path-advisor [your args]` -Dany_addl_cmake_vals=ON
    >>> cd objdir && make -j`getconf _NPROCESSORS_ONLN`
    >>> make install

Same for Linux/Mac/WSL. Substitute desired python version: 3.6, 3.7, 3.8, 3.9. Fine
to choose your own env name. Include ``-c psi4/label/dev`` to get dependencies to
build current master, as opposed to latest release.
Activate environment, ``conda activate
p4dev``.  Go to where you've cloned psi4. Execute ``psi4-path-advisor``.
It gives you a basic cmake command covering python, sphinx, link-time qc
addons, and run-time qc addons. There's a help menu -h that gives more
info. There's other options that will also pre-configure compilers. For
example, at GaTech ``psi4-path-advisor --intel`` works. On Macs with
XCode, ``psi4-path-advisor --clang`` works. Just read the help. For users
who want a minimal build, there's a ``--disable-addons``, but it is generally not
encouraged. It gives you a fully
functional cmake command, but those are just setting up CMake cache
|w---w| like the plugins you can always add your own CMake variables to
the command.

For run-time, you may also wish to install the optional runtime add-ons (*e.g.*, adcc)

.. code-block:: bash

    >>> conda install psi4-rt


.. _`sec:condadetails`:

What do the conda packages psi4 & psi4-dev and the installer psi4conda contain
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``psi4`` - has full-featured psi4 itself and necessarily all the link-time qc
addons (e.g., chemps2). It has python, pytest, numpy, and a few more python
modules for specialized functions. Of gcc-ness, it has minimal, run-time
libraries (*e.g.*, libgcc-ng) not compilers.
It doesn't have the run-time qc addons ``psi4-rt`` (*e.g.*, snsmp2) or build tools (*e.g.*, g++, sphinx, cmake).

``psi4-dev`` - does not have psi4 itself or the run-time addons ``psi4-rt`` or numpy (though fine to install them
alongside). Does have all the link-time addons. Does have
cmake and sphinx (and python). Of gcc-ness, has full packages, that is,
compilers as well as runtime packages.

Psi4conda installer - has full-featured ``psi4`` itself, all link-time qc addons, all
run-time qc addons, and minimal gcc runtime libraries. Developers should additionally install ``psi4-dev`` for build tools.

The :ref:`sec:psi4conda` uses a `conda constructor
<https://github.com/conda/constructor>`_ to package up Miniconda,
the |PSIfour| conda package, the |PSIfour| add-on conda packages, dependencies
thereof (possibly from particular channels), and the psi4 channel
as a default. This is very convenient for novice users and robust
against differing channel settings in ``~/.condarc``. But proficient
conda users may prefer to treat ``psi4`` as a normal conda package and
not have another large Miniconda installation (including the hefty MKL)
lying around just for |PSIfour|. Installing just the ``psi4`` package
itself will get you |PSIfour|, whatever add-ons require linking in to
|PSIfour| (*e.g.*, CheMPS2 and PCMSolver), and the correct versions of
packages. However, just the ``psi4`` package won't get you add-ons that
don't need linking (*e.g.*, adcc and v2rdm_casscf).

.. Conda Proficients
.. ^^^^^^^^^^^^^^^^^
..
.. or dependencies
.. from the "right" channels, which can be important for issues of fPIC
.. and libc++ vs. libstdc++. So ``conda create -c psi4 -n p4env psi4 dftd3
.. v2rdm_casscf`` *should* be equivalent to running the psi4conda installer,
.. but I wouldn't count on it. Instead, an `explicit environment spec
.. <http://conda.pydata.org/docs/using/envs.html#build-identical-conda-environments-with-urls>`_
.. will be available for download.
..
.. .. code-block:: bash
..
..     # Linux
..     >>> curl -o explicit-latest.sh "https://repo.continuum.io/miniconda/explicit2-latest-Linux-x86_64.txt"
..     >>> conda create --name p4env --file explicitenv2-latest-Linux-x86_64.txt
..     >>> conda activate p4env
..
.. .. code-block:: bash
..
..     # Mac
..     >>> curl -o explicit-latest.sh "https://repo.continuum.io/miniconda/explicit2-latest-MacOSX-x86_64.txt"
..     >>> conda create --name p4env --file explicitenv2-latest-MacOSX-x86_64.txt
..     >>> conda activate p4env

.. _`sec:quickconda`:

Quick Installation
^^^^^^^^^^^^^^^^^^

Sequence of commands to get you to a working |PSIfour|. Installs
Miniconda into ``$HOME/miniconda`` and the |PSIfour| executable into the
main conda environment at ``$HOME/miniconda/bin/psi4``.

.. code-block:: bash

    # Linux or Mac, Py2 or Py3 for main environment (immaterial to Py for Psi4): select between four lines
    # Windows: in Ubuntu shell, select either Linux line
    >>> curl -o Miniconda-latest.sh "https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh"
    >>> curl -o Miniconda-latest.sh "https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    >>> curl -o Miniconda-latest.sh "https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh"
    >>> curl -o Miniconda-latest.sh "https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh"

.. code-block:: bash

    >>> bash
    >>> bash Miniconda-latest.sh -b -p $HOME/miniconda  # agrees to conda's license terms
    >>> echo "export PATH=$HOME/miniconda/bin:\$PATH" >> ~/.bashrc  # Mac: use ~/.bash_profile
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

https://docs.conda.io/en/latest/miniconda.html

.. 0. Sanity check. If you already have Miniconda or Anaconda, skip to step 5. The whole installation takes ~5 min; reading this page takes far longer.
..
.. 1. Get ``bzip2``. You'll need this slightly exotic command so run ``which`` to test for availability, and install from ``yum``, source, *etc.* if unavailable. You'll also need an internet connection for downloading; computers behind a firewall or with restricted login domains are eligible. So long as you can ssh *into* the computer to an account with write permissions and can connect to the internet *from* the computer, all is well.
..
.. .. code-block:: bash
..
..     # check
..     >>> which bzip2
..     /usr/bin/bzip2
..     >>> curl -O "http://psicode.org/psi4manual/master/introduction.html"
..     >>> ls -1
..     introduction.html
..
.. 2. Get Miniconda installer script. Either issue the command below or download from http://conda.pydata.org/miniconda.html by clicking on the appropriate link for your OS. If you already have or would prefer to use Anaconda rather than Miniconda, that's fine. Locate or install Anaconda, check that ``conda`` is in your path, and skip to step 4.
..
.. .. code-block:: bash
..
..     >>> curl -O "http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh"
..     # check
..     >>> ls -1
..     Miniconda-latest-Linux-x86_64.sh
..
.. 3. Install Miniconda. Execute the script and answer its questions, particularly your choice of installation location. You may need to replace the filename below with the correct filename for the OS/version of installer you downloaded. Execute with ``bash`` regardless of ``csh``/``bash`` shell. If you're a ``bash`` user, it's convenient to agree to its offer to prepend ``conda`` commands to your :envvar:`PATH` in ``~/.bashrc``. If you're a ``csh``/``tcsh`` user, it's convenient to do the same by hand to your ``~/.tcshrc``: ``setenv PATH /path/to/miniconda/bin:${PATH}``. Further directions assume that the ``conda`` command is in your path; you may have to log out and log back in for ``which conda`` to return correctly.
..
.. .. code-block:: bash
..
..     >>> bash Miniconda-latest-Linux-x86_64.sh
..     # check
..     >>> which conda
..     /path/to/miniconda/bin/conda
..
.. 4. Update conda. This updates the package manager itself.
..
.. .. code-block:: bash
..
..     >>> conda update conda

.. _`sec:slowpsi4`:

Detailed Installation of |PSIfour|
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

5. Subscribe to |PSIfour|. Subscribe to our channel at https://anaconda.org/psi4 that contains the |PSIfour| package and several dependency packages. Make sure this shows up in your ``~/.condarc`` file.

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

Or, you can install into a `conda environment <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#>`_ that places |PSIfour| and its dependencies (including python) into a sandbox unaffected by any other software installed in Ana/Miniconda. This is recommended for developers to avoid interference between multiple versions (including github/conda versions) or to test python versions, *etc.*. In practical terms, installing into a conda environment means you can turn |PSIfours| availability on/off by switching conda environments without turning on/off the whole Ana/Miniconda installation. Below, |PSIfour| is installed into an environment called ``p4env``. Then the environment is activated, removing the main Ana/Miniconda ``bin`` and adding ``envs/p4env/bin`` to :envvar:`PATH`. The ``conda activate`` command (conda >=4.4; December 2017) works in all shells, but if you're using old ``source activate`` that only works for ``bash``; adjust as needed for ``csh``/``tcsh``.

.. code-block:: bash

    >>> conda create -n p4env psi4
    >>> conda activate p4env
    # check
    >>> which psi4
    /path/to/miniconda/envs/p4env/bin/psi4

.. The output for either of the installation commands above looks like the following. It checks what packages are needed, gets your approval for downloading them, fetches and installs them, prints out some useful information, and runs a |PSIfour| test case to check that all's well.
..
.. .. code-block:: bash
..
..     >>> conda install psi4
..     Using Anaconda Cloud api site https://api.anaconda.org
..     Fetching package metadata: ......
..     Solving package specifications: .........
..
..     Package plan for installation in environment /theoryfs2/ds/cdsgroup/miniconda/envs/tpsi4:
..
..     The following packages will be downloaded:
..
..         package                    |            build
..         ---------------------------|-----------------
..         psi4-0.4.322               |    py27_g84b3aa1        44.4 MB  http://conda.anaconda.org/psi4/linux-64/
..
..     The following NEW packages will be INSTALLED:
..
..         psi4: 0.4.322-py27_g84b3aa1 http://conda.anaconda.org/psi4/linux-64/
..
..     Proceed ([y]/n)? y
..
..     Fetching packages ...
..     psi4-0.4.322-p 100% |####################################################################################| Time: 0:00:08   5.77 MB/s
..     Extracting packages ...
..     [      COMPLETE      ]|#######################################################################################################| 100%
..     Linking packages ...
..
..
..       Thank you for installing psi4. Additional resources:
..         Website: www.psicode.org
..         Inputs:  /theoryfs2/ds/cdsgroup/miniconda/envs/tpsi4/share/psi4/samples
..         Manual:  http://psicode.org/psi4manual/master/index.html
..         GitHub:  https://github.com/psi4/psi4/wiki
..         Binary:  https://anaconda.org/psi4
..         Youtube: https://www.youtube.com/user/psitutorials
..
..       For csh/tcsh command-line use, add to shell or ~/.tcshrc file:
..         unsetenv PSIDATADIR
..         setenv PATH /theoryfs2/ds/cdsgroup/miniconda/envs/tpsi4/bin:$PATH
..         setenv PSI_SCRATCH /path/to/existing/writable/local-not-network/disk/for/scratch/files
..
..       For sh/bash command-line use, add to shell or ~/.bashrc file:
..         unset PSIDATADIR
..         export PATH=/theoryfs2/ds/cdsgroup/miniconda/envs/tpsi4/bin:$PATH
..         export PSI_SCRATCH=/path/to/existing/writable/local-not-network/disk/for/scratch/files
..
..       Report problems at http://forum.psicode.org/t/report-conda-update-psi4-oddities-here/32
..
..
..         Nuclear Repulsion Energy..........................................PASSED
..         SAPT0 Eelst.......................................................PASSED
..         SAPT0 Eexch.......................................................PASSED
..         SAPT0 Eind........................................................PASSED
..         SAPT0 Edisp.......................................................PASSED
..         SAPT0 Etotal......................................................PASSED
..
..     [      COMPLETE      ]|#######################################################################################################| 100%

7. Configure environment. Preceding steps have placed ``conda`` and ``psi4`` in your :envvar:`PATH`, either permanently through rc-files or temporarily in this terminal session. You can keep or undo these changes. For general psi4 use, you must enable the ``psi4`` executable to be found through any of:

  #. prepending to :envvar:`PATH` in shell, ``~/.bashrc``, ``~/.tcshrc``, or PBS ``cmd`` file
  #. activating the conda environment (p4env above) in shell, ``~/.bashrc``, or PBS ``cmd`` file
  #. supplying full path to executable (shell or PBS ``cmd`` file)

Similarly, the scratch directory (see :ref:`sec:Scratch`) must be specified through:

  #. defining :envvar:`PSI_SCRATCH` in shell, ``~/.bashrc``, ``~/.tcshrc``, or PBS ``cmd`` file

.. Suitable values for these variables have been printed to screen during installation (see last codeblock in step 6).

Useful Commands
^^^^^^^^^^^^^^^

* (A) Initially install |PSIfour| stable release

.. code-block:: console

   # equivalent
   >>> conda install psi4 -c psi4
   >>> conda install psi4 --channel psi4

* (B) Initially install |PSIfour| stable release with non-current python

.. code-block:: console

   >>> conda install psi4 python=3.8 -c psi4

* (C) Update to latest |PSIfour| stable release

.. code-block:: console

    >>> conda update psi4 -c psi4

* (D) Initially install stable release into a conda environment "p4env" instead of "root". This creates a sandbox with |PSIfour| and python (loaded as dependency).

.. code-block:: console

    >>> conda create -y -n p4env psi4 -c psi4
    >>> conda activate p4env

* (E) Install a particular |PSIfour| version

.. code-block:: console

    >>> conda install psi4=1.4 -c psi4

* (F) Uninstall |PSIfour| from current environment

.. code-block:: console

    >>> conda remove psi4

* (G) Initially install |PSIfour| nightly build

.. code-block:: console

   # equivalent
   >>> conda install psi4 -c psi4/label/dev
   >>> conda install psi4 --channel psi4/label/dev

* (H) Initially install |PSIfour| nightly build with non-current python

.. code-block:: console

   >>> conda install psi4 python=3.8 -c psi4/label/dev

* (I) Update to latest |PSIfour| nightly build

.. code-block:: console

    >>> conda update psi4 -c psi4/label/dev

* (J) Initially install nightly build into a conda environment "p4env" instead of "root". This creates a sandbox with |PSIfour| and python (loaded as dependency).

.. code-block:: console

    >>> conda create -y -n p4env psi4 -c psi4/label/dev
    >>> conda activate p4env

* (K) Install a particular |PSIfour| version

.. code-block:: console

    >>> conda install psi4=1.4 -c psi4/label/dev

.. Troubleshooting
.. ^^^^^^^^^^^^^^^
..
.. * If the target computer doesn't have libc >= 2.7 (released c.2007; for reference, 2.10 is newer than 2.7; unlike most libraries, libc generally not available in multiple versions on a computer), the |PSIfour| conda package won't work. ::
..
..     # unsuitable computer
..     >>> ldd --version
..     ldd (GNU libc) 2.5
..     # suitable computer
..     >>> ldd --version
..     ldd (GNU libc) 2.17
..
.. * It is of greatest importance that the |PSIfour| executable be linked against conda libpython.so *not* against any system libpython.so. This is arranged by setting ``RPATH`` to seek libraries relative to executable (thanks, conda binary relocation routine!). The conda |PSIfour| executable is not vulnerable to interference from your ``LD_LIBRARY_PATH`` settings. Below shows a well-linked executable.
..
..     * no libraries "not found"
..     * fundamental libraries like libc, ld-linux, pthreads found system libraries to link against
..     * libpython linked against conda python *not* system python
..     * libm is linked against conda *or* system
..     * blas, c++, and gcc libraries are absent because statically linked
..
.. .. code-block:: console
..
..       >>> conda install conda-build  # needed for next command
..       >>> conda inspect linkages psi4
..       python-2.7.9-2:
..         libpython2.7.so.1.0 (lib/libpython2.7.so.1.0)
..       system-5.8-1:
..         libm.so.6 (lib/libm.so.6)
..       system:
..         libc.so.6 (/lib64/libc.so.6)
..         libdl.so.2 (/lib64/libdl.so.2)
..         libpthread.so.0 (/lib64/libpthread.so.0)
..         librt.so.1 (/lib64/librt.so.1)
..         libutil.so.1 (/lib64/libutil.so.1)
..         linux-vdso.so.1 ()
..       not found:


.. comment find out about the current environment.
.. comment pythonhome should be empty
.. comment pythonpath should be empty or set to non-interfering packages (*e.g.*, qcdb)
.. comment ld_library_path shouldn't contain anything with a libpython
.. comment >>> conda info -a

