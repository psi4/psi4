
.. include:: autodoc_abbr_options_c.rst

.. _`sec:conda`:

Binary Distribution
===================

|PSIfour| is available as a pre-compiled binary for Linux architectures
(Mac coming soon) through `Continuum Analytics
<https://store.continuum.io/cshop/anaconda/>`_, the company that produces
`Anaconda Python <http://docs.continuum.io/anaconda/index.html>`_ (a
full-fledged scientific python environment with package manager `conda
<http://conda.pydata.org/index.html>`_) and, more particularly, `Miniconda
<http://conda.pydata.org/miniconda.html>`_ (a lightweight python
distribution with package manger ``conda``).

* no root, administrator, or sudo access required

* built with high-performance math (Intel MKL) libraries

* lightweight software stack (<100 MB w/o |PSIfour|; <500 MB including |PSIfour|)

* updated nightly so new features accessible

* standardizes python distribution so no need to find/install libpython packages

.. comment* easy to install programs interfaced to |PSIfour|

The |PSIfour| binary repository is at `Binstar <https://binstar.org/psi4>`_.

Quick Installation
^^^^^^^^^^^^^^^^^^

Sequence of commands to get you to a working psi4 on Linux. Installs
Miniconda into ``$HOME/miniconda`` and |PSIfour| into a conda environment
named "p4env" at ``$HOME/miniconda/envs/p4env``; feel free to change these
locations.

.. code-block:: bash

    >>> bash
    >>> curl -O "http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh"
    >>> bash Miniconda-latest-Linux-x86_64.sh -b -p $HOME/miniconda  # agrees to conda's license terms
    >>> echo "export PATH=$HOME/miniconda/bin:\$PATH" >> ~/.bashrc
    # log out, log back in so conda in path
    >>> conda update --yes --all
    >>> conda config --add channels http://conda.binstar.org/psi4
    >>> conda install --yes psi4
    >>> psi4 "$(dirname $(which psi4))"/../share/psi/samples/scf1/input.dat -o stdout  # works b/c PSI_SCRATCH defaults to /tmp

.. code-block:: bash

    >>> echo "export PSI_SCRATCH=/path/to/existing/writable/local-not-network/disk/for/scratch/files" >> ~/.bashrc
    # log out, log back in so variable takes effect

.. comment# works b/c scratch defaults to /tmp    
    >>> conda create --yes --name p4env psi4
    >>> source activate p4env
    >>> export PSI_SCRATCH=/path/to/existing/writable/local-not-network/disk/for/scratch/files

All done! For general psi4 use, you must enable the ``psi4`` executable to be found through any of:

  #. prepending to :envvar:`PATH` in shell, ``~/.bashrc``, ``~/.tcshrc``, or PBS ``cmd`` file
  #. activating the conda environment (p4env above) in shell, ``~/.bashrc``, or PBS ``cmd`` file
  #. supplying full path to executable (shell or PBC ``cmd`` file)

You must also specify a scratch directory.

.. comment  #. adding to :envvar:`PATH` in shell, (b) adding to :envvar:`PATH` in ``~/.bashrc`` or ``~/.tcshrc``, (c) always giving full path to executable, (d) activating conda environment in shell or rc-file, (e)



Detailed Installation of Miniconda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

0. If you alreay have Miniconda or Anaconda, skip to step 5. The whole installation takes ~5 min; reading this page takes far longer.

1. You'll need the slightly exotic command ``bzip2`` available to run the below. Run ``which`` to test for availability, and install from ``yum``, source, etc. if unavailable. You'll also need an internet connection for downloading; computers behind a firewall or with restricted login domains are eligible. So long as you can ssh *into* the computer to an account with write permissions and can connect to the internet *from* the computer, all is well.

.. code-block:: bash

    # check
    >>> which bzip2
    /usr/bin/bzip2
    >>> curl -O "http://sirius.chem.vt.edu/psi4manual/master/introduction.html"
    >>> ls -1
    introduction.html

2. Get the Miniconda installer script. Either issue the command below or download from http://conda.pydata.org/miniconda.html by clicking on the appropriate link for your OS. If you already have or would prefer to use Anaconda rather than Miniconda, that's fine. Locate or install Anaconda and skip the next step.

.. code-block:: bash

    >>> curl -O "http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh"
    # check
    >>> ls -1
    Miniconda-latest-Linux-x86_64.sh

3. Install Miniconda by executing the script and answering its questions, particularly your choice of installation location. You may need to replace the filename below with the correct filename for the OS/version of installer you downloaded. If you're a ``bash`` user, it's convenient to agree to its offer to prepend ``conda`` commands to your :envvar:`PATH` in ``~/.bashrc``. If you're a ``csh``/``tcsh`` user, it's convenient to do the same by hand to your ``~/.tcshrc``: ``setenv PATH /path/to/miniconda/bin:${PATH}``. Further directions assume that the ``conda`` command is in your path; you may have to log out and log back in for ``which conda`` to return correctly. ::

    >>> bash Miniconda-latest-Linux-x86_64.sh
    # check
    >>> which conda
    /path/to/miniconda/bin/conda

4. Update the package manager itself. ::

    >>> conda update conda

Detailed Installation of |PSIfour|
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

5. Subscribe to the |PSIfour| channel at http://binstar.org/psi4. Make sure this shows up in your ``~/.condarc`` file. ::

    >>> conda config --add channels http://conda.binstar.org/psi4
    # check
    >>> cat ~/.condarc
    channels:
      - http://conda.binstar.org/psi4
      - defaults

6. Make a |PSIfour| sandbox (here, called ``p4env``) with just |PSIfour| and python (loaded as dependency).

.. code-block:: bash

    >>> conda install psi4
    Fetching package metadata: ......
    Solving package specifications: .
    Package plan for installation in environment /path/to/miniconda:
    
    The following packages will be downloaded:
    
        package                    |            build
        ---------------------------|-----------------
        psi4-0.1.12                |    py27_g60a9afd        31.4 MB
    
    The following NEW packages will be INSTALLED:
    
        psi4: 0.1.12-py27_g60a9afd
    
    Fetching packages ...
    psi4-0.1.12-py 100% |################################################################################| Time: 0:00:01  32.49 MB/s
    Extracting packages ...
    [      COMPLETE      ] |##################################################################################################| 100%
    Linking packages ...
    
    
      Thank you for installing psi4. Additional resources:
        Website: www.psicode.org
        Inputs:  /path/to/miniconda/share/psi/samples
        Manual:  bit.ly/psi4manual
        GitHub:  https://github.com/psi4/psi4public/wiki
        Binstar: https://binstar.org/psi4
        Runtime Environment Diagnostic: /path/to/miniconda/share/psi/scripts/setenv.py
    
      For csh/tcsh command-line use, add to shell or ~/.tcshrc file:
        setenv PATH /path/to/miniconda/bin:$PATH
        setenv PSI_SCRATCH /path/to/existing/writable/local-not-network/disk/for/scratch/files
    
      For sh/bash command-line use, add to shell or ~/.bashrc file:
        export PATH=/path/to/miniconda/bin:$PATH
        export PSI_SCRATCH=/path/to/existing/writable/local-not-network/disk/for/scratch/files
    
    Nuclear Repulsion Energy..........................................PASSED SAPT0 Eelst.......................................................PASSED SAPT0 Eexch.......................................................PASSED SAPT0 Eind........................................................PASSED SAPT0 Edisp.......................................................PASSED SAPT0 Etotal......................................................PASSED
    
    [      COMPLETE      ] |##################################################################################################| 100%




Useful Commands
^^^^^^^^^^^^^^^

* Update to latest |PSIfour| version ::

    >>> conda update psi4

* Install into a conda environment "p4env" instead of "root". Second command only works on bash; for csh/tsch, ``setenv PATH /path/to/miniconda/envs/p4env/bin:$PATH`` instead. This creates a sandbox with just |PSIfour| and python (loaded as dependency). ::

    >>> conda create -n p4env psi4
    >>> source activate p4env

* Install a particular |PSIfour| version ::

    >>> conda install psi4=0.1.12

Troubleshooting
^^^^^^^^^^^^^^^

* If the target computer doesn't have libc >= 2.7 (released c.2007; for reference, 2.10 is newer than 2.7; unlike most libraries, libc generally not available in multiple versions on a computer), the |PSIfour| conda package won’t work. ::

    # unsuitable computer
    >>> ldd --version
    ldd (GNU libc) 2.5
    # suitable computer
    >>> ldd --version
    ldd (GNU libc) 2.17

* It is of greatest importance that the |PSIfour| executable be linked against conda libpython.so *not* against any system libpython.so. This is arranged by setting ``RPATH`` to seek libraries relative to executable (thanks, conda binary relocation routine!). The conda |PSIfour| executable should not be vulnerable to interference from your ``LD_LIBRARY_PATH`` settings. Below shows a well-linked executable.

    * no libraries "not found"
    * fundamental libraries like libc, ld-linux, pthreads found system libraries to link against
    * libpython linked against conda python not system python
    * libm is linked against conda *or* system

.. code-block:: bash

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


