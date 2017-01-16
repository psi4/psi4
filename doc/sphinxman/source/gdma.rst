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

.. index:: 
   DMA
   GDMA
   Distributed Multipole Analysis

.. _`sec:gdma`:

Interface to GDMA Distributed Multipole Analysis by A. J. Stone |w---w| :py:func:`~psi4.gdma`
=============================================================================================

.. codeauthor:: Anthony J. Stone, Andrew C. Simmonett
.. sectionauthor:: Andrew C. Simmonett

*Module:* :ref:`Keywords <apdx:gdma>`, :ref:`PSI Variables <apdx:gdma_psivar>`, :source:`GDMA_INTERFACE<psi4/src/psi4/gdma_interface>`

.. image:: https://img.shields.io/badge/home-gdma-5077AB.svg
   :target: https://github.com/psi4/gdma

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: http://www-stone.ch.cam.ac.uk/documentation/gdma/manual.pdf

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/gdma/badges/version.svg
     :target: https://anaconda.org/psi4/gdma

* GDMA is available as a conda package for Linux and macOS.

* If using the |PSIfour| binary, gdma has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  gdma can be obtained through ``conda install gdma``.
  Then enable it as a feature with :makevar:`ENABLE_gdma`,
  hint its location with :makevar:`CMAKE_PREFIX_PATH`,
  and rebuild |PSIfour| to detect gdma and activate dependent code.

* To remove a conda installation, ``conda remove gdma``.

**Source**

* .. image:: https://img.shields.io/github/tag/psi4/gdma.svg?maxAge=2592000
     :target: https://github.com/psi4/gdma

* If using |PSIfour| built from source and you want gdma built from
  from source also,
  enable it as a feature with :makevar:`ENABLE_gdma`,
  and let the build system fetch and build it and activate dependent code.

Input
~~~~~

The distributed multipole analysis (DMA) technique, developed by Anthony J.
Stone and implemented by him into the `GDMA package
<http://www-stone.ch.cam.ac.uk/programs.html>`_, is available in |PSIfour|.
The current implementation simply embeds Stone's GDMA code into the main
executable, and generates the appropriate input files automatically.  The
program takes as input a data file, and a Gaussian formatted checkpoint (see
Section :ref:`FCHK <sec:fchk>`) file.  The simplest usage of the GDMA code is
demonstrated below, along with a listing of the options supported; these
options correspond to the options described in the 
:download:`GDMA manual <gdma-2.2.06.pdf>`.

If more advanced usage is desired, which is not is permitted by the options
listed below, the user may provide their own data file containing keywords to
control the GDMA code.  Simply place the data file in the directory |PSIfour|
is called from, and provide the file name as the datafile argument to the
:py:func:`~psi4.gdma` routine.  For example, if GDMA data file is called
*control.dma*, the GDMA code is called as follows::

    grad, wfn = gradient('mp2', return_wfn=True)
    gdma(wfn, datafile='control.dma')

An FCHK file will be generated for the GDMA code to read; this file will have
the prefix given by |globals__writer_file_label| (if set), or else by the name
of the output file plus the name of the current molecule, and the suffix will
be '.fchk'.  This FCHK file name should be passed to the 'File' keyword in the
DGMA data file, to ensure that the GDMA code reads the correct wavefunction
information.

After running, two matrices of results can be accessed::

    dma_results = get_array_variable('DMA DISTRIBUTED MULTIPOLES')
    tot_results = get_array_variable('DMA TOTAL MULTIPOLES')

The first contains distributed multipoles, in units given by
|gdma__gdma_multipole_units|, with the row index corresponding to the site and
the column index referencing the multipole component.  Both indices are zero
based, and the :math:`Q^l_m` components of the multipoles are ordered as
:math:`Q^0_0, Q^1_0, Q^1_{1c}, Q^1_{1s}, Q^2_0, Q^2_{1c}, Q^2_{1s}, Q^2_{2c},
Q^2_{2s}, \ldots`  The second matrix returned has a single row, whose columns
are the total multipoles, translated to |gdma__gdma_origin|, and summed.


.. autofunction:: psi4.gdma(wfn)

Options
~~~~~~~

.. include:: autodir_options_c/gdma__gdma_limit.rst
.. include:: autodir_options_c/gdma__gdma_origin.rst
.. include:: autodir_options_c/gdma__gdma_multipole_units.rst
.. include:: autodir_options_c/gdma__gdma_radius.rst
.. include:: autodir_options_c/gdma__gdma_switch.rst

.. _`cmake:gdma`:

How to configure gdma for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, GDMA is a library that provides additional
  quantum chemical capabilities (multipole analysis).

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) gdma

* Upstream Dependencies |w---w| gdma |dr| Fortran

**CMake Variables**

* :makevar:`ENABLE_gdma` |w---w| CMake variable toggling whether Psi4 builds with gdma
* :makevar:`CMAKE_PREFIX_PATH` |w---w| CMake list variable to specify where pre-built dependencies can be found. For gdma, set to an installation directory containing ``include/GDMA/GDMA_MANGLE.h``
* :makevar:`gdma_DIR` |w---w| CMake variable to specify where pre-built gdma can be found. Set to installation directory containing ``share/cmake/gdma/gdmaConfig.cmake``
* :makevar:`CMAKE_DISABLE_FIND_PACKAGE_gdma` |w---w| CMake variable to force internal build of gdma instead of detecting pre-built

**Examples**

A. Build bundled

  .. code-block:: bash

    >>> cmake -DENABLE_gdma=ON

B. Build *without* gdma

  .. code-block:: bash

    >>> cmake

C. Link against pre-built

  .. code-block:: bash

    >>> cmake -DENABLE_gdma=ON -DCMAKE_PREFIX_PATH=/path/to/gdma/root

  .. code-block:: bash

    >>> cmake -DENABLE_gdma=ON -Dgdma_DIR=/path/to/gdma/configdir

D. Build bundled despite pre-built being detectable

  .. code-block:: bash

    >>> cmake -DENABLE_gdma=ON -DCMAKE_PREFIX_PATH=/path/to/unwanted/gdma/root/and/wanted/other/dependencies/root -DCMAKE_DISABLE_FIND_PACKAGE_gdma=ON

