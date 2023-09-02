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

.. index:: gau2grid, collocation

.. _`sec:gau2grid`:

Interface to gau2grid by D. G. A. Smith
=======================================

.. codeauthor:: D. G. A. Smith
.. sectionauthor:: Lori A. Burns

.. *Module:* :ref:`Keywords <apdx:efp>`, :ref:`PSI Variables <apdx:efp_psivar>`, :source:`LIBEFP <src/lib/libefp_solver>`

.. image:: https://img.shields.io/badge/home-gau2grid-5077AB.svg
   :target: https://github.com/dgasmith/gau2grid

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: https://github.com/dgasmith/gau2grid/blob/master/README.md

|PSIfour|, relies upon the gau2grid library for Gaussian collocations for DFT. gau2grid
requires no
additional licence, downloads, or configuration. Conversely, |Psifour|
cannot build *without* gau2grid.

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/gau2grid/badges/version.svg
     :target: https://anaconda.org/psi4/gau2grid

* gau2grid is available as a conda package for Linux and macOS (and Windows, through the Ubuntu shell).

* If using the |PSIfour| binary, gau2grid has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  gau2grid can be obtained through ``conda install gau2grid -c psi4``.
  Then, hint its location with :makevar:`CMAKE_PREFIX_PATH`,
  and rebuild |PSIfour| to detect gau2grid and activate dependent code.

.. * Previous bullet had details. To build |PSIfour| from source and use
..   gau2grid from conda without thinking, consult.

* To remove a conda installation, ``conda remove gau2grid``.

**Source**

* .. image:: https://img.shields.io/github/tag/dgasmith/gau2grid.svg?maxAge=2592000
     :target: https://github.com/dgasmith/gau2grid/tree/master

* If using |PSIfour| built from source and you want gau2grid built from
  from source also, let the build system fetch and build it and activate
  dependent code.


.. _`cmake:gau2grid`:

How to configure gau2gridfor building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, gau2grid is a library that provides essential
  grid operations for DFT.

* Downstream Dependencies |w---w| |PSIfour| |dr| gau2grid

* Upstream Dependencies |w---w| gau2grid |dr| None

**CMake Variables**

* :makevar:`CMAKE_PREFIX_PATH` |w---w| CMake list variable to specify where pre-built dependencies can be found. For gau2grid, set to an installation directory containing ``include/gau2grid/gau2grid.h``
* :makevar:`gau2grid_DIR` |w---w| CMake variable to specify where pre-built gau2grid can be found. Set to installation directory containing ``share/cmake/gau2grid/gau2gridConfig.cmake``
* :makevar:`CMAKE_DISABLE_FIND_PACKAGE_gau2grid` |w---w| CMake variable to force internal build of gau2grid instead of detecting pre-built
* :makevar:`CMAKE_INSIST_FIND_PACKAGE_gau2grid` |w---w| CMake variable to force detecting pre-built gau2grid and not falling back on internal build

**Examples**

A. Build bundled

  .. code-block:: bash

    >>> cmake

B. Link against pre-built

  .. code-block:: bash

    >>> cmake -DCMAKE_PREFIX_PATH=/path/to/gau2grid/root

  .. code-block:: bash

    >>> cmake -Dgau2grid_DIR=/path/to/gau2grid/configdir

C. Build bundled despite pre-built being detectable

  .. code-block:: bash

    >>> cmake -DCMAKE_PREFIX_PATH=/path/to/unwanted/gau2grid/root/and/wanted/other/dependencies/root -DCMAKE_DISABLE_FIND_PACKAGE_gau2grid=ON

