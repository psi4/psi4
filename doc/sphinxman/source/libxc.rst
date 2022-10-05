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

.. include:: autodoc_abbr_options_c.rst

.. index:: Libxc, DFT, functionals

.. _`sec:libxc`:

Interface to Libxc by M. A. L. Marques
======================================

.. codeauthor:: M. A. L. Marques and Micael Oliveira
.. sectionauthor:: Lori A. Burns

.. *Module:* :ref:`Keywords <apdx:efp>`, :ref:`PSI Variables <apdx:efp_psivar>`, :source:`LIBEFP <src/lib/libefp_solver>`

.. image:: https://img.shields.io/badge/home-libxc-5077AB.svg
   :target: https://gitlab.com/libxc/libxc

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: https://www.tddft.org/programs/libxc/manual/

|PSIfour|, relies upon the Libxc library for density functionals. Libxc
requires no
additional licence, downloads, or configuration. Conversely, |Psifour|
cannot build *without* Libxc.

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/libxc/badges/version.svg
     :target: https://anaconda.org/psi4/libxc

* Libxc is available as a conda package for Linux and macOS (and Windows, through the Ubuntu shell).

* If using the |PSIfour| binary, Libxc has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  Libxc can be obtained through ``conda install libxc -c psi4``.
  Then, hint its location with :makevar:`CMAKE_PREFIX_PATH`,
  and rebuild |PSIfour| to detect Libxc and activate dependent code.

.. * Previous bullet had details. To build |PSIfour| from source and use
..   Libxc from conda without thinking, consult.

* To remove a conda installation, ``conda remove libxc``.

**Source**

* .. image:: https://img.shields.io/github/tag/loriab/libxc.svg?maxAge=2592000
     :target: https://github.com/loriab/libxc/tree/libxc4retweaked

  Note that |PSIfour| has forked and slightly modified upstream Libxc from
  https://gitlab.com/libxc/libxc to regain functionality.

* If using |PSIfour| built from source and you want Libxc built from
  from source also,
  let the build system fetch and build it and activate dependent code.


.. _`cmake:libxc`:

How to configure Libxc for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, Libxc is a library that provides essential
  density functional definitions and compositions.

* Downstream Dependencies |w---w| |PSIfour| |dr| Libxc

* Upstream Dependencies |w---w| Libxc |dr| None

**CMake Variables**

* :makevar:`CMAKE_PREFIX_PATH` |w---w| CMake list variable to specify where pre-built dependencies can be found. For Libxc, set to an installation directory containing ``include/libxc/xc.h``
* :makevar:`Libxc_DIR` |w---w| CMake variable to specify where pre-built Libxc can be found. Set to installation directory containing ``share/cmake/Libxc/LibxcConfig.cmake``
* :makevar:`CMAKE_DISABLE_FIND_PACKAGE_Libxc` |w---w| CMake variable to force internal build of Libxc instead of detecting pre-built
* :makevar:`CMAKE_INSIST_FIND_PACKAGE_Libxc` |w---w| CMake variable to force detecting pre-built Libxc and not falling back on internal build

**Examples**

A. Build bundled

  .. code-block:: bash

    >>> cmake

B. Link against pre-built

  .. code-block:: bash

    >>> cmake -DCMAKE_PREFIX_PATH=/path/to/libxc/root

  .. code-block:: bash

    >>> cmake -DLibxc_DIR=/path/to/libxc/configdir

C. Build bundled despite pre-built being detectable

  .. code-block:: bash

    >>> cmake -DCMAKE_PREFIX_PATH=/path/to/unwanted/libxc/root/and/wanted/other/dependencies/root -DCMAKE_DISABLE_FIND_PACKAGE_Libxc=ON

