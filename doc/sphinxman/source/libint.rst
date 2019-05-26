.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2019 The Psi4 Developers.
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

.. index:: Libint, integrals

.. _`sec:libint`:

Interface to Libint by E. Valeev
================================

.. codeauthor:: Edward F. Valeev and Justin T. Fermann
.. sectionauthor:: Lori A. Burns

.. *Module:* :ref:`Keywords <apdx:efp>`, :ref:`PSI Variables <apdx:efp_psivar>`, :source:`LIBEFP <src/lib/libefp_solver>`

.. image:: https://img.shields.io/badge/home-libint-5077AB.svg
   :target: https://github.com/evaleev/libint

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: http://evaleev.github.io/libint/

|PSIfour|, particularly libmints utterly relies upon the Libint library
developed by E. Valeev from early roots by J. Fermann. Libint requires no
additional licence, downloads, or configuration. Conversely, |Psifour|
cannot build *without* Libint.

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/libint/badges/version.svg
     :target: https://anaconda.org/psi4/libint

* Libint is available as a conda package for Linux and macOS (and Windows, through the Ubuntu shell).

* If using the |PSIfour| binary, Libint has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  Libint can be obtained through ``conda install libint -c psi4``.
  Then, hint its location with :makevar:`CMAKE_PREFIX_PATH`,
  and rebuild |PSIfour| to detect Libint and activate dependent code.

* Previous bullet had details. To build |PSIfour| from source and use 
  Libint from conda without thinking, consult :ref:`sec:condapsi4dev`.

* To remove a conda installation, ``conda remove libint``.

**Source**

* .. image:: https://img.shields.io/github/tag/evaleev/libint.svg?maxAge=2592000
     :target: https://github.com/evaleev/libint/tree/v1

  Note that |PSIfour| uses v1.

* If using |PSIfour| built from source and you want Libint built from
  from source also,
  let the build system fetch and build it and activate dependent code.


.. _`cmake:libint`:

How to configure Libint for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, Libint is a library that provides essential
  two-body molecular integrals over Gaussian functions.

* Downstream Dependencies |w---w| |PSIfour| |dr| Libint

* Upstream Dependencies |w---w| Libint |dr| None

**CMake Variables**

* :makevar:`CMAKE_PREFIX_PATH` |w---w| CMake list variable to specify where pre-built dependencies can be found. For Libint, set to an installation directory containing ``include/libint/libint.h``
* :makevar:`Libint_DIR` |w---w| CMake variable to specify where pre-built Libint can be found. Set to installation directory containing ``share/cmake/Libint/LibintConfig.cmake``
* :makevar:`CMAKE_DISABLE_FIND_PACKAGE_Libint` |w---w| CMake variable to force internal build of Libint instead of detecting pre-built
* :makevar:`CMAKE_INSIST_FIND_PACKAGE_Libint` |w---w| CMake variable to force detecting pre-built Libint and not falling back on internal build
* :makevar:`MAX_AM_ERI` |w---w| CMake variable to specify minimum highest angular momentum built or detected

**Examples**

A. Build bundled

  .. code-block:: bash

    >>> cmake

B. Link against pre-built

  .. code-block:: bash

    >>> cmake -DCMAKE_PREFIX_PATH=/path/to/libint/root

  .. code-block:: bash

    >>> cmake -DLibint_DIR=/path/to/libint/configdir

C. Build bundled despite pre-built being detectable

  .. code-block:: bash

    >>> cmake -DCMAKE_PREFIX_PATH=/path/to/unwanted/libint/root/and/wanted/other/dependencies/root -DCMAKE_DISABLE_FIND_PACKAGE_Libint=ON

D. Build bundled shared library with AM=6

  .. code-block:: bash

    >>> cmake -DMAX_AM_ERI=6 -DBUILD_SHARED_LIBS=ON

