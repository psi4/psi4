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

.. index:: LIBINT, integrals

.. _`sec:libint`:

Interface to LIBINT by E. Valeev
================================

.. codeauthor:: Edward F. Valeev and Justin T. Fermann
.. sectionauthor:: Lori A. Burns

.. *Module:* :ref:`Keywords <apdx:efp>`, :ref:`PSI Variables <apdx:efp_psivar>`, :source:`LIBEFP <src/lib/libefp_solver>`

.. image:: https://img.shields.io/badge/home-libint-5077AB.svg
   :target: https://github.com/psi4/libint

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: http://evaleev.github.io/libint/

|PSIfour|, particularly libmints utterly relies upon the LIBINT library
developed by E. Valeev from early roots by J. Fermann. Libint requires no
additional licence, downloads, or configuration. Conversely, |Psifour|
cannot build *without* libint.

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/libint/badges/version.svg
     :target: https://anaconda.org/psi4/libint

* libint is available as a conda package for Linux and macOS.

* If using the |PSIfour| binary, libint has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  libint can be obtained through ``conda install libint``.
  Then, hint its location with :makevar:`CMAKE_PREFIX_PATH`,
  and rebuild |PSIfour| to detect libint and activate dependent code.

* To remove a conda installation, ``conda remove libint``.

**Source**

* .. image:: https://img.shields.io/github/tag/psi4/libint.svg?maxAge=2592000
     :target: https://github.com/psi4/libint

* If using |PSIfour| built from source and you want libint built from
  from source also,
  let the build system fetch and build it and activate dependent code.


.. _`cmake:libint`:

How to configure libint for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, libint is a library that provides essential
  two-body molecular integrals over Gaussian functions.

* Downstream Dependencies |w---w| |PSIfour| |dr| libint

* Upstream Dependencies |w---w| libint |dr| None

**CMake Variables**

* :makevar:`CMAKE_PREFIX_PATH` |w---w| CMake list variable to specify where pre-built dependencies can be found. For libint, set to an installation directory containing ``include/libint/libint.h``
* :makevar:`libint_DIR` |w---w| CMake variable to specify where pre-built libint can be found. Set to installation directory containing ``share/cmake/libint/libintConfig.cmake``
* :makevar:`CMAKE_DISABLE_FIND_PACKAGE_libint` |w---w| CMake variable to force internal build of libint instead of detecting pre-built

**Examples**

A. Build bundled

  .. code-block:: bash

    >>> cmake

B. Link against pre-built

  .. code-block:: bash

    >>> cmake -DCMAKE_PREFIX_PATH=/path/to/libint/root

  .. code-block:: bash

    >>> cmake -Dlibint_DIR=/path/to/libint/configdir

C. Build bundled despite pre-built being detectable

  .. code-block:: bash

    >>> cmake -DCMAKE_PREFIX_PATH=/path/to/unwanted/libint/root/and/wanted/other/dependencies/root -DCMAKE_DISABLE_FIND_PACKAGE_libint=ON

