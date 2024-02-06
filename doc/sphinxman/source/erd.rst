.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2024 The Psi4 Developers.
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

.. index:: ERD
.. _`sec:erd`:

Interface to ERD by N. Flocke and V. Lotrich
============================================

.. codeauthor:: Andrew C. Simmonett and Benjamin P. Pritchard
.. sectionauthor:: Lori A. Burns

.. *Module:* :ref:`Keywords <apdx:dkh>`, :ref:`Samples <apdx:testSuitedkh>`

.. image:: https://img.shields.io/badge/home-erd-5077AB.svg
   :target: https://github.com/psi4/erd

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: http://onlinelibrary.wiley.com/doi/10.1002/jcc.21018/abstract

.. _`sec:erdinstall`:

These are the AcesIII electron repulsion integrals that have
been partially interfaced into libmints. Enabling erd and adding
``set integral_package erd`` (do this in ``~/.psi4rc`` for universal
effect) runs libderiv from Libint for derivative integrals and erd for
non-derivative integrals.

.. warning:: The interface between erd and libderiv is not fully
   debugged. So analytic gradients, particularly density-fitted ones,
   are wrong, as are ESP calculations and some energies for long-range
   corrected ("omega") functionals. Insofar as faulty answers are
   anticipated with |globals__integral_package| ``erd``, |PSIfour| will
   throw an error if you try to execute that class of computation. But
   there may be more, so use with caution.

.. warning:: erd seems to be having some problems with Intel 2018 compilers. presently disabled in conda package.

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/erd/badges/version.svg
     :target: https://anaconda.org/psi4/erd

* ERD is available as a conda package for Linux and macOS (and Windows, through the Ubuntu shell).

.. * If using the |PSIfour| binary, erd has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  erd can be obtained through ``conda install erd -c psi4``.
  Then enable it as a feature with :makevar:`ENABLE_erd`,
  hint its location with :makevar:`CMAKE_PREFIX_PATH`,
  and rebuild |PSIfour| to detect erd and activate dependent code.

.. * Previous bullet had details. To build |PSIfour| from source and use
..   erd from conda without thinking, consult.

* To remove a conda installation, ``conda remove erd``.

**Source**

* .. image:: https://img.shields.io/github/tag/psi4/erd.svg?maxAge=2592000
     :target: https://github.com/psi4/erd

* If using |PSIfour| built from source and you want erd built from
  from source also,
  enable it as a feature with :makevar:`ENABLE_erd`,
  and let the build system fetch and build it and activate dependent code.

.. .. _`sec:erdinput`:

.. _`cmake:erd`:

How to configure erd for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, ERD is a library that provides alternate
  integrals.

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) erd

* Upstream Dependencies |w---w| erd |dr| Fortran

**CMake Variables**

* :makevar:`ENABLE_erd` |w---w| CMake variable toggling whether |PSIfour| builds with erd
* :makevar:`CMAKE_PREFIX_PATH` |w---w| CMake list variable to specify where pre-built dependencies can be found. For erd, set to an installation directory containing ``include/ERD/ERD_MANGLE.h``
* :makevar:`erd_DIR` |w---w| CMake variable to specify where pre-built erd can be found. Set to installation directory containing ``share/cmake/erd/erdConfig.cmake``
* :makevar:`CMAKE_DISABLE_FIND_PACKAGE_erd` |w---w| CMake variable to force internal build of erd instead of detecting pre-built
* :makevar:`CMAKE_INSIST_FIND_PACKAGE_erd` |w---w| CMake variable to force detecting pre-built erd and not falling back on internal build

**Examples**

A. Build bundled

  .. code-block:: bash

    >>> cmake -DENABLE_erd=ON

B. Build *without* erd

  .. code-block:: bash

    >>> cmake

C. Link against pre-built

  .. code-block:: bash

    >>> cmake -DENABLE_erd=ON -DCMAKE_PREFIX_PATH=/path/to/erd/root

  .. code-block:: bash

    >>> cmake -DENABLE_erd=ON -Derd_DIR=/path/to/erd/configdir

D. Build bundled despite pre-built being detectable

  .. code-block:: bash

    >>> cmake -DENABLE_erd=ON -DCMAKE_PREFIX_PATH=/path/to/unwanted/erd/root/and/wanted/other/dependencies/root -DCMAKE_DISABLE_FIND_PACKAGE_erd=ON

