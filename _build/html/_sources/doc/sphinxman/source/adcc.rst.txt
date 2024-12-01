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

.. index:: adcc, ADC

.. _`sec:adcc`:

Interface to adcc by M. F. Herbst and M. Scheurer
=================================================

.. codeauthor:: Michael F. Herbst
.. sectionauthor:: Michael F. Herbst

*Module:* :ref:`Keywords <apdx:adc>`, :ref:`PSI Variables <apdx:adc_psivar>`

.. image:: https://img.shields.io/badge/home-adcc-informational.svg
   :target: https://code.adc-connect.org

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: http://adc-connect.org/latest

|PSIfour| contains code to interface to the adcc python module developed
by M. F. Herbst *et. al.*. No additional licence or configuration
is required to use adcc. The module serves as the backend for
most algebraic-diagrammatic construction methods for correlated
excited states in |PSIfour|. For more details on ADC methods,
see :ref:`sec:adc`.

Installation
~~~~~~~~~~~~

For up to date information and more details,
see the `adcc installation documentation <https://adc-connect.org/latest/installation.html>`_.

**Binary**

* .. image:: https://anaconda.org/adcc/adcc/badges/version.svg
     :target: https://anaconda.org/adcc/adcc

* .. image:: https://img.shields.io/pypi/v/adcc
     :target: https://pypi.org/project/adcc

* adcc is available as a conda package for Linux and macOS
  and on pypi.

.. * If using the |PSIfour| binary, adcc has already been installed alongside.
..
.. * If using |PSIfour| built from source, and anaconda or miniconda has
..   already been installed (instructions at :ref:`sec:quickconda`),
..   adcc can be obtained through ``conda install adcc -c adcc``.
..   Then enable it as a feature with :makevar:`ENABLE_adcc`
..   and rebuild |PSIfour| to detect adcc and activate dependent code.
..
.. * Previous bullet had details. To build |PSIfour| from source and use
..   adcc from conda without thinking, consult :ref:`sec:condapsi4dev`.

* To remove a conda installation, ``conda remove adcc``.

**Source**

* .. image:: https://img.shields.io/github/tag-date/adc-connect/adcc.svg?maxAge=2592000
     :target: https://github.com/adc-connect/adcc

* If using |PSIfour| built from source and you want adcc installed as well,
  enable it as a feature with :makevar:`ENABLE_adcc`,
  and let the build system fetch and install it.


Keywords for adcc
~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/adc__cutoff_amps_print.rst
.. include:: autodir_options_c/adc__kind.rst
.. include:: autodir_options_c/adc__max_num_vecs.rst
.. include:: autodir_options_c/adc__maxiter.rst
.. include:: autodir_options_c/adc__num_core_orbitals.rst
.. include:: autodir_options_c/adc__num_guesses.rst
.. include:: autodir_options_c/adc__r_convergence.rst
.. include:: autodir_options_c/adc__reference.rst
.. include:: autodir_options_c/adc__roots_per_irrep.rst


.. _`cmake:adcc`:

How to configure adcc for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, adcc provides additional quantum-chemical methods
  (a wide range of ADC methods). In turn adcc can use |PSIfour| as the backend for
  self-consistent field calculations and required integrals.

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) adcc

* Upstream Dependencies |w---w| adcc (\ |dr| optional) |PSIfour|

**CMake Variables**

* :makevar:`ENABLE_adcc` |w---w| CMake variable toggling whether Psi4 automatically installs adcc

**Examples**

A. Build and install adcc if needed

  .. code-block:: bash

    >>> cmake -DENABLE_adcc=ON

B. Build *without* adcc

  .. code-block:: bash

    >>> cmake

