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

.. index::
   pair: plugin; v2rdm_casscf

.. _`sec:v2rdm_casscf`:

Plugin v2rdm_casscf by A. E. DePrince
=====================================

.. codeauthor:: A. E. DePrince
.. sectionauthor:: Lori A. Burns

.. *Module:* :ref:`Keywords <apdx:pcm>`, :ref:`PSI Variables <apdx:pcm_psivar>`, :source:`PCMSolver <src/lib/libpsipcm>`

.. image:: https://img.shields.io/badge/home-v2rdm_casscf-5077AB.svg
   :target: https://github.com/edeprince3/v2rdm_casscf

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: https://github.com/edeprince3/v2rdm_casscf/blob/master/README.md

A variational 2-RDM-driven CASSCF plugin to |PSIfour|

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/v2rdm_casscf/badges/version.svg
     :target: https://anaconda.org/psi4/v2rdm_casscf

* v2rdm_casscf is available as a conda package for Linux and macOS (and Windows, through the Ubuntu shell).

* If using the Psi4conda installer, v2rdm_casscf has already been installed alongside.

* If using the |PSIfour| conda package, the v2rdm_casscf conda package can
  be obtained through ``conda install v2rdm_casscf -c psi4`` or ``conda install
  psi4-rt -c psi4``.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  v2rdm_casscf can be obtained through ``conda install v2rdm_casscf -c psi4``.
  Then, hint its location with :envvar:`PYTHONPATH`.

* To remove a conda installation, ``conda remove v2rdm_casscf``.

**Source**

* .. image:: https://img.shields.io/github/tag/edeprince3/v2rdm_casscf.svg?maxAge=2592000
     :target: https://github.com/edeprince3/v2rdm_casscf

* If using |PSIfour| built from source and you want v2rdm_casscf built from
  from source also,
  build it, then hint its location with :envvar:`PYTHONPATH`.

