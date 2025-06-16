.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2025 The Psi4 Developers.
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

.. index:: gCP
.. _`sec:gcp`:

Interface to gCP by S. Grimme
=============================

.. codeauthor:: Lori A. Burns
.. sectionauthor:: Lori A. Burns

*Module:* :ref:`Samples <apdx:testSuitegcp>`

.. image:: https://img.shields.io/badge/home-gCP-5077AB.svg
   :target: https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/gcp

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/gcp/manGCP.pdf

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/gcp/badges/version.svg
     :target: https://anaconda.org/psi4/gcp

* There are two implementations of gCP; see :ref:`table:empdispimpl` . The newer
  "mctc" one is preferred, while the older "classic" one will work for the immediate future.
  |PSIfour| will automatically select whichever is available.
  Starting with v1.9, only "mctc-gcp" is supported, though the now untested
  "classic" continues to work for many applications.

* gCP is available as a conda package for Linux and macOS and Windows.

* If using the Psi4conda installer, gCP has already been installed alongside.

* If using the |PSIfour| conda package, the classic gcp conda package can
  be obtained through ``conda install gcp -c psi4`` or the newer implementation
  through ``conda install gcp-correction -c conda-forge``.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  the gcp executable can be obtained through ``conda install gcp -c psi4``
  or ``conda install gcp-correction -c conda-forge``.

* To remove a conda installation, ``conda remove gcp`` or ``conda remove gcp-correction``.

**Source**

* .. image:: https://img.shields.io/badge/home-gCP-5077AB.svg
     :target: https://www.chemie.uni-bonn.de/pctc/mulliken-center/software/gcp/getgcp

* If using |PSIfour| built from source and you want to build gCP from
  from source also, follow the instructions provided with the source
  (essentially, download the freely available tarball, unpack the source,
  edit the Makefile to select a Fortran compiler, and run make).

To be used by |PSIfour|, the program binary (``gcp`` or ``mctc-gcp``) must be
found in your :envvar:`PATH` so that QCEngine can detect it. Check if and where
found through ``qcengine info``. If
|PSIfour| is unable to execute the binary, an error will be reported.
To preferentially use a particular gcp compilation, simply adjust its
position in the path environment variables.


Running gCP
~~~~~~~~~~~

At present there is a limited interface to gCP that is used
only to implement the "HF-3c" [Sure:2013:1672]_, "PBEh-3c"
[Grimme:2015:054107]_, "B97-3c" [Brandenburg:2018:b973c]_, "r2SCAN-3c" [Grimme:2021:064103]_,
and "wB97X-3c" [Muller:2023:014103]_ methods (both energy and gradient).
The interface can use classic or mctc-gcp executables but only the latter implements "B97-3c" and "r2SCAN-3c".
The newest wB97X-3c method doesn't use a gcp correction (it does use ECPs down to first row elements)
but is listed here for completeness of the "3c" family.
A :ref:`DFTD3 <sec:dftd3>` executable, classic or simple-dftd3, must also be available for
the HF-3c, PBEh-3c, or B97-3c methods to run.
A :ref:`DFTD4 <sec:dftd3>` python module must also be available for
the r2SCAN-3c or wB97X-3c methods to run.
These method are defined with their own basis set and thus no basis set should be set by the user.
|PSIfour| will select the intended basis sets: HF-3c/MINIX, PBEh-3c/def2-mSVP, B97-3c/def2-mTZVP, r2SCAN-3c/def2-mTZVPP, wB97X-3c/vDZP.
If a basis has previously been set for another calculation, use the slash syntax to "empty" the basis
option for the scope of the current calculation, ``energy("hf3c/")``.

A few practical examples:

* HF-3c single point with default minix basis ::

   energy('hf3c')

* PBEh-3c optimization with default def2-mSVP basis ::

   optimize('pbeh3c')

* r2SCAN-3c with default basis after basis set ::

   set basis cc-pvdz
   energy('r2scan3c/')

If only BSSE/basis set corrections (rather than total energies) are of
interest, the ``gcp`` program can be run independently of the scf
through the python function :py:func:`~qcdb.Molecule.run_gcp`. (This function
is the same |PSIfour|/``gcp`` interface that is called during an scf job.)
This route is much faster than running a HF or DFT energy. ::

   molecule nene {
   Ne
   Ne 1 2.0
   }

   nene.update_geometry()

   >>> E, G = nene.run_gcp('hf3c')

   >>> E, G = nene.run_gcp(func='HF3c', verbose=True)

.. autofunction:: qcdb.Molecule.run_gcp

