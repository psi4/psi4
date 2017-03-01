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

.. index:: gCP
.. _`sec:gcp`:

Interface to gCP by S. Grimme
=============================

.. codeauthor:: Lori A. Burns
.. sectionauthor:: Lori A. Burns

*Module:* :ref:`Samples <apdx:testSuitegcp>`

.. image:: https://img.shields.io/badge/home-gCP-5077AB.svg
   :target: http://www.thch.uni-bonn.de/tc/downloads/gcp/index.html

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: http://www.thch.uni-bonn.de/tc/downloads/gcp/data/manGCP.pdf

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/gcp/badges/version.svg
     :target: https://anaconda.org/psi4/gcp

* gCP is available as a conda package for Linux and macOS.

* If using the |PSIfour| binary, gCP has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  the gcp executable can be obtained through ``conda install gcp``.

* To remove a conda installation, ``conda remove gcp``.

**Source**

* .. image:: https://img.shields.io/badge/home-gCP-5077AB.svg
     :target: http://www.thch.uni-bonn.de/tc/downloads/gcp/getgCP.html

* If using |PSIfour| built from source and you want to build gCP from
  from source also, follow the instructions provided with the source
  (essentially, download the freely available tarball, unpack the source,
  edit the Makefile to select a Fortran compiler, and run make).

To be used by |PSIfour|, the program binary (``gcp``) must be
found in your :envvar:`PSIPATH` or :envvar:`PATH` (in that order). If
|PSIfour| is unable to execute the binary, an error will be reported.
To preferentially use a particular gcp compilation, simply adjust its
position in the path environment variables.


Running gCP
~~~~~~~~~~~

At present there is a limited interface to gCP that is used
only to implement the "HF-3c" [Sure:2013:1672]_ and "PBEh-3c"
[Grimme:2015:054107]_ methods (both energy and gradient). The :ref:`DFTD3
<sec:dftd3>` executable must also be available for these methods to
run. Unlike every other method in |PSIfour|, if a basis set has not been
set, these will default to their intended basis sets: MINIX for HF-3c
and def2-mSVP for PBEh-3c. If a basis has previously been set, but you
want to use the default basis, use the slash syntax to "empty" the basis
option for the scope of the current calculation, ``energy("hf3c/")``.

A few practical examples:

* HF-3c single point with default minix basis ::

   energy('hf3c')

* PBEh-3c optimization with default def2-mSVP basis ::

   optimize('pbeh3c')

* HF-3c with non-standard basis ::

   set basis cc-pvdz
   energy('hf3c')

* PBEh-3c with default basis after basis set ::

   set basis cc-pvdz
   energy('pbeh3c/')

If only BSSE/basis set corrections (rather than total energies) are of
interest, the ``gcp`` program can be run independently of the scf
through the python function :py:func:`~qcdb.interface_gcp.run_gcp`. (This function
is the same |PSIfour|/``gcp`` interface that is called during an scf job.)
This route is much faster than running a HF or DFT energy. ::

   molecule nene {
   Ne
   Ne 1 2.0
   }
   
   nene.update_geometry()

   >>> E, G = nene.run_gcp('hf3c')

   >>> E, G = nene.run_gcp(func='HF3c', verbose=True)

.. autofunction:: qcdb.interface_gcp.run_gcp

