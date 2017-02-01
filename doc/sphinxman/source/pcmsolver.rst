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

.. index:: PCMSolver, PCM

.. _`sec:pcmsolver`:

Interface to PCMSolver by R. Di Remigio
=======================================

.. codeauthor:: Roberto Di Remigio, T. Daniel Crawford, Andrew C. Simmonett
.. sectionauthor:: Roberto Di Remigio

*Module:* :ref:`Keywords <apdx:pcm>`, :ref:`PSI Variables <apdx:pcm_psivar>`, :source:`PCMSolver <src/lib/libpsipcm>`

.. image:: https://img.shields.io/badge/home-PCMSolver-5077AB.svg
   :target: https://github.com/PCMSolver/pcmsolver

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: http://pcmsolver.readthedocs.io/en/latest/

|PSIfour| contains code to interface to the PCMSolver library developed
by R. Di Remigio and L. Frediani.
The PCMSolver library requires no additional licence, downloads, or
configuration. The library allows for calculations in solution with the
polarizable continuum model (PCM), a continuum solvation model.

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/pcmsolver/badges/version.svg
     :target: https://anaconda.org/psi4/pcmsolver

* PCMSolver is available as a conda package for Linux and macOS.

* If using the |PSIfour| binary, PCMSolver has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  PCMSolver can be obtained through ``conda install pcmsolver``.
  Then enable it as a feature with :makevar:`ENABLE_CheMPS2`,
  hint its location with :makevar:`CMAKE_PREFIX_PATH`,
  and rebuild |PSIfour| to detect PCMSolver and activate dependent code.

* To remove a conda installation, ``conda remove pcmsolver``.

**Source**

* .. image:: https://img.shields.io/github/tag/PCMSolver/pcmsolver.svg?maxAge=2592000
     :target: https://github.com/PCMSolver/pcmsolver

* If using |PSIfour| built from source and you want PCMSolver built from
  from source also,
  enable it as a feature with :makevar:`ENABLE_CheMPS2`,
  and let the build system fetch and build it and activate dependent code.

.. index:: PCM; Using PCM

.. _`sec:usingPCM`:

Using the polarizable continuum model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The inclusion of a PCM description of the solvent into your calculation
is achieved by setting |globals__pcm| ``true`` in your input file.
|Psifour| understands the additional option |globals__pcm_scf_type| with possible values ``total``
(the default) or ``separate``.
The latter forces the separate handling of nuclear and electronic electrostatic potentials and
polarization charges. It is mainly useful for debugging.

.. note:: At present PCM can only be used for energy calculations with SCF wavefunctions.
   Moreover, the PCMSolver library **cannot** exploit molecular point group symmetry.

The PCM model and molecular cavity are specified in a ``pcm`` section that has
to be explicitly typed in by the user. This additional section follows a syntax
that is slightly different from that of |Psifour| and is fully documented
`here <http://pcmsolver.readthedocs.org/en/latest/users/input.html>`_

A typical input for a Hartree--Fock calculation with PCM would look like the following: ::

    molecule NH3 {
    symmetry c1
    N     -0.0000000001    -0.1040380466      0.0000000000
    H     -0.9015844116     0.4818470201     -1.5615900098
    H     -0.9015844116     0.4818470201      1.5615900098
    H      1.8031688251     0.4818470204      0.0000000000
    units bohr
    no_reorient
    no_com
    }

    set {
      basis STO-3G
      scf_type pk
      pcm true
      pcm_scf_type total
    }

    pcm = {
       Units = Angstrom
       Medium {
       SolverType = IEFPCM
       Solvent = Water
       }

       Cavity {
       RadiiSet = UFF
       Type = GePol
       Scaling = False
       Area = 0.3
       Mode = Implicit
       }
    }

More examples can be found in the directories with PCM tests
:srcsample:`pcmsolver/pcm-scf`,
:srcsample:`pcmsolver/pcm-dft`, and
:srcsample:`pcmsolver/pcm-dipole`.

Keywords for PCMSolver
~~~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/globals__pcm.rst
.. include:: autodir_options_c/globals__pcm_scf_type.rst

.. _`cmake:pcmsolver`:

How to configure PCMSolver for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, PCMSolver is a library that provides additional
  quantum chemical capabilities (solvation modeling).

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) PCMSolver

* Upstream Dependencies |w---w| PCMSolver |dr| Fortran, ???

