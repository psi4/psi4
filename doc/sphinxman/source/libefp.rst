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

.. index:: LIBEFP, EFP

.. _`sec:libefp`:

Interface to LIBEFP by I. Kaliman
=================================

.. codeauthor:: Andrew C. Simmonett, A. Eugene DePrince III, Rollin A. King, and Lori A. Burns
.. sectionauthor:: Lori A. Burns

*Module:* :ref:`Keywords <apdx:efp>`, :ref:`PSI Variables <apdx:efp_psivar>`, :source:`LIBEFP <src/lib/libefp_solver>`

.. image:: https://img.shields.io/badge/home-libefp-5077AB.svg
   :target: https://github.com/ilyak/libefp

.. raw:: html

   <br>

.. image:: https://img.shields.io/badge/docs-latest-5077AB.svg
   :target: http://www.libefp.org/

|PSIfour| contains code to interface to the LIBEFP library developed
in L. Slipchenko's group by I. Kaliman. LIBEFP
requires no additional licence,
downloads, or configuration. Since February 2017, libefp is not required to build
|Psifour|.

Installation
~~~~~~~~~~~~

**Binary**

* .. image:: https://anaconda.org/psi4/libefp/badges/version.svg
     :target: https://anaconda.org/psi4/libefp

* libefp is available as a conda package for Linux and macOS.

* If using the |PSIfour| binary, libefp has already been installed alongside.

* If using |PSIfour| built from source, and anaconda or miniconda has
  already been installed (instructions at :ref:`sec:quickconda`),
  libefp can be obtained through ``conda install libefp``.
  Then enable it as a feature with :makevar:`ENABLE_libefp`,
  hint its location with :makevar:`CMAKE_PREFIX_PATH`,
  and rebuild |PSIfour| to detect libefp and activate dependent code.

* To remove a conda installation, ``conda remove libefp``.

**Source**

* .. image:: https://img.shields.io/github/tag/ilyak/libefp.svg?maxAge=2592000
     :target: https://github.com/ilyak/libefp

* If using |PSIfour| built from source and you want libefp built from
  from source also,
  enable it as a feature with :makevar:`ENABLE_libefp`,
  and let the build system fetch and build it and activate dependent code.

.. index:: EFP; library fragments
   pair: EFP; adding new

.. _`sec:findingEFPFragments`:

EFP Fragments
~~~~~~~~~~~~~

LIBEFP comes with a couple dozen ready-to-use fragments (water, benzene,
common solvents, etc.) listed :ref:`here <sec:availableEFPFragments>`
with source :source:`psi4/share/psi4/efpfrag`.  Any of these may be used directly in
a |PSIfour| input file as described :ref:`here <sec:usingEFPFragments>`.

.. comment .. note:: The built-in fragment library distributed with Q-Chem (as of version 4.0.1) is *not*
.. comment    equivalent to that distributed with LIBEFP. Although many of the same
.. comment    molecules are present and should perform similarly in computations,
.. comment    exact matches of fragment geometries or efp energies should not be
.. comment    expected. See files in test case directories :source:`qchem-efp-sp
.. comment    <tests/libefp/qchem-efp-sp>` and :source:`qchem-qmefp-sp
.. comment    <tests/libefp/qchem-qmefp-sp>` for equivalent Q-Chem and |PSIfour|
.. comment    EFP input files.

Creating new efp fragments requires the `GAMESS
<http://www.msg.ameslab.gov/gamess/>`_ quantum chemistry package.
Instructions on building new fragments are `here
<https://github.com/libefp/libefp#how-to-create-custom-efp-fragment-types>`_.
Once your new fragment is ready, make it accessible to |PSIfour| by
including the directory in which the ``.efp`` file is located to the colon
separated environment variable :envvar:`PSIPATH`. Fragments are searched
for first in the current directory, next in the paths of :envvar:`PSIPATH`, and
finally in built-in library. If |PSIfour| is unable to find the
fragment, an error will be reported.

.. note:: When constructing new fragment files, the name of the name of the
   file should be lowercase and have extension ``.efp``. The molecule name
   within the file, *e.g.*, ``$NH3`` must correspond to the name of the
   fragment file.

.. index:: molecule; EFP
   single: EFP; molecule specification

.. _`sec:usingEFPFragments`:

Molecule Specification
~~~~~~~~~~~~~~~~~~~~~~

EFP fragment geometries are specified alongside the quantum mechanical
(QM) molecule and make use of the ``--`` fragment separation scheme
described :ref:`here <sec:fragments>`. Each EFP fragment has its own
fragment section that includes the label ``efp``, the name of the file
*fragname* from which EFP parameters are to be read, and the position
specification for the fragment in one of two ways, XYZABC or POINTS. For
XYZABC, the fragment specification is all on one line: ``efp`` and
*fragname* are followed by two sets of three numbers: the coordinates
of the center of mass of the fragment and the three Euler angles that
specify orientation about the center of mass. This format is compact
but not readily generated from molecule viewing software. ::

    efp  nh3  0.0 0.0 5.0  5 2 8

More convenient is the POINTS fragment specification. This consists of
four lines, the first of which is ``efp`` and *fragname*. The next lines
are the coordinates (without element labels) of the first three atoms
in the fragment. Note that EFP fragment geometries are rigid, so the
first atom will be placed exactly where specified by the first point,
the second atom will be placed along the vector between the first and
second points, and the third atom will be placed in the plane formed
by the three points. ::

    efp ch3oh
    1.275    -2.447    -4.673
    0.709    -3.191    -3.592
    2.213    -1.978    -4.343

.. note:: At present, |PSIfour| has limited support for diatomic
   and monoatomic EFP fragments. Single points are allowed when the
   di-/mono-atomic fragments are specified in XYZABC format. Optimizations
   are not allowed.

:ref:`Just as for QM <sec:moleculeKeywords>`, the center of mass
coordinates in the XYZABC format and all coordinates in the POINTS format are
taken to be in Angstroms by default or in Bohr if ``units au`` is present.
Charge and multiplicity specifications are encoded in the fragment file
and so are not read from input.

Any combination of EFP and QM fragments can be placed in a molecule; even
the oddity below is legitimate. Note that symmetry and reorientation are
automatically turned off when EFP fragments are present (``symmetry c1``
and ``no_com`` and ``no_reorient`` are implied). ::

    molecule qmefp {
      efp nh3 0.0 0.0 5.0 5 2 8
      --
      C  0.0 0.0 0.0
      O  0.0 1.5 0.0
      O  0.0 -1.5 0.0
      --
      efp h2o 5.0 0.0 0.0 5 2 8
      --
      He  -3.0 4.0 4.0
      He  -4.0 5.0 4.0
      --
      efp ch3oh
      1.275    -2.447    -4.673
      0.709    -3.191    -3.592
      2.213    -1.978    -4.343
    }


Running EFP 
~~~~~~~~~~~~

EFP can be invoked in similar fashion as other theories provided in |PSIfour|.
For example, if you want to obtain the EFP interaction energy for benzene and two waters,
simply provide the following::

   molecule {
     efp c6h6  0.0 0.0 0.0   0.0 0.0 0.0
     --
     efp h2o   4.0 0.0 0.0   0.0 0.0 0.0
     --
     efp h2o  -4.0 0.0 0.0   0.0 0.0 0.0
   }
   
   energy('efp')

This computation involves purely EFP/EFP fragment interactions and is
performed entirely by the LIBEFP library.  |PSIfour| can also handle mixed
systems of quantum mechanical (QM) and EFP fragments through the native
:ref:`SCF <sec:scf>` code augmented by calls to the LIBEFP library. For
example, turning one of the waters in the example above into a QM
fragment is shown below. ::

   molecule {
     efp c6h6  0.0 0.0 0.0   0.0 0.0 0.0
     --
     O  4.0   0.0   0.0
     H  4.7   0.7   0.0
     H  3.3  -0.7   0.0 
     --
     efp h2o  -4.0 0.0 0.0   0.0 0.0 0.0
   }
   
   set basis 6-31g
   energy('scf')

Whenever an EFP fragment is present in the active molecule, the SCF energy
will include EFP contributions.

.. warning:: Although the EFP geometry is specified alongside the QM
   geometry in a ``molecule name {...}`` block, internally the handling
   of EFP is not so clean. In straightforward input files that involve
   any number of [molecule block, energy/opt/etc, clean()] portions,
   there should be no problem; the energy/opt computation will always
   be run on the molecule defined in the preceding block. For advanced
   users, unexpected difficulties may arise due to: (1) the EFP fragment
   from the last molecule block executed will always be active (and
   potentially interfering with SCF) and (2) recalling a molecule
   through ``activate(name)`` (where ``name`` was the python handle
   in the molecule block) will not load up any EFP portion of that
   molecule. This divergent treatment is a stopgap while we determine
   how best to handle molecules with different domains.

At this time, |PSIfour| is only able to perform pure-efp single-points and
geometry optimizations and mixed qm/efp SCF single-points.

.. _`table:libefpauto`:

    .. _`table:libefp_methods`:

    +-------------------------+----------------------+--------------------------------------------------------------------------+
    | name                    | molecule composition | calls method                                                             |
    +=========================+======================+==========================================================================+
    | efp                     | pure EFP             | EFP interaction energy (IE) on all frags                                 |
    +                         +----------------------+--------------------------------------------------------------------------+
    | efp                     | mixed QM/EFP         | EFP IE on EFP frags only                                                 |
    +                         +----------------------+--------------------------------------------------------------------------+
    | efp                     | pure QM              | *error*                                                                  |
    +-------------------------+----------------------+--------------------------------------------------------------------------+
    | scf                     | pure EFP             | *error*                                                                  |
    +                         +----------------------+--------------------------------------------------------------------------+
    | scf                     | mixed QM/EFP         | SCF energy on QM frags w/coupling to EFP frags, plus EFP IE on EFP frags | 
    +                         +----------------------+--------------------------------------------------------------------------+
    | scf                     | pure QM              | SCF energy on all frags (normal |Psifour| operation)                     |
    +-------------------------+----------------------+--------------------------------------------------------------------------+

.. index:: EFP; library fragments

.. _`sec:availableEFPFragments`:

Fragment Library
~~~~~~~~~~~~~~~~

Below are documented the EFP fragments available from the LIBEFP library.
These systems are accessible in ``molecule {...}`` blocks without
additional configuration.

----

.. comment This toctree directive only here to suppress warning at build time.
   include line below is doing the work.

.. toctree::
   :hidden:

   autodoc_available_efpfrag

.. include:: autodoc_available_efpfrag.rst


.. _`cmake:libefp`:

How to configure libefp for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, libefp is a library that provides additional
  molecular modeling capabilities (EFP).

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) libefp

* Upstream Dependencies |w---w| libefp |dr| None

**CMake Variables**

* :makevar:`ENABLE_libefp` |w---w| CMake variable toggling whether Psi4 builds with libefp
* :makevar:`CMAKE_PREFIX_PATH` |w---w| CMake list variable to specify where pre-built dependencies can be found. For libefp, set to an installation directory containing ``include/efp.h``
* :makevar:`libefp_DIR` |w---w| CMake variable to specify where pre-built libefp can be found. Set to installation directory containing ``share/cmake/libefp/libefpConfig.cmake``
* :makevar:`CMAKE_DISABLE_FIND_PACKAGE_libefp` |w---w| CMake variable to force internal build of libefp instead of detecting pre-built

**Examples**

A. Build bundled

  .. code-block:: bash

    >>> cmake -DENABLE_libefp=ON

B. Build *without* libefp

  .. code-block:: bash

    >>> cmake

C. Link against pre-built

  .. code-block:: bash

    >>> cmake -DENABLE_libefp=ON -DCMAKE_PREFIX_PATH=/path/to/libefp/root

  .. code-block:: bash

    >>> cmake -DENABLE_libefp=ON -Dlibefp_DIR=/path/to/libefp/configdir

D. Build bundled despite pre-built being detectable

  .. code-block:: bash

    >>> cmake -DENABLE_libefp=ON -DCMAKE_PREFIX_PATH=/path/to/unwanted/libefp/root/and/wanted/other/dependencies/root -DCMAKE_DISABLE_FIND_PACKAGE_libefp=ON

