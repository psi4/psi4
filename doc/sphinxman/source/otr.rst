.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2026 The Psi4 Developers.
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

.. index:: OTR, OpenTrustRegion

.. _`sec:otr`:

Interface to OpenTrustRegion by J. Greiner
==========================================

.. codeauthor:: Jonas Greiner
.. sectionauthor:: Lori A. Burns

.. image:: https://img.shields.io/badge/home-OpenTrustRegion-5077AB.svg
   :target: https://github.com/eriksen-lab/OpenTrustRegion

.. raw:: html

   <br>

OpenTrustRegion is a black-box second-order orbital optimizer developed by
J. Greiner in the Eriksen lab and interfaced to |PSIfour|. Rather than iterating the Fock matrix
to self-consistency, it minimizes the SCF energy directly with respect to the
orbital rotation parameters, using a trust-region method with the orbital
Hessian applied on the fly.
Enabling OTR and adding ``set orbital_optimizer_package opentrustregion`` or
``set orbital_optimizer_package otr``
replaces the SCF iterations with OTR. Set to ``internal`` to revoke.

Installation
~~~~~~~~~~~~

**Binary**

* A conda package for OpenTrustRegion is available. Obtain it through
  ``conda install opentrustregion -c conda-forge``, then enable it as a feature
  with :makevar:`ENABLE_OpenTrustRegion`, hint its location with
  :makevar:`CMAKE_PREFIX_PATH`, and rebuild |PSIfour| to detect OpenTrustRegion
  and activate dependent code.

**Source**

* .. image:: https://img.shields.io/github/tag/eriksen-lab/OpenTrustRegion.svg?maxAge=2592000
     :target: https://github.com/eriksen-lab/OpenTrustRegion

* If using |PSIfour| built from source and you want OpenTrustRegion built from
  source also,
  enable it as a feature with :makevar:`ENABLE_OpenTrustRegion=ON`,
  and let the build system fetch and build it and activate dependent code.
  Note that OpenTrustRegion is written in Fortran, so a Fortran compiler is
  required for the source build.


.. _`options:otr`:

OpenTrustRegion Options
~~~~~~~~~~~~~~~~~~~~~~~

OpenTrustRegion is ready to use for RHF, UHF, ROHF and the corresponding KS
references. Rather than taking dedicated keywords, it maps the |PSIfour| SCF
convergence keywords onto its own settings:

* |scf__d_convergence| sets the convergence threshold on the RMS orbital
  gradient. If |scf__soscf_conv| has been set, the tighter of the two is used.
* |scf__maxiter| caps the number of macro-iterations, and |scf__soscf_max_iter|,
  if set, caps the micro-iterations within each one.
* |globals__print| selects OpenTrustRegion's own verbosity level.

Because OpenTrustRegion replaces the whole SCF iteration loop, |PSIfour| falls
back on the internal optimizer whenever the computation needs something applied
per-iteration by the driver or something OpenTrustRegion cannot express:

* a ``CUHF`` reference, for which no orbital Hessian is implemented,
* meta-GGA and VV10 functionals,
* the LinK exchange build (|globals__scf_type| ``DFDIRJ+LINK``),
  which does not support the non-symmetric K matrices the Hessian needs,
* MOM (|scf__mom_start|) and fractional occupation (|scf__frac_start|), which
  change the occupation during the SCF,
* EFP, PCM, DDX and PE embedding, whose contributions are added to the Fock
  matrix by the Python driver on each iteration.

A note printed to the output file names whichever condition applied.

Some further differences from the internal optimizer are worth knowing about:

* OpenTrustRegion optimizes the orbitals at a *fixed* occupation, whereas the
  internal optimizer re-runs the aufbau assignment every iteration. |PSIfour|
  relaxes the occupation once against the guess Fock matrix before handing over,
  and reconverges if canonicalizing the result changes it again, but a poor guess
  can still settle on a different SCF solution than DIIS would find.
* Being a second-order method, it converges in far fewer macro-iterations, so
  ``SCF ITERATIONS`` is typically 3-6 rather than 8-15.
* With |globals__scf_type| ``DFDIRJ+COSX`` the final COSX grid is
  used throughout, since OpenTrustRegion has no notion of loosening screening
  for early iterations. As with OpenOrbitalOptimizer, this shifts the energy by
  a few times 1e-5 hartree relative to the internal optimizer's early-screened
  result.
* Stability analysis (|scf__stability_analysis|) runs through |PSIfour|'s own
  code as usual; OpenTrustRegion's internal stability check is left off so that
  both optimizers land on the same solution.

.. include:: autodir_options_c/globals__orbital_optimizer_package.rst
.. include:: autodir_options_c/scf__d_convergence.rst
.. include:: autodir_options_c/scf__maxiter.rst
.. include:: autodir_options_c/scf__soscf_conv.rst
.. include:: autodir_options_c/scf__soscf_max_iter.rst


.. _`cmake:otr`:

How to configure OpenTrustRegion for building Psi4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Role and Dependencies**

* Role |w---w| In |PSIfour|, OpenTrustRegion is a library that provides alternate
  orbital optimization.

* Downstream Dependencies |w---w| |PSIfour| (\ |dr| optional) OpenTrustRegion

* Upstream Dependencies |w---w| OpenTrustRegion |dr| LAPACK

**CMake Variables**

* :makevar:`ENABLE_OpenTrustRegion` |w---w| CMake variable toggling whether |PSIfour| builds with OpenTrustRegion
* :makevar:`CMAKE_PREFIX_PATH` |w---w| CMake list variable to specify where pre-built dependencies can be found. For OTR, set to an installation directory containing ``include/opentrustregion.h``
* :makevar:`OpenTrustRegion_DIR` |w---w| CMake variable to specify where pre-built OpenTrustRegion can be found. Set to installation directory containing ``lib/cmake/OpenTrustRegion/OpenTrustRegionConfig.cmake``
* :makevar:`CMAKE_DISABLE_FIND_PACKAGE_OpenTrustRegion` |w---w| CMake variable to force internal build of OpenTrustRegion instead of detecting pre-built
* :makevar:`CMAKE_INSIST_FIND_PACKAGE_OpenTrustRegion` |w---w| CMake variable to force detecting pre-built OpenTrustRegion and not falling back on internal build

**Examples**

A. Build bundled

  .. code-block:: bash

    >>> cmake -DENABLE_OpenTrustRegion=ON

B. Build *without* OpenTrustRegion

  .. code-block:: bash

    >>> cmake

C. Link against pre-built

  .. code-block:: bash

    >>> cmake -DENABLE_OpenTrustRegion=ON -DCMAKE_PREFIX_PATH=/path/to/OpenTrustRegion/root

  .. code-block:: bash

    >>> cmake -DENABLE_OpenTrustRegion=ON -DOpenTrustRegion_DIR=/path/to/otr/cmakeconfigdir

D. Build bundled despite pre-built being detectable

  .. code-block:: bash

    >>> cmake -DENABLE_OpenTrustRegion=ON -DCMAKE_PREFIX_PATH=/path/to/unwanted/OpenTrustRegion/root/and/wanted/other/dependencies/root -DCMAKE_DISABLE_FIND_PACKAGE_OpenTrustRegion=ON
