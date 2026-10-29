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

.. index:: cuEST
.. _`sec:cuest`:

Interface to the cuEST GPU library by NVIDIA
============================================

.. codeauthor:: NVIDIA Corporation
.. sectionauthor:: Lori A. Burns

|PSIfour| contains code to interface to the cuEST GPU library developed by
NVIDIA, which is available under a separate licence from NVIDIA. When
|PSIfour| is built with cuEST support, setting |globals__use_cuest| to ``true``
moves parts of an SCF or DFT computation onto the GPU.

.. caution:: The cuEST interface is under active development. Only the
   routes exercised by the test suite should be considered supported; others
   raise an error rather than fall back to the CPU.

.. _`sec:cuest_grid`:

The DFT quadrature grid
~~~~~~~~~~~~~~~~~~~~~~~

When |globals__use_cuest| is enabled, the DFT exchange-correlation grid is built
on and owned by the GPU. The division of labour is worth stating plainly,
because it determines which |PSIfour| grid keywords can reach the device at
all. |PSIfour| computes, on the host and with its own code, the radial
quadrature nodes and weights for each atom and the number of angular
(Lebedev) points for each radial shell. Those three arrays, together with the
atomic coordinates, are everything that is handed to cuEST. cuEST then
constructs the molecular grid on the device, applies its own nuclear
partitioning, and from that point on the grid exists only in GPU memory: it
has no host-side points or blocks, and the usual host post-processing
(blocking, weight screening, distant-point removal) never runs.

The consequence is a simple rule. A grid keyword is honoured under cuEST if,
and only if, it shapes the radial nodes, the radial weights, or the per-shell
angular counts. Keywords that feed the host post-processing have no route to
the device and are inert. In particular:

* Honoured: |scf__dft_radial_points|, |scf__dft_spherical_points|,
  |scf__dft_radial_scheme|, |scf__dft_bs_radius_alpha|,
  |scf__dft_pruning_scheme|, |scf__dft_pruning_alpha|, and
  |scf__dft_grid_name| (``SG0`` and ``SG1``). With these matched, a cuEST
  grid and the corresponding |PSIfour| grid agree to roughly :math:`10^{-10} E_h`.
* Inert: |scf__dft_block_scheme|, |scf__dft_block_max_points|,
  |scf__dft_block_min_points|, |scf__dft_block_max_radius|,
  |scf__dft_remove_distant_points|, |scf__dft_weights_tolerance|, and
  |scf__dft_basis_tolerance|. These are suppressed from the printed grid
  summary for a GPU-resident grid rather than reported as if they were in
  effect. The blocking keywords do not affect energies on the host either;
  the two tolerance keywords do, so a cuEST grid is not screened the way a
  host grid with a loosened tolerance would be.

Nuclear partitioning is the one case where a |PSIfour| keyword cannot simply
be passed through. cuEST applies Stratmann's partitioning scheme, fixed, with
no parameter through which another scheme could be requested. |PSIfour|
therefore treats |scf__dft_nuclear_scheme| as follows for a cuEST run:

* If it was left alone, it is reset to ``STRATMANN`` and the reset is
  announced, in the same spirit as resetting |globals__scf_type| to ``DF``
  for the density-fitted correlated methods::

      DFT Nuclear Scheme (re)set to STRATMANN for cuEST.

* If ``STRATMANN`` was requested explicitly, the run proceeds unremarkably.

* If any other scheme was requested explicitly, the computation raises a
  :py:exc:`ValidationError`. cuEST cannot deliver that scheme, and quietly
  substituting a different quadrature for the one that was asked for would be
  worse than refusing.

This matters when comparing GPU and CPU numbers: against the
|PSIfour| default of ``TREUTLER`` a cuEST energy differs by around
:math:`10^{-6} E_h` purely from the partitioning, and far more on a coarse
grid |w---w| an SG0 water computation differs by about :math:`10^{-2} E_h`.
Setting
|scf__dft_nuclear_scheme| to ``STRATMANN`` on the CPU side is the correct way
to isolate any other source of disagreement.

One further restriction comes from the angular grids. cuEST requires every
per-shell angular count to be a genuine Lebedev number, whereas |PSIfour|'s
own table begins with a degenerate one-point rule that it will select when a
pruning function decays far enough. The Gaussian pruning schemes
(``P_GAUSSIAN``, ``D_GAUSSIAN``, and ``LOG_GAUSSIAN`` with
|scf__dft_pruning_alpha| above roughly 0.2) reach that regime routinely. Such
shells are promoted to the smallest true Lebedev rule, of six points, and the
promotion is reported per atom in the output. The resulting grid is therefore
slightly denser than the host grid the same keywords would produce, so GPU
and CPU energies need not agree closely for those schemes; ``D_GAUSSIAN`` in
particular can differ by around :math:`10^{-4} E_h`. The region-based schemes
(``ROBUST``, ``TREUTLER``), the Slater-based schemes, and unpruned grids are
unaffected.

