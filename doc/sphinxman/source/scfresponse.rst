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

.. index:: SCF response, polarizability, hyperpolarizability, second-harmonic generation, SHG, optical rectification, Pockels effect
.. _`sec:scfresponse`:

SCF Response: Frequency-Dependent Polarizabilities and Hyperpolarizabilities
============================================================================

.. codeauthor:: edeprince3 (original scf_response plugin) and zw0930
.. sectionauthor:: zw0930

*Module:* :ref:`Keywords <apdx:scfresponse>`, :ref:`PSI Variables <apdx:scfresponse_psivar>`, :source:`SCFRESPONSE <psi4/src/psi4/scfresponse>`

The SCFRESPONSE module computes electric-dipole response properties of
Hartree--Fock and Kohn--Sham wavefunctions from the time-dependent SCF
(TDHF/TDDFT, random phase approximation) response equations:

* the frequency-dependent dipole polarizability :math:`\alpha(-\omega;\omega)`
  for HF and DFT, RHF/RKS and UHF/UKS;
* first hyperpolarizabilities for HF (RHF and UHF): static
  :math:`\beta(0;0,0)`, second-harmonic generation (SHG)
  :math:`\beta(-2\omega;\omega,\omega)`, optical rectification
  :math:`\beta(0;\omega,-\omega)`, and the Pockels effect
  :math:`\beta(-\omega;\omega,0)`, from the quadratic-response
  (:math:`2n+1`) expression using the first-order response vectors.

The module is based on the `scf_response <https://github.com/edeprince3/scf_response>`_
plugin. Static polarizabilities without frequencies continue to be computed by
the CPSCF solver described in :ref:`sec:oeprop`; the SCFRESPONSE module is used
when frequencies are given through |scfresponse__omega| or a hyperpolarizability
is requested.

Running a calculation
~~~~~~~~~~~~~~~~~~~~~

Properties are requested through :py:func:`~psi4.driver.properties` with an SCF
method (``scf``, ``hf``, or a density functional)::

    molecule {
    O  0.000000  0.000000  0.115720
    H  0.000000  0.748790 -0.462880
    H  0.000000 -0.748790 -0.462880
    symmetry c1
    }

    set basis aug-cc-pvdz
    set omega [589, 633, nm]

    # alpha at 589 and 633 nm, and SHG beta for an incident field at those wavelengths
    properties('scf', properties=['dipole_polarizabilities', 'scf_shg'])

    alpha_589 = variable("SCF DIPOLE POLARIZABILITY TENSOR @ 589NM")   # 3 x 3
    beta_shg = variable("SCF SHG HYPERPOLARIZABILITY TENSOR @ 589NM")  # 3 x 3 x 3

.. table:: Properties computed by the SCFRESPONSE module

   +------------------------------------+----------------------------------------------+-------------------------+
   | Name in ``properties=[...]``       | Quantity                                     | References              |
   +====================================+==============================================+=========================+
   | ``DIPOLE_POLARIZABILITIES``        | :math:`\alpha(-\omega;\omega)` at each       | RHF/UHF, RKS/UKS        |
   | (with |scfresponse__omega| set)    | frequency                                    |                         |
   +------------------------------------+----------------------------------------------+-------------------------+
   | ``SCF_DIPOLE_HYPERPOLARIZABILITY`` | static :math:`\beta(0;0,0)`                  | RHF/UHF                 |
   +------------------------------------+----------------------------------------------+-------------------------+
   | ``SCF_SHG``                        | :math:`\beta(-2\omega;\omega,\omega)`        | RHF/UHF                 |
   +------------------------------------+----------------------------------------------+-------------------------+
   | ``SCF_OPTICAL_RECTIFICATION``      | :math:`\beta(0;\omega,-\omega)`              | RHF/UHF                 |
   +------------------------------------+----------------------------------------------+-------------------------+
   | ``SCF_POCKELS``                    | :math:`\beta(-\omega;\omega,0)`              | RHF/UHF                 |
   +------------------------------------+----------------------------------------------+-------------------------+

Frequencies are given by |scfresponse__omega|, as for the coupled-cluster
response properties: a single value is in atomic units, and with more values
the last element names the unit (``au``, ``nm``, ``ev``, or ``hz``). Without
|scfresponse__omega|, the hyperpolarizabilities are evaluated at zero frequency.
``SCF_DIPOLE_HYPERPOLARIZABILITY`` is the static tensor and cannot be combined
with a nonzero frequency.

Results
~~~~~~~

Results are stored as PSI variables, both globally and on the returned
wavefunction, labeled by the wavelength of the incident field rounded to the
nearest nanometer (``INF NM`` for a static field):

* :psivar:`SCF DIPOLE POLARIZABILITY @ xNM` (isotropic) and
  :psivar:`SCF DIPOLE POLARIZABILITY TENSOR @ xNM` (3 x 3). Every
  polarizability a calculation solves for is stored; for example, an SHG
  calculation also stores :math:`\alpha` at :math:`\omega` and :math:`2\omega`.
* :psivar:`SCF DIPOLE HYPERPOLARIZABILITY TENSOR`,
  :psivar:`SCF SHG HYPERPOLARIZABILITY TENSOR @ xNM`,
  :psivar:`SCF OPTICAL RECTIFICATION HYPERPOLARIZABILITY TENSOR @ xNM`, and
  :psivar:`SCF POCKELS HYPERPOLARIZABILITY TENSOR @ xNM` (3 x 3 x 3). The first
  index is the induced dipole; for the static case
  :math:`\beta_{ijk} = \partial^2 \mu_i / \partial F_j \partial F_k`.

All values are in atomic units.

Capabilities and limitations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* References: RHF/RKS and UHF/UKS, closed and open shell. The restricted or
  unrestricted solver follows the wavefunction, so DFT with the default
  |scf__reference| works. ROHF, ROKS, and CUHF are not supported.
* Functionals: LDA, GGA, global hybrids, and range-separated hybrids (the
  long-range exact exchange is scaled by the functional's :math:`x_\beta`, as
  in the SCF). Meta-GGA functionals, VV10 nonlocal correlation, and
  |scf__wcombine| are not supported.
* Hyperpolarizabilities are available for Hartree--Fock only; DFT would need the
  third-order exchange-correlation kernel.
* Any |globals__scf_type|: the response uses the SCF's own JK object.
* Point-group symmetry is not used: run in C1 symmetry.
* Bases with near-linear dependencies removed (fewer orbitals than basis
  functions) are supported.
* The response equations are converged to |scfresponse__solver_convergence|
  (norm of the change in the response vectors, all three field directions) within
  |scfresponse__solver_maxiter| iterations. Near an excitation energy the
  response is large and convergence is slower.

Keywords
~~~~~~~~

.. include:: autodir_options_c/scfresponse__omega.rst
.. include:: autodir_options_c/scfresponse__solver_convergence.rst
.. include:: autodir_options_c/scfresponse__solver_maxiter.rst
