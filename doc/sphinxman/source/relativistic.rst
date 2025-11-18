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

.. index::
   relativistic

================================
Scalar Relativistic Hamiltonians
================================

.. _`sec:zora`:

Zeroth-order regular approximation (ZORA)
=========================================

.. codeauthor:: Nathan Gillispie and Daniel R. Nascimento
.. sectionauthor:: Nathan Gillispie and Daniel R. Nascimento

*Source*: :source:`zora.cc <psi4/src/psi4/libmints/zora.cc>` :source:`zora.h <psi4/src/psi4/libmints/zora.h>`

ZORA is a perturbative approximation of the full relativistic Hamiltonian for
DFT and wavefunction-based methods. It is commonly used to provide a more
accurate total energy in DFT calculations.

When ZORA is used in |PSIfour|, it creates a scalar relativistic kinetic
energy integral to be used in the SCF procedure. It has been tested
for HF and DFT references with the driver methods ``energy`` and ``optimize``.
Our implementation employs a grid-based scheme, therefore analytic energy
gradients are not supported.

Equations and implementation details are based on a paper by Pak, Dada and
Nascimento [Pak:2025:094110]_.

Usage
^^^^^

To use the ZORA Hamiltonian with default grid settings, use option
|globals__relativistic| with ``ZORA``. To change the number of grid points, use
the |globals__zora_radial_points| and |globals__zora_spherical_points| options. ::

    set {
        reference rhf
        scf_type pk
        relativistic zora
        zora_radial_points 160
        zora_spherical_points 1202
    }

.. note::
   The number of spherical points must be a Lebedev number.
   See :ref:`Grid Selection <sec:grid-selection>` for a list of all options.

It may be useful to compute the non-relativistic kinetic integral using the
grid points to compare against the analytic kinetic integrals. To do this, use
the |globals__zora_nr_debug| option. ::

    set {
        relativistic zora
        zora_nr_debug true
    }

See :ref:`sec:relativistic-keywords` for more options.

Theory
^^^^^^

In short, the FW-transformed Dirac Hamiltonian is perturbatively expanded with
respect to an expression (:math:`E/(2mc^2-V)`) involving the potential
:math:`V`. To the zeroth-order, we get the ZORA Hamiltonian. Below is the ZORA
Kohn-Sham (ZKS) Hamiltonian, however, results are analagous for other methods.

.. math::
   :label: zks
   
   \hat h^\text{ZKS}=\frac{\mathbf{p}^2}{2}+\mathbf{p}\left(\frac{\kappa -1}{2}\right)\mathbf{p}+\frac{\kappa^2}{4c^2}\mathbf{\sigma}\cdot\left(\nabla v^\text{KS}\times \mathbf{p}\right)+v^\text{KS}

:math:`\mathbf{p}` is the momentum operator, :math:`\mathbf{\sigma}` are
the Pauli matrices, :math:`\kappa=[1-v^\text{KS}/(2c^2)]^{-1}`, and
:math:`v^\text{KS}` is the usual Kohn-Sham potential. In this form, the terms
of the ZORA Hamiltonian are separated into the classical kinetic energy,
scalar relativistic correction to the kinetic energy, spin-orbit, and KS
potential terms, respectively. Were this implemented as is, the relativistic
terms would depend on the KS potential, which is impracticable due to
gauge-invariance and convergence issues.

The procedure given by Van W\ |u_dots|\ llen [vanWullen:1998:392]_ avoids these
problems by replacing :math:`v^\text{KS}` in the relativistic terms of
:eq:`zks` with an effective potential :math:`v_\text{eff}(\mathbf{r})`.
This is the potential experienced at a point due to nuclear attraction and
electron repulsion from a model potential that reproduces the nuclear cusp on 
each atom. The model potential comes from a model basis of Gaussian-type
orbitals, whose values were obtained from NWChem [NWChem:2020]_. This means that
geometry is the only thing affecting :math:`v_\text{eff}(\mathbf{r})`.

The ZORA scalar relativistic kinetic integral :math:`T^\text{SR}` is the
first two terms of :eq:`zks` in the atomic orbital basis. Once
:math:`v_\text{eff}(\mathbf{r})` is computed on a grid, :math:`T^\text{SR}` can
be calculated as

.. math:: T_{\mu\nu}^\text{SR}=\int d\mathbf{r}^3 \frac{c^2}{2c^2-v_\text{eff}(\mathbf{r})}\nabla\chi_\mu^\dagger(\mathbf{r}) \cdot\nabla\chi_\nu(\mathbf{r}),

given atomic orbitals :math:`\chi(\mathrm{r})`. Notice that as
:math:`v_\text{eff}\rightarrow 0`, :math:`T^\text{SR}` simply becomes the
non-relativistic kinetic energy. This is the basis behind the
|globals__zora_nr_debug| option. Because :math:`T^\text{SR}` is only evaluated
once before the SCF procedure, don't cheap out on the grid! In practice, the
time spent computing the ZORA Hamiltonian is relatively negligible.

Limitations
^^^^^^^^^^^

* Spin-orbit coupling effects are not available because they require a
  complex generalized SCF procedure.

* ZORA theory allows for adjustment of the molecular orbital energies.
  This provides an important correction for linear response calculations.
  Currently, this is not implemented in |PSIfour|.

* Analytic energy gradients are not available.

.. _`sec:relativistic-keywords`:

Keywords
^^^^^^^^

.. include:: autodir_options_c/globals__relativistic.rst
.. include:: autodir_options_c/globals__zora_radial_points.rst
.. include:: autodir_options_c/globals__zora_spherical_points.rst
.. include:: autodir_options_c/globals__zora_pruning_scheme.rst
.. include:: autodir_options_c/globals__zora_basis_tolerance.rst
.. include:: autodir_options_c/globals__zora_nr_debug.rst

.. _`sec:x2c`:

Exact two-component (X2C)
=========================

.. codeauthor:: Prakash Verma and Francesco A. Evangelista
.. sectionauthor:: Prakash Verma, Wallace D. Derricotte, and Francesco A. Evangelista

The X2C approach is a convenient way to introduce scalar
relativistic effects in DFT and wave function-based methods.
|PSIfour| implements the spin-free one-electron version of X2C, which produces
a modified one-electron Hamiltonian :math:`H_{\rm X2C}`:

.. math:: H_{\rm X2C} = T_{\rm X2C} + V_{\rm X2C}

that is a sum of a kinetic energy (:math:`T_{\rm X2C}`) and potential energy
(:math:`V_{\rm X2C}`) operator.
Our implementation is equivalent to the one reported by Cheng and Gauss [Cheng:084114]_.
X2C calculations require the use of special (alternatively fully uncontracted) basis sets designed for relativistic
calculations.  Common choices include the Dunning Douglass--Kroll basis sets
(cc-pVXZ-DK, cc-pCVXZ-DK, cc-pwCVXZ-DK) and Roos' ANO basis sets.

.. note:: See also :ref:`sec:DKH` for another relativistic Hamiltonian.

A First Example
^^^^^^^^^^^^^^^

The following is a simple input that will perform a Hartree--Fock calculation
using the X2C Hamiltonian. ::

    molecule {
      H
      F 1 0.92
    }

    set {
        scf_type pk
        basis cc-pvdz
        relativistic x2c
    }

    energy('hf')

This computation yields the following result::

  @RHF Final Energy:  -100.10007984692388

   => Energetics <=

    Nuclear Repulsion Energy =              5.1767335622934780
    One-Electron Energy =                -150.7611816259664579
    Two-Electron Energy =                  45.4843682167491039
    Total Energy =                       -100.1000798469238902

while a non-relativistic calculation yields the following energy::

  @RHF Final Energy:  -100.01928891411315

   => Energetics <=

    Nuclear Repulsion Energy =              5.1767335622934780
    One-Electron Energy =                -150.6645256529074572
    Two-Electron Energy =                  45.4685031765008461
    Total Energy =                       -100.0192889141131474

Basis sets options
^^^^^^^^^^^^^^^^^^

The X2C module in |PSIfour| supports different combinations of basis set.
By default, if the input file specifies only |mints__basis|, then the X2C
module will solve the modified Dirac equation in an uncontracted basis and then
recontract the X2C Hamiltonian in the original basis.
Alternatively, the user can use |globals__basis_relativistic| to specify a different
basis set to solve the modified Dirac equation. ::

    set {
        basis cc-pvdz-dk
        basis_relativistic cc-pvtz-dk
        relativistic x2c
    }

It is recommended that when employing the X2C relativistic Hamiltonian, that you use a fully
decontracted basis set. This can be done simply in the input by adding "-decon" to the 
name of the primary basis you want to use for the calculation as detailed in 
:ref:`Decontracted Basis Sets <sec:basisDecontracted>`. Publications resulting from the use 
of X2C should cite the following publication: [Verma:2015]_


Theory
^^^^^^

X2C is based on exact decoupling of 
positive-energy ( :math:`h^{FW}_{\rm ++}`
) and negative-energy (:math:`h^{FW}_{\rm --}` )
blocks of the Dirac Hamiltonian (:math:`h^{D}`). 

.. math:: 
   U^\dagger h^{\rm D} U = 
   U^\dagger
   \begin{pmatrix}
   h_{LL} & h_{LS} \\
   h_{SL} & h_{SS}
   \end{pmatrix}
   U
   =
  \begin{pmatrix}
  h^{\rm FW}_{++} & 0 \\
  0 & h^{\rm FW}_{--}
  \end{pmatrix}

The transformation ( :math:`U` ) is  obtained from the solutions of the Dirac equation in kinetically balanced basis [Kutzelnigg:1984]_ treatment. 
In the X2C treatment, the positive-energy block of the Hamiltonian ( :math:`h^{FW}_{\rm ++}` )
is given by the sum
of a transformed kinetic (:math:`T_{\rm X2C}`) and potential energy ( :math:`V_{\rm X2C}` ) contribution.
Relativistic kinetic energy ( :math:`T_{\rm X2C}` ) and nuclear-electron interaction potential ( :math:`V_{\rm X2C}` ) is given in terms of non-relativisitc kinetic (:math:`T=\hat{p}^2/2`) energy and nuclear-electron interaction potential (:math:`V`), coupling matrix ( :math:`X`) and renormalization matrix ( :math:`R`).  

.. math::
  T_{\rm X2C} = R^{\dagger} (TX +  {X}^{\dagger}T - {X}^{\dagger}TX ) R 

.. math::
  V_{\rm X2C} = R^{\dagger}(V + \frac{1}{4c^2} X^{\dagger}W^{\text{SF}}X) R

The coupling matrix ( :math:`{X} = C^{S} (C^{L})^{-1}` ) is obtained from the large (:math:`C^{\rm L}`) and small (:math:`C^{\rm S}`) components of the :math:`N` positive energy solutions of the Dirac equation.
The renormalization matrix 
:math:`{R}=S^{-1/2}(S^{-1/2}\tilde{S}S^{-1/2})^{-1/2}S^{1/2}`,
depends on the modified overlap matrix
:math:`\tilde{S}=S+\frac{1}{2c^2}X^{\dagger}TX`. The integrals :math:`W^{\rm SF}_{\mu\nu} = \langle {\chi_\mu} | \hat{p}\cdot (\hat{V}\hat{p}) |{\chi_\nu}\rangle` can be easily computed as derivatives of the nuclear-electron attraction integrals with respect to nuclear coordinates.
Existing nonrelativistic electronic structure code can be extended to include scalar relativistic effects
treated with the X2C method by replacing nonrelativistic kinetic and potential energy with the corresponding
X2C operators :math:`T_{X2C}` and :math:`V_{X2C}`. It is important to note that fully uncontracted basis in needed for the construction of X2C Hamiltonian as Foldy-Wouthuysen (FW [FW:1950]_) transformation is obtained in kinetically balance basis.

Keywords
^^^^^^^^

.. include:: autodir_options_c/globals__relativistic.rst
.. include:: autodir_options_c/globals__basis_relativistic.rst

