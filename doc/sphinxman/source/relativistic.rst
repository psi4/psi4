
.. include:: autodoc_abbr_options_c.rst

.. index::
   relativistic

.. _`sec:relativistic`:

Scalar relativistic Hamiltonians
================================

.. codeauthor:: Prakash Verma and Francesco A. Evangelista
.. sectionauthor:: Prakash Verma, Wallace D. Derricotte, and Francesco A. Evangelista

The exact-two-component (X2C) approch is a convenient way to introduce scalar
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
        basis cc-pvdz-dk
        relativistic x2c
    }

    energy('hf')


This computation yields the following result::

  @RHF Final Energy:  -100.10545426415609

   => Energetics <=

    Nuclear Repulsion Energy =              5.1767335622934780
    One-Electron Energy =                -150.7826788086396448
    Two-Electron Energy =                  45.5004909821901009
    Total Energy =                       -100.1054542641560516

while a non-relativistic calculation yields the following energy::

  @RHF Final Energy:  -100.01041683847258

   => Energetics <=

    Nuclear Repulsion Energy =              5.1767335622934780
    One-Electron Energy =                -150.6714586298456027
    Two-Electron Energy =                  45.4843082290795309
    Total Energy =                       -100.0104168384725796


Basis sets options
^^^^^^^^^^^^^^^^^^

The X2C module in |PSIfour| supports different combinations of basis set.
By default, if the input file specifies only the ``basis`` keyword, then the X2C
module will solve the modified Dirac equation in an uncontracted basis and then
recontract the X2C Hamiltonian in the original basis.
Alternatively, the user can use the ``rel_basis`` keyword to specify a different
basis set to solve the modified Dirac equation.::
    set {
        basis cc-pvdz-dk
        rel_basis cc-pvtz-dk
        relativistic x2c
    }
Based on our experience, we recommend the use of a fully uncontracted basis as the computational
basis throughout the computation. If this is not feasible, this recontraction technique should be helpful.
Publications resulting from the use of X2C should cite the following publication: [Verma:2015]_


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
Relativistic kinetic energy ( :math:`T_{\rm X2C}` ) and nuclear-electron interaction potential ( :math:`V_{\rm X2C}` ) is given interms of non-relativisitc kinetic (:math:`T=\hat{p}^2/2`) energy and nuclear-electron interaction potential (:math:`V`), coupling matrix ( :math:`X`) and renormalization matrix ( :math:`R`).  

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
