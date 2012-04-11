.. include:: autodoc_abbr_options_c.rst

.. index::
     single: multireference
     single: Mk-MRCC

PSIMRCC implementation of Mk-MRCC theory
==========================================================================
.. codeauthor:: Francesco A. Evangelista, and Andrew C. Simmonett 
.. sectionauthor:: Alexander E. Vaughn

State-specific Multireference coupled cluster theories provide highly accurate energies and properties of electronic states that require a multiconfigurational zeroth-order wavefunction.  The PSIMRCC module contained in |PSIfour| implements the state-specific multireference coupled-cluster approach of Mukherjee and co-workers (Mk-MRCC). This method is implemented and shown to be a powerful tool in `F. A. Evangelista, W. D. Allen, and H. F. Schaefer, III, J. Chem. Phys. 125 (2006)` and `F. A. Evangelista, A. C. Simmonett, W. D. Allen, H. F. Schaefer, III, and J. Gauss, J. Chem. Phys. 128, 124104 (2008)`.    
Mk-MRCC is based on the Jeziorski-Monkhorst ansatz [`B. Jeziorski and H. J. Monkhorst, Phys. Rev. A 24, 1668 (1981)`] for the wavefunction, :math: `\Psi` 
.. math:: \left| \Psi \right \rangle = \sum_\mu^d e^{\hat{T}^\mu} \left| \Phi_\mu \right\rangle c_\mu \, \text{,}
where $\Phi_\mu$ are the reference determinants, :math:`\hat{T}^\mu` are reference-specific excitation operators, and :math:`c_\mu` are expansion coefficients obtained through diagonalization of the Mk-MRCC effective Hamiltonian matrix that allows the various reference determinants to interact. As an example of how this works the Mk-MRCCSD excitation operators for each reference is contracted two-body terms
.. math:: \hat{T}^\mu = \hat{T}^\mu_1 + \hat{T}^\mu_2
where
.. math:: \hat{T}^\mu_1 = \sum_i^{\textrm{occ}(\mu)} \sum_a^{\textrm{vir}(\mu)} t_i^a (\mu) \hat{a}^\dagger_a \hat{a}_i
and
.. math:: \hat{T}^\mu_2 =\frac{1}{4} \sum_i^{\textrm{occ}(\mu)} \sum_a^{\textrm{vir}(\mu)} t_{ij}^{ab} (\mu) \hat{a}^\dagger_b \hat{a}_j \hat{a}^\dagger_a \hat{a}_i  
The Mk-MRCC energy is a chosen eigenvalue of the effective Hamiltonian, :math: `\textrm{H}^{eff}_{\mu \nu}`
.. math:: \sum_\nu \textrm{H}^{eff}_{\mu \nu} c_\nu =E c_\nu
where 
.. math:: \textrm{H}^{eff}_{\mu \nu} = \left \langle \Phi_\mu \right | \hat{H}e^{\hat{T}^\nu} \left | \Phi_\nu \right \rangle \, \textrm{.}
|PSIfour| currently has Mk-MRCC with singles and doubles [Mk-MRCCSD], single and doubles with perturbative triples[Mk-MRCCSD(T)]. A companion perturbation method(Mk-MRPT2) has been developed based on the Mukherjee formalisim as shown in `F. A. Evangelista, A. C. Simmonett, H. F. Schaefer, III, D. Mukherjee, and W. D. Allen, Phys. Chem. Chem. Phys. 11, 4728 (2009)`. 

A Simple Example
_________________ 
..
   molecule o2 {
      0 3
      O
      O 1 2.265122720724
      units au
   }
   set {
      basis cc-pvtz
   }
   set mcscf {
      reference       rohf
      docc            [3,0,0,0,0,2,1,1]      # Doubly occupied MOs
      socc            [0,0,1,1,0,0,0,0]      # Singly occupied MOs
   }
   set psimrcc {
      corr_wfn        ccsd                   # Do Mk-MRCCSD 
      frozen_docc     [1,0,0,0,0,1,0,0]      # Frozen MOs
      restricted_docc [2,0,0,0,0,1,1,1]      # Doubly occupied MOs
      active          [0,0,1,1,0,0,0,0]      # Active MOs
      frozen_uocc     [0,0,0,0,0,0,0,0]      # Frozen virtual MOs
      corr_multp      1                      # Select the Ms = 0 component
      root            1
      wfn_sym         B1g                    # Select the B1g state
   }
   energy('psimrcc')

The |psimrcc__corr_wfn| allows you to select one of three methods Mk-MRPT2 [``PT2``], Mk-MRCCSD [``CCSD``], or Mk-MRCCSD(T) [``CCSD_T``].   The |psimrcc__corr_multp| option allows you to select the Slater determinants with a particular :math:`M_s` value. The |psimrcc_wfn_sym| keyword is neccesary if you do not want to compute the energy of the all-symmetric state. The |psimrcc__root| option may be used to follow different roots of the effective Hamiltonian. A value of 1 instructs PSIMRCC to follow the solution with the lowest energy given a certain set of determinants.          

Orbital ordering and selection of the model space
_________________________________________________
The reference determinants :math:`\Phi_\mu` are specified in PSIMRCC via occupational numbers. PSIMRCC requires that four arrays be specified for this purpose.

- Frozen doubly occupied orbitals (|psimrcc__frozen-doccFROZEN-DOCC|) are doubly occupied in each reference determinant and are not correlated in the MRCC procedure.
- Doubly occupied orbitals(|psimrcc__retricted-docc|}) are doubly occupied in each reference determinant and are correlated in the MRCC procedure.
- Active orbitals(|psimrcc_active|) are partially occupied in each reference determinant.
- Frozen virtual orbitals (|psimrcc_frozen_uocc) are unoccupied in all reference determinants and are excluded from the correlated wave function.

The model space is selected by considering all possible occupations of the electrons among the orbitals in the active space that result in determinants with the correct symmetry (|psimrcc__wfn_sym|) and the correct :math:`\textrm{M}_s` value specified by the keyword MULTP. Note that this does not consider the multiplicity of the wavefunction. Thus, in order to obtain the wavefunction with a set of :math: `\textrm{M}_s = 0` reference determinants for an open-shell system you should request a MULTP of 1 within the PSIMRCC module, and select the root of the effective Hamiltonian that corresponds to the state of interest. In addition, the |psimrcc__wfn_sym| keyword needs to be specified otherwise the wavefunction belonging to the all-symmetric irrep will be selected. In addition, it should be noted that for an open-shell singlet based on two $\textrm{M}_s = 0$ determinants the eigenvector is [:math: `\frac{1}{\sqrt{2}}\text{,}\frac{1}{\sqrt{2}}`], which corresponds to a wavefunction of the following form:
.. math:: \frac{1}{\sqrt{2}} \left( \chi_1 \alpha (1) \chi_2 \beta (2) + \chi_2 \alpha(1) \chi_1 \beta (2) \right)

