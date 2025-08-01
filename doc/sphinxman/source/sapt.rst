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
   single: SAPT
   pair: SAPT; theory

.. _`sec:sapt`:

SAPT: Symmetry-Adapted Perturbation Theory
==========================================

.. codeauthor:: Edward G. Hohenstein, Rob M. Parrish, J\ |e_acute|\ r\ |o_circumflex|\ me F. Gonthier, and Daniel. G. A. Smith
.. sectionauthor:: Edward G. Hohenstein and J\ |e_acute|\ r\ |o_circumflex|\ me F. Gonthier

*Module:* :ref:`Keywords <apdx:sapt>`, :ref:`PSI Variables <apdx:sapt_psivar>`, :source:`LIBSAPT_SOLVER <psi4/src/psi4/libsapt_solver>`

.. warning:: In rare cases with systems having a high degree of symmetry, 
   |Psifour| gives (very obviously) wrong answers for SAPT computations 
   when the specification is in Z-matrix format. Use a Cartesian representation 
   to avoid this problem.

.. caution:: In early versions (notably |Psifour| alpha circa 2011
   and before), frozen core was implemented incompletely and for
   only selected terms. Comparisons with papers published using early
   |PSIfour| SAPT code may show discrepancies of 0.01-0.10 kcal/mol in
   individual terms, particularly :math:`E_{exch}^{(11)}` and :math:`E_{exch}^{(12)}`.

.. caution:: January 28th 2016, the default for all NAT_ORBS options
   was changed to true. Hence the code now by default uses natural
   orbital truncation to speed up the evaluation of energy terms
   wherever possible, according to literature recommendations.
   In early July 2016, some total SAPT energy psivars were renamed.

.. caution:: February 7, 2020, a missing term in :math:`E^{(30)}_{ind}` was added, causing
   possible discrepancies with prior versions of the code on the order of
   0.01 kcal/mol. See https://github.com/psi4/psi4/issues/1677

.. caution:: August 2021, the number of frozen core orbitals used in the dMP2 computations
   is now standardized. Specifically, we now rigorously enforce that the number of core orbitals 
   frozen in dimer computations is equal to the sum of frozen orbitals of each monomer. Prior to
   this, a discrepency between these values was possible when one of the monomers was (exclusively) 
   a charged alkali metal. 


.. _`sec:saptfitA`:


.. caution:: May 2022 c. v1.6, the default for |sapt__df_basis_elst|
   changed from the value of |sapt__df_basis_sapt| (which itself
   defaults to the RI of the orbital basis) to the JKFIT of the orbital
   basis. This affects SAPT0 and sSAPT0 computed with the :ref:`SAPT
   module<sec:sapt>` (the default code for ``energy("sapt0")`` that
   can also compute higher-order SAPT). Electostatics, exchange,
   and induction terms for SAPT0 and sSAPT0 accessed through
   ``energy("sapt0")`` or ``energy("ssapt0")`` change; the dispersion
   term does not change. The SAPT0 and sSAPT0 terms accessed as
   subsidiary calculations of higher-order SAPT do not change; that is,
   the :ref:`SAPT module<sec:sapt>` breaks the consistency of its SAPT0
   results. The reasoning and reward behind this change is that the JKFIT
   basis better describes the physics (see :ref:`fitting changes <sec:saptfitB>` ) and the
   default SAPT0 results from the :ref:`SAPT module<sec:sapt>` are now
   consistent with those from the :ref:`FISAPT module<sec:fisapt>` and
   the sapt(dft) module. See :srcsample:`sapt-compare` for an example.
   To reproduce former behavior, set |sapt__df_basis_elst| to the
   orbital basis set's RI auxiliary basis.

Symmetry-adapted perturbation theory (SAPT) provides a means of directly
computing the noncovalent interaction between two molecules, that is, the
interaction energy is determined without computing the total energy of the
monomers or dimer. In addition, SAPT provides a decomposition of the
interaction energy into physically meaningful components: *i.e.*,
electrostatic, exchange, induction, and dispersion terms. In SAPT, the 
Hamiltonian of the dimer is partitioned into contributions from each 
monomer and the interaction.

.. math:: H=F_A+W_A+F_B+W_B+V

Here, the Hamiltonian is written as a sum of the usual monomer Fock
operators, :math:`F`, the fluctuation potential of each monomer, :math:`W`, and the
interaction potential, :math:`V`. The monomer Fock operators, :math:`F_A+F_B`, are
treated as the zeroth-order Hamiltonian and the interaction energy is
evaluated through a perturbative expansion of :math:`V`, :math:`W_A`, and :math:`W_B`. 
Through first-order in :math:`V`, electrostatic and exchange interactions are
included; induction and dispersion first appear at second-order in :math:`V`. For
a complete description of SAPT, the reader is referred to the excellent
review by Jeziorski, Moszynski, and Szalewicz [Jeziorski:1994:1887]_.

Several truncations of the closed-shell SAPT expansion are available in the SAPT
module of |PSIfour|. The simplest truncation of SAPT is denoted SAPT0
and defined in Eq. :eq:`SAPT0`.

.. math:: E_{SAPT0} = E_{elst}^{(10)} + E_{exch}^{(10)} + E_{ind,resp}^{(20)} +
   E_{exch-ind,resp}^{(20)} + E_{disp}^{(20)} + E_{exch-disp}^{(20)} + \delta_{HF}^{(2)}
   :label: SAPT0

In this notation, :math:`E^{(vw)}` defines the order in :math:`V` and in :math:`W_A+W_B`; the
subscript, :math:`resp`, indicates that orbital relaxation effects are included.

.. math:: E_{SAPT2} = E_{SAPT0} + E_{elst,resp}^{(12)} + E_{exch}^{(11)} +
   E_{exch}^{(12)} + \; ^{t}\!E_{ind}^{(22)} + \; ^{t}\!E_{exch-ind}^{(22)}
   :label: SAPT2

.. math:: E_{SAPT2+} = E_{SAPT2} + E_{disp}^{(21)} + E_{disp}^{(22)}
   :label: SAPT2p

.. math:: E_{SAPT2+(3)} = E_{SAPT2+} + E_{elst,resp}^{(13)} + E_{disp}^{(30)}
   :label: SAPT2pparen3

.. math:: E_{SAPT2+3} = E_{SAPT2+(3)}
   + E_{exch-ind}^{(30)} + E_{ind,resp}^{(30)}
   + E_{exch-disp}^{(30)} + E_{ind-disp}^{(30)} + E_{exch-ind-disp}^{(30)}
   - \delta_{HF}^{(2)} + \delta_{HF}^{(3)}
   :label: SAPT2p3

For SAPT(DFT), the SAPT expansion is truncated at the same level of SAPT0, but the zeroth-order Hamiltonian is chosen to be :math:`K_A+K_B`, the monomer Kohn-Sham operators, instead of the Fock operators. The dispersion term needs to be computed with orbital relaxation for the result to be reasonable, and this is possible by computing dispersion energy through coupled frequency-dependent density susceptibility (FDDS). The exchange-dispersion term is estimated by scaling the uncoupled exchange-dispersion energy. 

.. math:: E_{SAPT(DFT)} =  E_{elst}^{(10)} + E_{exch}^{(10)} + E_{ind,resp}^{(20)} +
   E_{exch-ind,resp}^{(20)} + E_{disp,resp}^{(20)} + E_{exch-disp,resp}^{(20)}(est.) + \delta_{HF}^{(2)}
   :label: SAPTDFT
.. math:: E_{disp,resp}^{(20)} = -\frac{1}{2\pi}\int_0^\infty d\omega \int d\mathbf{r}_A d\mathbf{r}'_A d\mathbf{r}_B d\mathbf{r}'_B \frac{1}{\left|\mathbf{r}_A-\mathbf{r}_B\right|}\frac{1}{\left|\mathbf{r}'_A-\mathbf{r}'_B\right|}\chi^A_{coup}\left(\mathbf{r}_A,\mathbf{r}'_A|i\omega\right)\chi^B_{coup}\left(\mathbf{r}_B,\mathbf{r}'_B|i\omega\right)

.. math:: \chi_{coup}\left(\mathbf{r},\mathbf{r}'|\omega\right) = \chi_0\left(\mathbf{r},\mathbf{r}'|\omega\right) + \int d\omega d\mathbf{r}_1 d\mathbf{r}_2 \chi_0\left(\mathbf{r},\mathbf{r}_1|\omega\right)\left[\frac{1}{r_{12}} + f_{xc}\left(\mathbf{r}_1,\mathbf{r}_2,\omega\right)\right]\chi\left(\mathbf{r}_2,\mathbf{r}'|\omega\right)


The :math:`\delta_{HF}^{(2)}` and :math:`\delta_{HF}^{(3)}` terms take into
account higher-order induction effects and are included in the definition
of SAPT terms. They are computed from the Hartree--Fock supermolecular interaction energy
:math:`E_{int}^{HF}` and are only available in dimer-centered basis SAPT
computations, which is the default (see below for monomer-centered basis 
computations). They are defined by:

.. math:: \delta_{HF}^{(2)} = E_{int}^{HF} - (E_{elst}^{(10)} + E_{exch}^{(10)} 
          + E_{ind,resp}^{(20)} + E_{exch-ind,resp}^{(20)})
          :label: dHF2

.. math:: \delta_{HF}^{(3)} = \delta_{HF}^{(2)} - (E_{exch-ind}^{(30)} 
          + E_{ind,resp}^{(30)})
          :label: dHF3

Additionally, high-order coupling between induction and dispersion can be 
extracted from the supermolecular MP2 interaction energy:

.. math:: \delta_{MP2}^{(2)} = E_{int}^{MP2, corr} - (E_{elst}^{(12)} +
          E_{exch}^{(11)} + E_{exch}^{(12)} + \; ^{t}\!E_{ind}^{(22)}
          + \; ^{t}\!E_{exch-ind}^{(22)} + E_{disp}^{(20)} + E_{exch-disp}^{(20)})

.. math:: \delta_{MP2}^{(3)} = \delta_{MP2}^{(2)} - (E_{ind-disp}^{(30)} + E_{exch-ind-disp}^{(30)})

where :math:`E_{int}^{MP2, corr}` is the correlation part of the supermolecular MP2 
interaction energy. :math:`\delta_{MP2}^{(2)}` and :math:`\delta_{MP2}^{(3)}` also improve the 
description of electrostatically dominated complexes. :math:`\delta_{MP2}^{(2)}`
can be applied to SAPT2+ or SAPT2+(3) energies whereas :math:`\delta_{MP2}^{(3)}` 
should be applied to SAPT2+3 energies.

A thorough analysis of the performance of these truncations of closed-shell SAPT can be
found in a review by Hohenstein and Sherrill [Hohenstein:2012:WIREs]_,
and a systematic study of the accuracy of these truncations (with and 
without an improved CCD treatment of dispersion) using different basis sets
is reported in [Parker:2014:094106]_.

The closed-shell SAPT module relies entirely on the density-fitting approximation
of the two-electron integrals. The factorization of the SAPT energy
expressions, as implemented in |PSIfour|, assumes the use of density-fitted
two-electron integrals, therefore, the closed-shell SAPT module cannot be run with
exact integrals. In practice, we have found that the density-fitting
approximation introduces negligible errors into the SAPT energy 
(often less than 0.01 kcal/mol for small dimers) and greatly
improves efficiency. 

The latest addition to the SAPT code is the SAPT0 method for open-shell 
monomers [Gonthier:2016:134106]_. This code is available for both exact and density fitted integrals, 
except for the dispersion terms which implementation relies on a density fitting
factorization. Both ``UHF`` and ``ROHF`` |scf__reference| can be used, but coupled
induction computations are currently not supported with ``ROHF``. This means that orbital
relaxation is not included for ``ROHF`` and the uncoupled induction term is computed instead.
If both monomers are open-shell, their coupling is assumed to be high spin, i.e.
two doublets would interact to form a triplet.

The *S*\ :superscript:`2` approximation and scaling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All exchange terms in SAPT arise from the antisymmetrization
of the wavefunctions of monomers A and B. Taking into account exchange of all possible
electron pairs between the two monomers yields to complicated formulae.
For this reason, exchange terms are often evaluated in the :math:`S^{2}`
approximation, that can be interpreted as the exchange of a single electron 
pair between monomers.

The :math:`S^{2}` approximation is usually pretty good, but may 
break down for short intermolecular distance, particularly in high-order
terms. To compensate for these deviations, in 2014 Parker et al. [Parker:2014:094106]_ 
recommended to scale all :math:`S^{2}` approximated exchange terms by the ratio:

.. math:: p_{EX}(\alpha) = \left( \frac{E_{exch}^{(10)}}{E_{exch}^{(10)}(S^{2})} \right)^{\alpha}

and recommended using the ratio with exponent :math:`\alpha = 1`. To obtain SAPT energies with this scaling,
simply set the keyword ``exch_scale_alpha true``. Alternatively, another value for :math:`\alpha`
can be specified by setting |sapt__exch_scale_alpha| to a value. For example, ::

  set exch_scale_alpha 1.0

will set :math:`\alpha = 1.0` and scale exchange energies with :math:`p_{EX}(1.0)`.

However, as pointed out by Schaffer and Jansen [Schaffer:2013:2570]_
in the context of DFT-based SAPT, the ratios :math:`E_{exch}^{(1)}(S^{2})/E_{exch}^{(1)}`,
:math:`E_{\rm exch-ind}^{(2)}(S^2)/E_{\rm exch-ind}^{(2)}`, and 
:math:`E_{\rm exch-disp}^{(2)}(S^2)/E_{\rm exch-disp}^{(2)}` are not very similar to
each other.  Hence, in 2016 universal scaling of all :math:`S^{2}` approximated terms
was turned off by default.

Finally, in 2014 Parker et al. [Parker:2014:094106]_ empirically discovered that 
SAPT0 energies for van der Waals dimers with close contacts benefit from a slightly modified 
recipe that involves an empirically adjusted exponent :math:`\alpha = 3.0`. 
To distinguish it from its unscaled counterpart, this energy is denoted sSAPT0 (see [Parker:2014:094106]_).

.. math:: E_{sSAPT0} = E_{elst}^{(10)} + E_{exch}^{(10)} + E_{ind,resp}^{(20)} +
   p_{EX}(3.0) E_{exch-ind,resp}^{(20)} + E_{disp}^{(20)} + p_{EX}(3.0) E_{exch-disp}^{(20)} 
   + \delta_{HF}^{(2)}
   :label: sSAPT0

where :math:`\delta_{HF}^{(2)}` is computed *without* any scaling. Please note that 
sSAPT0 is thus not the same as requesting ``exch_scale_alpha 3.0``, and that the
scaling is automatically performed by requesting ``energy('ssapt0')``.


A First Example
^^^^^^^^^^^^^^^

The following is the simplest possible input that will perform all
available SAPT computations (normally, you would pick one of these methods,
not all of them). ::

	molecule water_dimer {
	     0 1
	     O  -1.551007  -0.114520   0.000000
	     H  -1.934259   0.762503   0.000000
	     H  -0.599677   0.040712   0.000000
	     --
	     0 1
	     O   1.350625   0.111469   0.000000
	     H   1.680398  -0.373741  -0.758561
	     H   1.680398  -0.373741   0.758561
	
	     units angstrom
	     no_reorient
	     symmetry c1
	}
	
	set basis aug-cc-pvdz
	
	energy('sapt0')
	energy('sapt2')
	energy('sapt2+')
	energy('sapt2+(3)')
	energy('sapt2+3')

The SAPT module uses the standard |PSIfour| partitioning of the dimer
into monomers. SAPT does not use spatial symmetry and needs the geometry
of the system to remain fixed throughout monomer and dimer calculations.
These requirements are imposed whenever a SAPT calculation is requested
but can also be set explicitly with the ``no_reorient`` and ``symmetry
c1`` molecule keywords, as in the example above. As a reminder, only
SAPT0 can handle the interaction of both closed- and open-shell monomers.
Higher-order SAPT and SAPT(DFT) are currently available only 
for computation of interactions between
closed-shell singlets.  The SAPT codes in |PSIfour| have been written
to utilize density fitting, which is much faster than using conventional
4-index electron repulsion integrals.  This happens automatically and
does not require any additional keywords to be specified (both the 
SAPT computations and the underlying Hartree-Fock computations will
utilize density fitting).


For SAPT(DFT), the user will need to manually specify the GRAC shift of 
both monomers, defined by the difference of ionization potential (IP) and
Kohn-Sham HOMO energy. The ionization potential data for many common molecules 
is available in NIST Chemistry Webbook. Alternatively, one can estimate the
ionization potential of molecule by computing the energy difference between
the molecule as given, and the molecule after one electron has been removed
(e.g., the energy difference between a neutral molecule and its cation).

The values of GRAC shifts should be given in Hartree. For the example above,
the GRAC shift value of both molecules are calculated to be 0.1307 (using experimental
IP from CCCBDB), and the user would specify them using the following keywords::


    set globals {
        sapt_dft_grac_shift_a 0.1307
        sapt_dft_grac_shift_b 0.1307 
    }

A complete, minimal example of a SAPT(DFT) computation is given below. ::

    molecule {
        0 1
        O  -1.551007  -0.114520   0.000000
        H  -1.934259   0.762503   0.000000
        H  -0.599677   0.040712   0.000000
        --
        0 1
        O   1.350625   0.111469   0.000000
        H   1.680398  -0.373741  -0.758561
        H   1.680398  -0.373741   0.758561

        units angstrom
    }
    
    set {
        basis                  aug-cc-pVDZ
        sapt_dft_grac_shift_a  0.1307
        sapt_dft_grac_shift_b  0.1307
    }
    
    energy('sapt(dft)')


As already mentioned above, the SAPT0 module for open-shell cases can also
use exact integrals for all terms except for dispersion. In practice,
density fitting is considerably faster and introduces negligible errors, thus
it is the preferred method for open-shell cases as well.
Below, you can find a minimum example of open-shell SAPT0 computation. ::

    molecule {
        0 1
        O 0.000000  0.000000  6.000000
        H 0.000000  1.431500  4.890600
        H 0.000000 -1.431500  4.890600
        --
        0 2
        O 0.000000  0.000000  0.000000
        O 0.000000  2.503900  0.000000
        H 0.000000 -0.424700 -1.839500
        units bohr
        symmetry c1
        no_reorient
        no_com
    }
    
    set {
        reference    uhf
        scf_type     df
        basis        aug-cc-pVDZ
    }
    
    energy('sapt0')
    
    

|scf__reference| needs to be ``UHF``  or ``ROHF`` for the open-shell computation to proceed.

Advanced example
^^^^^^^^^^^^^^^^

Open-shell computations can be difficult to converge in certain cases, thus you may
want to have more control over the SCF procedure. You have the option of 
doing the driver job in the input file, by performing the dimer and monomer computations
yourself. In the example below, we do a stability analysis for the open-shell monomer only ::

  molecule {
       0 2
       O 0.000000  0.000000  0.000000
       O 0.000000  2.503900  0.000000
       H 0.000000 -0.424700 -1.839500
       --
       0 1
       O 0.000000  0.000000  6.000000
       H 0.000000  1.431500  4.890600
       H 0.000000 -1.431500  4.890600
       units bohr
       symmetry c1
       no_reorient
       no_com
  }
  
  dimer = psi4.get_active_molecule()  
  
  set {
  reference uhf
  scf_type     df
  basis         cc-pVDZ
  df_basis_sapt cc-pVDZ-ri
  guess sad
  }
  
  dimer = psi4.get_active_molecule()
  
  set df_ints_io save
  psi4.IO.set_default_namespace('dimer')
  Edim, wfn_dimer = energy('scf',molecule=dimer,return_wfn=True)
  set df_ints_io load
  
  monomerA = dimer.extract_subsets(1,2)
  psi4.IO.change_file_namespace(97, 'dimer', 'monomerA')
  psi4.IO.set_default_namespace('monomerA')
  set {
  stability_analysis follow
  }
  EmonA, wfn_monA = energy('scf',molecule=monomerA,return_wfn=True)
  
  monomerB = dimer.extract_subsets(2,1)
  psi4.IO.change_file_namespace(97, 'monomerA', 'monomerB')
  psi4.IO.set_default_namespace('monomerB')
  set {
  stability_analysis none
  }
  EmonB, wfn_monB = energy('scf',molecule=monomerB,return_wfn=True)
  
  psi4.IO.change_file_namespace(97, 'monomerB', 'dimer')
  psi4.IO.set_default_namespace('dimer')
  
  aux_basis = psi4.core.BasisSet.build(wfn_dimer.molecule(), "DF_BASIS_SAPT",
                                psi4.core.get_global_option("DF_BASIS_SAPT"),
                                "RIFIT", psi4.core.get_global_option("BASIS"))
  wfn_dimer.set_basisset("DF_BASIS_SAPT", aux_basis)
  wfn_dimer.set_basisset("DF_BASIS_ELST", aux_basis)
  
  psi4.sapt(wfn_dimer,wfn_monA,wfn_monB)

In this way, any of the SCF options can be tweaked for individual fragments.
For optimal efficiency, the example uses ``set df_ints_io save`` to keep file 97,
which contains the three-index integrals for density fitting. ``set df_ints_io load``
then instructs the program to read these integrals from disk instead of recomputing
them. For each SCF computation, we use ``psi4.IO.set_default_namespace`` to uniquely
name scratch files. In the following SCF step, only file 97 is renamed using
``psi4.IO.change_file_namespace`` so that integrals can be read from it.
For more information on stability analysis, see the :ref:`stability <sec:scfstability_doc>`
documentation.

.. index:: SAPT; SAPT0

SAPT0
^^^^^

Generally speaking, SAPT0 should be applied to large systems or large data
sets. The performance of closed-shell SAPT0 relies entirely on error cancellation, which
seems to be optimal with a truncated aug-cc-pVDZ basis, namely,
jun-cc-pVDZ (which we have referred to in previous work as
aug-cc-pVDZ').  We do not recommend using SAPT0 with large basis sets
like aug-cc-pVTZ.  A systematic study of the accuracy of closed-shell SAPT0 and other SAPT 
truncations, using different basis sets, is reported in 
[Parker:2014:094106]_. In particular, an empirical recipe for scaled SAPT0
can yield improved performance and has been included in the output file as
the sSAPT0 interaction energy. sSAPT0 is a free by-product and is automatically 
computed when SAPT0 is requested (see above for more details).
The SAPT module has been used to perform SAPT0 computations with over
200 atoms and 2800 basis functions; this code should be scalable to 4000
basis functions. Publications resulting from the use of the SAPT0 code 
should cite the following publications: [Hohenstein:2010:184111]_ and 
[Hohenstein:2011:174107]_. If the open-shell SAPT0 code is used,
[Gonthier:2016:134106]_ should be additionally cited.

Basic SAPT0 Keywords
~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/sapt__sapt_level.rst
.. include:: autodir_options_c/sapt__basis.rst
.. include:: autodir_options_c/sapt__df_basis_sapt.rst
.. include:: autodir_options_c/sapt__df_basis_elst.rst
.. include:: autodir_options_c/sapt__freeze_core.rst
.. include:: autodir_options_c/sapt__cphf_r_convergence.rst
.. include:: autodir_options_c/sapt__maxiter.rst
.. include:: autodir_options_c/sapt__print.rst

Advanced SAPT0 Keywords
~~~~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/sapt__aio_cphf.rst
.. include:: autodir_options_c/sapt__aio_df_ints.rst
.. include:: autodir_options_c/sapt__coupled_induction.rst
.. include:: autodir_options_c/sapt__exch_scale_alpha.rst
.. include:: autodir_options_c/sapt__ints_tolerance.rst
.. include:: autodir_options_c/sapt__denominator_delta.rst
.. include:: autodir_options_c/sapt__denominator_algorithm.rst
.. include:: autodir_options_c/globals__debug.rst

Specific open-shell SAPT0 keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/sapt__sapt_mem_safety.rst
.. include:: autodir_options_c/sapt__coupled_induction.rst

.. index:: SAPT; SAPT(DFT)

.. _`sec:saptdft`:

SAPT(DFT)
^^^^^^^^^

In general, SAPT(DFT) should provide more accurate interaction energy 
components, and overall interaction energies, than SAPT0. The drawback 
is SAPT(DFT) method is more computationally demanding than SAPT0, 
SAPT(DFT) can still be applied to medium-sized or large systems. The 
SAPT(DFT) module was employed successfully in computations of systems 
with up to 2000 basis functions, and the code should be scalable to 
3000 basis functions. Like higher-order SAPT, SAPT(DFT) requires 
sufficient memory to hold :math:`2ovN_{\rm aux}` doubles.


SAPT(DFT) requires a few special treatments to obtain accurate
result. The DFT functionals used in SAPT(DFT) need to be asymptotically
corrected with Gradient Regulated Asymptotic Correction scheme (GRAC),
in order to recover the correct long-range asymptotic behavior
(approaching :math:`-1/r` as :math:`r` approaches infinity). The program
requires manual input of GRAC shift parameter for each monomer through
keywords |sapt__sapt_dft_grac_shift_a| and |sapt__sapt_dft_grac_shift_b|,
which should be equal to the difference of the actual ionization
potential and the corresponding Kohn-Sham HOMO energy. However,
|PSIfour| can automatically compute a GRAC shift for monomers A and B
if |sapt__sapt_dft_grac_compute| is set to ``SINGLE`` or ``ITERATIVE``. 


The dispersion term needs to be computed with orbital response for good
accuracy, and it is recommended to enable |sapt__sapt_dft_do_hybrid| (set to
``True`` by default). The coupled exchange-dispersion energy is usually
estimated by scaling from the uncoupled value either by a fitted fixed value
(suggested initially by [Hesselmann:2014:094107]_ for a local Hartree--Fock
(LHF) formulation and then revised by [Xie:2022:024801]_ for non-LHF) or by the
ratio of coupled and uncoupled dispersion energy (suggested by
[Podeszwa:2006:400]_ ). This can be controlled by keyword
|sapt__sapt_dft_exch_disp_scale_scheme|, with ``FIXED`` using the
Hesselmann/Xie approach (|PSIfours| default prior to Nov 2022), ``DISP`` using
the Podeszwa approach (|PSIfours| default after Nov 2022), or ``NONE`` for not
scaling and using the uncoupled exchange-dispersion energy directly.

.. warning:: Since Nov 2022, the defaults of options |sapt__sapt_dft_exch_disp_scale_scheme| and |sapt__sapt_dft_exch_disp_fixed_scale|
             have been changed. Before, the former defaulted to ``FIXED`` with Hesselmann value of 0.686 for the latter. Now, the former defaults to ``DISP`` and should you instead select ``FIXED``, the default for the latter is the Xie value of 0.770. This might cause
             an older version of |PSIfour| to produce a different value of
             exchange-dispersion energy from the latest version.

Basic Keywords for SAPT(DFT) 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/sapt__sapt_dft_grac_shift_a.rst
.. include:: autodir_options_c/sapt__sapt_dft_grac_shift_b.rst
.. include:: autodir_options_c/sapt__sapt_dft_do_dhf.rst
.. include:: autodir_options_c/sapt__sapt_dft_exch_disp_scale_scheme.rst
.. include:: autodir_options_c/sapt__sapt_dft_grac_compute.rst

Advanced Keywords for SAPT(DFT)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/sapt__sapt_dft_functional.rst
.. include:: autodir_options_c/sapt__sapt_dft_do_hybrid.rst
.. include:: autodir_options_c/sapt__sapt_dft_exch_disp_fixed_scale.rst
.. include:: autodir_options_c/sapt__sapt_dft_mp2_disp_alg.rst
.. include:: autodir_options_c/sapt__sapt_quiet.rst


.. index:: SAPT; higher-order

Higher-Order SAPT
^^^^^^^^^^^^^^^^^

For smaller systems (up to the size of a nucleic acid base pair), more
accurate interaction energies can be obtained through higher-order SAPT
computations. The SAPT module can perform density-fitted evaluations
of SAPT2, SAPT2+, SAPT2+(3), and SAPT2+3 energies for closed-shell systems 
only. Publications resulting from the use of the higher-order SAPT code 
should cite the following: 
[Hohenstein:2010:014101]_.

For methods SAPT2+ and above, one can replace the many-body treatment of
dispersion by an improved method based on coupled-cluster doubles (CCD).
This approach tends to give good improvements when dispersion effects
are very large, as in the PCCP dimer (see [Hohenstein:2011:2842]_).
As shown in [Parker:2014:094106]_, whether or not CCD dispersion offers
more accurate interaction energies tends to depend on the SAPT truncation
and basis set employed, due to cancellations of errors.  Thanks to
natural orbital methods [Parrish:2013:174102]_, the SAPT code
is able to include CCD dispersion with only a modest additional cost.
Computations employing CCD dispersion should cite [Parrish:2013:174102]_. 
To request CCD dispersion treatment in a SAPT computation, simply append
``(ccd)`` to the name of the method, as in the following examples ::

	energy('sapt2+(ccd)')
	energy('sapt2+(3)(ccd)')
	energy('sapt2+3(ccd)')

The :math:`\delta_{MP2}` corrections can also be computed automatically
by appending ``dmp2`` to the name of the method, with or without CCD dispersion ::

	energy('sapt2+dmp2')
	energy('sapt2+(3)dmp2')
	energy('sapt2+3dmp2')
	energy('sapt2+(ccd)dmp2')
	energy('sapt2+(3)(ccd)dmp2')
	energy('sapt2+3(ccd)dmp2')

A brief note on memory usage: the higher-order SAPT code assumes that
certain quantities can be held in core. This code requires sufficient
memory to hold :math:`3o^2v^2+v^2N_{aux}` arrays in core. With this
requirement computations on the adenine-thymine complex can be performed
with an aug-cc-pVTZ basis in less than 64GB of memory.

Higher-order SAPT is treated separately from the highly optimized SAPT0
code, therefore, higher-order SAPT uses a separate set of keywords. 
The following keywords are relevant for higher-order SAPT.

Basic Keywords for Higher-order SAPT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/sapt__basis.rst
.. include:: autodir_options_c/sapt__df_basis_sapt.rst
.. include:: autodir_options_c/globals__freeze_core.rst
.. include:: autodir_options_c/sapt__print.rst

Advanced Keywords for Higher-order SAPT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/sapt__do_ccd_disp.rst
.. include:: autodir_options_c/sapt__do_mbpt_disp.rst
.. include:: autodir_options_c/sapt__do_third_order.rst
.. include:: autodir_options_c/sapt__ints_tolerance.rst
.. include:: autodir_options_c/sapt__sapt_mem_check.rst
.. include:: autodir_options_c/globals__debug.rst

MP2 Natural Orbitals
^^^^^^^^^^^^^^^^^^^^

One of the unique features of the SAPT module is its ability to use
MP2 natural orbitals (NOs) to speed up the evaluation of the triples
contribution to dispersion. By transforming to the MP2 NO basis, we can
throw away virtual orbitals that are expected to contribute little to the
dispersion energy. Speedups in excess of :math:`50 \times` are possible. In
practice, this approximation is very good and should always be applied.
Publications resulting from the use of MP2 NO-based approximations should 
cite the following: [Hohenstein:2010:104107]_.

Basic Keywords Controlling MP2 NO Approximations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/sapt__nat_orbs_t2.rst
.. include:: autodir_options_c/sapt__nat_orbs_t3.rst
.. include:: autodir_options_c/sapt__nat_orbs_v4.rst
.. include:: autodir_options_c/sapt__occ_tolerance.rst

.. comment Advanced Keywords Controlling MP2 NO Approximations
.. comment ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. comment .. include:: autodir_options_c/sapt__nat_orbs_t2.rst

.. index:: SAPT; charge-transfer

.. _`sec:saptct`:

Charge-Transfer in SAPT
^^^^^^^^^^^^^^^^^^^^^^^

It is possible to obtain the stabilization energy of a complex due to
charge-transfer effects from a SAPT computation. The charge-transfer energy 
can be computed with the SAPT module as described by Stone
and Misquitta [Misquitta:2009:201]_.

Charge-transfer energies can be obtained from the following calls to the
energy function. ::

    energy('sapt0-ct')
    energy('sapt2-ct')
    energy('sapt2+-ct')
    energy('sapt2+(3)-ct')
    energy('sapt2+3-ct')
    energy('sapt2+(ccd)-ct')
    energy('sapt2+(3)(ccd)-ct')
    energy('sapt2+3(ccd)-ct')

For now, charge transfer computations are not available with open-shell SAPT0.

A SAPT charge-transfer analysis will perform 5 HF computations: the dimer
in the dimer basis, monomer A in the dimer basis, monomer B in the dimer
basis, monomer A in the monomer A basis, and monomer B in the monomer B
basis. Next, it performs two SAPT computations, one in the dimer basis and
one in the monomer basis. Finally, it will print a summary of the
charge-transfer results::

      SAPT Charge Transfer Analysis
    ------------------------------------------------------------------------------------------------
      SAPT Induction (Dimer Basis)       -2.0970 [mEh]      -1.3159 [kcal/mol]      -5.5057 [kJ/mol]
      SAPT Induction (Monomer Basis)     -1.1396 [mEh]      -0.7151 [kcal/mol]      -2.9920 [kJ/mol]
      SAPT Charge Transfer               -0.9574 [mEh]      -0.6008 [kcal/mol]      -2.5137 [kJ/mol]

These results are for the water dimer geometry shown above computed with 
SAPT0/aug-cc-pVDZ. 


.. index:: 
   pair: SAPT; output

Monomer-Centered Basis Computations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The charge-transfer analysis above is carried out by taking the
difference between SAPT induction as calculated in the dimer-centered
basis (*i.e.*, each monomer sees the basis functions on both monomers)
vs. the monomer-centered basis (*i.e.*, each monomer utilizes only its
own basis set).  It is also possible to run a closed-shell SAPT computation at any
level using only the monomer-centered basis. To do this, simply add
``sapt_basis='monomer'`` to the energy function, such as ::

    energy('sapt2',sapt_basis='monomer')

This procedure leads to faster compuations, but it converges more slowly
towards the complete basis set limit than the default procedure, which uses
the dimer-centered basis set.  Hence, monomer-centered basis SAPT
computations are not recommended. The open-shell SAPT0 code is not
compatible yet with monomer-centered computations.

Computations with Mid-bonds
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SAPT computations with midbonds can be accomplished by adding a third ghost
monomer to the computation.  For example ::

    molecule dimer {
        0 1
        He 0  0  5
        --
        0 1
        He 0  0 -5
        --
        0 1
        @He 0  0  0
    }


Here the functions of the third monomer will be added to the virtual space of
the entire computation. Note that an error will be thrown if each atom in the
third monomer is not a ghost to prevent confusion with three-body SAPT which is
not currently supported by Psi4.


Interpreting SAPT Results
^^^^^^^^^^^^^^^^^^^^^^^^^

We will examine the results of a SAPT2+3/aug-cc-pVDZ computation on the
water dimer. This computation can be performed with the following 
input::

    molecule water_dimer {
         0 1
         O  -1.551007  -0.114520   0.000000
         H  -1.934259   0.762503   0.000000
         H  -0.599677   0.040712   0.000000
         --
         0 1
         O   1.350625   0.111469   0.000000
         H   1.680398  -0.373741  -0.758561
         H   1.680398  -0.373741   0.758561
         units angstrom
    }
    
    set globals {
        basis          aug-cc-pvdz
        guess          sad
        scf_type       df
    }
    
    set sapt {
        print          1
        nat_orbs_t2    true
        freeze_core    true
    }
    
    energy('sapt2+3')

To reiterate some of the options mentioned above: the
|sapt__nat_orbs_t2| option will compute MP2 natural orbitals and use
them in the evaluation of the triples correction to dispersion, and the
|sapt__freeze_core| option will freeze the core throughout the SAPT
computation. This SAPT2+3/aug-cc-pVDZ computation produces the following
results::

    SAPT Results
  --------------------------------------------------------------------------------------------------------
    Electrostatics                -13.06509118 [mEh]      -8.19846883 [kcal/mol]     -34.30239689 [kJ/mol]
      Elst10,r                    -13.37542977 [mEh]      -8.39320925 [kcal/mol]     -35.11719087 [kJ/mol]
      Elst12,r                      0.04490350 [mEh]       0.02817737 [kcal/mol]       0.11789413 [kJ/mol]
      Elst13,r                      0.26543510 [mEh]       0.16656305 [kcal/mol]       0.69689985 [kJ/mol]

    Exchange                       13.41768202 [mEh]       8.41972294 [kcal/mol]      35.22812415 [kJ/mol]
      Exch10                       11.21822294 [mEh]       7.03954147 [kcal/mol]      29.45344432 [kJ/mol]
      Exch10(S^2)                  11.13802706 [mEh]       6.98921779 [kcal/mol]      29.24289005 [kJ/mol]
      Exch11(S^2)                   0.04558907 [mEh]       0.02860757 [kcal/mol]       0.11969410 [kJ/mol]
      Exch12(S^2)                   2.15387002 [mEh]       1.35157390 [kcal/mol]       5.65498573 [kJ/mol]

    Induction                      -3.91313050 [mEh]      -2.45552656 [kcal/mol]     -10.27392413 [kJ/mol]
      Ind20,r                      -4.57530818 [mEh]      -2.87104935 [kcal/mol]     -12.01247162 [kJ/mol]
      Ind30,r                      -4.91714746 [mEh]      -3.08555675 [kcal/mol]     -12.90997067 [kJ/mol]
      Ind22                        -0.83718642 [mEh]      -0.52534243 [kcal/mol]      -2.19803293 [kJ/mol]
      Exch-Ind20,r                  2.47828501 [mEh]       1.55514739 [kcal/mol]       6.50673730 [kJ/mol]
      Exch-Ind30,r                  4.33916119 [mEh]       2.72286487 [kcal/mol]      11.39246770 [kJ/mol]
      Exch-Ind22                    0.45347471 [mEh]       0.28455969 [kcal/mol]       1.19059785 [kJ/mol]
      delta HF,r (2)               -1.43239563 [mEh]      -0.89884187 [kcal/mol]      -3.76075473 [kJ/mol]
      delta HF,r (3)               -0.85440936 [mEh]      -0.53614999 [kcal/mol]      -2.24325177 [kJ/mol]

    Dispersion                     -3.62000698 [mEh]      -2.27158877 [kcal/mol]      -9.50432831 [kJ/mol]
      Disp20                       -3.54291925 [mEh]      -2.22321549 [kcal/mol]      -9.30193450 [kJ/mol]
      Disp30                        0.05959979 [mEh]       0.03739944 [kcal/mol]       0.15647926 [kJ/mol]
      Disp21                        0.11216169 [mEh]       0.07038252 [kcal/mol]       0.29448051 [kJ/mol]
      Disp22 (SDQ)                 -0.17892163 [mEh]      -0.11227502 [kcal/mol]      -0.46975875 [kJ/mol]
      Disp22 (T)                   -0.47692534 [mEh]      -0.29927518 [kcal/mol]      -1.25216749 [kJ/mol]
      Est. Disp22 (T)              -0.54385233 [mEh]      -0.34127251 [kcal/mol]      -1.42788430 [kJ/mol]
      Exch-Disp20                   0.64545587 [mEh]       0.40502969 [kcal/mol]       1.69464439 [kJ/mol]
      Exch-Disp30                  -0.01823410 [mEh]      -0.01144207 [kcal/mol]      -0.04787362 [kJ/mol]
      Ind-Disp30                   -0.91816882 [mEh]      -0.57615966 [kcal/mol]      -2.41065224 [kJ/mol]
      Exch-Ind-Disp30               0.76487181 [mEh]       0.47996433 [kcal/mol]       2.00817094 [kJ/mol]

  Total HF                         -5.68662563 [mEh]      -3.56841161 [kcal/mol]     -14.93023559 [kJ/mol]
  Total SAPT0                      -8.58408901 [mEh]      -5.38659740 [kcal/mol]     -22.53752571 [kJ/mol]
  Total SAPT2                      -6.72343814 [mEh]      -4.21902130 [kcal/mol]     -17.65238683 [kJ/mol]
  Total SAPT2+                     -7.33405042 [mEh]      -4.60218631 [kcal/mol]     -19.25554938 [kJ/mol]
  Total SAPT2+(3)                  -7.00901553 [mEh]      -4.39822383 [kcal/mol]     -18.40217026 [kJ/mol]
  Total SAPT2+3                    -7.18054663 [mEh]      -4.50586123 [kcal/mol]     -18.85252518 [kJ/mol]

  Special recipe for scaled SAPT0 (see Manual):
    Electrostatics sSAPT0         -13.37542977 [mEh]      -8.39320925 [kcal/mol]     -35.11719087 [kJ/mol]
    Exchange sSAPT0                11.21822294 [mEh]       7.03954147 [kcal/mol]      29.45344432 [kJ/mol]
    Induction sSAPT0               -3.47550008 [mEh]      -2.18090932 [kcal/mol]      -9.12492546 [kJ/mol]
    Dispersion sSAPT0              -2.88342055 [mEh]      -1.80937379 [kcal/mol]      -7.57042064 [kJ/mol]
  Total sSAPT0                     -8.51612746 [mEh]      -5.34395089 [kcal/mol]     -22.35909265 [kJ/mol]
  --------------------------------------------------------------------------------------------------------

At the bottom of this output are the total SAPT energies (defined above),
they are composed of subsets of the individual terms printed above. The
individual terms are grouped according to the component of the interaction
to which they contribute. The total component energies (*i.e.*,
electrostatics, exchange, induction, and dispersion) represent what we
regard as the best estimate available at a given level of SAPT computed
from a subset of the terms of that grouping. The groupings shown above are
not unique and are certainly not rigorously defined. We regard the groupings 
used in |PSIfour| as a "chemist's grouping" as opposed to a more
mathematically based grouping, which would group all exchange terms 
(*i.e.* :math:`E_{exch-ind,resp}^{(20)}`, :math:`E_{exch-disp}^{(20)}`, *etc.*) in
the exchange component. A final note is that both ``Disp22(T)``
and ``Est.Disp22(T)`` results appear if MP2 natural orbitals are 
used to evaluate the triples correction to dispersion. The ``Disp22(T)`` 
result is the triples correction as computed in the truncated NO basis;  
``Est.Disp22(T)`` is a scaled result that attempts to recover
the effect of the truncated virtual space and is our best estimate. The ``Est.Disp22(T)``
value is used in the SAPT energy and dispersion component (see [Hohenstein:2010:104107]_ 
for details). Finally, this part of the output file contains sSAPT0, a special scaling
scheme of the SAPT0 energy that can yield improved results and was described in more details
above. The corresponding scaled total component energies are printed as well.

As mentioned above, SAPT results with scaled exchange are also optionally available
by setting the |sapt__exch_scale_alpha| keyword. When activated, the unscaled results are
printed first as reported above, and then repeated with exchange scaling for all 
relevant terms: :: 

    SAPT Results ==> ALL S2 TERMS SCALED (see Manual) <==

    Scaling factor (Exch10/Exch10(S^2))^{Alpha} =     1.007200
    with Alpha =     1.000000
  --------------------------------------------------------------------------------------------------------
    Electrostatics                -13.06509118 [mEh]      -8.19846883 [kcal/mol]     -34.30239689 [kJ/mol]
      Elst10,r                    -13.37542977 [mEh]      -8.39320925 [kcal/mol]     -35.11719087 [kJ/mol]
      Elst12,r                      0.04490350 [mEh]       0.02817737 [kcal/mol]       0.11789413 [kJ/mol]
      Elst13,r                      0.26543510 [mEh]       0.16656305 [kcal/mol]       0.69689985 [kJ/mol]

    Exchange sc.                   13.43351854 [mEh]       8.42966050 [kcal/mol]      35.26970292 [kJ/mol]
      Exch10                       11.21822294 [mEh]       7.03954147 [kcal/mol]      29.45344432 [kJ/mol]
      Exch10(S^2)                  11.13802706 [mEh]       6.98921779 [kcal/mol]      29.24289005 [kJ/mol]
      Exch11(S^2) sc.               0.04591732 [mEh]       0.02881355 [kcal/mol]       0.12055592 [kJ/mol]
      Exch12(S^2) sc.               2.16937828 [mEh]       1.36130548 [kcal/mol]       5.69570268 [kJ/mol]

    Induction sc.                  -3.90986540 [mEh]      -2.45347768 [kcal/mol]     -10.26535160 [kJ/mol]
      Ind20,r                      -4.57530818 [mEh]      -2.87104935 [kcal/mol]     -12.01247162 [kJ/mol]
      Ind30,r                      -4.91714746 [mEh]      -3.08555675 [kcal/mol]     -12.90997067 [kJ/mol]
      Ind22                        -0.83718642 [mEh]      -0.52534243 [kcal/mol]      -2.19803293 [kJ/mol]
      Exch-Ind20,r sc.              2.49612913 [mEh]       1.56634474 [kcal/mol]       6.55358703 [kJ/mol]
      Exch-Ind30,r sc.              4.37040396 [mEh]       2.74247000 [kcal/mol]      11.47449560 [kJ/mol]
      Exch-Ind22 sc.                0.45673981 [mEh]       0.28660857 [kcal/mol]       1.19917038 [kJ/mol]
      delta HF,r (2) sc.           -1.45023975 [mEh]      -0.91003922 [kcal/mol]      -3.80760445 [kJ/mol]
      delta HF,r (3) sc.           -0.90349624 [mEh]      -0.56695248 [kcal/mol]      -2.37212939 [kJ/mol]

    Dispersion sc.                 -3.60998364 [mEh]      -2.26529903 [kcal/mol]      -9.47801205 [kJ/mol]
      Disp20                       -3.54291925 [mEh]      -2.22321549 [kcal/mol]      -9.30193450 [kJ/mol]
      Disp30                        0.05959979 [mEh]       0.03739944 [kcal/mol]       0.15647926 [kJ/mol]
      Disp21                        0.11216169 [mEh]       0.07038252 [kcal/mol]       0.29448051 [kJ/mol]
      Disp22 (SDQ)                 -0.17892163 [mEh]      -0.11227502 [kcal/mol]      -0.46975875 [kJ/mol]
      Disp22 (T)                   -0.47692534 [mEh]      -0.29927518 [kcal/mol]      -1.25216749 [kJ/mol]
      Est. Disp22 (T)              -0.54385233 [mEh]      -0.34127251 [kcal/mol]      -1.42788430 [kJ/mol]
      Exch-Disp20 sc.               0.65010327 [mEh]       0.40794598 [kcal/mol]       1.70684615 [kJ/mol]
      Exch-Disp30 sc.              -0.01836538 [mEh]      -0.01152445 [kcal/mol]      -0.04821832 [kJ/mol]
      Ind-Disp30                   -0.91816882 [mEh]      -0.57615966 [kcal/mol]      -2.41065224 [kJ/mol]
      Exch-Ind-Disp30 sc.           0.77037903 [mEh]       0.48342016 [kcal/mol]       2.02263015 [kJ/mol]

  Total HF                         -5.68662563 [mEh]      -3.56841161 [kcal/mol]     -14.93023559 [kJ/mol]
  Total SAPT0 sc.                  -8.57944161 [mEh]      -5.38368112 [kcal/mol]     -22.52532395 [kJ/mol]
  Total SAPT2 sc.                  -6.69968912 [mEh]      -4.20411857 [kcal/mol]     -17.59003378 [kJ/mol]
  Total SAPT2+ sc.                 -7.31030140 [mEh]      -4.58728357 [kcal/mol]     -19.19319632 [kJ/mol]
  Total SAPT2+(3) sc.              -6.98526650 [mEh]      -4.38332109 [kcal/mol]     -18.33981720 [kJ/mol]
  Total SAPT2+3 sc.                -7.15142168 [mEh]      -4.48758504 [kcal/mol]     -18.77605762 [kJ/mol]
  --------------------------------------------------------------------------------------------------------

The scaling factor is reported at the top (here ``1.0072``) together with the 
:math:`\alpha` parameter. All terms that are scaled are indicated by the ``sc.``
label. Note that if Exch10 is less than :math:`10^{-5}`, the scaling factor is
set to :math:`1.0`.


.. _`sec:saptfitB`:


.. caution:: To density fit the dispersion terms in SAPT, the RI auxiliary
   basis set (*e.g.*, aug-cc-pVDZ-RI) controlled through
   |sapt__df_basis_sapt| performs well. For Fock-type terms (*i.e.*,
   electrostatics, exchange, induction, and core Fock matrix elements in
   exchange-dispersion), the JKFIT density-fitting auxiliary basis
   (*e.g.*, aug-cc-pVDZ-JKFIT) is more appropriate. The :ref:`FISAPT
   module<sec:fisapt>` has always used JKFIT in this role. The
   :ref:`SAPT module<sec:sapt>` newly (see :ref:`fitting notes <sec:saptfitA>` ) uses
   JKFIT for computations targeting SAPT0 and sSAPT0 methods. But the
   :ref:`SAPT module<sec:sapt>` still uses the RI basis for higher-order
   SAPT. For heavier elements (*i.e.*, second-row and beyond), the RI
   auxiliary basis is unsound for this role (insufficiently flexible).
   For higher-order methods in :ref:`SAPT module<sec:sapt>`, there is
   no workaround; on-the-fly construction of an auxiliary basis through
   Cholesky decomposition (not implemented) is the long-term solution.

.. _`sec:sfsapt`:

Spin-Flip SAPT
^^^^^^^^^^^^^^

SAPT0 with two open-shell references will always yield a high-spin complex.  In
order to obtain a SAPT-based estimate of the splittings between different spin
states of a complex the first-order exchange energies for all multiplets can be
shown to be a linear combination of two matrix elements: a diagonal exchange
term that determines the spin-averaged effect and a spin-flip term responsible
for the splittings between the states. The numerical factors in this linear
combination are determined solely by the Clebsch-Gordan coefficients:
accordingly, the :math:`S^{2}` approximation implies a Heisenberg Hamiltonian
picture with a single coupling strength parameter determining all the
splittings. This method can be invoked with `energy("SF-SAPT")` and
publications resulting from the use of the SF-SAPT code should cite the
following publications: [Patkowski:2018:164110]_

.. _`sec:saptinf`:

Higher-Order Exchange Terms without Single-Exchange Approximation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Recently, several SAPT higher-order exchange terms have been derived without
the :math:`S^{2}` approximation: :math:`E_{exch-ind}^{(20)}` [Schaffer:2012:1235]_, 
:math:`E_{exch-disp}^{(20)}` [Schaffer:2013:2570]_, and :math:`E_{exch-ind}^{(30)}` 
[Waldrop:2021:024103]_. The second-order terms can be computed with the following
settings::

    set SAPT_DFT_FUNCTIONAL HF
    set DO_IND_EXCH_SINF true        # calculate Exch-Ind20 (S^inf) 
    set SAPT_DFT_MP2_DISP_ALG fisapt 
    set DO_DISP_EXCH_SINF true       # calculate Exch-Disp20 (S^inf)
    energy('sapt(dft)')
                       
and the third-order exchange-induction term is computed as follows::

    set DO_IND30_EXCH_SINF true        # calculate Exch-Ind30 (S^inf) 
    energy('sapt2+3')
                       
These calculations are performed with the atomic orbital and 
density-fitting scheme described in the Supplementary Material to
[Smith:2020:184108]_ for the second-order terms and in [Waldrop:2021:024103]_
for the third-order exchange induction. The coupled (response) version of the
exchange-induction corrections are also calculated, exactly for 
:math:`E_{exch-ind,resp}^{(20)}` and by scaling the uncoupled term for
:math:`E_{exch-ind,resp}^{(30)}`.

S^inf Keywords
~~~~~~~~~~~~~~

.. include:: autodir_options_c/sapt__do_ind_exch_sinf.rst
.. include:: autodir_options_c/sapt__do_disp_exch_sinf.rst
.. include:: autodir_options_c/sapt__do_ind30_exch_sinf.rst

.. _`sec:saptd`:

SAPT0-D
~~~~~~~

In SAPT0, the computation of :math:`E_{disp}^{(20)} + E_{exch-disp}^{(20)}` represents
the computational bottleneck. One can avoid this bottleneck by replacing these
dispersion terms with the empirical D3 corrections developed by Grimme.
  
:ref:`Grimme's dispersion corrections are discussed here. <sec:dftd3>`

The corresponding method, termed SAPT0-D, thus relies on empirically fit parameters
specific to SAPT0/jun-cc-pVDZ. While SAPT0-D can be used with any of the -D 
variants using default parameters optimized for Hartee--Fock interaction energies, 
we recommend using the refit parameters with Becke-Johnson damping, as described in
[Schriber:2021:234107]_. Again, use of SAPT0-D with a basis set other than
jun-cc-pVDZ is not tested and not guaranteed to give meaningful results without
refitting the dispersion parameters. 
A simple water dimer computation using SAPT0-D may look like::

	molecule water_dimer {
	     0 1
	     O  -1.551007  -0.114520   0.000000
	     H  -1.934259   0.762503   0.000000
	     H  -0.599677   0.040712   0.000000
	     --
	     0 1
	     O   1.350625   0.111469   0.000000
	     H   1.680398  -0.373741  -0.758561
	     H   1.680398  -0.373741   0.758561
	
	     units angstrom
	     no_reorient
	     symmetry c1
	}
	
	set basis jun-cc-pvdz

	energy('sapt0-d3mbj') # runs the recommended dispersion correction
    energy('sapt0-d3') # tests an alternative damping scheme/parameterization

Given the naturally pairwise-atomic nature of these empirical dispersion corrections,
integration with existing FSAPT functionality is also available simply by calling
``energy("fsapt0-d3mbj")``. See :ref:`FSAPT <sec:fisapt>` documentation for more details on using FSAPT
for functional group analyses.
