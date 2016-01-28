
.. include:: autodoc_abbr_options_c.rst

.. index::
   single: SAPT
   pair: SAPT; theory

.. _`sec:sapt`:

SAPT: Symmetry-Adapted Perturbation Theory
==========================================

.. codeauthor:: Edward G. Hohenstein and Rob M. Parrish
.. sectionauthor:: Edward G. Hohenstein

*Module:* :ref:`Keywords <apdx:sapt>`, :ref:`PSI Variables <apdx:sapt_psivar>`, :source:`LIBSAPT_SOLVER <src/lib/libsapt_solver>`

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

Symmetry-adapted perturbation theory (SAPT) provides a means of directly
computing the noncovalent interaction between two molecules, that is, the
interaction energy is determined without computing the total energy of the
monomers or dimer. In addition, SAPT provides a decomposition of the
interaction energy into physically meaningful components: *i.e.,*
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

Several truncations of the SAPT expansion are available in the SAPT
module of |PSIfour|. The simplest truncation of SAPT is denoted SAPT0
and defined in Eq. :eq:`SAPT0`.

.. math:: E_{SAPT0} = E_{elst}^{(10)} + E_{exch}^{(10)} + E_{ind,resp}^{(20)} +
   E_{exch-ind,resp}^{(20)} + E_{disp}^{(20)} + E_{exch-disp}^{(20)} + \delta_{HF}^{(2)}
   :label: SAPT0

In this notation, :math:`E^{(vw)}` defines the order in :math:`V` and in :math:`W_A+W_B`; the
subscript, :math:`resp`, indicates that orbital relaxation effects are included.

.. math:: E_{SAPT2} = E_{SAPT0} + E_{elst,resp}^{(12)} + E_{exch}^{(11)} +
   E_{exch}^{(12)} +\/ ^{t}\!E_{ind}^{(22)} +\/ ^{t}\!E_{exch-ind}^{(22)}
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

The :math:`\delta_{HF}^{(2)}` and :math:`\delta_{HF}^{(3)}` terms take into
account higher-order induction effects and are included in the definition
of SAPT terms. They are computed from the Hartree-Fock supermolecular interaction energy
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
          E_{exch}^{(11)} + E_{exch}^{(12)} +\/ ^{t}\!E_{ind}^{(22)}
          +\/ ^{t}\!E_{exch-ind}^{(22)} + E_{disp}^{(20)} + E_{exch-disp}^{(20)})

.. math:: \delta_{MP2}^{(3)} = \delta_{MP2}^{(2)} - (E_{ind-disp}^{(30)} + E_{exch-ind-disp}^{(30)})

where :math:`E_{int}^{MP2, corr}` is the correlation part of the supermolecular MP2 
interaction energy. :math:`\delta_{MP2}^{(2)}` and :math:`\delta_{MP2}^{(3)}` also improve the 
description of electrostatically dominated complexes. :math:`\delta_{MP2}^{(2)}`
can be applied to SAPT2+ or SAPT2+(3) energies whereas :math:`\delta_{MP2}^{(3)}` 
should be applied to SAPT2+3 energies.

A thorough analysis of the performance of these truncations of SAPT can be
found in a review by Hohenstein and Sherrill [Hohenstein:2012:WIREs]_,
and a systematic study of the accuracy of these truncations (with and 
without an improved CCD treatment of dispersion) using different basis sets
is reported in [Parker:2014:094106]_.

The SAPT module relies entirely on the density-fitting approximation
of the two-electron integrals. The factorization of the SAPT energy
expressions, as implemented in |PSIfour|, assumes the use of density-fitted
two-electron integrals, therefore, the SAPT module cannot be run with
exact integrals. In practice, we have found that the density-fitting
approximation introduces negligible errors into the SAPT energy 
(often less than 0.01 kcal/mol for small dimers) and greatly
improves efficiency. 

The S\ :superscript:`2` approximation and scaling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All exchange terms in SAPT arise from the antisymmetrization
of the wavefunctions of monomers A and B. Taking into account exchange of all possible
electron pairs between the two monomers yields to complicated formulae.
For this reason, exchange terms are often evaluated in the :math:`S^{2}`
approximation, that can be interpreted as the exchange of a single electron 
pair between monomers.

The :math:`S^{2}` approximation is usually pretty good, but may 
break down for short intermolecular distance, particularly in high-order
terms. To compensate these deviations, Parker et al. [Parker:2014:094106]_ 
recommend to scale all :math:`S^{2}` approximated exchange terms by the ratio:

.. math:: p_{EX}(\alpha) = \left( \frac{E_{exch}^{(10)}}{E_{exch}^{(10)}(S^{2})} \right)^{\alpha}

where the recommended exponent is :math:`\alpha = 1`. SAPT energies with and without
exchange scaling are reported in the output file.

In addition, the sSAPT0 method uses an empirically adjusted exponent :math:`\alpha = 3.0`, 
yielding improved results over regular SAPT0 in the jun-cc-pVDZ basis set (see [Parker:2014:094106]_).

.. math:: E_{sSAPT0} = E_{elst}^{(10)} + E_{exch}^{(10)} + E_{ind,resp}^{(20)} +
   p_{EX}(3.0) E_{exch-ind,resp}^{(20)} + E_{disp}^{(20)} + p_{EX}(3.0) E_{exch-disp}^{(20)} 
   + \delta_{HF}^{(2)}
   :label: SAPT0

where :math:`\delta_{HF}^{(2)}` is computed without any scaling.


A First Example
^^^^^^^^^^^^^^^

The following is the simplest possible input that will perform all
available SAPT computations (normally, you would pick one of these methods). ::

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
	
	set globals {
	    basis         aug-cc-pvdz
	}
	
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
c1`` molecule keywords, as in the example above. A final note is that the
SAPT module is only capable of performing SAPT computations for
interactions between closed-shell singlets.

The example input shown above would not be used in practice.
To exploit the efficiency of the density-fitted SAPT implementation in
|PSIfour|, the SCF computations should also be performed with density-fitted
(DF) integrals. ::

    set globals {
        basis         aug-cc-pvdz
        df_basis_scf  aug-cc-pvdz-jkfit
        df_basis_sapt aug-cc-pvdz-ri
        guess         sad
        scf_type      df
    }
    
    set sapt {
        print         1
    }

These options will perform the SAPT computation with DF-HF and a 
superposition-of-atomic-densities guess. This is the preferred method of 
running the SAPT module.

.. index:: SAPT; SAPT0

SAPT0
^^^^^

Generally speaking, SAPT0 should be applied to large systems or large data
sets. The performance of SAPT0 relies entirely on error cancellation, which
seems to be optimal with a truncated aug-cc-pVDZ basis, namely,
jun-cc-pVDZ (which we have referred to in previous work as
aug-cc-pVDZ').  We do not recommend using SAPT0 with large basis sets
like aug-cc-pVTZ.  A systematic study of the accuracy of SAPT0 and other SAPT 
truncations, using different basis sets, is reported in 
[Parker:2014:094106]_.
The SAPT module has been used to perform SAPT0 computations with over
200 atoms and 2800 basis functions; this code should be scalable to 4000
basis functions. Publications resulting from the use of the SAPT0 code 
should cite the following publications: [Hohenstein:2010:184111]_ and 
[Hohenstein:2011:174107]_.

Basic SAPT0 Keywords
~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/sapt__sapt_level.rst
.. include:: autodir_options_c/sapt__basis.rst
.. include:: autodir_options_c/sapt__df_basis_sapt.rst
.. include:: autodir_options_c/sapt__df_basis_elst.rst
.. include:: autodir_options_c/sapt__freeze_core.rst
.. include:: autodir_options_c/sapt__d_convergence.rst
.. include:: autodir_options_c/sapt__e_convergence.rst
.. include:: autodir_options_c/sapt__maxiter.rst
.. include:: autodir_options_c/sapt__print.rst

Advanced SAPT0 Keywords
~~~~~~~~~~~~~~~~~~~~~~~

.. include:: autodir_options_c/sapt__aio_cphf.rst
.. include:: autodir_options_c/sapt__aio_df_ints.rst
.. include:: autodir_options_c/sapt__no_response.rst
.. include:: autodir_options_c/sapt__ints_tolerance.rst
.. include:: autodir_options_c/sapt__denominator_delta.rst
.. include:: autodir_options_c/sapt__denominator_algorithm.rst
.. include:: autodir_options_c/sapt__sapt_os_scale.rst
.. include:: autodir_options_c/sapt__sapt_ss_scale.rst
.. include:: autodir_options_c/globals__debug.rst

.. index:: SAPT; higher-order

Higher-Order SAPT
^^^^^^^^^^^^^^^^^

For smaller systems (up to the size of a nucleic acid base pair), more
accurate interaction energies can be obtained through higher-order SAPT
computations. The SAPT module can perform density-fitted evaluations
of SAPT2, SAPT2+, SAPT2+(3), and SAPT2+3 energies. Publications resulting
from the use of the higher-order SAPT code should cite the following: 
[Hohenstein:2010:014101]_.

For methods SAPT2+ and above, one can replace the many-body treatment of
dispersion by an improved method based on coupled-cluster doubles (CCD).
This approach tends to give good improvements when dispersion effects
are very large, as in the PCCP dimer (see [Hohenstein:2011:2842]_).
As shown in [Parker:2014:094106]_, whether or not CCD dispersion offers
more accurate interaction energies tends to depend on the SAPT truncation
and basis set employed, due to cancellations of errors.  Thanks to
natural orbital methods [Parrish:2013:174102]_, the SAPT code in Psi
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

Higher-order SAPT is treated separately from the higly optimized SAPT0
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

A SAPT charge-transfer analysis will perform 5 HF computations: the dimer
in the dimer basis, monomer A in the dimer basis, monomer B in the dimer
basis, monomer A in the monomer A basis, and monomer B in the monomer B
basis. Next, it performs two SAPT computations, one in the dimer basis and
one in the monomer basis. Finally, it will print a summary of the
charge-transfer results::

      SAPT Charge Transfer Analysis
    -----------------------------------------------------------------------------
      SAPT Induction (Dimer Basis)         -2.0970 mH       -1.3159 kcal mol^-1
      SAPT Induction (Monomer Basis)       -1.1396 mH       -0.7151 kcal mol^-1
      SAPT Charge Transfer                 -0.9574 mH       -0.6008 kcal mol^-1

These results are for the water dimer geometry shown above computed with 
SAPT0/aug-cc-pVDZ. 


.. index:: 
   pair: SAPT; output

Monomer-Centered Basis Computations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The charge-transfer analysis above is carried out by taking the
difference between SAPT induction as calculated in the dimer-centered
basis (i.e., each monomer sees the basis functions on both monomers)
vs. the monomer-centered basis (i.e., each monomer utilizes only its
own basis set).  It is also possible to run a SAPT computation at any
level using only the monomer-centered basis.  To do this, simply add
``sapt_basis='monomer'`` to the energy function, such as ::

    energy('sapt2',sapt_basis='monomer')

This procedure leads to faster compuations, but it converges more slowly
towards the complete basis set limit than the default procedure, which uses
the dimer-centered basis set.  Hence, monomer-centered basis SAPT
computations are not recommended.


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

    SAPT Results ==> NO EXCHANGE SCALING APPLIED <==  
  --------------------------------------------------------------------------
    Electrostatics              -13.06509072 mH      -8.19846854 kcal mol^-1
      Elst10,r                  -13.37542914 mH      -8.39320885 kcal mol^-1
      Elst12,r                    0.04490321 mH       0.02817719 kcal mol^-1
      Elst13,r                    0.26543521 mH       0.16656311 kcal mol^-1

    Exchange                     13.41768624 mH       8.41972558 kcal mol^-1
      Exch10                     11.21822694 mH       7.03954398 kcal mol^-1
      Exch10(S^2)                11.13803101 mH       6.98922027 kcal mol^-1
      Exch11(S^2)                 0.04558910 mH       0.02860759 kcal mol^-1
      Exch12(S^2)                 2.15387020 mH       1.35157401 kcal mol^-1

    Induction                    -3.91313510 mH      -2.45552945 kcal mol^-1
      Ind20,r                    -4.57530912 mH      -2.87104994 kcal mol^-1
      Ind30,r                    -4.91715016 mH      -3.08555844 kcal mol^-1
      Ind22                      -0.83718660 mH      -0.52534254 kcal mol^-1
      Exch-Ind20,r                2.47828624 mH       1.55514816 kcal mol^-1
      Exch-Ind30,r                4.33916384 mH       2.72286653 kcal mol^-1
      Exch-Ind22                  0.45347494 mH       0.28455983 kcal mol^-1
      delta HF,r (2)             -1.43240056 mH      -0.89884496 kcal mol^-1
      delta HF,r (3)             -0.85441424 mH      -0.53615305 kcal mol^-1

    Dispersion                   -3.62000732 mH      -2.27158898 kcal mol^-1
      Disp20                     -3.54291985 mH      -2.22321586 kcal mol^-1
      Disp30                      0.05959980 mH       0.03739944 kcal mol^-1
      Disp21                      0.11216175 mH       0.07038256 kcal mol^-1
      Disp22 (SDQ)               -0.17892161 mH      -0.11227501 kcal mol^-1
      Disp22 (T)                 -0.47692539 mH      -0.29927522 kcal mol^-1
      Est. Disp22 (T)            -0.54385240 mH      -0.34127255 kcal mol^-1
      Exch-Disp20                 0.64545612 mH       0.40502984 kcal mol^-1
      Exch-Disp30                -0.01823410 mH      -0.01144207 kcal mol^-1
      Ind-Disp30                 -0.91816922 mH      -0.57615991 kcal mol^-1
      Exch-Ind-Disp30             0.76487219 mH       0.47996457 kcal mol^-1

  Total HF                           -5.68662563 mH      -3.56841161 kcal mol^-1
  Total SAPT0                        -8.58408936 mH      -5.38659762 kcal mol^-1
  Total SAPT2                        -6.72343851 mH      -4.21902154 kcal mol^-1
  Total SAPT2+                       -7.33405078 mH      -4.60218654 kcal mol^-1
  Total SAPT2+(3)                    -7.00901577 mH      -4.39822398 kcal mol^-1
  Total SAPT2+3                      -7.18054690 mH      -4.50586139 kcal mol^-1
  --------------------------------------------------------------------------


At the bottom of this output are the total SAPT energies (defined above),
they are composed of subsets of the individual terms printed above. The
individual terms are grouped according to the component of the interaction
to which they contribute. The total component energies (*i.e.,*
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
the effect of the truncated virtual space. The ``Est.Disp22(T)``
value is used in the SAPT energy and dispersion component (see [Hohenstein:2010:104107]_ 
for details). As indicated at the top of the result section, all results
are presented without exchange scaling. If the scaling factor :math:`p_{EX}` is 
significantly different from 1.0, results with exchange scaling are printed: ::

    SAPT Results ==> ALL S2 TERMS SCALED <== 

    Scaling factor:     1.007200  
  --------------------------------------------------------------------------
    Electrostatics              -13.06509072 mH      -8.19846854 kcal mol^-1
      Elst10,r                  -13.37542914 mH      -8.39320885 kcal mol^-1
      Elst12,r                    0.04490321 mH       0.02817719 kcal mol^-1
      Elst13,r                    0.26543521 mH       0.16656311 kcal mol^-1

    Exchange scal.               13.43352277 mH       8.42966315 kcal mol^-1
      Exch10                     11.21822694 mH       7.03954398 kcal mol^-1
      Exch10(S^2)                11.13803101 mH       6.98922027 kcal mol^-1
      Exch11(S^2) scal.           0.04591735 mH       0.02881357 kcal mol^-1
      Exch12(S^2) scal.           2.16937848 mH       1.36130560 kcal mol^-1

    Induction scal.              -3.90986999 mH      -2.45348056 kcal mol^-1
      Ind20,r                    -4.57530912 mH      -2.87104994 kcal mol^-1
      Ind30,r                    -4.91715016 mH      -3.08555844 kcal mol^-1
      Ind22                      -0.83718660 mH      -0.52534254 kcal mol^-1
      Exch-Ind20,r scal.          2.49613037 mH       1.56634552 kcal mol^-1
      Exch-Ind30,r scal.          4.37040664 mH       2.74247168 kcal mol^-1
      Exch-Ind22 scal.            0.45674004 mH       0.28660872 kcal mol^-1
      delta HF,r (2) scal.       -1.45024469 mH      -0.91004232 kcal mol^-1
      delta HF,r (3) scal.       -0.90350117 mH      -0.56695557 kcal mol^-1

    Dispersion scal.             -3.60998398 mH      -2.26529924 kcal mol^-1
      Disp20                     -3.54291985 mH      -2.22321586 kcal mol^-1
      Disp30                      0.05959980 mH       0.03739944 kcal mol^-1
      Disp21                      0.11216175 mH       0.07038256 kcal mol^-1
      Disp22 (SDQ)               -0.17892161 mH      -0.11227501 kcal mol^-1
      Disp22 (T)                 -0.47692539 mH      -0.29927522 kcal mol^-1
      Est. Disp22 (T)            -0.54385240 mH      -0.34127255 kcal mol^-1
      Exch-Disp20 scal.           0.65010352 mH       0.40794614 kcal mol^-1
      Exch-Disp30 scal.          -0.01836539 mH      -0.01152446 kcal mol^-1
      Ind-Disp30                 -0.91816922 mH      -0.57615991 kcal mol^-1
      Exch-Ind-Disp30 scal.       0.77037942 mH       0.48342040 kcal mol^-1

  Total HF                           -5.68662563 mH      -3.56841161 kcal mol^-1
  Total SAPT0 scal.                  -8.57944196 mH      -5.38368133 kcal mol^-1
  Total sSAPT0                       -8.51612775 mH      -5.34395107 kcal mol^-1
  Total SAPT2 scal.                  -6.69968948 mH      -4.20411879 kcal mol^-1
  Total SAPT2+ scal.                 -7.31030174 mH      -4.58728379 kcal mol^-1
  Total SAPT2+(3) scal.              -6.98526674 mH      -4.38332124 kcal mol^-1
  Total SAPT2+3 scal.                -7.15142193 mH      -4.48758520 kcal mol^-1
  --------------------------------------------------------------------------

Here, all previous results are repeated with all relevant exchange terms scaled. 
The scaling factor is reported at the top (here ``1.0072``) and all terms that
are scaled are indicated by the ``scal.`` keyword. Note that the sSAPT0 energy is 
reported here if the scaling factor is significantly different from 1.0, otherwise
it is reported in the unscaled results.
