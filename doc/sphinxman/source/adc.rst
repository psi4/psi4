.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2022 The Psi4 Developers.
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

.. index:: ADC
.. index:: Ab initio Polarization Propagator

.. index::
   pair: ADC; theory

.. _`sec:adc`:

ADC: Ab Initio Polarization Propagator
======================================
.. sectionauthor:: Michael F. Herbst

*Module:* :ref:`Keywords <apdx:adc>`, :ref:`PSI Variables <apdx:adc_psivar>`, :source:`ADC <psi4/src/psi4/adc>`

.. caution:: As of Spring 2022, v1.6, the built-in backend is deprecated and will only be active
   when adcc addon backend is unavailable. In v1.7, built-in ADC will be removed.

Algebraic-diagrammatic construction methods for the polarization propagator (ADC)
determine correlated excitation energies by investigating the pole structure
of said propagator. For this the propagator is expressed in a representation
constructed from so-called intermediate states, which in turn are based
upon a correlated |MollerPlesset| (MP) ground state. The original derivation
of the ADC scheme was purely diagrammatic [Schirmer:1982]_
and the connect to the intermediate states was developed only later [Trofimov:2006]_.
In general :math:`n`-th order ADC theory, ADC(:math:`n`),
is constructed upon an :math:`n`-th order MP ground state.
In this sense one can consider an ADC(:math:`n`) treatment of excited states
consistent to an MP(:math:`n`) perturbation expansion of the ground state.

In ADC methods the residue calculus of the propagator is translated into an eigenvalue
problem with respect to the so-called shifted Hamiltonian or ADC matrix.
Denoting this matrix as :math:`\mathbf{A}`, the eigenproblem can be written
in terms of several blocks

.. math:: \begin{pmatrix}
   \mathbf{A_{SS}} & \mathbf{A_{SD}}\\
   \mathbf{A_{DS}} & \mathbf{A_{DD}}
   \end{pmatrix}
   \begin{pmatrix}
   \mathbf{X_S}\\
   \mathbf{X_D}
   \end{pmatrix}
   =\omega
   \begin{pmatrix}
   \mathbf{X_S}\\
   \mathbf{X_D}
   \end{pmatrix}

where *S* refers to the single and *D* to the double excitation manifolds.
This matrix is typically sparse and thus may be diagonalised iteratively,
for example using Davidson's method [Dreuw:2014:82]_. An alternative viewpoint
has been addressed for example in [Haettig:2002]_, where ADC(2) is related
to other response theories such as CC2-LR, CIS(D) and CIS(D\ :math:`_n`).
In this sense one may consider the ADC matrix the correlated response matrix
to a response problem based on CIS
and apply the simultaneous expansion method (SEM),
in which the |sigma|-vectors (ADC matrix-vector products)
are constructed several times.

The structure and order of the blocks in the equation above
depend on the ADC level employed. With this also the computational cost changes.
The key computational step, namely the formation of the matrix-vector products
scales as :math:`{\cal O}(N^5)` for ADC(2) and :math:`{\cal O}(N^6)`
for ADC(2)-x and ADC(3). Several additional approximations,
such as frozen-core, frozen-virtual
may be applied to reduce the cost of the problem.
Using the core-valence separation (CVS) approximation
one may specifically target core-valence-excitations
at a substantial reduction in cost.
With the spin-flip modification few-reference ground states can
be tackled starting from a triplet reference by simultaneously
exciting an electron and flipping its spin.
A more detailed overview of such modifications gives [Dreuw:2014:82]_
and the `adcc theory documentation <https://adc-connect.org/latest/theory.html>`_.

Available ADC methods
---------------------
.. sectionauthor:: Michael F. Herbst

Several ADC methods are available in |PSIfour| for the computation of excited states,
see :ref:`table:adcsummary`.
The methods are implemented in two distinct codes,
one is part of |PSIfour| itself and will be referred to as the *built-in* implementation,
the other is available via an interface to the `adcc <https://adc-connect.org>`_ python module.
The two backends follow different approaches to compute ADC excited states
and as a result details and supported keywords differ.
After a more general introduction, specific aspects of the two implementations will be highlighted
in section :ref:`sec:interfaceadcc` and :ref:`sec:adcbuiltin`.

.. _`table:adcsummary`:

.. table:: ADC capabilities of Psi4

   +---------------+-----------+------------+---------------+-------+------------------------------------+
   | Method        | Backend   | References | Exc. Energies | Props | Supported values for kind keyword  |
   +===============+===========+============+===============+=======+====================================+
   | ADC(1)        | adcc      | RHF, UHF   | yes           | yes   |  any, singlet, triplet, spin_flip  |
   +---------------+-----------+------------+---------------+-------+------------------------------------+
   | ADC(2)        | adcc      | RHF, UHF   | yes           | yes   |  any, singlet, triplet, spin_flip  |
   +               +-----------+------------+---------------+-------+------------------------------------+
   |               | built-in  | RHF        | yes           | ---   |  singlet                           |
   +---------------+-----------+------------+---------------+-------+------------------------------------+
   | ADC(2)-x      | adcc      | RHF, UHF   | yes           | yes   |  any, singlet, triplet, spin_flip  |
   +---------------+-----------+------------+---------------+-------+------------------------------------+
   | ADC(3)        | adcc      | RHF, UHF   | yes           | yes   |  any, singlet, triplet, spin_flip  |
   +---------------+-----------+------------+---------------+--------------------------------------------+
   | CVS-ADC(1)    | adcc      | RHF, UHF   | yes           | yes   |  any, singlet, triplet             |
   +---------------+-----------+------------+---------------+-------+------------------------------------+
   | CVS-ADC(2)    | adcc      | RHF, UHF   | yes           | yes   |  any, singlet, triplet             |
   +---------------+-----------+------------+---------------+-------+------------------------------------+
   | CVS-ADC(2)-x  | adcc      | RHF, UHF   | yes           | yes   |  any, singlet, triplet             |
   +---------------+-----------+------------+---------------+-------+------------------------------------+
   | CVS-ADC(3)    | adcc      | RHF, UHF   | yes           | yes   |  any, singlet, triplet             |
   +---------------+-----------+------------+---------------+-------+------------------------------------+

The leftmost column of table :ref:`table:adcsummary` provides the supported ADC methods.
If only excitation energies are desired, one can simply pass one
of the listed method strings to the function :py:func:`~psi4.driver.energy`.
For example, ``energy('adc(2)-x')`` will compute
excitation energies at ADC(2)-x level.
Properties such as oscillator strengths, transition or state dipole moments
are available by calling the function :py:func:`~psi4.driver.properties`
with appropriate arguments.
Most commonly users will want to compute at least oscillator strengths
along with the excitation energies,
resulting in a call like ``properties('adc(2)', properties=["oscillator_strength"])``.

Running ADC calculations
------------------------
.. sectionauthor:: Michael F. Herbst

Running an ADC calculation with |PSIfour| requires
the call to :py:func:`~psi4.driver.properties` as discussed above
as well as one or more mandatory keyword arguments.

The most important keyword argument is |adc__roots_per_irrep|,
which is an array with the number of excited states desired
for each irreducible representation. Most ADC methods
are only supported at C1 symmetry at the moment, such that
this option should in most cases be set to an array with a single
element only. For example one can run an ADC(2) calculation for 10
(singlet) excited states using::

   set roots_per_irrep [10]
   properties('adc(2)', properties=["oscillator_strength"])

where the ``molecule`` section was dropped for brevity.

**Selecting the excitation manifold.**
To select between the possible excitation manifolds,
use the |adc__kind| keyword. For restricted references
by default only singlet excited states are computed,
corresponding to the keyword value ``'singlet'``.
To compute triplet states, select ``'triplet'``.
To compute both without making a spin distinction, select ``'any'``.
The latter is default for unrestricted references.

The special |adc__kind| value ``'spin_flip'`` selects
a spin-flip computation where a simultaneous flip of spin
and excitation is performed. This is only available
for unrestricted references and not for ``CVS-ADC(n)`` methods,
see table :ref:`table:adcsummary`.

**Using the core-valence separation.**
For tackling core-valence excitations using the ``CVS-ADC(n)``
methods, the keyword argument |adc__num_core_orbitals|
is additionally required. It is used to specify the number of
(spatial) orbitals to put into the core space and thus select
as target orbitals for a core-valence excitation process.
A value of ``2`` indicates, for example,
that the two lowest-energy :math:`\alpha` and the two
lowest-energy :math:`\beta` orbitals are placed in the core space.
Since the implemented ADC procedures tackle the
lowest-energy excitations, the value should be specified
such that the targeted core orbital is just inside the core space.

*Example:* Consider furane, :math:`C_4H_4O`. In order to tackle
the oxygen 1s edge, *i.e* simulate a O 1s XAS spectrum, one may
just set |adc__num_core_orbitals| to ``1``. This will select the
oxygen 1s orbital for the core space as it is energetically the lowest.
For C 1s core excitations the |adc__num_core_orbitals| value needs
to be set to ``5``, such that both the O 1s and all four C 1s orbitals
are part of the core space.

**Other keywords and examples.**
Apart from the mentioned keywords, the following are common:

.. include:: autodir_options_c/adc__reference.rst
.. include:: autodir_options_c/adc__r_convergence.rst
.. include:: autodir_options_c/adc__num_guesses.rst
.. include:: autodir_options_c/adc__cutoff_amps_print.rst

The full list is provided in appendix :ref:`apdx:adc_psivar`
and many more sample input files can be found in the adc and adcc
subfolders of :source:`samples`.
Note, that not all keywords are supported by all backends.

**Switching between ADC backends.**
Psi4 currently defaults to the built-in implementation for all ADC(2) energy calculations.
You can explicitly set the |globals__qc_module| option to ``'adcc'``
enforce using adcc also for this case.

.. _`sec:interfaceadcc`:

Interface to adcc
-----------------
.. codeauthor:: Michael F. Herbst
.. sectionauthor:: Michael F. Herbst

For most implemented ADC methods |PSIfour| relies
on an interface to the `adcc <https://adc-connect.org>`_ python package.
The approach of adcc is to directly diagonalise the
ADC matrix :math:`\mathbf{A}` in an iterative diagonalisation
procedure, usually a Jacobi-preconditioned Davidson. Expensive parts
of the ADC matrix-vector product are precomputed and stored
in memory. This approach is general in the sense
that it can be applied to a large range of ADC methods and variants.
So far levels up to ADC(3) and CVS-ADC(3) are available
and additional approximations such as
|globals__freeze_core| and |globals__num_frozen_uocc|
are supported with all ADC methods using the adcc backend.

Currently adcc is only capable of performing in-core calculations,
for which, however, permutational symmetry and spin symmetry is taken
into account for both tensor computations and tensor storage.
Inside adcc some heuristic checks for overly excessive memory requirements
are implemented, resulting in a warning in case a
successful execution is unlikely. There are no guarantees for the memory
to be sufficient in case such a warning is not displayed.

More detailed documentation about adcc and its features can be found
at `<https://adc-connect.org>`_,
especially the `theory section <https://adc-connect.org/latest/theory.html>`_.
If you are using adcc from |PSIfour| for your calculations,
please cite both |PSIfour| as well as adcc [Herbst2020]_
in your published work.

**The ADC wavefunction object.**
After running the ADC calculation in adcc, the interface code sets
a number of variables in the returned :py:class:`~psi4.core.Wavefunction`
in case they are computed.
In the following the ``<method>`` prefix refers to the ADC method (such as ``adc(1)``,
``adc(3)``, ``cvs-adc(2)-x``).


* Ground state energy terms like ``MP2 correlation energy``, ``MP3 correlation energy``,
  ``MP2 total energy``, ``MP3 total energy``, ``current correlation energy`` and ``current energy``.
* ``MP2 dipole X`` and the other components: Ground state dipole moments at MP(2) level.
* ``number of iterations``: The number of iterations the iterative solver required to converge.
* ``number of excited states``: The number of excited states, which were computed.
* More variables are summarized in :ref:`apdx:psivariables_alpha`.

The following attribute is set on returned wavefunctions:

* ``adcc_state``: The `adcc.ExcitedStates <https://adc-connect.org/q/excitedstates>`_
  object used by adcc to store the ADC(n) excitation energies and all precomputed data
  in the format used by adcc.
  Provides direct access to analysis and plotting capabilities from adcc.
  For example ``adcc_state.plot_spectrum()`` plots a broadened excited states spectrum
  in matplotlib. See the `adcc calculations documentation <https://adc-connect.org/latest/calculations.html>`_
  for details.

**Tips for convergence issues.**
If you encounter convergence issues inside adcc, the following parameters
are worth tweaking:

* |adc__max_num_vecs|: Specifies the maximal number of subspace vectors
  in the Jacobi-Davidson scheme before a restart occurs. The defaults are usually
  good, but do not be shy to increase this value if you encounter convergence problems.
* |adc__num_guesses|: By default adcc uses twice as many guess vectors as
  states to be computed. Sometimes increasing this value by a few vectors can be helpful.
  If you encounter a convergence to zero eigenvalues, than decreasing this parameter might
  solve the problems.

.. _`sec:adcbuiltin`:

Built-in ADC(2) code
--------------------
.. codeauthor:: Masaaki Saitow
.. sectionauthor:: Masaaki Saitow

.. warning:: The built-in ADC(2) method may give incorrect results if
             multiple roots are requested, due to an error in the Davidson solver,
             and is no longer maintained. It is slated for removal in Psi4 1.7.
             Use of the Psi interface to `adcc` instead is strongly recommended.
             To use this code regardless, either do not have `adcc` installed, or
             set `qc_module builtin`.

The ADC code built into |PSIfour| is capable of ADC(2) computations
of singlet excited states only.
It makes use of the libtrans library for efficient and flexible
integral-transformation and also the libdpd library to
utilize molecular symmetry in the tensorial manipulations in framework
of the direct-product decomposition algorithm. By this feature, the Ritz
space and intermediate tensors are blocked according to the irreducible
representations of the point group, and the excited states that belong
to different symmetry are sought separately.

In the output of ADC, the ADC(2) results may look as follows::

    ->  1 B1 state   :  0.2565095 (a.u.),  6.9799824 (eV)
    Non-iterative:  0.2565636 (a.u.),  6.9814532 (eV)
             Occ Vir        Coefficient
    ---------------------------------------------
              3   0        -0.9017047264
              3   2         0.3038332241
              3   1         0.2907567119
              3   5        -0.0790167706
              3   4        -0.0425829926
              
    Converged in   4 iteration.
    Squared norm of the S component:  0.9315336
    The S vector is rotated up to  8.102 (deg.)

in which the ADC(2) excitation energy is indicated with arrow symbol
and the pseudo-perturbative value, which is calculated in very similar
fashion to the CIS(D) energy, is also presented on the following line. In
this implementation, the ADC(2) secular matrix is treated effectively
by renormalization of the double excitation manifold into the single
excitation manifold. So, the effective secular equation is solved for
several times for the specific state due to the eigenvalue dependence of
the effective response matrix. Only the S component of the transition
amplitude is obtained explicitly and the squared norm of the S block
and the rotation angle from the corresponding CIS vector are given
below the element of the amplitude. The difference between the ADC(2)
value and its non-iterative counterpart is mostly negligible if the
mixture among the CIS excited states is small and the quasi-degeneracy
in the excited state is tolerably weak. But if there is a significant
discrepancy in these energies, or the rotation angle is visibly large,
special care may have to be taken for the strong effects caused by the
higher excited states.


**Partial Renormalization Scheme**

The built-in ADC code is capable of performing the partially-renormalized
ADC(2) computation, termed PR-ADC(2). In the perturbative treatment of
the singly-excited state, the doubly and triply excited configurations
are accounted for as in the case of CIS(D). In the language of 
CIS(D), the former is regarded to introduce the orbital relaxation (OR)
effect while the latter is argued to give rise to the differential
correlation (DC) correction to the excited state. In the PR-ADC(2)
scheme, the DC term is corrected according to the ground state
PR-MP2 correlation, in which the correlation between the electron pairs
is accounted for in size-consistent and unitary-invariant fashion by
modulating the MP1 amplitude. By utilizing the |adc__pr| scheme, substantial
resistance against quasi-degeneracy is readily granted as discussed
in Ref. [Saitow:2012]_.

**Theory of the built-in ADC(2) implementation**

For the built-in ADC(2) implementation
some very essential points shall be emphasized. In ADC(2) specifically
one may write the response equation as

.. math:: \begin{pmatrix}
   \mathbf{A_{SS}^{(2)}} & \mathbf{A_{SD}^{(1)}}\\
   \mathbf{A_{DS}^{(1)}} & \mathbf{A_{DD}^{(0)}}
   \end{pmatrix}
   \begin{pmatrix}
   \mathbf{X_S}\\
   \mathbf{X_D}
   \end{pmatrix}
   =\omega
   \begin{pmatrix}
   \mathbf{X_S}\\
   \mathbf{X_D}
   \end{pmatrix}

where the superscript on each matrix block indicates the order of
the fluctuation. Instead of solving the above equation explicitly,
the large D manifold is treated effectively as

.. math:: [\mathbf{A_{SS}^{(2)}}+
   \mathbf{A_{SD}^{(1)}}^{\dagger}(\omega-
   \mathbf{A_{DD}^{(0)}})^{-1}\mathbf{A_{DS}^{(1)}}]\mathbf{X_{S}}=
   \omega\mathbf{X_{S}}.

This form of the ADC(2) equation requires 7--10 iterations for
convergence on only one root. But thanks to Newton-Raphson
acceleration,

.. math:: \omega^{n+1}=\omega^{n}-
   \frac{\omega^n-\mathbf{X_{S}}(\omega^n)^{\dagger}
   [\mathbf{A_{SS}^{(2)}}+
   \mathbf{A_{SD}^{(1)}}^{\dagger}(\omega^n-\mathbf{A_{DD}^{(0)}})^{-1}
   \mathbf{A_{DS}^{(1)}}]\mathbf{X_{S}}(\omega^n)}{1+\mathbf{X_{S}}
   (\omega^n)^{\dagger}[\mathbf{A_{SD}^{(1)}}^{\dagger}
   (\omega^n-\mathbf{A_{DD}^{(0)}})^{-2}\mathbf{A_{DS}^{(1)}}]\mathbf{X_{S}}
   (\omega^n)}

the computational time reduces to shorter than half of the simple iterative
procedure. Construction of the denominator of the second term in the above
equation is less computationally expensive than construction of one :math:`\sigma`\ -vector
with respect to the effective response matrix. The non-iterative excitation energy stated
above is calculated as a diagonal element of the Davidson mini-Hamiltonian matrix in the SEM as,

.. math:: \omega^{Non-Iterative}=
   \mathbf{X_{CIS}}^{\dagger}[\mathbf{A_{SS}^{(2)}}+
   \mathbf{A_{SD}^{(1)}}^{\dagger}(\omega^{CIS}-\mathbf{A_{DD}^{(0)}})^{-1}
   \mathbf{A_{DS}^{(1)}}]\mathbf{X_{CIS}}

where :math:`\omega^{CIS}` and :math:`\mathbf{X_{CIS}}` denote the CIS 
excitation energy and wave function, respectively. The explicit form of the 
|sigma|-vector is provided in a note accompanying the source code,
in the file :source:`psi4/src/psi4/adc/sigma.pdf`.
