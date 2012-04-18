
.. include:: autodoc_abbr_options_c.rst

.. index:: ADC
.. index:: Ab initio Polarization Propagator

.. index::
   pair: ADC; theory

.. _`sec:adc`:

Ab Initio Polarization Propagator
=================================

.. codeauthor:: Masaaki Saitow
.. sectionauthor:: Masaaki Saitow

*Module:* :ref:`Keywords <apdx:adc>`, :ref:`PSI Variables <apdx:adc_psivar>`, :source:`ADC <src/bin/adc>`

The ADC code seeks the pole structure of the polarization
propagator, which is equivalent to the correlated excitation energy,
in terms of the second order algebraic-diagrammatic construction
(ADC(2)) scheme.  Originally, the ADC scheme was derived in purely
the diagrammatic language by Schirmer [Schirmer:1982]_ and later,
a sophisticated algebraic scheme was developed [Trofimov:2006]_
by Trofimov et al. In general by n-th order ADC theory, the
excited state is treated at completely equivalent level to the 
M\ |o_slash|\ ller--Plesset perturbation expansion of the same order. 
Hence the ADC(2)
can be described as MP2 theory for the excited state and the relation
to such other response theories as CC2-LR, CIS(D) and CIS(D:math:`_n`) has
been addressed [Haettig:2002]_ by Hattig et al.  In the ADC theory,
the residue calculus of the propagator is translated into the eigenvalue
problem with respect to the correlated response matrix, also known as the
shifted-Hamiltonian. The |sigma|-vectors (Hamiltonian-vector products)
are constructed several times in the simultaneous expansion method (SEM)
to solve the eigenvalue problem, and each |sigma|-vector construction
has a computational cost that scales as :math:`{\cal O}(N^5)`. In addition,
the tensorial form of the |sigma|-vector resembles to that of the
doubles correction in the CIS(D) energetic equation. As a consequence,
the pre-factor in the polynomial scaling becomes far larger than that
of the CIS(D) even though the quasi-degeneracy of the excited state is
properly accounted for in the ADC(2) model.

In |PSIfour| the quite efficient and flexible integral-transformation
library named libtrans is newly equipped and utilized in the
production level DCFT code. The ADC code is also based on
libtrans, and it is also based on libdpd, a library to
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

Partial Renormalization Scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ADC code is capable of performing the partially-renormalized
ADC(2) computation, termed PR-ADC(2). In the perturbative treatment of
the singly-excited state, the doubly and triply excited configurations
are accounted for as in the case of CIS(D). In the language of 
CIS(D), the former is regarded to introduce the orbital relaxation (OR)
effect while the latter is argued to give rise to the differential
correlation (DC) correction to the excited state. In the PR-ADC(2)
scheme, the the DC term is corrected according to the ground state
PR-MP2 correlation, in which the correlation between the electron pairs
is accounted for in size-consistent and unitary-invariant fashion by
modulating the MP1 amplitude. By utilizing the |adc__pr| scheme, substantial
resistance against quasi-degeneracy is readily granted as discussed
in Ref. [Saitow:2012]_.

Using the ADC(2) code
~~~~~~~~~~~~~~~~~~~~~

A complete list of keywords related to ADC(2) computations is provided in
Appendix :ref:`apdx:adc`.  Some sample inputs are provided in 
:source:`samples`, in directories starting with the name adc.  The most
important keyword is |adc__roots_per_irrep|, which is an array
giving the number of excited states desired for each irreducible
representation.

Implementation
~~~~~~~~~~~~~~
Some very essential points are emphasized for understanding of the
nature and the limitations of the theory. The ADC(2) response matrix,
denoted as :math:`\mathbf{A}`, is expanded in the single (S) and double (D)
excitation manifolds as

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

This form of the ADC(2) equation requires 7 -- 10 iterations for
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
equation is less computationally expensive than contruction of one $\sigma$-vector with respect to the effective response matrix. The non-iterative excitation energy stated above is calculated as a diagonal element of the Davidson mini-Hamiltonian matrix in the SEM as,

.. math:: \omega^{Non-Iterative}=
   \mathbf{X_{CIS}}^{\dagger}[\mathbf{A_{SS}^{(2)}}+
   \mathbf{A_{SD}^{(1)}}^{\dagger}(\omega^{CIS}-\mathbf{A_{DD}^{(0)}})^{-1}
   \mathbf{A_{DS}^{(1)}}]\mathbf{X_{CIS}}

where :math:`\omega^{CIS}` and :math:`\mathbf{X_{CIS}}` denote the CIS 
excitation energy and wave function, respectively. The explicit form of the 
|sigma|-vector is provided in a note accompanying the source code,
in the file :source:`src/bin/adc/sigma.pdf`. 

