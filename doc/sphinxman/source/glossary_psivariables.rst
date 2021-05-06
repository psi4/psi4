.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2021 The Psi4 Developers.
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

..  * NOTES (LAB 3-26-2012)
    * Any PSI variable added to the codebase should be added to this list
      (variables in the psi variable by module list will show up black
      and un-clickable if an entry isn't present here).
    * INCLUDE UNITS!
    * ALPHABETIZE!

.. include:: autodoc_abbr_options_c.rst

.. _`apdx:psivariables_alpha`:

PSI Variables by Alpha
======================

.. note:: Lowercase letters in PSI variable names represent portions of
   the variable name that vary by root number, calculation order, etc.
   See text for fuller description.

.. psivar:: [T] CORRECTION ENERGY

   The coupled-cluster bracket perturbative triples correction [Eh].

.. psivar:: (T) CORRECTION ENERGY

   The coupled-cluster perturbative triples correction [Eh].

.. psivar:: (AT) CORRECTION ENERGY

   The coupled-cluster asymmetric perturbative triples correction [Eh].

.. psivar:: AAA (T) CORRECTION ENERGY
   AAB (T) CORRECTION ENERGY
   ABB (T) CORRECTION ENERGY
   BBB (T) CORRECTION ENERGY

   Spin components of the UHF-based coupled-cluster perturbative triples correction [Eh].

.. psivar:: ACPF DIPOLE

   Dipole array [e a0] for the averaged coupled-pair functional level of theory, (3,).

.. psivar:: ACPF DIPOLE X
   ACPF DIPOLE Y
   ACPF DIPOLE Z

   The three components of the dipole [Debye] for the
   averaged coupled-pair functional level of theory.
   Deprecated in favor of :psivar:`ACPF DIPOLE`.

.. psivar:: ACPF QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the averaged coupled-pair functional level of theory, (3, 3).

.. psivar:: ACPF QUADRUPOLE XX
   ACPF QUADRUPOLE XY
   ACPF QUADRUPOLE XZ
   ACPF QUADRUPOLE YY
   ACPF QUADRUPOLE YZ
   ACPF QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the
   averaged coupled-pair functional level of theory.
   Deprecated in favor of :psivar:`ACPF QUADRUPOLE`.

.. psivar:: ACPF TOTAL ENERGY
   ACPF CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the averaged coupled-pair functional level of theory.

.. psivar:: AQCC DIPOLE

   Dipole array [e a0] for the averaged quadratic coupled-cluster level of theory, (3,).

.. psivar:: AQCC DIPOLE X
   AQCC DIPOLE Y
   AQCC DIPOLE Z

   The three components of the dipole [Debye] for the
   averaged quadratic coupled-cluster level of theory.
   Deprecated in favor of :psivar:`AQCC DIPOLE`.

.. psivar:: AQCC QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the averaged quadratic coupled-cluster level of theory, (3, 3).

.. psivar:: AQCC QUADRUPOLE XX
   AQCC QUADRUPOLE XY
   AQCC QUADRUPOLE XZ
   AQCC QUADRUPOLE YY
   AQCC QUADRUPOLE YZ
   AQCC QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the
   averaged quadratic coupled-cluster level of theory.
   Deprecated in favor of :psivar:`AQCC QUADRUPOLE`.

.. psivar:: AQCC TOTAL ENERGY
   AQCC CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the averaged quadratic coupled-cluster level of theory.
.. psivar:: BRUECKNER CONVERGED

   Value 1 (0) when the Brueckner orbitals have (have not) converged.

.. psivar:: CBS TOTAL ENERGY
   CBS CORRELATION ENERGY
   CBS REFERENCE ENERGY

   The total electronic energy [Eh] and its breakdown into reference total
   energy [Eh] and correlation correction components [Eh] for the compound
   method requested through cbs().

.. psivar:: CC ROOT n DIPOLE

   Dipole array [e a0] for the requested coupled cluster level of theory and root *n* (number starts at GS = 0), (3,).

.. psivar:: CC ROOT n DIPOLE X
   CC ROOT n DIPOLE Y
   CC ROOT n DIPOLE Z

   The three components of the dipole [Debye] for the requested
   coupled cluster level of theory and root *n* (number starts at GS = 0).
   Deprecated in favor of :psivar:`CC ROOT n DIPOLE`.

.. psivar:: CC ROOT n QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the requested coupled cluster level of theory and root *n* (number starts at GS = 0), (3, 3).

.. psivar:: CC ROOT n QUADRUPOLE XX
   CC ROOT n QUADRUPOLE XY
   CC ROOT n QUADRUPOLE XZ
   CC ROOT n QUADRUPOLE YY
   CC ROOT n QUADRUPOLE YZ
   CC ROOT n QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the requested
   coupled cluster level of theory and root *n* (numbering starts at GS = 0).
   Deprecated in favor of :psivar:`CC ROOT n QUADRUPOLE`.

.. psivar:: CC ROOT n TOTAL ENERGY
   CC ROOT n CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the requested coupled cluster level of theory and root
   *n* (numbering starts at GS = 0).

.. psivar:: CC TOTAL ENERGY
   CC CORRELATION ENERGY

.. psivar:: CC T1 DIAGNOSTIC
   CC D1 DIAGNOSTIC
   CC NEW D1 DIAGNOSTIC
   CC D2 DIAGNOSTIC

   Diagnostic of multireference character.

.. psivar:: CC2 TOTAL ENERGY
   CC2 CORRELATION ENERGY
   CC3 TOTAL ENERGY
   CC3 CORRELATION ENERGY
   CC4 TOTAL ENERGY
   CC4 CORRELATION ENERGY
   CCnn TOTAL ENERGY
   CCnn CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the requested approximate coupled-cluster (CC2, CC3, up to CC\ *nn*)
   level of theory.

.. psivar:: CC DIPOLE

   Dipole array [e a0] for the requested coupled cluster level of theory and root, (3,).

.. psivar:: CC DIPOLE X
   CC DIPOLE Y
   CC DIPOLE Z

   The three components of the dipole [Debye] for the requested
   coupled cluster level of theory and root.
   Deprecated in favor of :psivar:`CC DIPOLE`.

.. psivar:: CC2 DIPOLE POLARIZABILITY @ xNM

   The dipole polarizability [au] calculated at the CC2 level
   for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CC2 SPECIFIC ROTATION (LEN) @ xNM

   The specific rotation [deg/(dm (g/cm^3))] calculated at the CC2 level in the
   length gauge for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CC2 SPECIFIC ROTATION (VEL) @ xNM

   The specific rotation [deg/(dm (g/cm^3))] calculated at the CC2 level in the
   velocity gauge for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CC2 SPECIFIC ROTATION (MVG) @ xNM

   The specific rotation [deg/(dm (g/cm^3))] calculated at the CC2 level in the
   modified velocity gauge for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CC QUADRUPOLE XX
   CC QUADRUPOLE XY
   CC QUADRUPOLE XZ
   CC QUADRUPOLE YY
   CC QUADRUPOLE YZ
   CC QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the requested
   coupled cluster level of theory and root.

.. psivar:: CCD TOTAL ENERGY
   CCD CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the coupled-cluster doubles level of theory.

.. psivar:: CCSD PAIR ENERGIES

   The restricted-reference pair energies for coupled-cluster singles and doubles
   level of theory. Size number of active doubly occupied orbitals, square.

.. psivar:: CCSD TOTAL ENERGY
   CCSD CORRELATION ENERGY
   CCSDT TOTAL ENERGY
   CCSDT CORRELATION ENERGY
   CCSDTQ TOTAL ENERGY
   CCSDTQ CORRELATION ENERGY
   CCn TOTAL ENERGY
   CCn CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the requested full coupled-cluster (CCSD, CCSDT, up to CC\ *n*)
   level of theory.

.. psivar:: CCSD(T) TOTAL ENERGY
   CCSD(T) CORRELATION ENERGY
   CCSD(AT) TOTAL ENERGY
   CCSD(AT) CORRELATION ENERGY
   CCSDT(Q) TOTAL ENERGY
   CCSDT(Q) CORRELATION ENERGY
   CC(n-1)(n) TOTAL ENERGY
   CC(n-1)(n) CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the perturbatively corrected coupled-cluster (CCSD(T), CCSD(AT), CCSDT(Q),
   up to CC(\ *n*\ -1)(\ *n*\ ) level of theory.

.. psivar:: CCSDT-1a TOTAL ENERGY
   CCSDT-1a CORRELATION ENERGY
   CCSDTQ-1a TOTAL ENERGY
   CCSDTQ-1a CORRELATION ENERGY
   CCn-1a TOTAL ENERGY
   CCn-1a CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the approximate coupled-cluster (CCSD(T)-1a, CCSDT(Q)-1a,
   up to CC\ *n*\ -1a) level of theory.

.. psivar:: CCSDT-1b TOTAL ENERGY
   CCSDT-1b CORRELATION ENERGY
   CCSDTQ-1b TOTAL ENERGY
   CCSDTQ-1b CORRELATION ENERGY
   CCn-1b TOTAL ENERGY
   CCn-1b CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the approximate coupled-cluster (CCSD(T)-1b, CCSDT(Q)-1b,
   up to CC\ *n*\ -1b) level of theory.

.. psivar:: CCSDT-3 TOTAL ENERGY
   CCSDT-3 CORRELATION ENERGY
   CCSDTQ-3 TOTAL ENERGY
   CCSDTQ-3 CORRELATION ENERGY
   CCn-3 TOTAL ENERGY
   CCn-3 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the approximate coupled-cluster (CCSD(T)-3, CCSDT(Q)-3,
   up to CC\ *n*\ -3) level of theory.

.. psivar:: CCSD(T)_L TOTAL ENERGY
   CCSD(T)_L CORRELATION ENERGY
   CCSDT(Q)_L TOTAL ENERGY
   CCSDT(Q)_L CORRELATION ENERGY
   CC(n-1)(n)_L TOTAL ENERGY
   CC(n-1)(n)_L CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the approximate coupled-cluster (CCSD(T)_L, CCSDT(Q)_L,
   up to CC(\ *n*\ -1)(\ *n*\ )L level of theory.

.. psivar:: CCSDT(Q)/A TOTAL ENERGY
   CCSDT(Q)/A CORRELATION ENERGY
   CCSDT(Q)/B TOTAL ENERGY
   CCSDT(Q)/B CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the modified CCSDT(Q) level of theory.

.. psivar:: CCSD DIPOLE POLARIZABILITY @ xNM

   The dipole polarizability [au] calculated at the CCSD level
   for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CCSD SPECIFIC ROTATION (LEN) @ xNM

   The specific rotation [deg/(dm (g/cm^3))] calculated at the CCSD level in the
   length gauge for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CCSD SPECIFIC ROTATION (VEL) @ xNM

   The specific rotation [deg/(dm (g/cm^3))] calculated at the CCSD level in the
   velocity gauge for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CCSD SPECIFIC ROTATION (MVG) @ xNM

   The specific rotation [deg/(dm (g/cm^3))] calculated at the CCSD level in the
   modified velocity gauge for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CEPA(0) DIPOLE

   Dipole array [e a0] for the coupled electron pair approximation variant 0 level of theory, (3,).

.. psivar:: CEPA(0) DIPOLE X
   CEPA(0) DIPOLE Y
   CEPA(0) DIPOLE Z

   The three components of the dipole [Debye] for the
   coupled electron pair approximation variant 0 level of theory.
   Deprecated in favor of :psivar:`CEPA(0) DIPOLE`.

.. psivar:: CEPA(0) QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the coupled electron pair approximation variant 0 level of theory, (3, 3).

.. psivar:: CEPA(0) QUADRUPOLE XX
   CEPA(0) QUADRUPOLE XY
   CEPA(0) QUADRUPOLE XZ
   CEPA(0) QUADRUPOLE YY
   CEPA(0) QUADRUPOLE YZ
   CEPA(0) QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the
   coupled electron pair approximation variant 0 level of theory.
   Deprecated in favor of :psivar:`CEPA(0) QUADRUPOLE`.

.. psivar:: CEPA(0) TOTAL ENERGY
   CEPA(0) CORRELATION ENERGY
   CEPA(1) TOTAL ENERGY
   CEPA(1) CORRELATION ENERGY
   CEPA(2) TOTAL ENERGY
   CEPA(2) CORRELATION ENERGY
   CEPA(3) TOTAL ENERGY
   CEPA(3) CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the requested variant of coupled electron pair approximation level of theory.

.. psivar:: CFOUR ERROR CODE

   The non-zero return value from a Cfour execution.

.. psivar:: CI DIPOLE

   Dipole array [e a0] for the requested configuration interaction level of theory, (3,).

.. psivar:: CI DIPOLE X
   CI DIPOLE Y
   CI DIPOLE Z

   The three components of the dipole [Debye] for the requested
   configuration interaction level of theory and root.
   Deprecated in favor of :psivar:`CI DIPOLE`.

.. psivar:: CI QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the requested configuration interaction level of theory, (3, 3).

.. psivar:: CI QUADRUPOLE XX
   CI QUADRUPOLE XY
   CI QUADRUPOLE XZ
   CI QUADRUPOLE YY
   CI QUADRUPOLE YZ
   CI QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the requested
   configuration interaction level of theory and root.
   Deprecated in favor of :psivar:`CI QUADRUPOLE`.

.. psivar:: CI ROOT n -> ROOT m DIPOLE

   Transition dipole array [e a0] between roots *n* and *m* for the requested configuration interaction level of theory, (3,).

.. psivar:: CI ROOT n -> ROOT m DIPOLE X
   CI ROOT n -> ROOT m DIPOLE Y
   CI ROOT n -> ROOT m DIPOLE Z

   The three components of the transition dipole [Debye] between roots *n*
   and *m* for the requested configuration interaction level of theory.
   Deprecated in favor of :psivar:`CI ROOT n -> ROOT m DIPOLE`.

.. psivar:: CI ROOT n -> ROOT m QUADRUPOLE

   Redundant transition quadrupole array [e a0^2] between roots *n* and *m* for the requested configuration interaction level of theory, (3, 3).

.. psivar:: CI ROOT n -> ROOT m QUADRUPOLE XX
   CI ROOT n -> ROOT m QUADRUPOLE XY
   CI ROOT n -> ROOT m QUADRUPOLE XZ
   CI ROOT n -> ROOT m QUADRUPOLE YY
   CI ROOT n -> ROOT m QUADRUPOLE YZ
   CI ROOT n -> ROOT m QUADRUPOLE ZZ

   The three components of the transition quadrupole [Debye Ang] between
   roots *n* and *m* for the requested configuration interaction level of
   theory.
   Deprecated in favor of :psivar:`CI ROOT n -> ROOT m QUADRUPOLE`.

.. psivar:: CI ROOT n DIPOLE

   Dipole array [e a0] for the requested configuration interaction level of theory and root *n*, (3,).

.. psivar:: CI ROOT n DIPOLE X
   CI ROOT n DIPOLE Y
   CI ROOT n DIPOLE Z

   The three components of the dipole [Debye] for the requested
   configuration interaction level of theory and root *n*.
   Deprecated in favor of :psivar:`CI ROOT n DIPOLE`.

.. psivar:: CI ROOT n QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the requested configuration interaction level of theory and root *n*, (3, 3).

.. psivar:: CI ROOT n QUADRUPOLE XX
   CI ROOT n QUADRUPOLE XY
   CI ROOT n QUADRUPOLE XZ
   CI ROOT n QUADRUPOLE YY
   CI ROOT n QUADRUPOLE YZ
   CI ROOT n QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the requested
   configuration interaction level of theory and root *n*.
   Deprecated in favor of :psivar:`CI ROOT n QUADRUPOLE`.

.. psivar:: CI ROOT n TOTAL ENERGY
   CI ROOT n CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the requested configuration interaction level of theory and root
   *n* (numbering starts at 0).

.. psivar:: CI STATE-AVERAGED TOTAL ENERGY
   CI STATE-AVERAGED CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for state-averaged CI/CASSCF levels of theory.

.. psivar:: CI TOTAL ENERGY
   CI CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the requested configuration interaction level of theory and root.

.. psivar:: CISD DIPOLE

   Dipole array [e a0] for the configuration interaction singles and doubles level of theory, (3,).

.. psivar:: CISD DIPOLE X
   CISD DIPOLE Y
   CISD DIPOLE Z

   The three components of the dipole [Debye] for the
   configuration interaction singles and doubles level of theory and root.
   Deprecated in favor of :psivar:`CISD DIPOLE`.

.. psivar:: CISD QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the configuration interaction singles and doubles level of theory, (3, 3).

.. psivar:: CISD QUADRUPOLE XX
   CISD QUADRUPOLE XY
   CISD QUADRUPOLE XZ
   CISD QUADRUPOLE YY
   CISD QUADRUPOLE YZ
   CISD QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the
   configuration interaction singles and doubles level of theory and root.
   Deprecated in favor of :psivar:`CISD QUADRUPOLE`.

.. psivar:: CISD TOTAL ENERGY
   CISD CORRELATION ENERGY
   CISDT TOTAL ENERGY
   CISDT CORRELATION ENERGY
   CISDTQ CORRELATION ENERGY
   CISDTQ TOTAL ENERGY
   CIn CORRELATION ENERGY
   CIn TOTAL ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the labeled configuration interaction level of theory and root.
   *n* is CI order for *n* > 4.

.. psivar:: CP-CORRECTED 2-BODY INTERACTION ENERGY

   The interaction energy [Eh] considering only two-body interactions,
   computed with counterpoise correction.
   Related variable :psivar:`UNCP-CORRECTED 2-BODY INTERACTION ENERGY`.

   .. math:: E_{\text{IE}} = E_{dimer} - \sum_{monomer}^{n}{E_{monomer}^{\text{CP}}}

.. psivar:: CURRENT CORRELATION ENERGY

   The correlation energy [Eh] corresponding to the :psivar:`CURRENT ENERGY` variable.

.. psivar:: CURRENT ENERGY

   The total electronic energy [Eh] of the most recent stage of a
   calculation (frequently overwritten). This is the quantity tracked by
   the geometry optimizer.

.. psivar:: CURRENT REFERENCE ENERGY

   The total electronic energy [Eh] of the reference stage corresponding to
   the :psivar:`CURRENT ENERGY` variable.

.. psivar:: CURRENT DIPOLE

   The total dipole [e a0] of the most recent stage of a calculation (frequently overwritten), (3,).

.. psivar:: CURRENT GRADIENT

   The total electronic gradient [E_h/a0] of the most recent stage of a
   calculation (frequently overwritten). This is the quantity tracked by
   the geometry optimizer, ({nat}, 3).

.. psivar:: CURRENT DIPOLE GRADIENT

   The derivative of the dipole with respect to nuclear perturbations [E_h a0/u] = [(e a0/a0)^2/u]
   as a degree-of-freedom by dipole component array, (3 * {nat}, 3).

.. psivar:: CURRENT HESSIAN

   The total electronic Hessian [E_h/a0/a0] of the most recent stage of a
   calculation, (3 * {nat}, 3 * {nat}).

.. psivar:: CUSTOM SCS-MP2 TOTAL ENERGY
   CUSTOM SCS-MP2 CORRELATION ENERGY

   Changeable quantities based on options.
   The total electronic energy [Eh] and correlation energy component [Eh]
   for the MP2-like method formed by any reweighting of :psivar:`MP2 DOUBLES ENERGY`
   for opposite-spin and same-spin contributions, with
   any singles carried along.
   Depending on weights, may equal any of MP2, SCS-MP2, SCS(N)-MP2, etc. quantities.
   Contrast with :psivar:`SCS-MP2 TOTAL ENERGY`.

.. psivar:: CUSTOM SCS-MP2.5 TOTAL ENERGY
   CUSTOM SCS-MP2.5 CORRELATION ENERGY
   CUSTOM SCS-MP3 TOTAL ENERGY
   CUSTOM SCS-MP3 CORRELATION ENERGY
   CUSTOM SCS-LCCD TOTAL ENERGY
   CUSTOM SCS-LCCD CORRELATION ENERGY
   CUSTOM SCS-OMP2 TOTAL ENERGY
   CUSTOM SCS-OMP2 CORRELATION ENERGY
   CUSTOM SCS-OMP2.5 TOTAL ENERGY
   CUSTOM SCS-OMP2.5 CORRELATION ENERGY
   CUSTOM SCS-OMP3 TOTAL ENERGY
   CUSTOM SCS-OMP3 CORRELATION ENERGY
   CUSTOM SCS-OLCCD TOTAL ENERGY
   CUSTOM SCS-OLCCD CORRELATION ENERGY

   Changeable quantities based on options.
   The total electronic energy [Eh] and correlation energy component [Eh]
   for the method formed by any reweighting of the named :samp:`{method} DOUBLES ENERGY`
   for opposite-spin and same-spin contributions, with
   any singles carried along.
   Contrast with :samp`SCS-{method} TOTAL ENERGY`.

.. psivar:: db_name DATABASE MEAN ABSOLUTE DEVIATION

   The mean absolute deviation [\ |kcalpermol|\ ] of the requested method
   *name* from the stored reference values for the requested reactions in
   database *db_name*. If no reference is available, this will be a large
   and nonsensical value.

   .. math:: \frac{1}{n}\sum_{rxn}^{n}{| \textsf{\textsl{name}}_{rxn}-\text{REF}_{rxn} | }

.. psivar:: db_name DATABASE MEAN SIGNED DEVIATION

   The mean deviation [\ |kcalpermol|\ ] of the requested method *name*
   from the stored reference values for the requested reactions in
   database *db_name*. If no reference is available, this will be a large
   and nonsensical value.

   .. math:: \frac{1}{n}\sum_{rxn}^{n}{\textsf{\textsl{name}}_{rxn}-\text{REF}_{rxn}}

.. psivar:: db_name DATABASE ROOT-MEAN-SQUARE DEVIATION

   The rms deviation [\ |kcalpermol|\ ] of the requested method *name*
   from the stored reference values for the requested reactions in
   database *db_name*. If no reference is available, this will be a large
   and nonsensical value.

   .. math:: \sqrt{\frac{1}{n}\sum_{rxn}^{n}{(\textsf{\textsl{name}}_{rxn}-\text{REF}_{rxn})^2}}

.. psivar:: DCT LAMBDA ENERGY

   An energy term in density cumulant theory [Eh]. This term is the
   2-electron cumulant's contribution contribution to the reduced
   density matrix energy expression. Not recommended for interpretative
   use except by reduced density matrix specialists.

.. psivar:: DCT SCF ENERGY

   An energy term in density cumulant theory [Eh]. This term is the
   1-electron reduced density matrix (1RDM) contribution to the reduced
   density matrix energy expression, plus the contribution of the
   antisymmetrized product of 1RDMs. Not recommended for interpretative
   use except by reduced density matrix specialists.

.. psivar:: DCT THREE-PARTICLE ENERGY

   The three-particle correlation energy correction [Eh] in density cumulant
   theory, akin to :psivar:`(T) CORRECTION ENERGY` in coupled-cluster.

.. psivar:: DCT TOTAL ENERGY

   Total energy [Eh] in density cumulant theory. Sum of :psivar:`DCT SCF ENERGY`,
   :psivar:`DCT LAMBDA ENERGY`, and :psivar:`DCT THREE-PARTICLE ENERGY` when present.

.. psivar:: DETCI AVG DVEC NORM

   A measure of configuration interaction convergence.

.. psivar:: DFT FUNCTIONAL TOTAL ENERGY

   The total electronic energy [Eh] for the underlying functional of the
   requested DFT method, without any dispersion correction; the first four
   terms in Eq. :eq:`SCFterms` or :eq:`DFTterms`. Quantity
   :math:`E_{\text{FCTL}}` in Eqs.  :eq:`SCFterms` and :eq:`DFTterms`.
   Unless the method includes a dispersion correction, this quantity is
   equal to :psivar:`SCF TOTAL ENERGY`.

.. psivar:: DFT TOTAL ENERGY

   The total electronic energy [Eh] for the requested DFT method,
   :math:`E_{\text{DFT}}` in Eq. :eq:`DFTterms`.

   .. math::
      :nowrap:
      :label: DFTterms

         \begin{align*}
            E_{\text{DFT}} & = E_{NN} + E_{1e^-} + E_{2e^-} + E_{xc} + E_{\text{-D}} + E_{\text{DH}} \\
                           & = E_{\text{FCTL}} + E_{\text{-D}} + E_{\text{DH}} \\
                           & = E_{\text{SCF}} + E_{\text{DH}}
         \end{align*}

   Unless the method is a DFT double-hybrid, this quantity is equal to
   :psivar:`SCF TOTAL ENERGY`. If the method is neither a
   double-hybrid, nor dispersion corrected, this quantity is equal to
   :psivar:`DFT FUNCTIONAL TOTAL ENERGY`.

.. psivar:: DFT TOTAL GRADIENT

   The total electronic gradient [E_h/a0] of the requested DFT method, ({nat}, 3).

.. psivar:: DFT DIPOLE GRADIENT

   The derivative of the requested DFT method dipole [E_h a0/u] = [(e a0/a0)^2/u] with respect to nuclear perturbations
   as a degree-of-freedom by dipole component array, (3 * {nat}, 3).

.. psivar:: DFT TOTAL HESSIAN

   The total electronic second derivative [Eh/a0/a0] for the requested DFT method, (3 * {nat}, 3 * {nat}).

.. psivar:: DFT XC ENERGY

   The functional energy contribution [Eh] to the total SCF energy (DFT only).
   Quantity :math:`E_{xc}` in Eqs. :eq:`SCFterms` and :eq:`DFTterms`.

.. psivar:: DFT VV10 ENERGY

   The VV10 nonlocal contribution [Eh] to the total SCF energy (DFT only).
   Included in :psivar:`DFT FUNCTIONAL TOTAL ENERGY`.

.. psivar:: DISPERSION CORRECTION ENERGY
   fctl DISPERSION CORRECTION ENERGY

   The dispersion correction [Eh] appended to an underlying functional
   when a DFT-D method is requested. Quantity :math:`E_{\text{-D}}`
   in Eqs. :eq:`SCFterms` and :eq:`DFTterms`.
   When dispersion parameters are untweaked for a functional and dispersion
   level, labeled QCVariable also defined.

.. psivar:: DOUBLE-HYBRID CORRECTION ENERGY

   The scaled MP2 correlation energy correction [Eh] appended to an
   underlying functional when a DH-DFT method is requested.
   Quantity :math:`E_{\text{DH}}` in Eq. :eq:`DFTterms`.

.. psivar:: DMA DISTRIBUTED MULTIPOLES

   Distributed multipoles in units given by |gdma__gdma_multipole_units|
   with the row index corresponding to the site and the column index
   referencing the multipole component. Both indices are zero based,
   and the Qlm components of the multipoles are ordered as Q00, Q10,
   Q11c, Q11s, Q20, Q21c, Q21s, Q22c, Q22s, etc.

.. psivar:: DMA TOTAL MULTIPOLES

   Distributed multipoles as a single row, whose columns are the total
   multipoles, translated to |gdma__gdma_origin|, and summed.

.. psivar:: DMRG-SCF TOTAL ENERGY

   The total DMRG total electonic energy [Eh]. Not unique because oribital spaces vary.

.. psivar:: DMRG-CASPT2 TOTAL ENERGY

   The total DMRG plus CASPT2 total electonic energy [Eh] . Not unique because orbital spaces vary.

.. psivar:: EFP DISP ENERGY
   EFP ELST ENERGY
   EFP EXCH ENERGY
   EFP IND ENERGY

   Respectively, the dispersion, electrostatics, exchange, and induction
   components of the total electronic interaction energy [Eh] for EFP/EFP
   computations. The sum of these four components yields
   :psivar:`EFP TOTAL ENERGY`.

.. psivar:: EFP TOTAL ENERGY

   The total electronic interaction energy [Eh] for EFP/EFP computations.

.. psivar:: EFP TORQUE

   The torque, not gradient for EFP/EFP computations.

.. psivar:: ENTHALPY

   Total enthalpy H [Eh] at given temperature.

.. psivar:: ENTHALPY CORRECTION

   Sum of electronic, translational, rotational, and vibrational corrections [Eh] to the enthalpy at given temperature.

.. psivar:: ESP AT CENTER n

   Property of electrostatic potential [Eh / e] at location, usually atom center, n.

.. psivar:: FCI TOTAL ENERGY
   FCI CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the full configuration interaction level of theory.

.. psivar:: GIBBS FREE ENERGY

   Total Gibbs free energy [Eh], free enthalpy at given temperature.

.. psivar:: GIBBS FREE ENERGY CORRECTION

   Sum of electronic, translational, rotational, and vibrational corrections [Eh] to the free enthalpy at given temperature.

.. psivar:: GRID ELECTRONS TOTAL
   GRID ELECTRONS ALPHA
   GRID ELECTRONS BETA

   The number of electrons integrated by the xc quadrature grid.

.. psivar:: HF TOTAL ENERGY

   The total electronic energy [Eh] for the Hartree--Fock method, without
   any dispersion correction; the first three (or four, since
   :math:`E_{xc} = 0`) terms in Eq. :eq:`SCFterms`. Quantity :math:`E_{\text{HF}}`
   in Eq. :eq:`SCFterms`.

.. psivar:: HF TOTAL GRADIENT

   The total electronic gradient [E_h/a0] of the Hartree--Fock method, ({nat}, 3).

.. psivar:: HF DIPOLE GRADIENT

   The derivative of the Hartree--Fock method dipole [E_h a0/u] = [(e a0/a0)^2/u] with respect to nuclear perturbations
   as a degree-of-freedom by dipole component array, (3 * {nat}, 3).

.. psivar:: HF TOTAL HESSIAN

   The total electronic second derivative [Eh/a0/a0] for the Hartree-Fock method, (3 * {nat}, 3 * {nat}).

.. psivar:: LCCD TOTAL ENERGY
   LCCD CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the linearized coupled cluster doubles level of theory.

.. psivar:: LCCSD TOTAL ENERGY
   LCCSD CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the linearized coupled cluster singles and doubles level of theory.

.. psivar:: LCC2 (+LMP2) TOTAL ENERGY

   The total electronic energy [Eh] for the local CC2 level of theory.

.. psivar:: LCCSD (+LMP2) TOTAL ENERGY

   The total electronic energy [Eh] for the local CCSD level of theory.

.. psivar:: LOWDIN CHARGES

   Property of partial atomic charges [e] by the method of L\ |o_dots|\ wdin, (nat,).

.. psivar:: MAYER INDICES

   Property of Mayer bond indices, (nat, nat).

.. psivar:: MBIS CHARGES
   MBIS DIPOLES
   MBIS OCTUPOLES
   MBIS QUADRUPOLES

   Per-atom charges [e], dipoles [e a0], quadrupoles [e a0^2], and octupoles [e a0^3]
   resulting from partitioning the total electron density through the Minimal Basis
   Iterative Stockholder (MBIS) Charge Partitioning Scheme.

.. psivar:: MBIS FREE ATOM n VOLUME

   Free-atom volume [a0^3] for atom n, computed using the MBIS charge
   partitioning scheme. Free atom densities are computed at the same
   level of theory as the parent MBIS calculation, with UHF turned on
   as needed.

.. psivar:: MBIS RADIAL MOMENTS <R^3>

   Per-atom expectation value of r^3 [a0^3], equivalent to the volume
   of the MBIS-partitioned density.

.. psivar:: MBIS VALENCE WIDTHS

   Per-atom density width [a0] of the associated valence charge computed
   from an MBIS partitioned density. Equivalent to the inverse of the
   linear decay rate of the atomic density.

.. psivar:: MBIS VOLUME RATIOS

   Per-atom ratio between the atomic volume (<R^3>) and the free-atomic
   volume, unitless.

.. psivar:: MCSCF TOTAL ENERGY

   Multiconfigurational self-consistent-field energy [Eh] in the course of
   a configuration interaction computation. May be single-root or state-averaged.

.. psivar:: mtd DIPOLE

   Dipole array [e a0] for the named method, (3,).

.. psivar:: mtd DIPOLE X
   mtd DIPOLE Y
   mtd DIPOLE Z

   The three components of the named method dipole [Debye].
   Deprecated in favor of :psivar:`mtd DIPOLE`.

.. psivar:: mtd QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the named method, (3, 3).

.. psivar:: mtd QUADRUPOLE XX
   mtd QUADRUPOLE XY
   mtd QUADRUPOLE XZ
   mtd QUADRUPOLE YY
   mtd QUADRUPOLE YZ
   mtd QUADRUPOLE ZZ

   The six components of the named method quadrupole [Debye Ang].
   Deprecated in favor of :psivar:`mtd QUADRUPOLE`.

.. psivar:: mtd OCTUPOLE

   Redundant octupole array [e a0^3] for the named method, (3, 3, 3).

.. psivar:: mtd OCTUPOLE XXX
   mtd OCTUPOLE XXY
   mtd OCTUPOLE XXZ
   mtd OCTUPOLE XYY
   mtd OCTUPOLE XYZ
   mtd OCTUPOLE XZZ
   mtd OCTUPOLE YYY
   mtd OCTUPOLE YYZ
   mtd OCTUPOLE YZZ
   mtd OCTUPOLE ZZZ

   The ten components of the named method octupole [Debye Ang^2].
   Deprecated in favor of :psivar:`mtd OCTUPOLE`.

.. psivar:: mtd HEXADECAPOLE

   Redundant hexadecapole array [e a0^4] for the named method, (3, 3, 3, 3).

.. psivar:: mtd HEXADECAPOLE XXXX
   mtd HEXADECAPOLE XXXY
   mtd HEXADECAPOLE ZZZZ

   The 15 components of the named method hexadecapole [Debye Ang^3].
   Deprecated in favor of :psivar:`mtd HEXADECAPOLE`.

.. psivar:: mtd 32-POLE

   Redundant 32-pole array [e a0^5] for the named method, (3, 3, 3, 3, 3).

.. psivar:: mtd 32-POLE XXXXX
   mtd 32-POLE XXXXY
   mtd 32-POLE ZZZZZ

   The 21 components of the named method 32-pole [Debye Ang^4].
   Deprecated in favor of :psivar:`mtd 32-POLE`.

.. psivar:: mtd 64-POLE

   Redundant 64-pole array [e a0^6] for the named method, (3, 3, 3, 3, 3, 3).

.. psivar:: mtd 64-POLE XXXXXX
   mtd 64-POLE XXXXXY
   mtd 64-POLE ZZZZZZ

   The 28 components of the named method 64-pole [Debye Ang^5].
   Deprecated in favor of :psivar:`mtd 64-POLE`.

.. psivar:: mtd 128-POLE

   Redundant 128-pole array [e a0^7] for the named method, (3, 3, 3, 3, 3, 3, 3).

.. psivar:: mtd 128-POLE XXXXXXX
   mtd 128-POLE XXXXXXY
   mtd 128-POLE ZZZZZZZ

   The 36 components of the named method 128-pole [Debye Ang^6].
   Deprecated in favor of :psivar:`mtd 128-POLE`.

.. psivar:: MP2 TOTAL ENERGY
   MP2 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the MP2 level of theory.

.. psivar:: MP2 TOTAL GRADIENT
   The total electronic gradient [E_h/a0] of the MP2 level of theory, ({nat}, 3).

.. psivar:: MP2 DIPOLE GRADIENT

   The derivative of the MP2 level of theory dipole [E_h a0/u] = [(e a0/a0)^2/u] with respect to nuclear perturbations
   as a degree-of-freedom by dipole component array, (3 * {nat}, 3).

.. psivar:: MP2 TOTAL HESSIAN

   The total electronic second derivative [Eh/a0/a0] for the MP2 level of theory, (3 * {nat}, 3 * {nat}).

.. psivar:: MP2.5 TOTAL ENERGY
   MP2.5 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the MP2.5 level of theory.

.. psivar:: MP3 TOTAL ENERGY
   MP3 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the MP3 level of theory.

.. psivar:: MP4(T) CORRECTION ENERGY

   The MP4 triples component [Eh]. Quantity is second right-hand term in
   Eq. :eq:`MP4terms`.

.. psivar:: MP4(SDQ) TOTAL ENERGY
   MP4(SDQ) CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the MP4 singles, doubles, quadruples level of theory.  Quantity
   :psivar:`MP4(SDQ) CORRELATION ENERGY` is
   first right-hand term in Eq. :eq:`MP4terms`.

.. psivar:: MP4 TOTAL ENERGY
   MP4 CORRELATION ENERGY
   MP4(SDTQ) TOTAL ENERGY
   MP4(SDTQ) CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the full MP4 level of theory. Quantity :psivar:`MP4 CORRELATION
   ENERGY` / :psivar:`MP4(SDTQ) CORRELATION ENERGY`
   is left-hand term in Eq. :eq:`MP4terms`.

   .. math:: E_{\text{MP4}} = E_{\text{MP4(SDQ)}} + E_{\text{MP4(T)}}
      :label: MP4terms

.. psivar:: MPn TOTAL ENERGY
   MPn CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the labeled |MollerPlesset| perturbation theory level.
   *n* is MP perturbation order.

.. psivar:: MP2 DOUBLES ENERGY
   MP2.5 DOUBLES ENERGY
   MP3 DOUBLES ENERGY
   CEPA(0) DOUBLES ENERGY
   CEPA(1) DOUBLES ENERGY
   CEPA(2) DOUBLES ENERGY
   CEPA(3) DOUBLES ENERGY
   CISD DOUBLES ENERGY
   QCISD DOUBLES ENERGY
   LCCD DOUBLES ENERGY
   CCD DOUBLES ENERGY
   LCCSD DOUBLES ENERGY
   CCSD DOUBLES ENERGY
   OMP2 DOUBLES ENERGY
   OMP2.5 DOUBLES ENERGY
   OMP3 DOUBLES ENERGY
   OLCCD DOUBLES ENERGY

   The doubles portion [Eh] of the named correlation energy
   including same-spin and opposite-spin correlations.

.. psivar:: MP2 SINGLES ENERGY
   MP2.5 SINGLES ENERGY
   MP3 SINGLES ENERGY
   CEPA(0) SINGLES ENERGY
   CEPA(1) SINGLES ENERGY
   CEPA(2) SINGLES ENERGY
   CEPA(3) SINGLES ENERGY
   CISD SINGLES ENERGY
   QCISD SINGLES ENERGY
   LCCD SINGLES ENERGY
   CCD SINGLES ENERGY
   LCCSD SINGLES ENERGY
   CCSD SINGLES ENERGY
   OLCCD SINGLES ENERGY

   The singles portion [Eh] of the named correlation energy.
   Zero except in ROHF.

.. psivar:: MP2 SAME-SPIN CORRELATION ENERGY
   MP2.5 SAME-SPIN CORRELATION ENERGY
   MP3 SAME-SPIN CORRELATION ENERGY
   CEPA(0) SAME-SPIN CORRELATION ENERGY
   CEPA(1) SAME-SPIN CORRELATION ENERGY
   CEPA(2) SAME-SPIN CORRELATION ENERGY
   CEPA(3) SAME-SPIN CORRELATION ENERGY
   CISD SAME-SPIN CORRELATION ENERGY
   QCISD SAME-SPIN CORRELATION ENERGY
   ACPF SAME-SPIN CORRELATION ENERGY
   AQCC SAME-SPIN CORRELATION ENERGY
   LCCD SAME-SPIN CORRELATION ENERGY
   CCD SAME-SPIN CORRELATION ENERGY
   LCCSD SAME-SPIN CORRELATION ENERGY
   CCSD SAME-SPIN CORRELATION ENERGY
   OLCCD SAME-SPIN CORRELATION ENERGY

   The unscaled portion [Eh] of the named correlation energy
   from same-spin or triplet doubles correlations.

.. psivar:: MP2 OPPOSITE-SPIN CORRELATION ENERGY
   MP2.5 OPPOSITE-SPIN CORRELATION ENERGY
   MP3 OPPOSITE-SPIN CORRELATION ENERGY
   CEPA(0) OPPOSITE-SPIN CORRELATION ENERGY
   CEPA(1) OPPOSITE-SPIN CORRELATION ENERGY
   CEPA(2) OPPOSITE-SPIN CORRELATION ENERGY
   CEPA(3) OPPOSITE-SPIN CORRELATION ENERGY
   CISD OPPOSITE-SPIN CORRELATION ENERGY
   QCISD OPPOSITE-SPIN CORRELATION ENERGY
   ACPF OPPOSITE-SPIN CORRELATION ENERGY
   AQCC OPPOSITE-SPIN CORRELATION ENERGY
   LCCD OPPOSITE-SPIN CORRELATION ENERGY
   CCD OPPOSITE-SPIN CORRELATION ENERGY
   LCCSD OPPOSITE-SPIN CORRELATION ENERGY
   CCSD OPPOSITE-SPIN CORRELATION ENERGY
   OLCCD OPPOSITE-SPIN CORRELATION ENERGY

   The unscaled portion [Eh] of the named correlation energy
   from opposite-spin or singlet doubles correlations.

.. psivar:: MRPT TOTAL ENERGY
   MP2-CCSD TOTAL ENERGY
   MRCC TOTAL ENERGY

   Energies [Eh] from correlated multi-reference theories.

.. psivar:: MULLIKEN CHARGES

   Property of partial atomic charges [e] by the method of Mulliken, (nat,).

.. psivar:: NAUX (SCF)
   NAUX (CC)

   Convenience storage of number of functions [] in the auxiliary basis
   set for named stage of the calculation.

.. psivar:: NBODY (i, j, ..., k)@(a, b, ..., c) TOTAL ENERGY

   The total energy [Eh] of a component of the requested N-Body energy.
   The first parenthetical list over *i*, *j*, ..., *k* enumerates
   molecular fragments included in the computation in 1-indexed,
   input-file order, while the second enumerates list over *a*, *b*,
   ..., *c* enumerates which fragments contribute basis functions to the
   computation.  For example, ``(1, 2)@(1, 2, 3, 4)`` indicates that the
   fragments 1 and 2 are explicitly included in the energy computation,
   with basis functions from each of fragments 1, 2, 3, & 4 included in
   the basis set.  Therefore, the basis functions from fragments 3 and 4
   are included as ghost functions within the energy computation.

.. psivar:: NUCLEAR REPULSION ENERGY

   The nuclear repulsion energy contribution [Eh] to the total SCF energy.
   Quantity :math:`E_{NN}` in Eq. :eq:`SCFterms`.

   .. math:: E_{NN} = \sum_{i, j<i}^{N_{atom}}\frac{Z_i Z_j}{|\mathbf{R}_i - \mathbf{R}_j|}
      :label: ENN

.. psivar:: OCEPA(0) TOTAL ENERGY
   OCEPA(0) CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the orbital-optimized CEPA(0) level of theory.

.. psivar:: OLCCD TOTAL ENERGY
   OLCCD CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the orbital-optimized linearized coupled cluster doubles level of theory.

.. psivar:: OLCCD REFERENCE CORRECTION ENERGY

   The additional correction to the SCF reference energy [Eh]
   for the orbital-optimized linearized coupled cluster doubles level of theory.

.. psivar:: OMP2 TOTAL ENERGY
   OMP2 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the orbital-optimized MP2 level of theory.

.. psivar:: OMP2.5 TOTAL ENERGY
   OMP2.5 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the orbital-optimized MP2.5 level of theory.

.. psivar:: OMP3 TOTAL ENERGY
   OMP3 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the orbital-optimized MP3 level of theory.

.. psivar:: ONE-ELECTRON ENERGY

   The one-electron energy contribution [Eh] to the total SCF energy.
   Quantity :math:`E_{1e^-}` in Eq. :eq:`SCFterms`.

.. psivar:: PCM POLARIZATION ENERGY

   The energy contribution [Eh] from the polarizable continuum model for solvation.

.. psivar:: PE ENERGY

   The energy contribution [Eh] from the polarizable embedding model for solvation.

.. psivar:: QCISD TOTAL ENERGY
   QCISD CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the quadratic configuration interaction singles and doubles level
   of theory.

.. psivar:: QCISD(T) TOTAL ENERGY
   QCISD(T) CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the quadratic configuration interaction singles and doubles with
   perturbative triples correction level of theory.

.. psivar:: SAPT DISP ENERGY
   SAPT ELST ENERGY
   SAPT EXCH ENERGY
   SAPT IND ENERGY

   Respectively, the dispersion, electrostatics, exchange, and induction
   components of the total electronic interaction energy [Eh] for the
   requested SAPT level of theory. The sum of these four components yields
   :psivar:`SAPT TOTAL ENERGY`.

.. psivar:: SAPT TOTAL ENERGY
   SAPT ENERGY

   The total electronic interaction energy [Eh] for the requested SAPT
   level of theory.

.. psivar:: SAPT ELST10,R ENERGY

   An electrostatics-classified SAPT term energy [Eh] implemented for SAPT0.

.. psivar:: SAPT ELST EXTERN-EXTERN ENERGY

   Electrostatic interaction [Eh] between the point charges in fragments
   A and B in F/I-SAPT.

.. psivar:: SAPT EXCH10 ENERGY

   An exchange-classified SAPT term energy [Eh] implemented for SAPT0.

.. psivar:: SAPT EXCH10(S^2) ENERGY

   An exchange-classified SAPT term energy [Eh] implemented for SAPT0.

.. psivar:: SAPT IND20,R ENERGY
   SAPT EXCH-IND20,R ENERGY
   SAPT IND20,U ENERGY
   SAPT EXCH-IND20,U ENERGY

   An induction-classified SAPT term energy [Eh] implemented for SAPT0.

.. psivar:: SAPT DISP20 ENERGY
   SAPT EXCH-DISP20 ENERGY

   A dispersion-classified SAPT term energy [Eh] implemented for SAPT0.

.. psivar:: SAPT EXCH-DISP20(S^INF) ENERGY

   A dispersion-classified SAPT term energy [Eh] implemented for SAPT0. See :ref:`sec:saptinf`.

.. psivar:: SAPT SAME-SPIN DISP20 ENERGY
   SAPT SAME-SPIN EXCH-DISP20 ENERGY

   The portion of :psivar:`SAPT DISP20 ENERGY` or
   :psivar:`SAPT EXCH-DISP20 ENERGY` resulting from
   from same-spin or triplet doubles correlations.

.. psivar:: SAPT HF(2) ENERGY ABC(HF)

   The total Hartree--Fock energy [Eh] of the supersystem implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY AC(0)

   The Hartree--Fock energy [Eh] of subsystems A and C implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY BC(0)

   The Hartree--Fock energy [Eh] of subsystems B and C implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY A(0)

   The Hartree--Fock energy [Eh] of subsystem A implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY B(0)

   The Hartree--Fock energy [Eh] of subsystem B implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY AC(HF)

   The Hartree--Fock localized energy [Eh] of subsystems A and C implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY BC(HF)

   The Hartree--Fock localized energy [Eh] of subsystems B and C implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY AB(HF)

   The Hartree--Fock localized energy [Eh] of subsystems A and B implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY A(HF)

   The Hartree--Fock localized energy [Eh] of subsystem A implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY B(HF)

   The Hartree--Fock localized energy [Eh] of subsystem B implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY C

   The Hartree--Fock energy [Eh] of subsystem C implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY HF

   The FI-SAPT Hartree--Fock interaction energy [Eh] implemented for F/I-SAPT.

.. psivar:: SAPT ELST12,R ENERGY

   An electrostatics-classified SAPT term energy [Eh] implemented for SAPT2.

.. psivar:: SAPT EXCH11(S^2) ENERGY
   SAPT EXCH12(S^2) ENERGY

   An exchange-classified SAPT term energy [Eh] implemented for SAPT2.

.. psivar:: SAPT IND22 ENERGY
   SAPT EXCH-IND22 ENERGY

   An induction-classified SAPT term energy [Eh] implemented for SAPT2.

.. .. psivar:: SAPT HF TOTAL ENERGY
.. .. psivar:: SAPT CT ENERGY

.. psivar:: SAPT DISP21 ENERGY

   A dispersion-classified SAPT term energy [Eh] implemented for SAPT2+.

.. psivar:: SAPT DISP22(SDQ) ENERGY
   SAPT DISP22(T) ENERGY
   SAPT EST.DISP22(T) ENERGY

   Dispersion-classified MBPT-based SAPT term energy [Eh] implemented for SAPT2+.

.. psivar:: SAPT DISP2(CCD) ENERGY
   SAPT DISP22(S)(CCD) ENERGY
   SAPT DISP22(T)(CCD) ENERGY
   SAPT EST.DISP22(T)(CCD) ENERGY

   Dispersion-classified coupled-cluster-based SAPT term energy [Eh] implemented for SAPT2+.

.. psivar:: SAPT ELST13,R ENERGY

   An electrostatics-classified SAPT term energy [Eh] implemented for SAPT2+(3).

.. psivar:: SAPT IND30,R ENERGY
   SAPT IND-DISP30 ENERGY
   SAPT EXCH-IND30,R ENERGY

   A induction-classified SAPT term energy [Eh] implemented for SAPT2+3.

.. psivar:: SAPT DISP30 ENERGY
   SAPT EXCH-DISP30 ENERGY
   SAPT EXCH-IND-DISP30 ENERGY

   A dispersion-classified SAPT term energy [Eh] implemented for SAPT2+3.

.. psivar:: SAPT ALPHA

   SAPT exchange-scaling alpha.

.. psivar:: SAPT CT ENERGY

   SAPT charge-transfer energy.

.. psivar:: SAPT HF TOTAL ENERGY

   An induction-classified correction from HF implemented for SAPT0.
   Value varies by SAPT level.

.. psivar:: SAPT MP2 CORRELATION ENERGY

   An induction-classified correction from MP2 implemented for SAPT2.
   Value varies by SAPT level.

.. psivar:: SAPT0 DISP ENERGY
   SAPT0 ELST ENERGY
   SAPT0 EXCH ENERGY
   SAPT0 IND ENERGY
   SSAPT0 DISP ENERGY
   SSAPT0 ELST ENERGY
   SSAPT0 EXCH ENERGY
   SSAPT0 IND ENERGY
   SAPT2 DISP ENERGY
   SAPT2 ELST ENERGY
   SAPT2 EXCH ENERGY
   SAPT2 IND ENERGY
   SAPT2+ DISP ENERGY
   SAPT2+ ELST ENERGY
   SAPT2+ EXCH ENERGY
   SAPT2+ IND ENERGY
   SAPT2+(3) DISP ENERGY
   SAPT2+(3) ELST ENERGY
   SAPT2+(3) EXCH ENERGY
   SAPT2+(3) IND ENERGY
   SAPT2+3 DISP ENERGY
   SAPT2+3 ELST ENERGY
   SAPT2+3 EXCH ENERGY
   SAPT2+3 IND ENERGY

   Respectively, the dispersion, electrostatics, exchange, and induction
   components of the total electronic interaction energy [Eh] for the
   given SAPT level of theory. The sum of these four components yields
   the :samp:`{SAPT Level} TOTAL ENERGY`

.. psivar:: SAPT0 TOTAL ENERGY
   SSAPT0 TOTAL ENERGY
   SAPT2 TOTAL ENERGY
   SAPT2+ TOTAL ENERGY
   SAPT2+(3) TOTAL ENERGY
   SAPT2+3 TOTAL ENERGY

   The total electronic interaction energy [Eh] for the labeled SAPT level
   of theory.

.. psivar:: SAPT2+(CCD) DISP ENERGY
   SAPT2+(CCD) ELST ENERGY
   SAPT2+(CCD) EXCH ENERGY
   SAPT2+(CCD) IND ENERGY
   SAPT2+(3)(CCD) DISP ENERGY
   SAPT2+(3)(CCD) ELST ENERGY
   SAPT2+(3)(CCD) EXCH ENERGY
   SAPT2+(3)(CCD) IND ENERGY
   SAPT2+3(CCD) DISP ENERGY
   SAPT2+3(CCD) ELST ENERGY
   SAPT2+3(CCD) EXCH ENERGY
   SAPT2+3(CCD) IND ENERGY

   Respectively, the dispersion, electrostatics, exchange, and induction
   components of the total electronic interaction energy [Eh] for the
   given SAPT level of theory that incorporates coupled-cluster dispersion.
   The sum of these four components yields the :samp:`{SAPT Level} TOTAL ENERGY`

.. psivar:: SAPT2+(CCD) TOTAL ENERGY
   SAPT2+(3)(CCD) TOTAL ENERGY
   SAPT2+3(CCD) TOTAL ENERGY

   The total electronic interaction energy [Eh] for the labeled SAPT level
   of theory that incorporates coupled-cluster dispersion.

.. psivar:: SAPT2+DMP2 DISP ENERGY
   SAPT2+DMP2 ELST ENERGY
   SAPT2+DMP2 EXCH ENERGY
   SAPT2+DMP2 IND ENERGY
   SAPT2+(3)DMP2 DISP ENERGY
   SAPT2+(3)DMP2 ELST ENERGY
   SAPT2+(3)DMP2 EXCH ENERGY
   SAPT2+(3)DMP2 IND ENERGY
   SAPT2+3DMP2 DISP ENERGY
   SAPT2+3DMP2 ELST ENERGY
   SAPT2+3DMP2 EXCH ENERGY
   SAPT2+3DMP2 IND ENERGY
   SAPT2+(CCD)DMP2 DISP ENERGY
   SAPT2+(CCD)DMP2 ELST ENERGY
   SAPT2+(CCD)DMP2 EXCH ENERGY
   SAPT2+(CCD)DMP2 IND ENERGY
   SAPT2+(3)(CCD)DMP2 DISP ENERGY
   SAPT2+(3)(CCD)DMP2 ELST ENERGY
   SAPT2+(3)(CCD)DMP2 EXCH ENERGY
   SAPT2+(3)(CCD)DMP2 IND ENERGY
   SAPT2+3(CCD)DMP2 DISP ENERGY
   SAPT2+3(CCD)DMP2 ELST ENERGY
   SAPT2+3(CCD)DMP2 EXCH ENERGY
   SAPT2+3(CCD)DMP2 IND ENERGY

   Respectively, the dispersion, electrostatics, exchange, and induction
   components of the total electronic interaction energy [Eh] for the
   given SAPT level of theory that incorporates MP2 induction correction.
   The sum of these four components yields the :samp:`{SAPT Level} TOTAL ENERGY`

.. psivar:: SAPT2+DMP2 TOTAL ENERGY
   SAPT2+(3)DMP2 TOTAL ENERGY
   SAPT2+3DMP2 TOTAL ENERGY
   SAPT2+(CCD)DMP2 TOTAL ENERGY
   SAPT2+(3)(CCD)DMP2 TOTAL ENERGY
   SAPT2+3(CCD)DMP2 TOTAL ENERGY

   The total electronic interaction energy [Eh] for the labeled SAPT level
   of theory that incorporates MP2 induction correction.

.. psivar:: SCF ITERATIONS
   ADC ITERATIONS
   CCSD ITERATIONS
   OPTIMIZATION ITERATIONS

   Number of iterations [] in the named iterative method or optimization procedure.

.. psivar:: SCF DIPOLE

   Dipole array [e a0] for the SCF stage, (3,).

.. psivar:: SCF DIPOLE X
   SCF DIPOLE Y
   SCF DIPOLE Z

   The three components of the SCF dipole [Debye].
   Deprecated in favor of :psivar:`SCF DIPOLE`.

.. psivar:: SCF QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the SCF stage, (3, 3).

.. psivar:: SCF QUADRUPOLE XX
   SCF QUADRUPOLE XY
   SCF QUADRUPOLE XZ
   SCF QUADRUPOLE YY
   SCF QUADRUPOLE YZ
   SCF QUADRUPOLE ZZ

   The six components of the SCF quadrupole [Debye Ang].
   Deprecated in favor of :psivar:`SCF QUADRUPOLE`.

.. psivar:: SCF TOTAL ENERGY
   SCF TOTAL ENERGY (CHKPT)

   The total electronic energy [Eh] of the SCF stage of the calculation.
   The :samp:`{method} CORRELATION ENERGY` variables from subsequent stages of a
   calculation are often the corresponding :samp:`{method} TOTAL ENERGY`
   variables less this quantity. Constructed from Eq. :eq:`SCFterms`,
   where this quantity is :math:`E_{\text{SCF}}`.

   .. math::
      :nowrap:
      :label: SCFterms

         \begin{align*}
            E_{\text{SCF}} & = E_{NN} + E_{1e^-} + E_{2e^-} + E_{xc} + E_{\text{-D}} \\
                           & = E_{\text{FCTL/HF}} + E_{\text{-D}}
         \end{align*}

   Unless the method includes a dispersion correction, this quantity is
   equal to :psivar:`HF TOTAL ENERGY` (for HF) or
   :psivar:`DFT FUNCTIONAL TOTAL ENERGY` (for
   DFT). Unless the method is a DFT double-hybrid, this quantity is equal
   to :psivar:`DFT TOTAL ENERGY`.

.. psivar:: SCF TOTAL GRADIENT

   The total electronic gradient [E_h/a0] of the SCF stage of the calculation, ({nat}, 3).

.. psivar:: SCF DIPOLE GRADIENT

   The derivative of the SCF stage dipole [E_h a0/u] = [(e a0/a0)^2/u] with respect to nuclear perturbations
   as a degree-of-freedom by dipole component array, (3 * {nat}, 3).

.. psivar:: SCF TOTAL HESSIAN

   The total electronic second derivative [Eh/a0/a0] for the SCF stage, (3 * {nat}, 3 * {nat}).

.. psivar:: SCF STABILITY EIGENVALUES

   Array of eigenvalues from UHF or ROHF stability analysis.

.. psivar:: SCS-CCSD TOTAL ENERGY
   SCS-CCSD CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the CCSD-like method formed by reweighting :psivar:`CCSD DOUBLES ENERGY`
   by 1.27 opposite-spin and 1.13 same-spin contributions, with
   any singles carried along.

.. psivar:: SCS-MP2 TOTAL ENERGY
   SCS-MP2 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the MP2-like method formed by reweighting :psivar:`MP2 DOUBLES ENERGY`
   by 6/5 opposite-spin and 1/3 same-spin contributions, with
   any singles carried along.

.. psivar:: SCS-MP2-VDW TOTAL ENERGY
   SCS-MP2-VDW CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the MP2-like method formed by reweighting :psivar:`MP2 DOUBLES ENERGY`
   by 1.28 opposite-spin and 0.50 same-spin contributions, with
   any singles carried along. DOI: 10.1080/00268970802641242

.. psivar:: SCS(N)-MP2 TOTAL ENERGY
   SCS(N)-MP2 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the MP2-like method formed by reweighting :psivar:`MP2 DOUBLES ENERGY`
   by 0 opposite-spin and 1.76 same-spin contributions, with
   any singles carried along. doi: 10.1021/ct6002737

.. psivar:: SCS(N)-OMP2 CORRELATION ENERGY
   SCS(N)-OMP2 TOTAL ENERGY
   SCSN-OMP2 CORRELATION ENERGY
   SCSN-OMP2 TOTAL ENERGY

   Two spellings of a discontinued QCVariable that may still appear
   because the code is frozen pending an update.

.. psivar:: SCS-OMP2 TOTAL ENERGY
   SCS-OMP2 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the OMP2-like method formed by reweighting :psivar:`OMP2 DOUBLES ENERGY`
   by 6/5 opposite-spin and 1/3 same-spin contributions, with
   any singles carried along.

.. psivar:: SCS-MP3 TOTAL ENERGY
   SCS-MP3 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the MP3-like method formed by reweighting the difference between
   :psivar:`MP3 DOUBLES ENERGY` and :psivar:`MP2 DOUBLES ENERGY`
   by 0.25, atop the SCS-MP2 energy, with any singles carried along.

.. psivar:: SCS-OMP3 TOTAL ENERGY
   SCS-OMP3 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the OMP3-like method formed by reweighting the difference between
   :psivar:`OMP3 DOUBLES ENERGY` and :psivar:`OMP2 DOUBLES ENERGY`
   by 0.25, atop the SCS-OMP2 energy, with any singles carried along.

.. psivar:: SOS-MP2 TOTAL ENERGY
   SOS-MP2 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the MP2-like method formed by reweighting :psivar:`MP2 DOUBLES ENERGY`
   by 1.3 opposite-spin and 0 same-spin contributions, with
   any singles carried along.

.. psivar:: SOS-OMP2 TOTAL ENERGY
   SOS-OMP2 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the OMP2-like method formed by reweighting :psivar:`OMP2 DOUBLES ENERGY`
   by 1.2 opposite-spin and 0 same-spin contributions, with
   any singles carried along.

.. psivar:: SOS-OMP3 TOTAL ENERGY
   SOS-OMP3 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the OMP3-like method formed by reweighting the difference between
   :psivar:`OMP3 DOUBLES ENERGY` and :psivar:`OMP2 DOUBLES ENERGY`
   by 0.25, atop the SOS-OMP2
   energy using non-canonical weighting, with any singles carried along.

.. psivar:: SOS-PI-MP2 TOTAL ENERGY
   SOS-PI-MP2 CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the MP2-like method formed by reweighting :psivar:`MP2 DOUBLES ENERGY`
   by 1.4 opposite-spin and 0 same-spin contributions, with
   any singles carried along.

.. psivar:: TDDFT ROOT 0 -> ROOT m EXCITATION ENERGY - h SYMMETRY
   TD-fctl ROOT 0 -> ROOT m EXCITATION ENERGY - h SYMMETRY
   ADC ROOT 0 -> ROOT m EXCITATION ENERGY - h SYMMETRY
   EOM-CCSD ROOT 0 -> ROOT m EXCITATION ENERGY - h SYMMETRY

   The excitation energy of given method from ground state to root m
   in h symmetry (if available). DFT functional labeled if canonical.

.. psivar:: TDDFT ROOT n TOTAL ENERGY - h SYMMETRY
   TD-fctl ROOT n TOTAL ENERGY - h SYMMETRY
   ADC ROOT n TOTAL ENERGY - h SYMMETRY
   EOM-CCSD ROOT n TOTAL ENERGY - h SYMMETRY

   The total energy of given method from ground state to root m in h symmetry.

.. psivar:: ADC ROOT 0 -> ROOT m CORRELATION ENERGY - h SYMMETRY
   EOM-CCSD ROOT 0 -> ROOT m CORRELATION ENERGY - h SYMMETRY

   The correlation energy of given method from ground state reference energy to root m in h symmetry.

.. psivar:: TD-fctl ROOT 0 -> ROOT m OSCILLATOR STRENGTH (LEN) - h SYMMETRY
   TD-fctl ROOT 0 -> ROOT m OSCILLATOR STRENGTH (VEL) - h SYMMETRY

   The oscillator strength in length or velocity gauge of named method
   from ground state to root m in h symmetry (if available). DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 -> ROOT m ROTATORY STRENGTH (LEN) - h SYMMETRY
   TD-fctl ROOT 0 -> ROOT m ROTATORY STRENGTH (VEL) - h SYMMETRY

   The rotatory strength in length or velocity gauge of named method
   from ground state to root m in h symmetry (if available). DFT
   functional labeled if canonical.

.. psivar:: THERMAL ENERGY

   Total thermal energy E [Eh] at given temperature.

.. psivar:: THERMAL ENERGY CORRECTION

   Sum of electronic, translational, rotational, and vibrational corrections [Eh] to the thermal energy at given temperature.

.. psivar:: TWO-ELECTRON ENERGY

   The two-electron energy contribution [Eh] to the total SCF energy.
   Quantity :math:`E_{2e^-}` in Eq. :eq:`SCFterms`.

.. psivar:: UNCP-CORRECTED 2-BODY INTERACTION ENERGY

   The interaction energy [Eh] considering only two-body interactions,
   computed without counterpoise correction.
   Related variable :psivar:`CP-CORRECTED 2-BODY INTERACTION ENERGY`.

   .. math:: E_{\text{IE}} = E_{dimer} - \sum_{monomer}^{n}{E_{monomer}^{\text{unCP}}}

.. psivar:: WIBERG LOWDIN INDICES

   Property of Wiberg bond indices using orthogonal L\ |o_dots|\ wdin orbitals, (nat, nat).

.. psivar:: ZAPTn TOTAL ENERGY
   ZAPTn CORRELATION ENERGY

   The total electronic energy [Eh] and correlation energy component [Eh]
   for the labeled Z-averaged perturbation theory level.
   *n* is ZAPT perturbation order.

.. psivar:: ZERO K ENTHALPY

   Total electronic and zero-point energy [Eh] at 0 [K].

.. psivar:: ZPVE

   Vibrational zero-point energy [Eh] at 0 [K].

