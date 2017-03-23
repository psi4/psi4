.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2017 The Psi4 Developers.
.. #
.. # The copyrights for code used from other parties are included in
.. # the corresponding files.
.. #
.. # This program is free software; you can redistribute it and/or modify
.. # it under the terms of the GNU General Public License as published by
.. # the Free Software Foundation; either version 2 of the License, or
.. # (at your option) any later version.
.. #
.. # This program is distributed in the hope that it will be useful,
.. # but WITHOUT ANY WARRANTY; without even the implied warranty of
.. # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
.. # GNU General Public License for more details.
.. #
.. # You should have received a copy of the GNU General Public License along
.. # with this program; if not, write to the Free Software Foundation, Inc.,
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

.. psivar:: (T) CORRECTION ENERGY

   The coupled-cluster perturbative triples correction [H].

.. psivar:: A-(T) CORRECTION ENERGY

   The coupled-cluster asymmetric perturbative triples correction [H].

.. psivar:: MP4(T) CORRECTION ENERGY

   The MP4 triples component [H]. Quantity is second right-hand term in
   Eq. :eq:`MP4terms`.

.. psivar:: AAA (T) CORRECTION ENERGY
   AAB (T) CORRECTION ENERGY
   ABB (T) CORRECTION ENERGY
   BBB (T) CORRECTION ENERGY

   Spin components of the UHF-based coupled-cluster perturbative triples correction [H].

.. psivar:: ACPF DIPOLE X
   ACPF DIPOLE Y
   ACPF DIPOLE Z

   The three components of the dipole [Debye] for the 
   averaged coupled-pair functional level of theory.

.. psivar:: ACPF QUADRUPOLE XX
   ACPF QUADRUPOLE XY
   ACPF QUADRUPOLE XZ
   ACPF QUADRUPOLE YY
   ACPF QUADRUPOLE YZ
   ACPF QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the 
   averaged coupled-pair functional level of theory.

.. psivar:: ACPF TOTAL ENERGY
   ACPF CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the averaged coupled-pair functional level of theory.

.. psivar:: AQCC DIPOLE X
   AQCC DIPOLE Y
   AQCC DIPOLE Z

   The three components of the dipole [Debye] for the 
   averaged quadratic coupled-cluster level of theory.

.. psivar:: AQCC QUADRUPOLE XX
   AQCC QUADRUPOLE XY
   AQCC QUADRUPOLE XZ
   AQCC QUADRUPOLE YY
   AQCC QUADRUPOLE YZ
   AQCC QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the 
   averaged quadratic coupled-cluster level of theory.

.. psivar:: AQCC TOTAL ENERGY
   AQCC CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the averaged quadratic coupled-cluster level of theory.
.. psivar:: BRUECKNER CONVERGED

   Value 1 (0) when the Brueckner orbitals have (have not) converged.

.. psivar:: CBS TOTAL ENERGY
   CBS CORRELATION ENERGY
   CBS REFERENCE ENERGY

   The total electronic energy [H] and its breakdown into reference total
   energy [H] and correlation correction components [H] for the compound
   method requested through cbs().

.. psivar:: CC ROOT n DIPOLE X
   CC ROOT n DIPOLE Y 
   CC ROOT n DIPOLE Z

   The three components of the dipole [Debye] for the requested
   coupled cluster level of theory and root *n* (number starts at GS = 0).

.. psivar:: CC ROOT n QUADRUPOLE XX
   CC ROOT n QUADRUPOLE XY
   CC ROOT n QUADRUPOLE XZ
   CC ROOT n QUADRUPOLE YY
   CC ROOT n QUADRUPOLE YZ
   CC ROOT n QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the requested
   coupled cluster level of theory and root *n* (numbering starts at GS = 0).

.. psivar:: CC ROOT n TOTAL ENERGY
   CC ROOT n CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the requested coupled cluster level of theory and root 
   *n* (numbering starts at GS = 0).

.. psivar:: CC TOTAL ENERGY
   CC CORRELATION ENERGY

.. psivar:: CC2 TOTAL ENERGY
   CC2 CORRELATION ENERGY
   CC3 TOTAL ENERGY
   CC3 CORRELATION ENERGY
   CC4 TOTAL ENERGY
   CC4 CORRELATION ENERGY
   CCnn TOTAL ENERGY
   CCnn CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the requested approximate coupled-cluster (CC2, CC3, up to CC\ *nn*)
   level of theory.

.. psivar:: CC DIPOLE X
   CC DIPOLE Y
   CC DIPOLE Z

   The three components of the dipole [Debye] for the requested
   coupled cluster level of theory and root.

.. psivar:: CC QUADRUPOLE XX
   CC QUADRUPOLE XY
   CC QUADRUPOLE XZ
   CC QUADRUPOLE YY
   CC QUADRUPOLE YZ
   CC QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the requested
   coupled cluster level of theory and root.

.. psivar:: CCSD TOTAL ENERGY
   CCSD CORRELATION ENERGY
   CCSDT TOTAL ENERGY
   CCSDT CORRELATION ENERGY
   CCSDTQ TOTAL ENERGY
   CCSDTQ CORRELATION ENERGY
   CCn TOTAL ENERGY
   CCn CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the requested full coupled-cluster (CCSD, CCSDT, up to CC\ *n*) 
   level of theory.

.. psivar:: CCSD(T) TOTAL ENERGY
   CCSD(T) CORRELATION ENERGY
   A-CCSD(T) TOTAL ENERGY
   A-CCSD(T) CORRELATION ENERGY
   CCSDT(Q) TOTAL ENERGY
   CCSDT(Q) CORRELATION ENERGY
   CC(n-1)(n) TOTAL ENERGY
   CC(n-1)(n) CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the perturbatively corrected coupled-cluster (CCSD(T), a-CCSD(T), CCSDT(Q), 
   up to CC(\ *n*\ -1)(\ *n*\ ) level of theory.

.. psivar:: CCSDT-1a TOTAL ENERGY
   CCSDT-1a CORRELATION ENERGY
   CCSDTQ-1a TOTAL ENERGY
   CCSDTQ-1a CORRELATION ENERGY
   CCn-1a TOTAL ENERGY
   CCn-1a CORRELATION ENERGY
   
   The total electronic energy [H] and correlation energy component [H]
   for the approximate coupled-cluster (CCSD(T)-1a, CCSDT(Q)-1a, 
   up to CC\ *n*\ -1a) level of theory.

.. psivar:: CCSDT-1b TOTAL ENERGY
   CCSDT-1b CORRELATION ENERGY
   CCSDTQ-1b TOTAL ENERGY
   CCSDTQ-1b CORRELATION ENERGY
   CCn-1b TOTAL ENERGY
   CCn-1b CORRELATION ENERGY
   
   The total electronic energy [H] and correlation energy component [H]
   for the approximate coupled-cluster (CCSD(T)-1b, CCSDT(Q)-1b, 
   up to CC\ *n*\ -1b) level of theory.

.. psivar:: CCSDT-3 TOTAL ENERGY
   CCSDT-3 CORRELATION ENERGY
   CCSDTQ-3 TOTAL ENERGY
   CCSDTQ-3 CORRELATION ENERGY
   CCn-3 TOTAL ENERGY
   CCn-3 CORRELATION ENERGY
   
   The total electronic energy [H] and correlation energy component [H]
   for the approximate coupled-cluster (CCSD(T)-3, CCSDT(Q)-3, 
   up to CC\ *n*\ -3) level of theory.

.. psivar:: CCSD(T)_L TOTAL ENERGY
   CCSD(T)_L CORRELATION ENERGY
   CCSDT(Q)_L TOTAL ENERGY
   CCSDT(Q)_L CORRELATION ENERGY
   CC(n-1)(n)_L TOTAL ENERGY
   CC(n-1)(n)_L CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the approximate coupled-cluster (CCSD(T)_L, CCSDT(Q)_L, 
   up to CC(\ *n*\ -1)(\ *n*\ )L level of theory.

.. psivar:: CEPA(0) DIPOLE X
   CEPA(0) DIPOLE Y
   CEPA(0) DIPOLE Z

   The three components of the dipole [Debye] for the 
   coupled electron pair approximation variant 0 level of theory.

.. psivar:: CEPA(0) QUADRUPOLE XX
   CEPA(0) QUADRUPOLE XY
   CEPA(0) QUADRUPOLE XZ
   CEPA(0) QUADRUPOLE YY
   CEPA(0) QUADRUPOLE YZ
   CEPA(0) QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the 
   coupled electron pair approximation variant 0 level of theory.

.. psivar:: CEPA(0) TOTAL ENERGY
   CEPA(0) CORRELATION ENERGY
   CEPA(1) TOTAL ENERGY
   CEPA(1) CORRELATION ENERGY
   CEPA(2) TOTAL ENERGY
   CEPA(2) CORRELATION ENERGY
   CEPA(3) TOTAL ENERGY
   CEPA(3) CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the requested variant of coupled electron pair approximation level of theory.

.. psivar:: CFOUR ERROR CODE

   The non-zero return value from a Cfour execution.

.. psivar:: CI DIPOLE X
   CI DIPOLE Y
   CI DIPOLE Z

   The three components of the dipole [Debye] for the requested
   configuration interaction level of theory and root.

.. psivar:: CI QUADRUPOLE XX
   CI QUADRUPOLE XY
   CI QUADRUPOLE XZ
   CI QUADRUPOLE YY
   CI QUADRUPOLE YZ
   CI QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the requested
   configuration interaction level of theory and root.

.. psivar:: CI ROOT n -> ROOT m DIPOLE X
   CI ROOT n -> ROOT m DIPOLE Y
   CI ROOT n -> ROOT m DIPOLE Z

   The three components of the transition dipole [Debye] between roots *n*
   and *m* for the requested configuration interaction level of theory.

.. psivar:: CI ROOT n -> ROOT m QUADRUPOLE XX
   CI ROOT n -> ROOT m QUADRUPOLE XY
   CI ROOT n -> ROOT m QUADRUPOLE XZ
   CI ROOT n -> ROOT m QUADRUPOLE YY
   CI ROOT n -> ROOT m QUADRUPOLE YZ
   CI ROOT n -> ROOT m QUADRUPOLE ZZ

   The three components of the transition quadrupole [Debye Ang] between
   roots *n* and *m* for the requested configuration interaction level of
   theory.

.. psivar:: CI ROOT n DIPOLE X
   CI ROOT n DIPOLE Y 
   CI ROOT n DIPOLE Z

   The three components of the dipole [Debye] for the requested
   configuration interaction level of theory and root *n*.

.. psivar:: CI ROOT n QUADRUPOLE XX
   CI ROOT n QUADRUPOLE XY
   CI ROOT n QUADRUPOLE XZ
   CI ROOT n QUADRUPOLE YY
   CI ROOT n QUADRUPOLE YZ
   CI ROOT n QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the requested
   configuration interaction level of theory and root *n*.

.. psivar:: CI ROOT n TOTAL ENERGY
   CI ROOT n CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the requested configuration interaction level of theory and root 
   *n* (numbering starts at 0).

.. psivar:: CI STATE-AVERAGED TOTAL ENERGY
   CI STATE-AVERAGED CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for state-averaged CI/CASSCF levels of theory.
   
.. psivar:: CI TOTAL ENERGY
   CI CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the requested configuration interaction level of theory and root.

.. psivar:: CISD DIPOLE X
   CISD DIPOLE Y
   CISD DIPOLE Z

   The three components of the dipole [Debye] for the 
   configuration interaction singles and doubles level of theory and root.

.. psivar:: CISD QUADRUPOLE XX
   CISD QUADRUPOLE XY
   CISD QUADRUPOLE XZ
   CISD QUADRUPOLE YY
   CISD QUADRUPOLE YZ
   CISD QUADRUPOLE ZZ

   The six components of the quadrupole [Debye Ang] for the 
   configuration interaction singles and doubles level of theory and root.

.. psivar:: CISD TOTAL ENERGY
   CISD CORRELATION ENERGY
   CISDT TOTAL ENERGY
   CISDT CORRELATION ENERGY
   CISDTQ CORRELATION ENERGY
   CISDTQ TOTAL ENERGY
   CIn CORRELATION ENERGY
   CIn TOTAL ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the labeled configuration interaction level of theory and root.
   *n* is CI order for *n* > 4.

.. psivar:: CP-CORRECTED 2-BODY INTERACTION ENERGY

   The interaction energy [H] considering only two-body interactions,
   computed with counterpoise correction.
   Related variable :psivar:`UNCP-CORRECTED 2-BODY INTERACTION ENERGY <UNCP-CORRECTED2-BODYINTERACTIONENERGY>`.

   .. math:: E_{\text{IE}} = E_{dimer} - \sum_{monomer}^{n}{E_{monomer}^{\text{CP}}}

.. psivar:: CURRENT CORRELATION ENERGY

   The correlation energy [H] corresponding to the :psivar:`CURRENT ENERGY <CURRENTENERGY>` variable.

.. psivar:: CURRENT ENERGY

   The total electronic energy [H] of the most recent stage of a
   calculation (frequently overwritten). This is the quantity tracked by
   the geometry optimizer.

.. psivar:: CURRENT REFERENCE ENERGY

   The total electronic energy [H] of the reference stage corresponding to
   the :psivar:`CURRENT ENERGY <CURRENTENERGY>` variable.

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

.. psivar:: db_name DATABASE ROOT-MEAN-SQUARE SIGNED DEVIATION

   The rms deviation [\ |kcalpermol|\ ] of the requested method *name*
   from the stored reference values for the requested reactions in
   database *db_name*. If no reference is available, this will be a large
   and nonsensical value.

   .. math:: \sqrt{\frac{1}{n}\sum_{rxn}^{n}{(\textsf{\textsl{name}}_{rxn}-\text{REF}_{rxn})^2}}

.. psivar:: DFT FUNCTIONAL TOTAL ENERGY

   The total electronic energy [H] for the underlying functional of the
   requested DFT method, without any dispersion correction; the first four
   terms in Eq. :eq:`SCFterms` or :eq:`DFTterms`. Quantity
   :math:`E_{\text{FCTL}}` in Eqs.  :eq:`SCFterms` and :eq:`DFTterms`.
   Unless the method includes a dispersion correction, this quantity is
   equal to :psivar:`SCF TOTAL ENERGY <SCFTOTALENERGY>`.

.. psivar:: DFT TOTAL ENERGY

   The total electronic energy [H] for the requested DFT method, 
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
   :psivar:`SCF TOTAL ENERGY <SCFTOTALENERGY>`. If the method is neither a
   double-hybrid, nor dispersion corrected, this quantity is equal to
   :psivar:`DFT FUNCTIONAL TOTAL ENERGY <DFTFUNCTIONALTOTALENERGY>`.

.. psivar:: DFT XC ENERGY

   The functional energy contribution [H] to the total SCF energy (DFT only).
   Quantity :math:`E_{xc}` in Eqs. :eq:`SCFterms` and :eq:`DFTterms`.

.. psivar:: DISPERSION CORRECTION ENERGY

   The dispersion correction [H] appended to an underlying functional
   when a DFT-D method is requested. Quantity :math:`E_{\text{-D}}`
   in Eqs. :eq:`SCFterms` and :eq:`DFTterms`.

.. psivar:: DOUBLE-HYBRID CORRECTION ENERGY

   The scaled MP2 correlation energy correction [H] appended to an 
   underlying functional when a DH-DFT method is requested.
   Quantity :math:`E_{\text{DH}}` in Eq. :eq:`DFTterms`.

.. psivar:: FCI TOTAL ENERGY
   FCI CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the full configuration interaction level of theory.

.. psivar:: HF TOTAL ENERGY

   The total electronic energy [H] for the Hartree--Fock method, without
   any dispersion correction; the first three (or four, since 
   :math:`E_{xc} = 0`) terms in Eq. :eq:`SCFterms`. Quantity :math:`E_{\text{HF}}`
   in Eq. :eq:`SCFterms`.
..   Unless the method includes a dispersion correction, this quantity is
   equal to :psivar:`SCF TOTAL ENERGY <SCFTOTALENERGY>`.

.. psivar:: LCC2 (+LMP2) TOTAL ENERGY

   The total electronic energy [H] for the local CC2 level of theory.

.. psivar:: LCCSD (+LMP2) TOTAL ENERGY

   The total electronic energy [H] for the local CCSD level of theory.

.. psivar:: MP2 TOTAL ENERGY
   MP2 CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the MP2 level of theory.

.. psivar:: MP2.5 TOTAL ENERGY
   MP2.5 CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the MP2.5 level of theory.

.. psivar:: MP3 TOTAL ENERGY
   MP3 CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the MP3 level of theory.

.. psivar:: MP4(SDQ) TOTAL ENERGY
   MP4(SDQ) CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the MP4 singles, doubles, quadruples level of theory.  Quantity
   :psivar:`MP4(SDQ) CORRELATION ENERGY <MP4(SDQ)CORRELATIONENERGY>` is
   first right-hand term in Eq. :eq:`MP4terms`.

.. psivar:: MP4 TOTAL ENERGY
   MP4 CORRELATION ENERGY
   MP4(SDTQ) TOTAL ENERGY
   MP4(SDTQ) CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the full MP4 level of theory. Quantity :psivar:`MP4 CORRELATION
   ENERGY <MP4CORRELATIONENERGY>` / :psivar:`MP4(SDTQ) CORRELATION ENERGY
   <MP4(SDTQ)CORRELATIONENERGY>` is left-hand term in Eq. :eq:`MP4terms`.

   .. math:: E_{\text{MP4}} = E_{\text{MP4(SDQ)}} + E_{\text{MP4(T)}}
      :label: MP4terms

.. psivar:: MPn TOTAL ENERGY
   MPn CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the labeled |MollerPlesset| perturbation theory level.
   *n* is MP perturbation order.

.. psivar:: NUCLEAR REPULSION ENERGY

   The nuclear repulsion energy contribution [H] to the total SCF energy.
   Quantity :math:`E_{NN}` in Eq. :eq:`SCFterms`.

   .. math:: E_{NN} = \sum_{i, j<i}^{N_{atom}}\frac{Z_i Z_j}{|\mathbf{R}_i - \mathbf{R}_j|}
      :label: ENN

.. psivar:: OCEPA(0) TOTAL ENERGY
   OCEPA(0) CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the orbital-optimized CEPA(0) level of theory.

.. psivar:: OMP2 TOTAL ENERGY
   OMP2 CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the orbital-optimized MP2 level of theory.

.. psivar:: OMP3 TOTAL ENERGY
   OMP3 CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the orbital-optimized MP3 level of theory.

.. psivar:: ONE-ELECTRON ENERGY

   The one-electron energy contribution [H] to the total SCF energy.
   Quantity :math:`E_{1e^-}` in Eq. :eq:`SCFterms`.

.. psivar:: QCISD TOTAL ENERGY
   QCISD CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the quadratic configuration interaction singles and doubles level
   of theory.

.. psivar:: QCISD(T) TOTAL ENERGY
   QCISD(T) CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the quadratic configuration interaction singles and doubles with
   perturbative triples correction level of theory.

.. psivar:: SAPT DISP ENERGY
   SAPT ELST ENERGY
   SAPT EXCH ENERGY
   SAPT IND ENERGY

   Respectively, the dispersion, electrostatics, exchange, and induction
   components of the total electronic interaction energy [H] for the the
   requested SAPT level of theory. The sum of these four components yields
   :psivar:`SAPT TOTAL ENERGY <SAPTTOTALENERGY>`.

.. psivar:: SAPT TOTAL ENERGY

   The total electronic interaction energy [H] for the requested SAPT
   level of theory.

.. psivar:: SAPT0 TOTAL ENERGY
   SSAPT0 TOTAL ENERGY
   SAPT2 TOTAL ENERGY
   SAPT2+ TOTAL ENERGY
   SAPT2+(3) TOTAL ENERGY
   SAPT2+3 TOTAL ENERGY

   The total electronic interaction energy [H] for the labeled SAPT level
   of theory.

.. psivar:: SAPT2+(CCD) TOTAL ENERGY
   SAPT2+(3)(CCD) TOTAL ENERGY
   SAPT2+3(CCD) TOTAL ENERGY

   The total electronic interaction energy [H] for the labeled SAPT level
   of theory that incorporates coupled-cluster dispersion.

.. psivar:: SAPT2+DMP2 TOTAL ENERGY
   SAPT2+(3)DMP2 TOTAL ENERGY
   SAPT2+3DMP2 TOTAL ENERGY
   SAPT2+(CCD)DMP2 TOTAL ENERGY
   SAPT2+(3)(CCD)DMP2 TOTAL ENERGY
   SAPT2+3(CCD)DMP2 TOTAL ENERGY

   The total electronic interaction energy [H] for the labeled SAPT level
   of theory that incorporates MP2 induction correction.

.. psivar:: SCF DIPOLE X
   SCF DIPOLE Y
   SCF DIPOLE Z

   The three components of the SCF dipole [Debye].

.. psivar:: SCF QUADRUPOLE XX
   SCF QUADRUPOLE XY
   SCF QUADRUPOLE XZ
   SCF QUADRUPOLE YY
   SCF QUADRUPOLE YZ
   SCF QUADRUPOLE ZZ

   The six components of the SCF quadrupole [Debye Ang].

.. psivar:: SCF TOTAL ENERGY

   The total electronic energy [H] of the SCF stage of the calculation.
   The :psivar:`CORRELATION ENERGY` variables from subsequent stages of a
   calculation are often the corresponding :psivar:`TOTAL ENERGY`
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
   equal to :psivar:`HF TOTAL ENERGY <HFTOTALENERGY>` (for HF) or
   :psivar:`DFT FUNCTIONAL TOTAL ENERGY <DFTFUNCTIONALTOTALENERGY>` (for
   DFT). Unless the method is a DFT double-hybrid, this quantity is equal
   to :psivar:`DFT TOTAL ENERGY <DFTTOTALENERGY>`.

.. psivar:: TWO-ELECTRON ENERGY

   The two-electron energy contribution [H] to the total SCF energy.
   Quantity :math:`E_{2e^-}` in Eq. :eq:`SCFterms`.

.. psivar:: UNCP-CORRECTED 2-BODY INTERACTION ENERGY

   The interaction energy [H] considering only two-body interactions,
   computed without counterpoise correction.
   Related variable :psivar:`CP-CORRECTED 2-BODY INTERACTION ENERGY <CP-CORRECTED2-BODYINTERACTIONENERGY>`.

   .. math:: E_{\text{IE}} = E_{dimer} - \sum_{monomer}^{n}{E_{monomer}^{\text{unCP}}}

.. psivar:: ZAPTn TOTAL ENERGY
   ZAPTn CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the labeled Z-averaged perturbation theory level.
   *n* is ZAPT perturbation order.
   
