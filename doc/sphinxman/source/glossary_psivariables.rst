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

.. psivar:: AAA (T) CORRECTION ENERGY
   AAB (T) CORRECTION ENERGY
   ABB (T) CORRECTION ENERGY
   BBB (T) CORRECTION ENERGY

   Components of the coupled-cluster perturbative triples correction [H].

.. psivar:: BRUECKNER CONVERGED

   Value 1 (0) when the Brueckner orbitals have (have not) converged.

.. psivar:: CBS TOTAL ENERGY
   CBS CORRELATION ENERGY
   CBS REFERENCE ENERGY

   The total electronic energy [H] and its breakdown into reference total
   energy [H] and correlation correction components [H] for the compound
   method requested through cbs().

.. psivar:: CC ROOT n TOTAL ENERGY

   The total electronic energy [H]
   for the requested coupled cluster level of theory and root 
   *n* (numbering starts at GS = 0).

.. psivar:: CC TOTAL ENERGY
   CC CORRELATION ENERGY

.. psivar:: CC2 TOTAL ENERGY
   CC2 CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the CC2 level of theory.

.. psivar:: CC3 TOTAL ENERGY
   CC3 CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the CC3 level of theory.

.. psivar:: CCSD TOTAL ENERGY
   CCSD CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the CCSD level of theory.

.. psivar:: CCSD(T) TOTAL ENERGY
   CCSD(T) CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the CCSD(T) level of theory.

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
   *n* (numbering starts at 1).

.. psivar:: CI STATE-AVERAGED TOTAL ENERGY
   CI STATE-AVERAGED CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for state-averaged CI/CASSCF levels of theory.
   
.. psivar:: CI TOTAL ENERGY
   CI CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the requested configuration interaction level of theory and root.

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

.. psivar:: DF-MP2 TOTAL ENERGY
   DF-MP2 CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the density-fitted MP2 level of theory.

.. psivar:: DFT FUNCTIONAL ENERGY

   The functional energy contribution [H] to the total SCF energy (DFT only).
   Quantity :math:`E_{xc}` in Eq. :eq:`SCFterms`.

.. psivar:: DFT FUNCTIONAL TOTAL ENERGY

   The total electronic energy [H] for the underlying functional
   of the requested DFT method, without any dispersion correction,
   the first four terms in Eq. :eq:`SCFterms`.
   When the requested method includes a dispersion correction, this 
   quantity is :math:`E_{\text{DFT}}` in Eq. :eq:`DFTDterms`.
   Otherwise, quantity equal to :psivar:`DFT TOTAL ENERGY <DFTTOTALENERGY>`
   and :psivar:`SCF TOTAL ENERGY <SCFTOTALENERGY>`.

.. psivar:: DFT TOTAL ENERGY

   The total electronic energy [H] for the requested DFT method, 
   :math:`E_{\text{SCF}}` in Eq. :eq:`SCFterms`.
   When the method includes a dispersion correction, this quantity
   is :math:`E_{\text{DFT-D}}` in Eq. :eq:`DFTDterms`.

.. psivar:: DISPERSION CORRECTION ENERGY

   The dispersion correction [H] appended to an underlying functional
   when a DFT-D method is requested. Quantity :math:`E_{\text{-D}}`
   in Eqs. :eq:`DFTDterms` and :eq:`SCFterms`.

   .. math:: E_{\text{DFT-D}} = E_{\text{DFT}} + E_{\text{-D}}
      :label: DFTDterms

.. psivar:: FCI TOTAL ENERGY
   FCI CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the full configuration interaction level of theory.

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

.. psivar:: MPn TOTAL ENERGY
   MPn CORRELATION ENERGY

   The total electronic energy [H] and correlation energy component [H]
   for the labeled M\ |o_dots|\ ller–Plesset perturbation theory level.
   *n* is MP perturbation order.

.. psivar:: NUCLEAR REPULSION ENERGY

   The nuclear repulsion energy contribution [H] to the total SCF energy.
   Quantity :math:`E_{NN}` in Eq. :eq:`SCFterms`.

.. psivar:: ONE-ELECTRON ENERGY

   The one-electron energy contribution [H] to the total SCF energy.
   Quantity :math:`E_{1e^-}` in Eq. :eq:`SCFterms`.

.. psivar:: SAPT DISP ENERGY
   SAPT ELST ENERGY
   SAPT EXCH ENERGY
   SAPT IND ENERGY

   Respectively, the dispersion, electrostatics, exchange, and induction
   components of the total electronic interaction energy [H] for the the
   requested SAPT level of theory. The sum of these four components yields
   :psivar:`SAPT ENERGY <SAPTENERGY>`.

.. psivar:: SAPT ENERGY

   The total electronic interaction energy [H] for the requested SAPT
   level of theory.

.. psivar:: SAPT SAPT0 ENERGY
   SAPT SAPT2 ENERGY
   SAPT SAPT2+ ENERGY
   SAPT SAPT2+(3) ENERGY
   SAPT SAPT2+3 ENERGY

   The total electronic interaction energy [H] for the labeled SAPT level
   of theory.

.. psivar:: SCF TOTAL ENERGY

   The total electronic energy [H] of the SCF stage of the calculation.
   The :psivar:`CORRELATION ENERGY` variables from subsequent stages of a
   calculation are often the corresponding :psivar:`TOTAL ENERGY` variables
   less this quantity. Constructed from Eq. :eq:`SCFterms`.

   .. math:: E_{\text{SCF}} = E_{NN} + E_{1e^-} + E_{2e^-} + E_{xc} + E_{\text{-D}}
      :label: SCFterms

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
   
