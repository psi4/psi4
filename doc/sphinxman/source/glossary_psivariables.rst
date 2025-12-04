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

   The coupled-cluster bracket perturbative triples correction [E_h].

.. psivar:: (T) CORRECTION ENERGY

   The coupled-cluster perturbative triples correction [E_h].

.. psivar:: (AT) CORRECTION ENERGY
   A-(T) CORRECTION ENERGY

   The coupled-cluster asymmetric perturbative triples correction [E_h].

.. psivar:: AAA (T) CORRECTION ENERGY
   AAB (T) CORRECTION ENERGY
   ABB (T) CORRECTION ENERGY
   BBB (T) CORRECTION ENERGY

   Spin components of the UHF-based coupled-cluster perturbative triples correction [E_h].

.. psivar:: ACPF DIPOLE

   Dipole array [e a0] for the averaged coupled-pair functional level of theory, (3,).

.. psivar:: ACPF QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the averaged coupled-pair functional level of theory, (3, 3).

.. psivar:: ACPF TOTAL ENERGY
   ACPF CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the averaged coupled-pair functional level of theory.

.. psivar:: ADC ROOT 0 -> ROOT n EXCITATION ENERGY
   TD-fctl ROOT 0 -> ROOT n EXCITATION ENERGY

   The excitation energy [E_h] from ground state to root *n*.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 (IN h) -> ROOT n (IN i) EXCITATION ENERGY
   TD-fctl ROOT 0 (IN h) -> ROOT n (IN i) EXCITATION ENERGY

   The excitation energy [E_h] from the ground state (which is of irrep *h*)
   to root *n* within irrep *i*.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 (h) -> ROOT n (i) EXCITATION ENERGY
   TD-fctl ROOT 0 (h) -> ROOT n (i) EXCITATION ENERGY

   The excitation energy [E_h] from the ground state (which is of irrep *h*)
   to root *n* (which is of irrep *i*).
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 -> ROOT n EXCITATION ENERGY - h TRANSITION
   TD-fctl ROOT 0 -> ROOT n EXCITATION ENERGY - h TRANSITION

   The excitation energy [E_h] from the ground state to root *n*, and the
   transition is of irrep *h*.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 -> ROOT n ELECTRIC TRANSITION DIPOLE MOMENT (LEN)
   TD-fctl ROOT 0 -> ROOT n ELECTRIC TRANSITION DIPOLE MOMENT (LEN)

   The electric transition dipole moment [e a0] in length gauge, for the transition
   from the ground state to root *n*.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 (IN h) -> ROOT n (IN i) ELECTRIC TRANSITION DIPOLE MOMENT (LEN)
   TD-fctl ROOT 0 (IN h) -> ROOT n (IN i) ELECTRIC TRANSITION DIPOLE MOMENT (LEN)

   The electric transition dipole moment [e a0] in length gauge, for the transition
   from the ground state, which is of irrep *h*, to root *n* within irrep *i*.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 (h) -> ROOT n (i) ELECTRIC TRANSITION DIPOLE MOMENT (LEN)
   TD-fctl ROOT 0 (h) -> ROOT n (i) ELECTRIC TRANSITION DIPOLE MOMENT (LEN)

   The electric transition dipole moment [e a0] in length gauge, for the transition
   from the ground state, which is of irrep *h*, to root *n*, which is of irrep *i*.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 -> ROOT n ELECTRIC TRANSITION DIPOLE MOMENT (LEN) - h TRANSITION
   TD-fctl ROOT 0 -> ROOT n ELECTRIC TRANSITION DIPOLE MOMENT (LEN) - h TRANSITION

   The electric transition dipole moment [e a0] in length gauge, for the transition
   from the ground state to root *n*, and the transition is of *h* symmetry.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 -> ROOT n OSCILLATOR STRENGTH (LEN)
   CCname ROOT m -> ROOT n OSCILLATOR STRENGTH (LEN)
   TD-fctl ROOT 0 -> ROOT n OSCILLATOR STRENGTH (LEN)

   The length-gauge oscillator strength of the transition from root *m* to root *n*.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 (IN h) -> ROOT n (IN i) OSCILLATOR STRENGTH (LEN)
   CCname ROOT m (IN h) -> ROOT n (IN i) OSCILLATOR STRENGTH (LEN)
   TD-fctl ROOT 0 (IN h) -> ROOT n (IN i) OSCILLATOR STRENGTH (LEN)

   The length-gauge oscillator strength of the transition from root *m* within irrep *h*
   to root *n* within irrep *i*.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 (h) -> ROOT n (i) OSCILLATOR STRENGTH (LEN)
   CCname ROOT m (h) -> ROOT n (i) OSCILLATOR STRENGTH (LEN)
   TD-fctl ROOT 0 (h) -> ROOT n (i) OSCILLATOR STRENGTH (LEN)

   The length-gauge oscillator strength of the transition from root *m* to root *n*,
   which are in irreps *h* and *i*, respectively..
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 -> ROOT n OSCILLATOR STRENGTH (LEN) - h TRANSITION
   CCname ROOT m -> ROOT n OSCILLATOR STRENGTH (LEN) - h TRANSITION
   TD-fctl ROOT 0 -> ROOT n OSCILLATOR STRENGTH (LEN) - h TRANSITION

   The length-gauge oscillator strength of the transition from root *m* to root *n*,
   and the transition is of irrep *h*.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 -> ROOT n OSCILLATOR STRENGTH (VEL)
   TD-fctl ROOT 0 -> ROOT n OSCILLATOR STRENGTH (VEL)

   The velocity-gauge oscillator strength of the transition from the ground state to root *n*.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 (IN h) -> ROOT n (IN i) OSCILLATOR STRENGTH (VEL)
   TD-fctl ROOT 0 (IN h) -> ROOT n (IN i) OSCILLATOR STRENGTH (VEL)

   The velocity-gauge oscillator strength of the transition from the ground state within irrep *h*
   to root *n* within irrep *i*.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 (h) -> ROOT n (i) OSCILLATOR STRENGTH (VEL)
   TD-fctl ROOT 0 (h) -> ROOT n (i) OSCILLATOR STRENGTH (VEL)

   The velocity-gauge oscillator strength of the transition from the ground state to root *n*,
   which are in irreps *h* and *i*, respectively..
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 -> ROOT n OSCILLATOR STRENGTH (VEL) - h TRANSITION
   TD-fctl ROOT 0 -> ROOT n OSCILLATOR STRENGTH (VEL) - h TRANSITION

   The velocity-gauge oscillator strength of the transition from the ground state to root *n*,
   and the transition is of irrep *h*.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 -> ROOT n ROTATORY STRENGTH (VEL)
   CCname ROOT m -> ROOT n ROTATORY STRENGTH (VEL)
   TD-fctl ROOT 0 -> ROOT n ROTATORY STRENGTH (VEL)

   The velocity-gauge oscillator strength of the transition from root *m* to root *n*.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 (IN h) -> ROOT n (IN i) ROTATORY STRENGTH (VEL)
   CCname ROOT m (IN h) -> ROOT n (IN i) ROTATORY STRENGTH (VEL)
   TD-fctl ROOT 0 (IN h) -> ROOT n (IN i) ROTATORY STRENGTH (VEL)

   The velocity-gauge oscillator strength of the transition from root *m* within irrep *h*
   to root *n* within irrep *i*.
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 (h) -> ROOT n (i) ROTATORY STRENGTH (VEL)
   CCname ROOT m (h) -> ROOT n (i) ROTATORY STRENGTH (VEL)
   TD-fctl ROOT 0 (h) -> ROOT n (i) ROTATORY STRENGTH (VEL)

   The velocity-gauge oscillator strength of the transition from root *m* to root *n*,
   which are in irreps *h* and *i*, respectively..
   DFT functional labeled if canonical.

.. psivar:: ADC ROOT 0 -> ROOT n ROTATORY STRENGTH (VEL) - h TRANSITION
   CCname ROOT m -> ROOT n ROTATORY STRENGTH (VEL) - h TRANSITION
   TD-fctl ROOT 0 -> ROOT n ROTATORY STRENGTH (VEL) - h TRANSITION

   The velocity-gauge oscillator strength of the transition from root *m* to root *n*,
   and the transition is of irrep *h*.
   DFT functional labeled if canonical.

.. psivar:: AQCC DIPOLE

   Dipole array [e a0] for the averaged quadratic coupled-cluster level of theory, (3,).

.. psivar:: AQCC QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the averaged quadratic coupled-cluster level of theory, (3, 3).

.. psivar:: AQCC TOTAL ENERGY
   AQCC CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the averaged quadratic coupled-cluster level of theory.
.. psivar:: BRUECKNER CONVERGED

   Value 1 (0) when the Brueckner orbitals have (have not) converged.

.. psivar:: CBS NUMBER
   NBODY NUMBER
   FINDIF NUMBER

   Number of tasks [] the named procedure performs. These are immediate
   tasks, so if procedures are nested, the total number of tasks is
   the product.

.. psivar:: CBS TOTAL ENERGY
   CBS CORRELATION ENERGY
   CBS REFERENCE ENERGY

   The total electronic energy [E_h] and its breakdown into reference total
   energy [E_h] and correlation correction components [E_h] for the compound
   method requested through cbs().

.. psivar:: CCname ROOT n TOTAL ENERGY
   TD-fctl ROOT n TOTAL ENERGY

   The total electronic energy [E_h] for the requested theory and root *n* (*n* starts at 0).
   DFT functional labeled if canonical.

.. psivar:: CCname ROOT n (IN h) TOTAL ENERGY
   TD-fctl ROOT n (IN h) TOTAL ENERGY

   The total electronic energy [E_h] for the requested theory and root *n* within irrep *h* (*n* starts at 0).
   DFT functional labeled if canonical.

.. psivar:: CCname ROOT n (h) TOTAL ENERGY
   TD-fctl ROOT n (h) TOTAL ENERGY

   The total electronic energy [E_h] for the requested theory and root *n*, which is of irrep *h* (*n* starts at 0).
   DFT functional labeled if canonical.

.. psivar:: CCname ROOT n TOTAL ENERGY - h TRANSITION
   TD-fctl ROOT n TOTAL ENERGY - h TRANSITION

   The total electronic energy [E_h] for the requested theory and root *n*, and the transition is of irrep *h*, (*n* starts at 0).

.. psivar:: CCname ROOT n CORRELATION ENERGY

   The correlation energy [E_h] for the requested coupled cluster level of theory and root *n* (*n* starts at 0).
   DFT functional labeled if canonical.

.. psivar:: CCname ROOT n (IN h) CORRELATION ENERGY

   The correlation energy [E_h] for the requested coupled cluster level of theory and root *n* within irrep *h* (*n* starts at 0).

.. psivar:: CCname ROOT n (h) CORRELATION ENERGY

   The correlation energy [E_h] for the requested coupled cluster level of theory and root *n*, which is of irrep *h* (*n* starts at 0).

.. psivar:: CCname ROOT n CORRELATION ENERGY - h TRANSITION

   The correlation energy [E_h] for the requested coupled cluster level of theory and root *n*, and the transition is of irrep *h*, (*n* starts at 0).

.. psivar:: CCname ROOT n DIPOLE

   Dipole array [e a0] for the requested coupled cluster level of theory and root *n* (*n* starts at 0), (3,).

.. psivar:: CCname ROOT n (IN h) DIPOLE

   Dipole array [e a0] for the requested coupled cluster level of theory and root *n* within irrep *h* (*n* starts at 0), (3,).

.. psivar:: CCname ROOT n (h) DIPOLE

   Dipole array [e a0] for the requested coupled cluster level of theory and root *n*, which is of irrep *h* (*n* starts at 0), (3,).

.. psivar:: CCname ROOT n DIPOLE - h TRANSITION

   Dipole array [e a0] for the requested coupled cluster level of theory and root *n*, and the transition is of irrep *h*, (*n* starts at 0), (3,).

.. psivar:: CCname ROOT n QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the requested coupled cluster level of theory and root *n* (*n* starts at 0), (3,3).

.. psivar:: CCname ROOT n (IN h) QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the requested coupled cluster level of theory and root *n* within irrep *h* (*n* starts at 0), (3,3).

.. psivar:: CCname ROOT n (h) QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the requested coupled cluster level of theory and root *n*, which is of irrep *h* (*n* starts at 0), (3,3).

.. psivar:: CCname ROOT n QUADRUPOLE - h TRANSITION

   Redundant quadrupole array [e a0^2] for the requested coupled cluster level of theory and root *n*, and the transition is of irrep *h*, (*n* starts at 0), (3,3).

.. psivar:: CCname ROOT m -> ROOT n EINSTEIN A (LEN)

   The Einstein A coefficient, the spontaneous emission 'probability.'
   Units are in [1/s].
   Describes the transition between roots *m* and *n*.

.. psivar:: CCname ROOT m (IN h) -> ROOT n (IN i) EINSTEIN A (LEN)

   The Einstein A coefficient, the spontaneous emission 'probability.'
   Units are in [1/s].
   Describes the transition between root *m* within irrep *h* and root *n* which irrep *i*.

.. psivar:: CCname ROOT m (h) -> ROOT n (i) EINSTEIN A (LEN)

   The Einstein A coefficient, the spontaneous emission 'probability.'
   Units are in [1/s].
   Describes the transition between roots *m* and *n*, which are in irreps *h* and *i*, respectively..

.. psivar:: CCname ROOT m -> ROOT n EINSTEIN A (LEN) - h TRANSITION

   The Einstein A coefficient, the spontaneous emission 'probability.'
   Units are in [1/s].
   Describes the irrep *h* transition between roots *m* and *n*.

.. psivar:: CCname ROOT m -> ROOT n EINSTEIN B (LEN)

   The Einstein B coefficient, the stimulated emission 'probability'
   in terms of energy density. Units are in [m^3 / J / s^2].
   Describes the transition between roots *m* and *n*.

.. psivar:: CCname ROOT m (IN h) -> ROOT n (IN i) EINSTEIN B (LEN)

   The Einstein B coefficient, the stimulated emission 'probability'
   in terms of energy density. Units are in [m^3 / J / s^2].
   Describes the transition between root *m* within irrep *h* and root *n* which irrep *i*.

.. psivar:: CCname ROOT m (h) -> ROOT n (i) EINSTEIN B (LEN)

   The Einstein B coefficient, the stimulated emission 'probability'
   in terms of energy density. Units are in [m^3 / J / s^2].
   Describes the transition between roots *m* and *n*, which are in irreps *h* and *i*, respectively..

.. psivar:: CCname ROOT m -> ROOT n EINSTEIN B (LEN) - h TRANSITION

   The Einstein B coefficient, the stimulated emission 'probability'
   in terms of energy density. Units are in [m^3 / J / s^2].
   Describes the irrep *h* transition between roots *m* and *n*.

.. psivar:: CCname ROOT m -> ROOT n ROTATORY STRENGTH (LEN)
   TD-fctl ROOT 0 -> ROOT n ROTATORY STRENGTH (LEN)

   The length-gauge rotatory strength of the transition from root *m* to root *n*.
   DFT functional labeled if canonical.

.. psivar:: CCname ROOT m (IN h) -> ROOT n (IN i) ROTATORY STRENGTH (LEN)
   TD-fctl ROOT 0 (IN h) -> ROOT n (IN i) ROTATORY STRENGTH (LEN)

   The length-gauge oscillator strength of the transition from root *m* within irrep *h*
   to root *n* within irrep *i*.
   DFT functional labeled if canonical.

.. psivar:: CCname ROOT m (h) -> ROOT n (i) ROTATORY STRENGTH (LEN)
   TD-fctl ROOT 0 (h) -> ROOT n (i) ROTATORY STRENGTH (LEN)

   The length-gauge oscillator strength of the transition from root *m* to root *n*,
   which are in irreps *h* and *i*, respectively..
   DFT functional labeled if canonical.

.. psivar:: CCname ROOT m -> ROOT n ROTATORY STRENGTH (LEN) - h TRANSITION
   TD-fctl ROOT 0 -> ROOT n ROTATORY STRENGTH (LEN) - h TRANSITION

   The length-gauge oscillator strength of the transition from root *m* to root *n*,
   and the transition is of irrep *h*.
   DFT functional labeled if canonical.

.. psivar:: CC TOTAL ENERGY
   CC CORRELATION ENERGY

.. psivar:: CC CORRELATION KINETIC ENERGY

   The correlation correction to the kinetic energy [E_h], as computed by a coupled cluster method.

.. psivar:: CC CORRELATION POTENTIAL ENERGY

   The correlation correction to the potential energy [E_h], as computed by a coupled cluster method.

.. psivar:: CC CORRELATION VIRIAL RATIO

   The correlation virial ratio, as defined in https://doi/org/10.1063/1.1535440 for basis set completeness analysis. Computed using coupled cluster.

.. psivar:: CC VIRIAL RATIO

   The virial ratio, as computed by a coupled cluster method. Only defined for a fully quantum mechanical computation, i.e., not QM/MM or solvated.

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

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the requested approximate coupled-cluster (CC2, CC3, up to CC\ *nn*)
   level of theory.

.. psivar:: CC DIPOLE

   Dipole array [e a0] for the requested coupled cluster level of theory and root, (3,).

.. psivar:: CC2 DIPOLE POLARIZABILITY @ xNM
   CCSD DIPOLE POLARIZABILITY @ xNM

   The dipole polarizability in atomic units [(e^2 a0^2)/E_h] calculated at
   the CC level for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CC2 DIPOLE POLARIZABILITY TENSOR @ xNM
   CCSD DIPOLE POLARIZABILITY TENSOR @ xNM

   The dipole polarizability tensor in atomic units [(e^2 a0^2)/E_h] calculated at
   the CC level for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CC2 QUADRUPOLE POLARIZABILITY @ xNM
   CCSD QUADRUPOLE POLARIZABILITY @ xNM

   The quadrupole polarizability in atomic units [(e^2 a0^3)/E_h] calculated at
   the CC level for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CC2 QUADRUPOLE POLARIZABILITY TENSOR @ xNM
   CCSD QUADRUPOLE POLARIZABILITY TENSOR @ xNM

   The quadrupole polarizability in atomic units [(e^2 a0^3)/E_h] calculated at
   the CC level for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CC2 SPECIFIC ROTATION (LEN) @ xNM
   CCSD SPECIFIC ROTATION (LEN) @ xNM

   The specific rotation [deg/(dm (g/cm^3))] calculated at the CC level in the
   length gauge for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CC2 SPECIFIC ROTATION (VEL) @ xNM
   CCSD SPECIFIC ROTATION (VEL) @ xNM

   The specific rotation [deg/(dm (g/cm^3))] calculated at the CC level in the
   velocity gauge for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CC2 SPECIFIC ROTATION (MVG) @ xNM
   CCSD SPECIFIC ROTATION (MVG) @ xNM

   The specific rotation [deg/(dm (g/cm^3))] calculated at the CC level in the
   modified velocity gauge for a given (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CC2 ROTATION (LEN) ORIGIN-DEPENDENCE @ xNM
   CCSD ROTATION (LEN) ORIGIN-DEPENDENCE @ xNM

   The origin-dependence of the CC specific rotation in deg/[dm (g/cm^3)]/bohr and the
   length gauge, computed at (x) wavelength, (x) rounded to nearest integer.

.. psivar:: CCD TOTAL ENERGY
   CCD CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the coupled-cluster doubles level of theory.

.. psivar:: CC ALPHA-ALPHA PAIR ENERGIES
   CCSD ALPHA-ALPHA PAIR ENERGIES
   CC2 ALPHA-ALPHA PAIR ENERGIES
   CC3 ALPHA-ALPHA PAIR ENERGIES
   MP2 ALPHA-ALPHA PAIR ENERGIES

   Restricted-reference same-spin pair energies for coupled-cluster theories.
   Size number of active doubly occupied orbitals, square.

.. psivar:: CC ALPHA-BETA PAIR ENERGIES
   CCSD ALPHA-BETA PAIR ENERGIES
   CC2 ALPHA-BETA PAIR ENERGIES
   CC3 ALPHA-BETA PAIR ENERGIES
   MP2 ALPHA-BETA PAIR ENERGIES

   Restricted-reference opposite-spin (alpha first) pair energies for coupled-cluster
   theories. Size number of active doubly occupied orbitals, square.

.. psivar:: CC SINGLET PAIR ENERGIES
   CCSD SINGLET PAIR ENERGIES
   CC2 SINGLET PAIR ENERGIES
   CC3 SINGLET PAIR ENERGIES
   MP2 SINGLET PAIR ENERGIES

   Restricted-reference singlet-adapted pair energies for coupled-cluster theories.
   Size number of active doubly occupied orbitals, square.

.. psivar:: CC TRIPLET PAIR ENERGIES
   CCSD TRIPLET PAIR ENERGIES
   CC2 TRIPLET PAIR ENERGIES
   CC3 TRIPLET PAIR ENERGIES
   MP2 TRIPLET PAIR ENERGIES

   Restricted-reference triplet-adapted pair energies for coupled-cluster theories.
   Size number of active doubly occupied orbitals, square.

.. psivar:: CCSD TOTAL ENERGY
   CCSD CORRELATION ENERGY
   CCSDT TOTAL ENERGY
   CCSDT CORRELATION ENERGY
   CCSDTQ TOTAL ENERGY
   CCSDTQ CORRELATION ENERGY
   CCn TOTAL ENERGY
   CCn CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the requested full coupled-cluster (CCSD, CCSDT, up to CC\ *n*)
   level of theory.

.. psivar:: CCSD(T) TOTAL ENERGY
   CCSD(T) CORRELATION ENERGY
   CCSD(AT) TOTAL ENERGY
   CCSD(AT) CORRELATION ENERGY
   A-CCSD(T) TOTAL ENERGY
   A-CCSD(T) CORRELATION ENERGY
   CCSDT(Q) TOTAL ENERGY
   CCSDT(Q) CORRELATION ENERGY
   CC(n-1)(n) TOTAL ENERGY
   CC(n-1)(n) CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the perturbatively corrected coupled-cluster (CCSD(T), A-CCSD(T) = CCSD(AT), CCSDT(Q),
   up to CC(\ *n*\ -1)(\ *n*\ ) level of theory.

.. psivar:: CCSDT-1a TOTAL ENERGY
   CCSDT-1a CORRELATION ENERGY
   CCSDTQ-1a TOTAL ENERGY
   CCSDTQ-1a CORRELATION ENERGY
   CCn-1a TOTAL ENERGY
   CCn-1a CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the approximate coupled-cluster (CCSD(T)-1a, CCSDT(Q)-1a,
   up to CC\ *n*\ -1a) level of theory.

.. psivar:: CCSDT-1b TOTAL ENERGY
   CCSDT-1b CORRELATION ENERGY
   CCSDTQ-1b TOTAL ENERGY
   CCSDTQ-1b CORRELATION ENERGY
   CCn-1b TOTAL ENERGY
   CCn-1b CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the approximate coupled-cluster (CCSD(T)-1b, CCSDT(Q)-1b,
   up to CC\ *n*\ -1b) level of theory.

.. psivar:: CCSDT-3 TOTAL ENERGY
   CCSDT-3 CORRELATION ENERGY
   CCSDTQ-3 TOTAL ENERGY
   CCSDTQ-3 CORRELATION ENERGY
   CCn-3 TOTAL ENERGY
   CCn-3 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the approximate coupled-cluster (CCSD(T)-3, CCSDT(Q)-3,
   up to CC\ *n*\ -3) level of theory.

.. psivar:: CCSD(T)_L TOTAL ENERGY
   CCSD(T)_L CORRELATION ENERGY
   CCSDT(Q)_L TOTAL ENERGY
   CCSDT(Q)_L CORRELATION ENERGY
   CC(n-1)(n)_L TOTAL ENERGY
   CC(n-1)(n)_L CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the approximate coupled-cluster (CCSD(T)_L, CCSDT(Q)_L,
   up to CC(\ *n*\ -1)(\ *n*\ )L level of theory.

.. psivar:: CCSDT(Q)/A TOTAL ENERGY
   CCSDT(Q)/A CORRELATION ENERGY
   CCSDT(Q)/B TOTAL ENERGY
   CCSDT(Q)/B CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the modified CCSDT(Q) level of theory.

.. psivar:: CEPA(0) DIPOLE

   Dipole array [e a0] for the coupled electron pair approximation variant 0 level of theory, (3,).

.. psivar:: CEPA(0) QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the coupled electron pair approximation variant 0 level of theory, (3, 3).

.. psivar:: CEPA(0) TOTAL ENERGY
   CEPA(0) CORRELATION ENERGY
   CEPA(1) TOTAL ENERGY
   CEPA(1) CORRELATION ENERGY
   CEPA(2) TOTAL ENERGY
   CEPA(2) CORRELATION ENERGY
   CEPA(3) TOTAL ENERGY
   CEPA(3) CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the requested variant of coupled electron pair approximation level of theory.

.. psivar:: CFOUR ERROR CODE

   The non-zero return value from a Cfour execution.

.. psivar:: CI DIPOLE

   Dipole array [e a0] for the requested configuration interaction level of theory, (3,).

.. psivar:: CI QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the requested configuration interaction level of theory, (3, 3).

.. psivar:: CI ROOT n -> ROOT m DIPOLE

   Transition dipole array [e a0] between roots *n* and *m* for the requested configuration interaction level of theory, (3,).

.. psivar:: CI ROOT n -> ROOT m QUADRUPOLE

   Redundant transition quadrupole array [e a0^2] between roots *n* and *m* for the requested configuration interaction level of theory, (3, 3).

.. psivar:: CI ROOT n DIPOLE

   Dipole array [e a0] for the requested configuration interaction level of theory and root *n*, (3,).

.. psivar:: CI ROOT n QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the requested configuration interaction level of theory and root *n*, (3, 3).

.. psivar:: CI ROOT n TOTAL ENERGY
   CI ROOT n CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the requested configuration interaction level of theory and root
   *n* (numbering starts at 0).

.. psivar:: CI STATE-AVERAGED TOTAL ENERGY
   CI STATE-AVERAGED CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for state-averaged CI/CASSCF levels of theory.

.. psivar:: CI TOTAL ENERGY
   CI CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the requested configuration interaction level of theory and root.

.. psivar:: CISD DIPOLE

   Dipole array [e a0] for the configuration interaction singles and doubles level of theory, (3,).

.. psivar:: CISD QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the configuration interaction singles and doubles level of theory, (3, 3).

.. psivar:: CISD TOTAL ENERGY
   CISD CORRELATION ENERGY
   CISDT TOTAL ENERGY
   CISDT CORRELATION ENERGY
   CISDTQ CORRELATION ENERGY
   CISDTQ TOTAL ENERGY
   CIn CORRELATION ENERGY
   CIn TOTAL ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the labeled configuration interaction level of theory and root.
   *n* is CI order for *n* > 4.

.. psivar:: CP-CORRECTED 2-BODY INTERACTION ENERGY

   The interaction energy [E_h] considering only two-body interactions,
   computed with counterpoise correction.
   Related variable :psivar:`UNCP-CORRECTED 2-BODY INTERACTION ENERGY`.

   .. math:: E_{\text{IE}} = E_{dimer} - \sum_{monomer}^{n}{E_{monomer}^{\text{CP}}}

.. psivar:: CURRENT CORRELATION ENERGY

   The correlation energy [E_h] corresponding to the :psivar:`CURRENT ENERGY` variable.

.. psivar:: CURRENT ENERGY

   The total electronic energy [E_h] of the most recent stage of a
   calculation (frequently overwritten). This is the quantity tracked by
   the geometry optimizer.

.. psivar:: CURRENT REFERENCE ENERGY

   The total electronic energy [E_h] of the reference stage corresponding to
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
   The total electronic energy [E_h] and correlation energy component [E_h]
   for the MP2-like method formed by any reweighting of :psivar:`MP2 DOUBLES ENERGY`
   for opposite-spin and same-spin contributions, with
   any singles carried along.
   Depending on weights, may equal any of MP2, SCS-MP2, SCS(N)-MP2, etc. quantities.
   Contrast with :psivar:`SCS-MP2 TOTAL ENERGY`.

.. psivar:: CUSTOM SCS-MP2.5 TOTAL ENERGY
   CUSTOM SCS-MP2.5 CORRELATION ENERGY
   CUSTOM SCS-MP3 TOTAL ENERGY
   CUSTOM SCS-MP3 CORRELATION ENERGY
   CUSTOM SCS-REMP2 TOTAL ENERGY
   CUSTOM SCS-REMP2 CORRELATION ENERGY
   CUSTOM SCS-LCCD TOTAL ENERGY
   CUSTOM SCS-LCCD CORRELATION ENERGY
   CUSTOM SCS-OMP2 TOTAL ENERGY
   CUSTOM SCS-OMP2 CORRELATION ENERGY
   CUSTOM SCS-OMP2.5 TOTAL ENERGY
   CUSTOM SCS-OMP2.5 CORRELATION ENERGY
   CUSTOM SCS-OMP3 TOTAL ENERGY
   CUSTOM SCS-OMP3 CORRELATION ENERGY
   CUSTOM SCS-OREMP2 TOTAL ENERGY
   CUSTOM SCS-OREMP2 CORRELATION ENERGY
   CUSTOM SCS-OLCCD TOTAL ENERGY
   CUSTOM SCS-OLCCD CORRELATION ENERGY

   Changeable quantities based on options.
   The total electronic energy [E_h] and correlation energy component [E_h]
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

   An energy term in density cumulant theory [E_h]. This term is the
   2-electron cumulant's contribution contribution to the reduced
   density matrix energy expression. Not recommended for interpretative
   use except by reduced density matrix specialists.

.. psivar:: DCT SCF ENERGY

   An energy term in density cumulant theory [E_h]. This term is the
   1-electron reduced density matrix (1RDM) contribution to the reduced
   density matrix energy expression, plus the contribution of the
   antisymmetrized product of 1RDMs. Not recommended for interpretative
   use except by reduced density matrix specialists.

.. psivar:: DCT THREE-PARTICLE ENERGY

   The three-particle correlation energy correction [E_h] in density cumulant
   theory, akin to :psivar:`(T) CORRECTION ENERGY` in coupled-cluster.

.. psivar:: DCT TOTAL ENERGY

   Total energy [E_h] in density cumulant theory. Sum of :psivar:`DCT SCF ENERGY`,
   :psivar:`DCT LAMBDA ENERGY`, and :psivar:`DCT THREE-PARTICLE ENERGY` when present.

.. psivar:: DETCI AVG DVEC NORM

   A measure of configuration interaction convergence.

.. psivar:: DFT FUNCTIONAL TOTAL ENERGY

   The total electronic energy [E_h] for the underlying functional of the
   requested DFT method, without any dispersion correction; the first four
   terms in Eq. :eq:`SCFterms` or :eq:`DFTterms`. Quantity
   :math:`E_{\text{FCTL}}` in Eqs.  :eq:`SCFterms` and :eq:`DFTterms`.
   Unless the method includes a dispersion correction, this quantity is
   equal to :psivar:`SCF TOTAL ENERGY`.

.. psivar:: DFT TOTAL ENERGY

   The total electronic energy [E_h] for the requested DFT method,
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

   The total electronic second derivative [E_h/a0/a0] for the requested DFT method, (3 * {nat}, 3 * {nat}).

.. psivar:: DFT XC ENERGY

   The functional energy contribution [E_h] to the total SCF energy (DFT only).
   Quantity :math:`E_{xc}` in Eqs. :eq:`SCFterms` and :eq:`DFTterms`.

.. psivar:: DFT VV10 ENERGY

   The VV10 nonlocal contribution [E_h] to the total SCF energy (DFT only).
   Included in :psivar:`DFT FUNCTIONAL TOTAL ENERGY`.

.. psivar:: DISPERSION CORRECTION ENERGY
   fctl DISPERSION CORRECTION ENERGY

   The dispersion correction [E_h] appended to an underlying functional
   when a DFT-D method is requested. Quantity :math:`E_{\text{-D}}`
   in Eqs. :eq:`SCFterms` and :eq:`DFTterms`.
   When dispersion parameters are untweaked for a functional and dispersion
   level, labeled QCVariable also defined.

.. psivar:: DLPNO LMP2 WEAK PAIR ENERGY
   DLPNO SC-LMP2 PAIR ENERGY
   DLPNO DIPOLE ENERGY
   DLPNO PNO TRUNCATION ERROR

   Various corrections in the overall DLPNO-CCSD correlation energy

.. psivar:: DLPNO SEMICANONICAL (T0) ENERGY
   DLPNO SCREENED TRIPLETS ENERGY

   Various components to the overall DLPNO-(T) correlation energy

.. psivar:: DOUBLE-HYBRID CORRECTION ENERGY

   The scaled MP2 correlation energy correction [E_h] appended to an
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

   The total DMRG total electonic energy [E_h]. Not unique because oribital spaces vary.

.. psivar:: DMRG-CASPT2 TOTAL ENERGY

   The total DMRG plus CASPT2 total electonic energy [E_h] . Not unique because orbital spaces vary.

.. psivar:: EFP DISP ENERGY
   EFP ELST ENERGY
   EFP EXCH ENERGY
   EFP IND ENERGY

   Respectively, the dispersion, electrostatics, exchange, and induction
   components of the total electronic interaction energy [E_h] for EFP/EFP
   computations. The sum of these four components yields
   :psivar:`EFP TOTAL ENERGY`.

.. psivar:: EFP TOTAL ENERGY

   The total electronic interaction energy [E_h] for EFP/EFP computations.

.. psivar:: EFP TORQUE

   The torque, not gradient for EFP/EFP computations.

.. psivar:: ENTHALPY

   Total enthalpy H [E_h] at given temperature.

.. psivar:: ENTHALPY CORRECTION

   Sum of electronic, translational, rotational, and vibrational corrections [E_h] to the enthalpy at given temperature.

.. psivar:: ESP AT CENTER n

   Property of electrostatic potential [E_h / e] at location, usually atom center, n.

.. psivar:: FCI TOTAL ENERGY
   FCI CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the full configuration interaction level of theory.

.. psivar:: GIBBS FREE ENERGY

   Total Gibbs free energy [E_h], free enthalpy at given temperature.

.. psivar:: GIBBS FREE ENERGY CORRECTION

   Sum of electronic, translational, rotational, and vibrational corrections [E_h] to the free enthalpy at given temperature.

.. psivar:: GRID ELECTRONS TOTAL
   GRID ELECTRONS ALPHA
   GRID ELECTRONS BETA

   The number of electrons integrated by the xc quadrature grid.

.. psivar:: HF TOTAL ENERGY

   The total electronic energy [E_h] for the Hartree--Fock method, without
   any dispersion correction; the first three (or four, since
   :math:`E_{xc} = 0`) terms in Eq. :eq:`SCFterms`. Quantity :math:`E_{\text{HF}}`
   in Eq. :eq:`SCFterms`.

.. psivar:: HF KINETIC ENERGY

   The total kinetic energy [E_h] of the Hartree--Fock method.

.. psivar:: HF POTENTIAL ENERGY

   The total potential energy [E_h] of the Hartree--Fock method.

.. psivar:: HF VIRIAL RATIO

   The virial ratio of the Hartree--Fock method. Only defined for a fully quantum mechanical computation, i.e., not QM/MM.

.. psivar:: HF TOTAL GRADIENT

   The total electronic gradient [E_h/a0] of the Hartree--Fock method, ({nat}, 3).

.. psivar:: HF DIPOLE GRADIENT

   The derivative of the Hartree--Fock method dipole [E_h a0/u] = [(e a0/a0)^2/u] with respect to nuclear perturbations
   as a degree-of-freedom by dipole component array, (3 * {nat}, 3).

.. psivar:: HF TOTAL HESSIAN

   The total electronic second derivative [E_h/a0/a0] for the Hartree-Fock method, (3 * {nat}, 3 * {nat}).

.. psivar:: HF-CABS TOTAL ENERGY
   F12 CABS CORRECTION ENERGY

   The correction and total electronic energy [E_h] for the Hartree--Fock method
   including the complementary auxiliary basis set (CABS) correction. Defined for
   explicity correlated methods.

.. psivar:: LCCD TOTAL ENERGY
   LCCD CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the linearized coupled cluster doubles level of theory.

.. psivar:: LCCSD TOTAL ENERGY
   LCCSD CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the linearized coupled cluster singles and doubles level of theory.

.. psivar:: LCC2 (+LMP2) TOTAL ENERGY

   The total electronic energy [E_h] for the local CC2 level of theory.

.. psivar:: LCCSD (+LMP2) TOTAL ENERGY

   The total electronic energy [E_h] for the local CCSD level of theory.

.. psivar:: LEFT-RIGHT CC2 EIGENVECTOR OVERLAP
   LEFT-RIGHT CC3 EIGENVECTOR OVERLAP
   LEFT-RIGHT CCSD EIGENVECTOR OVERLAP
   LEFT-RIGHT CCSD(T) EIGENVECTOR OVERLAP

   The overlap between the right-hand coupled coupled cluster eigenvector and the
   left-hand eigenvector from the coupled cluster lambda (response) equations.

.. psivar:: LOWDIN CHARGES

   Property of partial atomic charges [e] by the method of L\ |o_dots|\ wdin, (nat,).

.. psivar:: LOWDIN SPINS

   Property of partial atomic spin population by the method of L\ |o_dots|\ wdin, (nat,).
   Consider this the fractional number of unpaired electrons.

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

.. psivar:: MBIS VALENCE CHARGES

   Per-atom valence charges [e] computed from an MBIS partitioned density.

.. psivar:: MBIS VALENCE WIDTHS

   Per-atom density width [a0] of the associated valence charge computed
   from an MBIS partitioned density. Equivalent to the inverse of the
   linear decay rate of the atomic density.

.. psivar:: MBIS VOLUME RATIOS

   Per-atom ratio between the atomic volume (<R^3>) and the free-atomic
   volume, unitless.

.. psivar:: MCSCF TOTAL ENERGY

   Multiconfigurational self-consistent-field energy [E_h] in the course of
   a configuration interaction computation. May be single-root or state-averaged.

.. psivar:: mtd DIPOLE

   Dipole array [e a0] for the named method, (3,).

.. psivar:: mtd QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the named method, (3, 3).

.. psivar:: mtd OCTUPOLE

   Redundant octupole array [e a0^3] for the named method, (3, 3, 3).

.. psivar:: mtd HEXADECAPOLE

   Redundant hexadecapole array [e a0^4] for the named method, (3, 3, 3, 3).

.. psivar:: mtd 32-POLE

   Redundant 32-pole array [e a0^5] for the named method, (3, 3, 3, 3, 3).

.. psivar:: mtd 64-POLE

   Redundant 64-pole array [e a0^6] for the named method, (3, 3, 3, 3, 3, 3).

.. psivar:: mtd 128-POLE

   Redundant 128-pole array [e a0^7] for the named method, (3, 3, 3, 3, 3, 3, 3).

.. psivar:: MP2 TOTAL ENERGY
   MP2 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the MP2 level of theory.

.. psivar:: MP2 TOTAL GRADIENT
   The total electronic gradient [E_h/a0] of the MP2 level of theory, ({nat}, 3).

.. psivar:: MP2 DIPOLE GRADIENT

   The derivative of the MP2 level of theory dipole [E_h a0/u] = [(e a0/a0)^2/u] with respect to nuclear perturbations
   as a degree-of-freedom by dipole component array, (3 * {nat}, 3).

.. psivar:: MP2 TOTAL HESSIAN

   The total electronic second derivative [E_h/a0/a0] for the MP2 level of theory, (3 * {nat}, 3 * {nat}).

.. psivar:: MP2.5 TOTAL ENERGY
   MP2.5 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the MP2.5 level of theory.

.. psivar:: MP2-F12 CORRECTION ENERGY

   The component [E_h] correcting :psivar:`MP2 CORRELATION ENERGY` for
   explicit correlation for the MP2-F12/3C(FIX) level of theory.

.. psivar:: MP2-F12 TOTAL ENERGY
   MP2-F12 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   atop the :psivar:`HF-CABS TOTAL ENERGY` ,
   for the MP2-F12/3C(FIX) level of theory.

.. psivar:: MP3 TOTAL ENERGY
   MP3 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the MP3 level of theory.

.. psivar:: MP4(T) CORRECTION ENERGY

   The MP4 triples component [E_h]. Quantity is second right-hand term in
   Eq. :eq:`MP4terms`.

.. psivar:: MP4(SDQ) TOTAL ENERGY
   MP4(SDQ) CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the MP4 singles, doubles, quadruples level of theory.  Quantity
   :psivar:`MP4(SDQ) CORRELATION ENERGY` is
   first right-hand term in Eq. :eq:`MP4terms`.

.. psivar:: MP4 TOTAL ENERGY
   MP4 CORRELATION ENERGY
   MP4(SDTQ) TOTAL ENERGY
   MP4(SDTQ) CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the full MP4 level of theory. Quantity :psivar:`MP4 CORRELATION
   ENERGY` / :psivar:`MP4(SDTQ) CORRELATION ENERGY`
   is left-hand term in Eq. :eq:`MP4terms`.

   .. math:: E_{\text{MP4}} = E_{\text{MP4(SDQ)}} + E_{\text{MP4(T)}}
      :label: MP4terms

.. psivar:: MPn TOTAL ENERGY
   MPn CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the labeled |MollerPlesset| perturbation theory level.
   *n* is MP perturbation order.

.. psivar:: MP2 DOUBLES ENERGY
   MP2.5 DOUBLES ENERGY
   MP3 DOUBLES ENERGY
   CEPA(0) DOUBLES ENERGY
   CEPA(1) DOUBLES ENERGY
   CEPA(2) DOUBLES ENERGY
   CEPA(3) DOUBLES ENERGY
   ACPF DOUBLES ENERGY
   AQCC DOUBLES ENERGY
   CISD DOUBLES ENERGY
   QCISD DOUBLES ENERGY
   REMP2 DOUBLES ENERGY
   MP2-F12 DOUBLES ENERGY
   LCCD DOUBLES ENERGY
   CCD DOUBLES ENERGY
   LCCSD DOUBLES ENERGY
   CCSD DOUBLES ENERGY
   OMP2 DOUBLES ENERGY
   OMP2.5 DOUBLES ENERGY
   OMP3 DOUBLES ENERGY
   OREMP2 DOUBLES ENERGY
   OLCCD DOUBLES ENERGY

   The doubles portion [E_h] of the named correlation energy
   including same-spin and opposite-spin correlations.
   Explicitly correlated methods are with respect to the CABS-corrected reference.

.. psivar:: MP2 SINGLES ENERGY
   MP2.5 SINGLES ENERGY
   MP3 SINGLES ENERGY
   CEPA(0) SINGLES ENERGY
   CEPA(1) SINGLES ENERGY
   CEPA(2) SINGLES ENERGY
   CEPA(3) SINGLES ENERGY
   ACPF SINGLES ENERGY
   AQCC SINGLES ENERGY
   CISD SINGLES ENERGY
   QCISD SINGLES ENERGY
   REMP2 SINGLES ENERGY
   MP2-F12 SINGLES ENERGY
   LCCD SINGLES ENERGY
   CCD SINGLES ENERGY
   LCCSD SINGLES ENERGY
   CCSD SINGLES ENERGY
   OREMP2 SINGLES ENERGY
   OLCCD SINGLES ENERGY

   The singles portion [E_h] of the named correlation energy.
   Zero except in ROHF.
   Explicitly correlated methods are with respect to the CABS-corrected reference.

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
   REMP2 SAME-SPIN CORRELATION ENERGY
   MP2-F12 SAME-SPIN CORRELATION ENERGY
   LCCD SAME-SPIN CORRELATION ENERGY
   CCD SAME-SPIN CORRELATION ENERGY
   LCCSD SAME-SPIN CORRELATION ENERGY
   CCSD SAME-SPIN CORRELATION ENERGY
   OMP2 SAME-SPIN CORRELATION ENERGY
   OMP2.5 SAME-SPIN CORRELATION ENERGY
   OMP3 SAME-SPIN CORRELATION ENERGY
   OREMP2 SAME-SPIN CORRELATION ENERGY
   OLCCD SAME-SPIN CORRELATION ENERGY

   The unscaled portion [E_h] of the named correlation energy
   from same-spin or triplet doubles correlations.
   Explicitly correlated methods are with respect to the CABS-corrected reference.

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
   REMP2 OPPOSITE-SPIN CORRELATION ENERGY
   MP2-F12 OPPOSITE-SPIN CORRELATION ENERGY
   LCCD OPPOSITE-SPIN CORRELATION ENERGY
   CCD OPPOSITE-SPIN CORRELATION ENERGY
   LCCSD OPPOSITE-SPIN CORRELATION ENERGY
   CCSD OPPOSITE-SPIN CORRELATION ENERGY
   OMP2 OPPOSITE-SPIN CORRELATION ENERGY
   OMP2.5 OPPOSITE-SPIN CORRELATION ENERGY
   OMP3 OPPOSITE-SPIN CORRELATION ENERGY
   OREMP2 OPPOSITE-SPIN CORRELATION ENERGY
   OLCCD OPPOSITE-SPIN CORRELATION ENERGY

   The unscaled portion [E_h] of the named correlation energy
   from opposite-spin or singlet doubles correlations.
   Explicitly correlated methods are with respect to the CABS-corrected reference.

.. psivar:: MRPT TOTAL ENERGY
   MP2-CCSD TOTAL ENERGY
   MRCC TOTAL ENERGY

   Energies [E_h] from correlated multi-reference theories.

.. psivar:: MULLIKEN CHARGES

   Property of partial atomic charges [e] by the method of Mulliken, (nat,).

.. psivar:: NAUX (SCF)
   NAUX (CC)

   Convenience storage of number of functions [] in the auxiliary basis
   set for named stage of the calculation.

.. psivar:: NBODY (i, j, ..., k)@(a, b, ..., c) TOTAL ENERGY

   The total energy [E_h] of a component of the requested N-Body energy.
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

   The nuclear repulsion energy contribution [E_h] to the total SCF energy.
   Quantity :math:`E_{NN}` in Eq. :eq:`SCFterms`.

   .. math:: E_{NN} = \sum_{i, j<i}^{N_{atom}}\frac{Z_i Z_j}{|\mathbf{R}_i - \mathbf{R}_j|}
      :label: ENN

.. psivar:: OCEPA(0) TOTAL ENERGY
   OCEPA(0) CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the orbital-optimized CEPA(0) level of theory.

.. psivar:: OLCCD TOTAL ENERGY
   OLCCD CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the orbital-optimized linearized coupled cluster doubles level of theory.

.. psivar:: OLCCD REFERENCE CORRECTION ENERGY

   The difference [E_h] between the single-determinant energy of the final and
   initial orbitals for the orbital-optimized linearized coupled cluster
   doubles level of theory.

.. psivar:: OMP2 TOTAL ENERGY
   OMP2 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the orbital-optimized MP2 level of theory.

.. psivar:: OMP2 REFERENCE CORRECTION ENERGY

   The difference [E_h] between the single-determinant energy of the final and
   initial orbitals for the orbital-optimized MP2 level of theory.

.. psivar:: OMP2.5 TOTAL ENERGY
   OMP2.5 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the orbital-optimized MP2.5 level of theory.

.. psivar:: OMP2.5 REFERENCE CORRECTION ENERGY

   The difference [E_h] between the single-determinant energy of the final and
   initial orbitals for the orbital-optimized MP2.5 level of theory.

.. psivar:: OMP3 TOTAL ENERGY
   OMP3 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the orbital-optimized MP3 level of theory.

.. psivar:: OMP3 REFERENCE CORRECTION ENERGY

   The difference [E_h] between the single-determinant energy of the final and
   initial orbitals for the orbital-optimized MP3 level of theory.

.. psivar:: OREMP2 TOTAL ENERGY
   OREMP2 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the orbital-optimized retaining-the-excitation-degree |MollerPlesset|
   hybrid perturbation theory level.

.. psivar:: OREMP2 REFERENCE CORRECTION ENERGY

   The difference [E_h] between the single-determinant energy of the final and
   initial orbitals for the orbital-optimized retaining-the-excitation-degree
   |MollerPlesset| hybrid perturbation theory level.

.. psivar:: ONE-ELECTRON ENERGY

   The one-electron energy contribution [E_h] to the total SCF energy.
   Quantity :math:`E_{1e^-}` in Eq. :eq:`SCFterms`.

.. psivar:: PCM POLARIZATION ENERGY

   The energy contribution [E_h] from the polarizable continuum model for solvation.

.. psivar:: DD SOLVATION ENERGY

   The energy contribution [Eh] from continuum solvation models based on a
   domain-decomposition ansatz.

.. psivar:: PE ENERGY

   The energy contribution [E_h] from the polarizable embedding model for solvation.

.. psivar:: QCISD TOTAL ENERGY
   QCISD CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the quadratic configuration interaction singles and doubles level
   of theory.

.. psivar:: QCISD(T) TOTAL ENERGY
   QCISD(T) CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the quadratic configuration interaction singles and doubles with
   perturbative triples correction level of theory.

.. psivar:: QCISD(T) CORRECTION ENERGY

   The quadratic configuration interaction singles and doubles perturbative
   triples correction [E_h].

.. psivar:: REMP2 TOTAL ENERGY
   REMP2 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the retaining-the-excitation-degree |MollerPlesset| hybrid perturbation
   theory level.

.. psivar:: SAPT DISP ENERGY
   SAPT ELST ENERGY
   SAPT EXCH ENERGY
   SAPT IND ENERGY

   Respectively, the dispersion, electrostatics, exchange, and induction
   components of the total electronic interaction energy [E_h] for the
   requested SAPT level of theory. The sum of these four components yields
   :psivar:`SAPT TOTAL ENERGY`.

.. psivar:: SAPT TOTAL ENERGY
   SAPT ENERGY

   The total electronic interaction energy [E_h] for the requested SAPT
   level of theory.

.. psivar:: SAPT ELST10,R ENERGY

   An electrostatics-classified SAPT term energy [E_h] implemented for SAPT0.

.. psivar:: SAPT ELST EXTERN-EXTERN ENERGY

   Electrostatic interaction [E_h] between the point charges in fragments
   A and B in F/I-SAPT.

.. psivar:: SAPT EXCH10 ENERGY

   An exchange-classified SAPT term energy [E_h] implemented for SAPT0.

.. psivar:: SAPT EXCH10(S^2) ENERGY

   An exchange-classified SAPT term energy [E_h] implemented for SAPT0.

.. psivar:: SAPT IND20,R ENERGY
   SAPT EXCH-IND20,R ENERGY
   SAPT IND20,U ENERGY
   SAPT EXCH-IND20,U ENERGY

   An induction-classified SAPT term energy [E_h] implemented for SAPT0.

.. psivar:: SAPT DISP20 ENERGY
   SAPT EXCH-DISP20 ENERGY

   A dispersion-classified SAPT term energy [E_h] implemented for SAPT0.

.. psivar:: SAPT EXCH-DISP20(S^INF) ENERGY

   A dispersion-classified SAPT term energy [E_h] implemented for SAPT0. See :ref:`sec:saptinf`.

.. psivar:: SAPT SAME-SPIN DISP20 ENERGY
   SAPT SAME-SPIN EXCH-DISP20 ENERGY

   The portion of :psivar:`SAPT DISP20 ENERGY` or
   :psivar:`SAPT EXCH-DISP20 ENERGY` resulting from
   from same-spin or triplet doubles correlations.

.. psivar:: SAPT HF(2) ENERGY ABC(HF)

   The total Hartree--Fock energy [E_h] of the supersystem implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY AC(0)

   The Hartree--Fock energy [E_h] of subsystems A and C implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY BC(0)

   The Hartree--Fock energy [E_h] of subsystems B and C implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY A(0)

   The Hartree--Fock energy [E_h] of subsystem A implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY B(0)

   The Hartree--Fock energy [E_h] of subsystem B implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY AC(HF)

   The Hartree--Fock localized energy [E_h] of subsystems A and C implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY BC(HF)

   The Hartree--Fock localized energy [E_h] of subsystems B and C implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY AB(HF)

   The Hartree--Fock localized energy [E_h] of subsystems A and B implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY A(HF)

   The Hartree--Fock localized energy [E_h] of subsystem A implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY B(HF)

   The Hartree--Fock localized energy [E_h] of subsystem B implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY C

   The Hartree--Fock energy [E_h] of subsystem C implemented for F/I-SAPT.

.. psivar:: SAPT HF(2) ENERGY HF

   The FI-SAPT Hartree--Fock interaction energy [E_h] implemented for F/I-SAPT.

.. psivar:: SAPT ELST12,R ENERGY

   An electrostatics-classified SAPT term energy [E_h] implemented for SAPT2.

.. psivar:: SAPT EXCH11(S^2) ENERGY
   SAPT EXCH12(S^2) ENERGY

   An exchange-classified SAPT term energy [E_h] implemented for SAPT2.

.. psivar:: SAPT IND22 ENERGY
   SAPT EXCH-IND22 ENERGY

   An induction-classified SAPT term energy [E_h] implemented for SAPT2.

.. .. psivar:: SAPT HF TOTAL ENERGY
.. .. psivar:: SAPT CT ENERGY

.. psivar:: SAPT DISP21 ENERGY

   A dispersion-classified SAPT term energy [E_h] implemented for SAPT2+.

.. psivar:: SAPT DISP22(SDQ) ENERGY
   SAPT DISP22(T) ENERGY
   SAPT EST.DISP22(T) ENERGY

   Dispersion-classified MBPT-based SAPT term energy [E_h] implemented for SAPT2+.

.. psivar:: SAPT DISP2(CCD) ENERGY
   SAPT DISP22(S)(CCD) ENERGY
   SAPT DISP22(T)(CCD) ENERGY
   SAPT EST.DISP22(T)(CCD) ENERGY

   Dispersion-classified coupled-cluster-based SAPT term energy [E_h] implemented for SAPT2+.

.. psivar:: SAPT ELST13,R ENERGY

   An electrostatics-classified SAPT term energy [E_h] implemented for SAPT2+(3).

.. psivar:: SAPT IND30,R ENERGY
   SAPT IND-DISP30 ENERGY
   SAPT EXCH-IND30,R ENERGY

   A induction-classified SAPT term energy [E_h] implemented for SAPT2+3.

.. psivar:: SAPT EXCH-IND30(S^INF) ENERGY
   SAPT EXCH-IND30,R(S^INF) ENERGY

   A induction-classified SAPT term energy [E_h] implemented for SAPT2+3. See :ref:`sec:saptinf`.

.. psivar:: SAPT DISP30 ENERGY
   SAPT EXCH-DISP30 ENERGY
   SAPT EXCH-IND-DISP30 ENERGY

   A dispersion-classified SAPT term energy [E_h] implemented for SAPT2+3.

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
   components of the total electronic interaction energy [E_h] for the
   given SAPT level of theory. The sum of these four components yields
   the :samp:`{SAPT Level} TOTAL ENERGY`

.. psivar:: SAPT0 TOTAL ENERGY
   SSAPT0 TOTAL ENERGY
   SAPT2 TOTAL ENERGY
   SAPT2+ TOTAL ENERGY
   SAPT2+(3) TOTAL ENERGY
   SAPT2+3 TOTAL ENERGY

   The total electronic interaction energy [E_h] for the labeled SAPT level
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
   components of the total electronic interaction energy [E_h] for the
   given SAPT level of theory that incorporates coupled-cluster dispersion.
   The sum of these four components yields the :samp:`{SAPT Level} TOTAL ENERGY`

.. psivar:: SAPT2+(CCD) TOTAL ENERGY
   SAPT2+(3)(CCD) TOTAL ENERGY
   SAPT2+3(CCD) TOTAL ENERGY

   The total electronic interaction energy [E_h] for the labeled SAPT level
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
   components of the total electronic interaction energy [E_h] for the
   given SAPT level of theory that incorporates MP2 induction correction.
   The sum of these four components yields the :samp:`{SAPT Level} TOTAL ENERGY`

.. psivar:: SAPT2+DMP2 TOTAL ENERGY
   SAPT2+(3)DMP2 TOTAL ENERGY
   SAPT2+3DMP2 TOTAL ENERGY
   SAPT2+(CCD)DMP2 TOTAL ENERGY
   SAPT2+(3)(CCD)DMP2 TOTAL ENERGY
   SAPT2+3(CCD)DMP2 TOTAL ENERGY

   The total electronic interaction energy [E_h] for the labeled SAPT level
   of theory that incorporates MP2 induction correction.

.. psivar:: SAPT DFT GRAC SHIFT A
   SAPT DFT GRAC SHIFT B

   The gradient-regulation asymptotic correction (GRAC) [E_h] 
   used in a SAPT(DFT) computation for monomer A or B, respectively,
   to improve the accuracy by correctly describing the electron density
   at long-range.


.. psivar:: SCF ITERATIONS
   ADC ITERATIONS
   CCSD ITERATIONS
   OPTIMIZATION ITERATIONS

   Number of iterations [] in the named iterative method or optimization procedure.

.. psivar:: SCF DIPOLE

   Dipole array [e a0] for the SCF stage, (3,).

.. psivar:: SCF QUADRUPOLE

   Redundant quadrupole array [e a0^2] for the SCF stage, (3, 3).

.. psivar:: SCF TOTAL ENERGY

   The total electronic energy [E_h] of the SCF stage of the calculation.
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

   The total electronic second derivative [E_h/a0/a0] for the SCF stage, (3 * {nat}, 3 * {nat}).

.. psivar:: SCF STABILITY EIGENVALUES

   Array of eigenvalues from UHF or ROHF stability analysis.

.. psivar:: SCS-CCSD TOTAL ENERGY
   SCS-CCSD CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the CCSD-like method formed by reweighting :psivar:`CCSD DOUBLES ENERGY`
   by 1.27 opposite-spin and 1.13 same-spin contributions, with
   any singles carried along.

.. psivar:: SCS-MP2 TOTAL ENERGY
   SCS-MP2 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the MP2-like method formed by reweighting :psivar:`MP2 DOUBLES ENERGY`
   by 6/5 opposite-spin and 1/3 same-spin contributions, with
   any singles carried along.

.. psivar:: SCS-MP2-VDW TOTAL ENERGY
   SCS-MP2-VDW CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the MP2-like method formed by reweighting :psivar:`MP2 DOUBLES ENERGY`
   by 1.28 opposite-spin and 0.50 same-spin contributions, with
   any singles carried along. DOI: 10.1080/00268970802641242

.. psivar:: SCS(N)-MP2 TOTAL ENERGY
   SCS(N)-MP2 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
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

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the OMP2-like method formed by reweighting :psivar:`OMP2 DOUBLES ENERGY`
   by 6/5 opposite-spin and 1/3 same-spin contributions, with
   any singles carried along.

.. psivar:: SCS-MP3 TOTAL ENERGY
   SCS-MP3 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the MP3-like method formed by reweighting the difference between
   :psivar:`MP3 DOUBLES ENERGY` and :psivar:`MP2 DOUBLES ENERGY`
   by 0.25, atop the SCS-MP2 energy, with any singles carried along.

.. psivar:: SCS-OMP3 TOTAL ENERGY
   SCS-OMP3 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the OMP3-like method formed by reweighting the difference between
   :psivar:`OMP3 DOUBLES ENERGY` and :psivar:`OMP2 DOUBLES ENERGY`
   by 0.25, atop the SCS-OMP2 energy, with any singles carried along.

.. psivar:: SOS-MP2 TOTAL ENERGY
   SOS-MP2 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the MP2-like method formed by reweighting :psivar:`MP2 DOUBLES ENERGY`
   by 1.3 opposite-spin and 0 same-spin contributions, with
   any singles carried along.

.. psivar:: SOS-OMP2 TOTAL ENERGY
   SOS-OMP2 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the OMP2-like method formed by reweighting :psivar:`OMP2 DOUBLES ENERGY`
   by 1.2 opposite-spin and 0 same-spin contributions, with
   any singles carried along.

.. psivar:: SOS-OMP3 TOTAL ENERGY
   SOS-OMP3 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the OMP3-like method formed by reweighting the difference between
   :psivar:`OMP3 DOUBLES ENERGY` and :psivar:`OMP2 DOUBLES ENERGY`
   by 0.25, atop the SOS-OMP2
   energy using non-canonical weighting, with any singles carried along.

.. psivar:: SOS-PI-MP2 TOTAL ENERGY
   SOS-PI-MP2 CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the MP2-like method formed by reweighting :psivar:`MP2 DOUBLES ENERGY`
   by 1.4 opposite-spin and 0 same-spin contributions, with
   any singles carried along.

.. psivar:: TD-fctl ROOT 0 -> ROOT n ELECTRIC TRANSITION DIPOLE MOMENT (VEL)

   The electric transition dipole moment [e a0] in velocity gauge, for the transition
   from the ground state to root *m*.
   DFT functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 (IN h) -> ROOT n (IN i) ELECTRIC TRANSITION DIPOLE MOMENT (VEL)

   The electric transition dipole moment [e a0] in velocity gauge, for the transition
   from the ground state, which is of irrep *h*, to root *n* within irrep *i*.
   DFT functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 (h) -> ROOT n (i) ELECTRIC TRANSITION DIPOLE MOMENT (VEL)

   The electric transition dipole moment [e a0] in velocity gauge, for the transition
   from the ground state, which is of irrep *h*, to root *n*, which is of irrep *i*.
   DFT functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 -> ROOT n ELECTRIC TRANSITION DIPOLE MOMENT (VEL) - h TRANSITION

   The electric transition dipole moment [e a0] in velocity gauge, for the transition
   from the ground state to root *n*, and the transition is of *h* symmetry.
   DFT functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 -> ROOT n MAGNETIC TRANSITION DIPOLE MOMENT

   The magnetic transition dipole moment, for the transition
   from the ground state to root *n*.
   DFT functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 (IN h) -> ROOT n (IN i) MAGNETIC TRANSITION DIPOLE MOMENT

   The magnetic transition dipole moment, for the transition
   from the ground state, which is of irrep *h*, to root *n* within irrep *i*.
   DFT functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 (h) -> ROOT n (i) MAGNETIC TRANSITION DIPOLE MOMENT

   The magnetic transition dipole moment, for the transition
   from the ground state, which is of irrep *h*, to root *n*, which is of irrep *i*.
   DFT functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 -> ROOT n MAGNETIC TRANSITION DIPOLE MOMENT - h TRANSITION

   The magnetic transition dipole moment, for the transition
   from the ground state to root *n*, and the transition is of *h* symmetry.
   DFT functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 -> ROOT n LEFT EIGENVECTOR ALPHA

   The left alpha spin eigenvectors of the named method
   from ground state to root *n*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 (IN h) -> ROOT n (IN i) LEFT EIGENVECTOR ALPHA

   The left alpha spin eigenvectors of the named method
   from ground state, which is in irrep *h*, to root *n* within irrep *i*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 (h) -> ROOT n (i) LEFT EIGENVECTOR ALPHA

   The left alpha spin eigenvectors of the named method
   from ground state, which is in irrep *h*, to root *n*, which is in irrep *i*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 -> ROOT n LEFT EIGENVECTOR ALPHA - h TRANSITION

   The left alpha spin eigenvectors of the named method
   from ground state to root *n*, and the transition is of irrep *h*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 -> ROOT n LEFT EIGENVECTOR BETA

   The left beta spin eigenvectors of the named method
   from ground state to root *n*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 (IN h) -> ROOT n (IN i) LEFT EIGENVECTOR BETA

   The left beta spin eigenvectors of the named method
   from ground state, which is in irrep *h*, to root *n* within irrep *i*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 (h) -> ROOT n (i) LEFT EIGENVECTOR BETA

   The left beta spin eigenvectors of the named method
   from ground state, which is in irrep *h*, to root *n*, which is in irrep *i*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 -> ROOT n LEFT EIGENVECTOR BETA - h TRANSITION

   The left beta spin eigenvectors of the named method
   from ground state to root *n*, and the transition is of irrep *h*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 -> ROOT n RIGHT EIGENVECTOR ALPHA

   The right alpha spin eigenvectors of the named method
   from ground state to root *n*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 (IN h) -> ROOT n (IN i) RIGHT EIGENVECTOR ALPHA

   The right alpha spin eigenvectors of the named method
   from ground state, which is in irrep *h*, to root *n* within irrep *i*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 (h) -> ROOT n (i) RIGHT EIGENVECTOR ALPHA

   The right alpha spin eigenvectors of the named method
   from ground state, which is in irrep *h*, to root *n*, which is in irrep *i*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 -> ROOT n RIGHT EIGENVECTOR ALPHA - h TRANSITION

   The right alpha spin eigenvectors of the named method
   from ground state to root *n*, and the transition is of irrep *h*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 -> ROOT n RIGHT EIGENVECTOR BETA

   The right beta spin eigenvectors of the named method
   from ground state to root *n*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 (IN h) -> ROOT n (IN i) RIGHT EIGENVECTOR BETA

   The right beta spin eigenvectors of the named method
   from ground state, which is in irrep *h*, to root *n* within irrep *i*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 (h) -> ROOT n (i) RIGHT EIGENVECTOR BETA

   The right beta spin eigenvectors of the named method
   from ground state, which is in irrep *h*, to root *n*, which is in irrep *i*. DFT
   functional labeled if canonical.

.. psivar:: TD-fctl ROOT 0 -> ROOT n RIGHT EIGENVECTOR BETA - h TRANSITION

   The right alpha and beta spin eigenvectors of the named method
   from ground state to root *n*, and the transition is of irrep *h*. DFT
   functional labeled if canonical.

.. psivar:: THERMAL ENERGY

   Total thermal energy E [E_h] at given temperature.

.. psivar:: THERMAL ENERGY CORRECTION

   Sum of electronic, translational, rotational, and vibrational corrections [E_h] to the thermal energy at given temperature.

.. psivar:: TWO-ELECTRON ENERGY

   The two-electron energy contribution [E_h] to the total SCF energy.
   Quantity :math:`E_{2e^-}` in Eq. :eq:`SCFterms`.

.. psivar:: UNCP-CORRECTED 2-BODY INTERACTION ENERGY

   The interaction energy [E_h] considering only two-body interactions,
   computed without counterpoise correction.
   Related variable :psivar:`CP-CORRECTED 2-BODY INTERACTION ENERGY`.

   .. math:: E_{\text{IE}} = E_{dimer} - \sum_{monomer}^{n}{E_{monomer}^{\text{unCP}}}

.. psivar:: WIBERG LOWDIN INDICES

   Property of Wiberg bond indices using orthogonal L\ |o_dots|\ wdin orbitals, (nat, nat).

.. psivar:: ZAPTn TOTAL ENERGY
   ZAPTn CORRELATION ENERGY

   The total electronic energy [E_h] and correlation energy component [E_h]
   for the labeled Z-averaged perturbation theory level.
   *n* is ZAPT perturbation order.

.. psivar:: ZERO K ENTHALPY

   Total electronic and zero-point energy [E_h] at 0 [K].

.. psivar:: ZPVE

   Vibrational zero-point energy [E_h] at 0 [K].

.. psivar:: 2-BODY PAIRWISE DISPERSION CORRECTION ANALYSIS

   The interatomic contributions to the dispersion correction [E_h].
   Sums to the dispersion energy.

