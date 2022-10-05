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

.. _`sec:prog_integrals`:

======================
Integrals in |PSIfour|
======================

Introduction
------------

|PSIfour| has a number of backends available to compute integrals. In order to
accomodate these options, while providing a clean interface to the programmer,
an abstraction layer is implemented within Libmints.  A recent upgrade to the
primary integral engine has seen some important changes to the way this
interface layer is used; this document is designed to aid new developers as
well as those familiar with the older calling conventions to ensure that the
most efficient calling conventions are applied.

The older style
---------------

A very simple loop that does not use permutational symmetry might look
something like this in the old scheme:

.. code-block:: cpp

    auto sieve = std::make_shared<ERISieve>(basisset, cutoff);
    auto factory= std::make_shared<IntegralFactory>(basisset);
    int deriv_level = 0;
    bool use_shell_pairs = true;
    auto eri = factory->eri(deriv_level, use_shell_pairs);
    const double* buffer = eri_->buffer();
    for (int P = 0; P < basisset->nshell(); ++P) {
        const auto& Pshell = basisset->shell(P);
        for (int Q = 0; Q < basisset->nshell(); ++Q) {
            const auto& Qshell = basisset->shell(Q);
            for (int R = 0; R < basisset->nshell(); ++R) {
                const auto& Rshell = basisset->shell(R);
                for (int S = 0; S < basisset->nshell(); ++S) {
                    const auto& Sshell = basisset->shell(S);
                    if(sieve->shell_significant(P, Q, R, S) {
                        eri->compute_shell(P, Q, R, S);
                        // results are in buffer, do something with them..
                    }
                }
            }
        }
    }

An integral factory is used, which can then produce integral object for various
operator types and derivative levels.  A sieve is also constructed; this allows
a quick determination of whether an integral shell quartet will be significant
in magnitude or not, potentially saving a lot of work.  This simple scheme is
clean and easy to understand, and is still supported in the latest version of
|PSIfour| with only a small change to the sieve syntax and handling of buffer
addresses, noted below.

The new syntax
--------------

The newer integral engines being interfaced to |PSIfour| may or may not require
a group of similar integrals to be computed together in a block using
vectorized instructions.  To accomodate this possibility, a new syntax has been
introduced in Libmints:

.. code-block:: cpp

    auto blocksPQ = ints[0]->get_blocks12();
    auto blocksRS = ints[0]->get_blocks34();

    auto factory= std::make_shared<IntegralFactory>(basisset);
    int deriv_level = 0;
    bool use_shell_pairs = true;
    bool needs_exchange = true;
    auto eri = factory->eri(deriv_level, use_shell_pairs, needs_exchange);
    const auto &buffers = eri->buffers();

    eri->update_density(D);
    bool use_batching = eri->maximum_block_size() > 1;

    // loop over all the blocks of (P>=Q|
    for (size_t blockPQ_idx = 0; blockPQ_idx < blocksPQ.size(); blockPQ_idx++) {
        const auto& blockPQ = blocksPQ[blockPQ_idx];
        // loop over all the blocks of |R>=S)
        size_t start = eri->first_RS_shell_block(blockPQ_idx);
        for (int blockRS_idx = loop_start; blockRS_idx < blocksRS.size(); ++blockRS_idx) {
            const auto& blockRS = blocksRS[blockRS_idx];

            if (!eri->shell_block_significant(blockPQ_idx, blockRS_idx)) continue;

            eri->compute_shell_blocks(blockPQ_idx, blockRS_idx);
            const auto* block_start = buffers[0];

            // Loop over all of the P,Q,R,S shells within the blocks.  We have P>=Q, R>=S and PQ<=RS.
            for (const auto& pairPQ : blockPQ) {
                const auto &P = pairPQ.first;
                const auto &Q = pairPQ.second;
                const auto& Pshell = basisset->shell(P);
                const auto& Qshell = basisset->shell(Q);
                const auto Pam = Pshell.am();
                const auto Qam = Qshell.am();
                for (const auto& pairRS : blockRS) {
                    const auto &R = pairRS.first;
                    const auto &S = pairRS.second;
                    const auto& Rshell = basisset->shell(R);
                    const auto& Sshell = basisset->shell(S);
                    const auto Ram = Rshell.am();
                    const auto Sam = Sshell.am();

                    size_t block_size = Psize * Qsize * Rsize * Ssize;
                    // When there are chunks of shellpairs in RS, we need to make sure
                    // we filter out redundant combinations.
                    if (use_batching && Pam == Ram && Qam == Sam && ((P > R) || (P == R && Q > S))) {
                        block_start += block_size;
                        continue;
                    }
                    const double* int_ptr = block_start;
                    // Query P,Q,R,S shells for metadata and loop over that quartet
                    // as usual, getting the integrals from the int_ptr buffer.
                    block_start += block_size;
                }
            }
        }
    }

Although this looks more complex, it's essentially doing the same thing.  There
are a number of differences that we'll highlight now.

Sieving
.......

This is one of two breaking changes to the old style syntax.  Instead of
constructing a sieve object, the integral object should be queried directly
using the exact same syntax.  Requests for whether a shell is significant or a
shell block is significant are both supported.  A sieve object is created if
matching basis sets are found in either the bra or the ket.  For a density
fitting integral (PQ|0A) where 0 is the null basis set and A is an auxiliary
basis set the (PQ| pair will be used to construct all of the sieving data.

Buffer address
..............

The old code copied integrals into a buffer owned by the integral object, whose
address remained constant and could be retrieved by the ``buffer()`` member
function.  To avoid unnecessary copies, the new code instead uses the integrals
directly from the underlying integral engine's memory, which may change with
each call to compute integrals.  The integral engine provides a
``std::vector<const double*>`` containing the pointers to the start of each
"chunk" of integrals.  For first derivatives there are 12 such "chunks", which
are ordered Px,Py,Pz,Qx,Qy,Qz,Rx,Ry,Rz,Sx,Sy,Sz, where the Px refers to the x
derivative with respect to the basis functions in shell P.  Note that all
integral derivatives are provided by the new integral code, unlike the previous
version where only 9 of 12 were provided and the user was responsible for using
translation invariance relationships to fill in the rest.  The addresses for
each chunk are updated in the vector after each call to compute integrals, so
the user should keep a const reference to that object, and query that for the
address of interest.

Density Screening
.................

The old code looked only at the integral to determine whether terms can be
avoided *a priori*.  However, if the integral is to be contracted with a
density or a density-like quantity, the screening can be performed on the
product, which yields more sparsity.  To enable this, simply call the integral
object's ``update_density`` member, passing it a SharedMatrix holding the
current density (remember that it changes during each iteration of the SCF) and
the product will be considered during screening.  If only coulomb-like terms
are to be computed, the ``needs_exchange`` argument to the integral object
constructor should be set to false, otherwise it should be true to correcly
account for products of the density and integrals that contribute to
exchange-like terms.

Shell blocking
..............

Each underlying integral engine knows whether it will use blocks, and will set up
the metadata automatically.  Instead of looping over individual shells, the
user should loop over blocks supplied by the integral object; these blocks will
be just a single shell quartet combination for the case where blocking is not
used. It is simple to loop over pairs within each block using C++11 syntax, as
demonstrated in the code snippet above.  Only shell pairs with significant
overlap are included in the shell block information, making this an efficient
way to loop over non-negligible terms.

Permutational symmetry
......................

The pairs within each block are optimized for efficiency.  First, they are
screened during the integral object's creation to ensure that only terms with
appreciable overlap are stored.  Second, only P,Q combinations that are
permutationally unique are stored, ordered with the higher angular momentum
first.  Therefore care must be taken to ensure that the missing permutations
are correctly accounted for when processing the integrals within the loop.  See
the DirectJK code in libfock for an example of using this scheme for a Fock
matrix build.

Using bra-ket symmetry
......................

In cases where there is no batching performed, bra-ket symmetry can be
trivially enforced by ensuring that one of the block indices is greater than or
equal to the other.  When batching is used, the situation is trickier; some ket
batches may contain a mixture of integrals that are bra-ket unique and those
that are not.  To handle this we must do a coarse check at the top of the loop
to see if *any* integrals in the batch are needed, which is implemented by
asking the integral engine where to start looping in the ket via the call to
``eri->first_RS_shell_block(PQpair_idx)``.  This is followed by a more fine
grained check within the loops to filter individual integrals in the case where
bra and ket have the same angular momentum and there's a possibility of a
handful of integrals coming from the ket that are redundant.  Note that the bra
is not batched in any of our engines currently: only the ket is.  For this
reason, density fitting integrals should be written as (A0|PQ) rather than
(PQ|A0) where possible, because we want the ket to contain more functions than
the bra for efficient blocking.

Instantiating integral objects
..............................

With sieving being introduced in the new integral objects, the cost of their
construction has increased.  Although significantly cheaper than computing
integrals themselves, construction of integral objects can be non-negligible,
especially if many threads are used.  For example, this pattern can be found in
old versions of the code:

.. code-block:: cpp

    std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
    ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri()));
    for (int thread = 1; thread < num_threads; thread++) {
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri()));
    }

This builds many objects and the cost can add up.  With the new scheme,
integral objects are forced to implement a `clone()` member that can be used as
follows:

.. code-block:: cpp

    std::vector<std::shared_ptr<TwoBodyAOInt>> ints;
    ints.push_back(std::shared_ptr<TwoBodyAOInt>(factory->eri()));
    for (int thread = 1; thread < num_threads; thread++) {
        ints.push_back(std::shared_ptr<TwoBodyAOInt>(ints[0]->clone()));
    }

This method only incurs the cost of creating a single integral object, and
performs much cheaper cloning operations to create the other objects for each
thread.  Moreover, if integral objects are created only in the initialization
of each code that uses them, and stored persistently, the cost of integral
object creation is further reduced.


One Electron Integrals in |PSIfour|
-----------------------------------

After version 1.5, we started transitioning the one electron integral code over to
use Libint2 instead of the old handwritten Obara-Saika code.  There are a
number of reasons motivating this switch.  For methods requiring
potentials and fields evaluated at many external sites, such as PCM and
polarizable embedding, the efficiency of the one electron integrals can be rate
limiting.  We also started to introduce integral screening, and it is important
to balance the screening used for one- and two-electron terms carefully, so this
is a good opportunity to re-evaluate the code.  Finally, given the complexity
of the OS recursion code, the switch to an external library leaves a more
compact codebase to maintain.  The one electron integrals which are not provided by Libint2
are now handled by a new implementation of the McMurchie-Davidson (M-D) algorithm,
leading to removal of the OS code in version 1.6. An overview of the one electron integrals
is shown in table :ref:`table:oei_impl_summary`, together with the implementation they use.
The tips below serve as a guide to what changed,
why it changed, and how to interface with |PSIfour|'s one-electron integral
machinery now.

Calling ``compute_shell(int P, int Q)``
.......................................

The hand-implemented OS recursion code also took care of the Cartesian->pure
transformation (if required by the basis set).  The mechanism for handling this
was to provide a public facing ``compute_shell(int P, int Q)`` method for the
caller; this then looked up the appropriate ``GaussianShell`` objects that were
passed into the corresponding (private) ``compute_pair(GaussianShell &s1,
GaussianShell &s2)`` function that computed the integrals and transformed them
to the spherical harmonic basis, if needed.  The switch to Libint2 integrals
preserves this mechanism, but the ``compute_shell(int P, int Q)`` simply looks
up the appropriate Libint2-compatible shells and hands them off to the
re-written, private ``compute_pair()`` routines, which call Libint2 directly.
Therefore, any calls to shell-pair level integral computations should look the
same as before the introduction of Libint2, however access to the integrals has
changed, as described below.

Accessing integrals
...................

Before the Libint2 transition, one electron integrals were computed in a flat
array, internally called `buffer_`, which was accessed through the integral
object's ``buffer()`` method.  For integrals with multiple operators, e.g.,
dipole operators that have three distinct components, the buffer was simply
elongated by the appropriate amount and the caller was responsible for striding
through each resulting batch correctly.  The Libint2 engines instead return a
list of pointers into each operator's batch of integrals, the ordering of which
are detailed on the Libint2 wiki.  For this reason, the call to ``buffer()``
that returns a single buffer must be replaced with a call to ``buffer()`` to
get a list of pointers; we recommend that be assigned the type ``const auto
&``.  For simple integrals, such as overlap or kinetic, only the buffer
corresponding to the zeroth element of this array contains integrals.

Derivative Integrals
....................

The old one electron integral code used translational invariance relations to
minimze the number of integrals to be computed, leaving the caller with some
bookkeeping to do to compute all terms.  For example, consider an overlap
integral: its value depends only on the relative separation of the two centers
and not their absolute positions in space.  Therefore, the derivative with
respect to center A is the negative of the same derivative with respect to
center B, so one is trivially gleaned from the other.  Extending this to second
derivatives, the same principle leads to the fact that double derivatives with
respect to center A are equal to double derivatives with respect to center B,
which are also equal to the negative of the mixed double derivatives with
respect to both center A and B.  The old code only provided the double
derivative with respect to center A, leaving the caller to determine the other
values.  The Libint2 engine instead provides all integrals, so the caller
simply needs to loop over all of the buffers provided in the appropriate order.

Changes to External Potential Engines
.....................................

Benchmarking showed that early versions of the old code spent a non-negligible
amount of time performing the Cartesian to spherical harmonic transformation of
the integrals, which is needed for most modern basis sets.  To improve
performance, we instead backtransformed the density to the Cartesian
representation (denoted "CartAO") and computed / contracted all integrals in
this Cartesian basis, eliminating the need to transform to spherical harmonics
as the integrals are computed.  This bottleneck no longer exists, so these
extra transformation steps have been removed as part of the switch to Libint2,
and the affected codes (PCM and CPPE interfaces) now compute the potential and
field integrals in the representation required by the basis set.

Also, note that the way external point charges are specified has changed.
Previously, a set of N external point charges would be specified by passing a
matrix with dimensions N rows and 4 columns -- corresponding to charge, x, y, z
-- to the ``set_charge_field()`` member of the potential integral class.  The
same information is now passed using the more verbose
``std::vector<std::pair<double, std::array<double, 3>>>`` type instead, to be
consistent with Libint2's convention.

New Operators Available
.......................

Libint2 provides a range of integrals that were previously not available in
|PSIfour|, such as the Erfc attenuated nuclear potential integrals needed for
Ewald methods.  If new integrals are added to Libint2 but are not yet
interfaced to |PSIfour|, please open an issue on the |PSIfour| GitHub page to
alert the developers, who will be able to add the appropriate code.
Available integrals classes and parameters currently documented at 
`Libint2 C++11 Interface Wiki <https://github.com/evaleev/libint/wiki/using-modern-CPlusPlus-API#create-an-integral-engine>`_


Shell Pairs
...........

To ensure consistency between one- and two-electron terms when screening, and
for efficiency reasons, shell pair lists should be used to iterate over pairs
of Gaussian shells.  These lists contain integer pair numbers, corresponding to
the pairs of shells that have sufficient overlap to survive the screening
process.  Iterating over these lists is simple:

.. code-block:: cpp

    const auto& shell_pairs = Vint->shellpairs();
    size_t n_pairs = shell_pairs.size();
    for (size_t p = 0; p < n_pairs; ++p) {
         auto P = shell_pairs[p].first;
         auto Q = shell_pairs[p].second;
         // do something with shells P and Q
    }

Note that list considers all P,Q pairs if the two basis sets differ, but only
P>=Q if the basis sets are the same; the caller should account for this
restricted summation in the latter case.


One Electron Integral Algorithm Overview
........................................

The following table summarizes which implementation is used
for each type of one electron integral in |PSIfour|.

.. _`table:oei_impl_summary`:

.. table:: Algorithms used for One Electron Integrals

    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | Integral                   | Class                      | Implementation  | Comment                                                            |
    +============================+============================+=================+====================================================================+
    | Three-Center Overlap       | ``ThreeCenterOverlapInt``  | Libint2         | using ``libint2::Operator::delta`` for 4-center integrals          |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | Angular Momentum           | ``AngularMomentumInt``     | M-D             |                                                                    |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | Dipole                     | ``DipoleInt``              | Libint2         | no derivatives supported                                           |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | Electric Field             | ``ElectricFieldInt``       | Libint2         | using first derivative of ``libint2::Operator::nuclear``           |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | Coulomb Potential          | ``ElectrostaticInt``       | Libint2         | evaluated for a single origin and unity charge                     |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | Kinetic                    | ``KineticInt``             | Libint2         |                                                                    |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | Multipole Potential        | ``MultipolePotentialInt``  | M-D             | arbitrary order derivative of 1/R supported                        |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | Multipole Moments          | ``MultipoleInt``           | M-D             | arbitrary order multipoles supported, including nuclear gradients  |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | Nabla Operator             | ``NablaInt``               | Libint2         | using first derivative of ``libint2::Operator::overlap``           |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | Overlap                    | ``OverlapInt``             | Libint2         |                                                                    |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | Nuclear Coulomb Potential  | ``PotentialInt``           | Libint2         | assumes nuclear centers/charges as the potential                   |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | PCM Potential              | ``PCMPotentialInt``        | Libint2         | parallelized over charge points                                    |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | Quadrupole                 | ``QuadrupoleInt``          | Libint2         |                                                                    |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | Traceless Quadrupole       | ``TracelessQuadrupoleInt`` | Libint2         |                                                                    |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+
    | Relativistic Potential     | ``RelPotentialInt``        | Libint2         |                                                                    |
    +----------------------------+----------------------------+-----------------+--------------------------------------------------------------------+

