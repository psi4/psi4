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

.. _`sec:blaslapack`:

===========================
Linear Algebra in |PSIfour|
===========================

.. _`faq:blaswrappers`:

How to call BLAS & LAPACK in |PSIfour|
--------------------------------------

Computational chemistry is essentially linear algebra on molecular
systems, so using stable, portable, scalable, and efficient numerical
linear algebra methods in |PSIfour| is critical. To that end, we use BLAS1
(vector-vector operations, like dot products), BLAS2 (matrix-vector
operations, like rank-1 update), BLAS3 (matrix-matrix operations, like
matrix multiplication), and LAPACK (advanced matrix decompositions and
solutions). The methods provided by BLAS and LAPACK are standard, but the
performance of actual implementations differ greatly from one version to
another. Moreover, the standard interfaces to the libraries are Fortran,
so |PSIfour| provides a common set of wrappers in :source:`psi4/src/psi4/libqt/qt.h` .

.. warning:: Although block_matrix, init_array, and print_mat are still
   around, their use is discouraged in favor of operations on
   `psi4.core.Matrix` itself. The advice in these docs will catch up
   shortly.

BLAS Wrappers
^^^^^^^^^^^^^

BLAS wrappers are currently fully supported at double precision.

BLAS commands involving matrices are wrapped so as to be conventional
C-style "row-major" indexing, meaning that the column is the fast index
like normal.

* The calls to BLAS1 routines are wrapped so as to allow for operations on
  vectors with more than 2^{31} elements (~16 GB, getting to be a problem).
  So passing a signed or unsigned long works, though the stride arguments
  must be integers.

* All routines are declared in ``qt.h``. Each routine is prefixed with a
  ``C_``, followed by the standard Fortran name of the routine, in capital
  letters. Input parameters of single primitives (``int``, ``double``,
  ``unsigned long int``, ``char``, ...) are passed by value. Arrays,
  including multidimensional arrays, are required to be in contiguous
  memory (as provided by block_matrix, for example), and are passed by
  providing a pointer to the first double or int element of the data (this
  is array[0] if array is ``double**``). BLAS1 routines occasionally
  return values (DDOT for instance), BLAS2 and BLAS3 always return void.
  For char arguments, case is insensitive. A few examples are provided::

    // BLAS/LAPACK
    #include "psi4/libqt/qt.h"
    // block_matrix, init_array
    #include "psi4/libciomr/libciomr.h"
    
    using namespace psi;
    ...
    
    // Allocate a,b vectors
    int n = 100;
    double* a = init_array(n);
    double* b = init_array(n);
    
    // Allocate A matrix;
    double** A = block_matrix(n,n);
    double** B = block_matrix(n,n);
    double** C = block_matrix(n,n);
    
    // Call the BLAS1 dot product between a and b
    // n can be a ULI with the BLAS1 wrappers,
    // All strides must be ints though
    double dot = C_DDOT(n, a, 1, b, 1);
    
    // Call the BLAS2 GEMV without transposition
    // Note this works in row-major order
    C_DGEMV('N', n, n, 1.0, A[0], n, a, 1, 0.0, b, 1);
    
    // Call the BLAS3 GEMM without transposition
    // Note this works in row-major order
    C_DGEMM('N','N', n, n, n, 1.0, A[0], n, B[0], n, 0.0, C[0], n);
    
    // Array's init'd with init_array must be free'd, not delete[]'d
    free(a);
    free(b);
    
    // Block matrix should be free_blocked
    free_block(A);
    free_block(B);
    free_block(C);

Important BLAS Routines
^^^^^^^^^^^^^^^^^^^^^^^

* BLAS1

  * DDOT: dot product
  * DCOPY: efficient memory copy (with variable stride)
  * DAXPY: y = y + alpha*x
  * DROT: Givens Rotation
  * DNRM2: Vector norm square

* BLAS2

  * DGEMV: General Matrix-Vector product
  * DTRMV: Triangular Matrix-Vector product (2x faster, not wrapped yet)
  * DTRSM: Triangular Matrix-Vector solution via back substitution (just as fast as DTRMV)
  * DGER: Rank-1 update (not wrapped yet)

* BLAS3

  * DGEMM: General Matrix-Matrix product
  * DTRMM: General Triangular Matrix-General Matrix product
  * DTRSM: Triangular Matrix-General Matrix solution via back substitution (just as fast as DTRMM)
  * DSYMM/DSYMV calls are not appreciably faster than DGEMM calls, and should only be used in expert situations (like using the other half of the matrix for some form of other transformation).
  * DTRMM/DTRMV calls are 2x faster than DGEMM, and should be used where possible.

LAPACK Wrappers
^^^^^^^^^^^^^^^

All standard LAPACK 3.2 double precision routines are provided.

LAPACK commands remain in Fortran's "column-major" indexing, so all the
results will be transposed, and leading dimensions may have to be fiddled
with (using ``lda = n`` in both directions for square matrices is highly
recommended). An example of the former problem is a Cholesky
Decomposition: you expect to get back a lower triangular matrix L such
that ``L L^T = A``, but this is returned in column-major order, so the actual
recovery of the matrix A with the row-major BLAS wrappers effectively
involves ``L^T L = A``. On of the biggest consequences is in linear equations:
The input/output forcing/solution vector must be explicitly formed in
column-major indexing (each vector is placed in a C++ row, with its
entries along the C++ column). This is visualized in C++ as the transpose
of the forcing/solution vector.  All routines are declared in qt.h. Each
routine is prefixed with a ``C_``, followed by the standard Fortran name of
the routine, in capital letters. Input parameters of single primitives
(int, double, unsigned long int, char, ...) are passed by value. Arrays,
including multidimensional arrays, are required to be in contiguous memory
(as provided by block_matrix, for example), and are passed by providing a
pointer to the first double or int element of the data (this is array[0]
if array is ``double**``). All routines return an int INFO with error and
calculation information specific to the routine, In Fortran, this is the
last argument in all LAPACK calls, but should not be provided as an
argument here. For char arguments, case is insensitive. A Cholesky
transform example is shown::

    // BLAS/LAPACK
    #include "psi4/libqt/qt.h"
    // block_matrix, init_array
    #include "psi4/libciomr/libciomr.h"
    
    using namespace psi;
    ...
    int n = 100;
    
    // Allocate A matrix;
    double** A = block_matrix(n,n);
    
    // Call the LAPACK DPOTRF to get the Cholesky factor
    // Note this works in column-major order
    // The result fills like:
    //   * * * *
    //     * * *
    //       * *
    //         *
    // instead of the expected:
    //   *
    //   * *
    //   * * *
    //   * * * *
    //
    int info = C_DPOTRF('L', n, A[0], n);
    
    // A bit painful, see below
    fprintf(outfile, "A:\n");
    print_mat(A,n,n,outfile);

    // Block matrix should be free_blocked
    free_block(A);

Important Lapack Routines
^^^^^^^^^^^^^^^^^^^^^^^^^

* DSYEV: Eigenvalues and, optionally eigenvectors of a symmetric matrix. Eigenvectors take up to 10x longer than eigenvalues.
* DGEEV: Eigenvalues and, optionally eigenvectors of a general matrix. Up to 10x slower than DSYEV.
* DGESV: General solver (uses LU decomposition).
* DGESVD: General singular value decomposition.
* DGETRF: LU decomposition.
* DPOTRF: Cholesky decomposition (much more stable/faster)
* DGETRS: Solver, given LU decomposition by DGETRF
* DPOTRS: Solver, given Cholesky decomposition by DPOTRF
* DGETRI: Inverse, given LU decomposition by DGETRF (Warning: it's faster and more stable just to solve with DGETRS)
* DPOTRI: Inverse, given Cholesky decomposition by DPOTRF (Warning: it's faster and more stable just to solve with DPOTRS)


.. _`faq:blasmatrix`:

How to use low-level BLAS/LAPACK with ``psi4.core.Matrix``
----------------------------------------------------------

Jet's awesome new Matrix object has a lot of simple BLAS/LAPACK built in,
but you can just as easily use the ``double***`` array underneath if you are
careful (the outer index is the submatrix for each irrep). Here's an
example:

.. code-block:: cpp

    // BLAS/LAPACK
    #include "psi4/libqt/qt.h"
    // Matrix
    #include "psi4/libmints/matrix.h"
    
    using namespace psi;
    ...
    int n = 100;
    
    // Allocate A Matrix (new C1 convenience constructor);
    shared_ptr<Matrix> A(new Matrix("Still A, but way cooler", n,n));
    // Get the pointer to the 0 irrep (C1 for now, it errors if you ask for too high of an index)
    double** A_pointer = A->get_pointer(0);
    
    // Call the LAPACK DPOTRF to get the Cholesky factor
    // Note this works in column-major order
    // The result fills like:
    //   * * * *
    //     * * *
    //       * *
    //         *
    // instead of the expected:
    //   *
    //   * *
    //   * * *
    //   * * * *
    //
    int info = C_DPOTRF('L', n, A_pointer[0], n);
    
    // Wow that's a lot easier
    A->print();
    
    // Don't free, it's shared_ptr!


.. _`faq:labas`:

How to name orbital bases (e.g., AO & SO)
-----------------------------------------

Many different working bases (the internal linear algebraic basis, not the
name of the Gaussian basis) are used within |PSIfour|, each with a unique
and important purpose. It is critical to keep them all distinct to prevent
weird results from occurring.

* ``AO`` (Atomic Orbitals): Cartesian Gaussians (6D, 10F, etc.),
  ``(L + 1)(L + 2)/2`` functions per shell of angular momentum L. The
  ordering of Cartesian exponents for a given L is given by the standard
  ordering below (MATLAB code)::

    ncart = (L + 1) * (L + 2) / 2;
    exps = zeros(ncart,3);
    index = 1;
    for i = 0:L
        for j = 0:i
            lx = L - i;
            ly = i - j;
            lz = j;
            exps(index,:) = [lx ly lz];
          index = index + 1;
        end
    end

* ``SO`` (Spherical Atomic Orbitals): Pure Gaussians (5D, 7F, etc.) or
  Cartesian Gaussians, as determined by the user. This is typically the
  first layer encountered, Libmints handles the transform from AO to SO
  automatically. If Cartesian functions are used, the number of functions
  per shell remains ``(L + 1)(L + 2)/2``, and the ordering remains the same
  as above. Note that the individual functions are not normalized for
  angular momentum as in most codes: the self-overlap of a |PSIfour| Cartesian D
  or higher function with more than one nonzero Cartesian exponent (e.g., lx
  = 1, ly = 1, lz = 0) will be less than one. If Spherical Harmonics are
  used, 2L + 1 real combinations of the spherical harmonics are built from
  the ``(L+1)(L+2)/2`` Cartesian Gaussians, according to H. Schlegel and M.
  Frish, IJQC, 54, 83-87, 1995. Unlike Cartesian functions these functions
  are all strictly normalized. Note that in |PSIfour|, the real combinations of
  spherical harmonic functions (see the paragraph below Eq. 15 in the
  Schlegel paper) are ordered as: 0, 1+, 1-, 2+, 2-, ....

* ``USO`` (Unique Symmetry-Adapted Orbitals): Spatial symmetry-adapted
  combinations of SOs, blocked according to irrep. The total number of USOs
  is the same as the number of SOs, but the number of USOs within each irrep
  is usually much smaller, which can lead to significant performance
  improvements. Note that this basis is sometimes unfortunately referred to
  as the SO basis, so it's a bit context specific.

* ``OSO`` (Orthogonal Symmetry-Adapted Orbitals): USOs orthogonalized by
  Symmetric or Canonical Orthogonalization. The number of OSOs may be
  slightly smaller than the total number of USOs, due to removal of linear
  dependencies via Canonical Orthogonalization. The OSOs are rarely
  encountered, as usually we go straight from USOs to MOs.

* ``MO`` (Molecular Orbitals): The combination of OSOs that diagonalizes
  the Fock Matrix, so each basis function is a Hartree-Fock (or Kohn-Sham)
  molecular orbital. The number of OSOs and MOs is always the same. MOs are
  orthonormal.

* ``LO`` (Localized Orbitals): Localized occupied orbitals, a different
  combination of the occupied molecular orbitals which enhances spatial
  locality. LOs do not diagonalize the occ-occ block of the Fock Matrix, but
  remain orthonormal to each other and the virtual space.


.. _`faq:orbdims`:

How to name orbital dimensions
------------------------------

There are a number of different names used to refer to the basis set size.
These may seem redundant, but they have subtly different meanings, as
detailed below.

A calculation can use either pure (5D, 7F, 9G, etc.) basis functions or
Cartesian (6D, 10F, 15G, etc.), as dictated by the input file / basis set
specification. Also, the basis can be represented in terms of atomic
orbitals (AO) or symmetry-adapted orbitals (SO). Further complications
come from the fact that a nearly linearly-dependent basis set will have
functions removed from it to prevent redundancies. With all of these
factors in mind, here are the conventions used internally:

* nao |w---w| The number of atomic orbitals in Cartesian representation.
* nso |w---w| The number of atomic orbitals but in the pure representation if the current basis uses pure functions, number of Cartesian AOs otherwise.
* nbf |w---w| The number of basis functions, which is the same as nso.
* nmo |w---w| The number of basis functions, after projecting out redundancies in the basis.

When molecular symmetry is utilized, a small array of sizes per irrep is
usually allocated on the stack, and is named by augmenting the name above
with a pi (per-irrep), e.g. nmopi. Note that the number of irreps is
always the singular nirrep, and that the index variable h is always used
in a for-loop traverse of irreps.


.. _`faq:orbspaces`:

How to name orbital spaces (e.g., docc)
---------------------------------------

As with basis sets, a number of names are used to refer to refer to the
quantity of electrons, virtuals, and active sub-quantities of a |PSIfour|
calculation. All of these can be defined per irrep as above. Some common
conventions are:

* nelec |w---w| The number of electrons, rarely used due to specialization of alphas and betas or soccs and doccs.
* nalpha |w---w| The number of alpha electrons.
* nbeta |w---w| The number of beta electrons
* docc |w---w| The number of doubly-occupied orbitals
* socc |w---w| The number of singly-occupied orbitals (Almost always alpha, we don't like open-shell singlets much).
* nvir |w---w| The number of virtual orbitals

Multireference Dimensions
^^^^^^^^^^^^^^^^^^^^^^^^^

A orbital diagram of the nomenclature used for CI and MCSCF calculations.

Diagrammatically::

    -----------------------------------------------
           CI       |      RAS      |     CAS
    -----------------------------------------------
                    | frozen_uocc   | frozen_uocc
    dropped_uocc    | rstr_uocc     | rstr_uocc
    -----------------------------------------------
                    | RAS IV        |
                    | RAS III       |
    active          |               | active
                    | RAS II        |
                    | RAS I         |
    -----------------------------------------------
    dropped_docc    | rstr_docc     | rstr_docc
                    | frozen_docc   | frozen_dcc
    -----------------------------------------------

Notation:

* uocc |w---w| Unoccupied orbitals.
* active |w---w| Variable occupation orbitals.
* socc |w---w| Singly occupied orbitals.
* docc |w---w| Doubly occupied orbitals.

Orbital spaces:

* frozen_uocc |w---w| Absolutely frozen virtual orbital.
* rstr_uocc |w---w| Can have rotations, no excitations into.
* dropped_uocc |w---w| rstr_uocc + frozen_uocc

----- end CI active -----

* RAS IV |w---w| uocc, limited number of excitations into.
* RAS III |w---w| uocc, limited number of excitations into.
* RAS II |w---w| docc/socc/uocc, equivalent to active in CAS.
* RAS I |w---w| docc/socc/uocc, limited excitations out of.

----- start CI active -----

* dropped_docc |w---w| rstr_docc + frozen_docc
* rstr_docc |w---w| Can have rotations, no excitations from.
* frozen_docc |w---w| Absolutely frozen core orbital.

Orbitals are sorted by space, irrep, eigenvalue.

