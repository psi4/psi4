.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2024 The Psi4 Developers.
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

.. _`sec:prog_ways_to_add`:

======================================================
Ways to Add Code: Psi4NumPy, Plugins, Full Integration 
======================================================

Easier and more rapid development
---------------------------------

Fully-featured electronic structure programs are large and complex.  However,
the |PSIfour| ecosystem provides a path for easier and more rapid development
of new features.  The earliest versions of |PSIfour| were written by merging
individual executables that performed specific tasks into a unified C++
executable.  By linking this C++ executable against the Python interpreter, the
individual modules could be called in any sequence, enabling a very diverse
range of tasks to be accomplished within a given input.  Although Python-driven
model allowed for great flexibility from a user's perspective, programming was
still relatively difficult because it required modifications to be made in C++
code.

Since those early days, the code has undergone some important structural
changes that have greatly simplified the development workflow.  These changes
were motivated by the realization that only a few bottlenecks exist in a typical
calculation; by focusing on optimized C++ implementations of these bottlenecks
and making these C++ functions available in Python, most of the code to implement
the overall calculation can be written in simpler Python code.  Python is far
better suited to management tasks such as directory navigation and retrieval,
making it a natural choice for overall calculation layout than C++.  With the
emergence of `NumPy <https://numpy.org/>`_ as a standard tool for executing almost any
mathematical technique efficiently in Python, the transitioning of code from
C++ to Python has facilitated a much simpler work flow for prototyping and
developing methods: this is detailed in the next section.

.. _`sec:prog_psi4numpy`:

Rapid initial development using Psi4NumPy
-----------------------------------------

The `Psi4NumPy <https://github.com/psi4/psi4numpy>`_ project [Smith:2018:3504]_ is the recommended
mechanism for developing and prototyping new methods in Psi4.  Because
`NumPy <https://numpy.org/>`_ provides such a rich set of features for efficient linear
algebra, Fourier transforms, and general tensor manipulations, a massive number
of methods can be easily implemented very easily using that library.  To
facilitate this workflow, |PSIfour| exports key quantities such as integrals,
densities and molecular orbitals in NumPy format.  From this point, the
programmer can simply call the appropriate |PSIfour| functions to compute the
desired input quantities, retrieve them in NumPy format, and then write the
remaining code using standard Python and/or NumPy syntax.  This approach does
not require any recompilation of code, resulting in a particularly facile
development workflow.  Detailed examples and tutorials are available in the
`Psi4NumPy <https://github.com/psi4/psi4numpy>`_ repository.

.. _`sec:prog_plugins`:

Avoiding the need to modify Psi4, using plugins
-----------------------------------------------

In the early days when |PSIfour| was still primarily a C++ code, development
was very cumbersome due to a lengthy build process.  To expedite development, a
plugin system was developed.  This plugin machinery allows developers to access
the classes defined in the innards of |PSIfour|, with only the small plugin
code requiring recompilation during development.  The resulting lightweight
code can be maintained and distributed independently of |PSIfour|, making this
a good strategy for development, especially in cases where tighter integration
of the new code with existing |PSIfour| machinery is required than that
afforded by the Numpy based strategy outlined in the :ref:`sec:prog_psi4numpy`
section.  For details about how to write these plugins, see the
:ref:`sec:plugins` section.

.. _`sec:prog_fullintegration`:

Incorporating code into |PSIfour|
---------------------------------

For features to be incorporated fully into the |PSIfour| ecosystem, changes to
the core routines are inevitable.  However, the programmer should think very
carefully about the most appropriate language for the task in hand.  Let's
consider a new feature that downloads some data from an external source and
then performs some kind of expensive matrix operation on those data.  Because
Python has a rich set of tools for obtaining data from external sources,
writing this tool in the Python layer is a natural choice.  If we know that the
matrix will always be small enough to fit in memory, we can simply rely on the
routines present in NumPy to do the heavy lifting and the code is easy to
implement entirely in the Python layer.  In the case where the matrix operation
is non-standard and requires some specialized code to handle disk-based
storage, the decision to write in Python is less clear cut.  It is certainly
possible to write these out-of-core routines using Numpy primitives, but there
are a number of tools in |PSIfour| already to perform tasks like these that are
required, *e.g.*, for cluster.  In this case, a good design would be to write a
simple piece of code in the C++ layer that performs the matrix operation on a
given input, using the I/O routines available in |PSIfour| and the parallelism
afforded by OpenMP, and to make that code available to the front end as
described in :ref:`sec:prog_tour-exposing`.  The Python layer could then be
responsible for obtaining the input data and calling this C++ code to do the
manipulations, allowing each language layer to handle the subset of the work
that caters to their individual strengths.

A number of concrete examples of this workflow exist in the code already.  For
finite difference computations of energy derivatives, the logic to determine
the type of stencil and which displacements are needed is not going to be rate
limiting for any reasonable quantum mechanical energy function.  Therefore,
doing that work in the Python layer is a good idea, as it allows the many
Python tools for farming out *embarrassingly parallel* workloads to be used,
while the C++ layer can be used to implement the energy function to be
differentiated.

In SCF, we have a number of sources of external embedding potentials that could
enter the calculation.  Allowing Python to handle only the details of driving
the SCF iterations, such as external potentials and convergence acceleration
methods, but deferring to C++ to do the heavy lifting for building and
diagonalizing the Fock matrix also takes advantage of the two languages'
strengths and improves maintainability of the code.
