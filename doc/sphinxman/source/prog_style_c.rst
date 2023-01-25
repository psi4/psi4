.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2023 The Psi4 Developers.
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

.. _`sec:style_c`:

C++ Style
=========


.. _`faq:nullptr`:

Prefer ``nullptr`` to ``0`` or ``NULL``
---------------------------------------

``0`` is an ``int`` not a pointer. Almost the same goes for ``NULL``,
though implementations of the language can differ in the details. If you
want to overload on pointer types and/or use pointer types with templates,
use ``nullptr`` to signal the null pointer. The correct overload/template
parameter will then be deduced. Using ``nullptr`` also makes the code more
readable, especially if ``auto`` is used consistently throughout.

*Reference:* Item 8 in `[Effective Modern C++] <https://isbnsearch.org/isbn/9781491903995>`_


.. _`faq:automakeshared`:

Prefer ``std::make_shared`` to direct use of ``new``
----------------------------------------------------

Using ``std::make_shared``:

1. Reduces code verbosity, especially when coupled with ``auto``:

.. code-block:: cpp

    // Type information given 3 TIMES!!!
    std::shared_ptr<Matrix> F = std::shared_ptr<Matrix>(new Matrix("Fock matrix", nso, nso));

    // So much typing...
    std::shared_ptr<Matrix> F = std::make_shared<Matrix>("Fock matrix", nso, nso);

    // Much better!!!!
    auto F = std::make_shared<Matrix>("Fock matrix", nso, nso);

2. Ensures exception safety and prevents resource leaks.

3. Improves efficiency:

.. code-block:: cpp

    // Performs TWO allocations
    std::shared_ptr<Matrix> F = std::shared_ptr<Matrix>(new Matrix("Fock matrix", nso, nso));

    // Performs ONE allocation
    auto F = std::make_shared<Matrix>("Fock matrix", nso, nso);

*Reference:* Item 21 in `[Effective Modern C++] <https://isbnsearch.org/isbn/9781491903995>`_


.. _`faq:autodecl`:

Prefer ``auto`` to explicit type declarations
---------------------------------------------

Using ``auto`` reduces and/or avoids:

1. Verbosity in variable declarations:

.. code-block:: cpp

    std::shared_ptr<Matrix> F = std::make_shared<Matrix>("Fock matrix", nso, nso);  // So much typing...
    auto F = std::make_shared<Matrix>("Fock matrix", nso, nso);  // Much better!

2. Problems with uninitialized variables. auto works like template type
   deduction, hence the right-hand side of the declaration needs to have an
   initializer:

.. code-block:: cpp

    int x1;  // fine, but uninitialized :(
    auto x2;  // WON'T COMPILE!!!
    auto x3 = 1;  // fine and initialized

3. Problems with unintended type casts and type mismatches that are hard
   to impossible to catch:

.. code-block:: cpp

    std::vector<int> v;
    // !!! The size of a vector is of type std::vector<int>::size_type and is compiler- AND architecture-DEPENDENT
    unsigned sz = v.size();  // might not be correct on some compiler/machines
    auto size = v.size();  // size is ALWAYS of the correct type

*Reference:* Items 2 and 5 in `[Effective Modern C++] <https://isbnsearch.org/isbn/9781491903995>`_

Mark virtual functions in derived classes with override
-------------------------------------------------------

The ``override`` keyword introduced in C++11 is used to mark a function in a
derived class and guarantee that it is overloading a function *with the same
signature* in the base class.  This behavior is `checked at compile time
<https://en.cppreference.com/w/cpp/language/override>`_.

.. _`faq:printmem`:

Prefer `GiB` for memory printing
--------------------------------

As memory sizes get larger, we should work in giga (requires decimal printing to not round to zero) rather than mega units.
As it's what we're computing anyways, we should work in 1024-based (mebi, gibi, etc. https://en.wikipedia.org/wiki/Gibibyte) rather than 1000-based units.
As it's a unit, put it in brackets.
Note that users can supply MiB, GB, bytes, or whatever; this guideline is for output printing. ::

        outfile->Printf("  DFHelper Memory: AOs need %.3f [GiB]; user supplied %.3f [GiB]. ",
                        (required *  8 / (1024 * 1024 * 1024.0)),

