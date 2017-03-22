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

.. include:: autodoc_abbr_options_c.rst

.. index:: 
   NumPy

.. _`sec:numpy`:

Interface to NumPy
==================

.. codeauthor:: Daniel G. A. Smith
.. sectionauthor:: Daniel G. A. Smith

*Module:* :ref:`Keywords <apdx:numpy>`, :source:`psi4/driver/p4util/numpy_helper.py`

Basics
~~~~~~

Converting between the |PSIfour| Data classes and a NumPy array is easy through
various helper functions as detailed in this section.  A quick overview NumPy
functionality can be found `here
<https://docs.scipy.org/doc/numpy-dev/user/quickstart.html>`_.  In addition,
numerous example of hybrid NumPy and Psi4 can be found at the `Psi4Numpy
project <https://github.com/dgasmith/psi4numpy`_.  Currently only the Matrix
and Vector objects support NumPy interfacing. Let us begin with a simple
conversion from these objects to a NumPy array::

    >>> import psi4
    >>> import numpy as np

    # Build the Psi4 data objects
    >>> mat = psi4.Matrix(3, 3) 
    >>> vec = psi4.Vector(3)

    # Convert to a NumPy array
    >>> numpy_mat = np.array(mat)
    >>> numpy_vec = np.array(vec)

Here the data is copied into new NumPy arrays. NumPy arrays can be converted
back to |PSIfour| objects using the ``from_array`` interface::

    >>> new_mat = psi4.Matrix.from_array(mat)
    >>> new_vec = psi4.Vector.from_array(vec)


NumPy Views
~~~~~~~~~~~

Copying the data between NumPy and Psi4 objects can lead to excessive data
movement and convoluted code. Here we introduce the idea of "Views" where the
same data can be viewed by multiple objects. However, this can lead to very
subtle errors if used incorrectly and care needs to be taken when using these
views.  Views can be created in two ways::

    >>> numpy_mat_view = np.asarray(mat)

    # Access the NumPy object and set all values to 1 through broadcasting
    >>> numpy_mat_view[:] = 1
    
    >>> print(np.array(mat))
    [[ 1.  1.  1.]
     [ 1.  1.  1.]
     [ 1.  1.  1.]]


Secondly, these objects have a ``.np`` attribute for easy access to the underlying data::

    >>> mat.np[:] = 1

this operation is identical to the above.

.. warning:: The following will lead to reference errors: ``view =
   psi4.Matrix(3, 3).np``. Here, the Python garbage collection deletes the Matrix
   object, the view then points to deleted data resulting in the view effectively
   reading random data. As a general rule, never assign the ``.nph`` or ``.np``
   accessors.


|PSIfour| Data Objects with Irreps
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

|PSIfour| data objects natively support multiple irreducible representations
which is quite useful for Quantum Chemistry. However, this is not fundamental
to NumPy and some work around are required to natively support these
operations. Take the following irreped Matrix::

    >>> dim = psi4.Dimension.from_list([1, 2, 3])
    >>> irreped_mat = psi4.Matrix("New Matrix", dim, dim)

    # Create a list of Psi4 arrays
    >>> list_of_arrays = irreped_mat.to_array()

    # Or, use the .nph irreped accessor
    >>> irreped_mat.nph[0][:] = 1

Where ``.nph`` is the irreped accessor form. If ``.np`` or ``np.array`` are
called on irreped Matrices or Vectors an error will be thrown; however, the
irreped form is always valid for non-irreped matrices.

Array to Matrix
~~~~~~~~~~~~~~~
A general function that converts |PSIfour| data objects to NumPy arrays.

.. autofunction:: psi4.driver.p4util.numpy_helper.array_to_matrix

Matrix to Array
~~~~~~~~~~~~~~~
A general function that converts NumPy arrays to |PSIfour| data objects.

.. autofunction:: psi4.driver.p4util.numpy_helper._to_array



