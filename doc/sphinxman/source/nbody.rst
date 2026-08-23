.. #
.. # @BEGIN LICENSE
.. #
.. # Psi4: an open-source quantum chemistry software package
.. #
.. # Copyright (c) 2007-2026 The Psi4 Developers.
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

.. index:: BSSE

.. index::
   triple: setting; keywords; cp
   triple: setting; keywords; vmfc
   single: counterpoise correction
   triple: setting; keywords; mbe
   triple: setting; keywords; ssfc
   single: QCManyBody

.. _`sec:cp()`:

Basis Set Superposition Corrections
===================================

.. codeauthor:: Daniel G. A. Smith

.. autofunction:: psi4.driver.driver_nbody.nbody(func, method_string [, molecule, bsse_type, max_nbody, ptype, return_total_data, supersystem_ie_only, external_potentials])


The nbody function computes counterpoise-corrected (CP), non-CP (noCP), and Valiron-Mayer Function Counterpoise (VMFC) interaction energies for complexes composed of arbitrary numbers of monomers.

.. caution:: July 2026, v1.12 QCManyBody is no longer required but may need to be
   separately installed. August 2025, v1.10 many-body computations are no longer
   using internal driver code but have been offloaded to QCManyBody
   software.

External potentials in many-body computations
----------------------------------------------

For a calculation activated by ``bsse_type``, ``external_potentials`` may be a
mapping from 1-based molecular-fragment indices to lists of external potentials.
Each potential list uses the conventional ``[[charge, [x, y, z]], ...]`` form,
with charges and coordinates in atomic units. For example::

    potentials_a = [[q1, [x1, y1, z1]], ...]
    potentials_c = [[q2, [x2, y2, z2]], ...]

    energy(
        'scf',
        bsse_type='cp',
        external_potentials={1: potentials_a, 3: potentials_c},
    )

The mapping associates a potential with the *real* atoms of a fragment. A
component calculation receives the potentials belonging to all of its real
fragments and no others. This rule is the same for ``bsse_type='nocp'``,
``'cp'``, ``'vmfc'``, or a list containing several treatments. Consequently:

* In a noCP component, the basis and real-fragment sets coincide, and only the
  potentials assigned to those fragments are present.
* In a CP component, a fragment represented only by ghost atoms contributes
  basis functions but does **not** activate its external potentials. For a dimer
  with ``external_potentials={2: potentials_b}``, the monomer-A-in-dimer-basis
  calculation has no external potential, whereas the monomer-B-in-dimer-basis
  and dimer calculations include ``potentials_b``.
* VMFC components follow the same real-fragment rule: ghost fragments do not
  activate potentials.

Thus the fragment mapping follows the physical subsystem in each term of the
many-body expansion, rather than applying every potential to every generated
calculation. Keys must be integer fragment indices between 1 and the number of
fragments. Potentials from multiple real fragments are combined for that
component.

Fragment scoping applies to point charges only. Each value must be the flat
``[[charge, [x, y, z]], ...]`` list (equivalently ``[[charge, x, y, z], ...]``
or the corresponding NumPy array); a positional ``[points, diffuse, matrix]``
list, a ``{'points': ..., 'diffuse': ...}`` dictionary, or a potential matrix
raises a :py:exc:`~psi4.driver.p4util.exceptions.ValidationError` rather than
being combined with another fragment's charges. Use one of the whole-molecule
forms below for diffuse charges and potential matrices. A fragment mapping is
also only interpreted when ``bsse_type`` is given; supplying integer keys
without it raises a ``ValidationError``.

Passing a flat list instead of a fragment mapping retains the conventional
single-calculation behavior: the complete list is applied to **every** generated
component, including CP calculations where some molecular fragments are ghosts.
Use it only when the potential represents a fixed background that should be
unchanged across all components. The whole-molecule dictionary spellings
described at :ref:`External potentials and QM/MM <sec:scfqmmm>` -- keys among
``points``, ``diffuse``, and ``matrix``, or the fragment-interaction keys ``A``,
``B``, and ``C`` -- are likewise treated as a fixed background rather than a
fragment mapping. Psi4 emits a warning for any such unscoped form in an MBE
calculation.

When ``embedding_charges`` is also active, the point charges generated for the
fragments outside a component are *added to* the potentials that component
would otherwise receive, rather than replacing them. For the ``A``/``B``/``C``
form the embedding charges join the whole-environment ``C`` partition, leaving
the ``A`` and ``B`` partitions untouched, so the charges enter the total
potential exactly once.

For the general point-charge format and further external-potential forms, see
:ref:`External potentials and QM/MM <sec:scfqmmm>`.

**Examples :** ::

    set {
      basis def2-svp
    }

    # Counterpoise corrected CCSD(T) energies for the Helium dimer
    molecule mol {
      He
      --
      He 1 3
    }
    # Calculate interaction energies only (skips monomers in monomer basis):
    energy('CCSD(T)', bsse_type='cp')
    # Calculate interaction and total energies, return interaction energies:
    energy('CCSD(T)', bsse_type=['cp','nocp'])
    # Calculate and return counterpoise-corrected gradient
    # Useful for e.g. CP-corrected geometry optimization
    gradient('CCSD(T)', bsse_type='cp', return_total_data=True)


    # noCP, VMFC, and CP energy for a helium cluster, limited at 3 bodies
    molecule mol {
      He 0 0 0
      --
      He 0 0 4
      --
      He 0 4 0
      --
      He 4 0 0
    }

    # Returns the nocp energy as its first in the list
    energy('CCSD(T)', bsse_type=['nocp', 'cp', 'vmfc'], max_nbody=3)
    # Calculate cp geometry optimization skipping the MBE intermediate levels
    optimize("ccsd(t)/cc-pv[dt]z", bsse_type="cp", supersystem_ie_only=True)

API
---

.. autoclass:: psi4.driver.driver_nbody.BsseEnum
   :members:
   :undoc-members:

.. autopydantic_model:: psi4.driver.driver_nbody.ManyBodyComputer
   :members:
   :undoc-members:
   :inherited-members: BaseModel, ProtoModel

.. .. autopydantic_model:: qcmanybody.ManyBodyComputer
..    :members:
..    :undoc-members:
..    :noindex:

