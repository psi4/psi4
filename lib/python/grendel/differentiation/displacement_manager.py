from __future__ import print_function
from collections import Iterable
from copy import copy
import sys
from grendel import type_checking_enabled, sanity_checking_enabled
from grendel.chemistry.molecule import Molecule
from grendel.chemistry.derivative_properties import RepresentationDependentProperty
from grendel.chemistry.molecular_properties import MolecularProperty
from grendel.coordinates.coordinate import Coordinate
from grendel.differentiation.finite_difference import FiniteDifferenceFunction, Differentiable, FiniteDifferenceVariable
from grendel.interface.computation_details import ComputationDetails
from grendel.interface.legacy_xml import LegacyXMLResultGetter
from grendel.representations.representation import Representation
from grendel.util.decorators import typechecked, SequenceOf
from grendel.util.metaprogramming import ReadOnlyAttribute
from grendel.util.strings import indented
from grendel.util.units import hasunits, compatible_units, IncompatibleUnitsError
from grendel.util.units.unit import DistanceUnit

__author__ = 'dhollman'

class DisplacementManager(FiniteDifferenceFunction):
    """ A function (in the `finite_difference.FiniteDifferenceFunction` sense) that can get values for displacements of
    some base molecule.
    """

    ##############
    # Attributes #
    ##############

    representation = ReadOnlyAttribute('representation',
        doc=""" The representation in which to carry out the displacements."""
    )
    differentiable = ReadOnlyAttribute('differentiable',
        doc=""" The subclass of Differentiable that `value_for_displacements` should return an instance of."""
    )

    variables = ReadOnlyAttribute('variables',
        doc="""The list of `FiniteDifferenceVariable` instances on which the function depends."""
    )

    output_type = ReadOnlyAttribute('output_type',
        doc="""The subclass of `Differentiable` that the function generates."""
    )

    deltas = None

    displacements = None

    details = None

    ##################
    # Initialization #
    ##################

    @typechecked(
        base_molecule_or_representation=('Molecule', 'Representation'),
        differentiable=type,
        deltas=(SequenceOf(float), None),
        details=(ComputationDetails, None)
    )
    def __init__(self, base_molecule_or_representation, differentiable, deltas=None, details=None):
        if type_checking_enabled:
            if not isinstance(base_molecule_or_representation, (Molecule, Representation)):
                raise TypeError
            if not isinstance(differentiable, type) and issubclass(differentiable, Differentiable):
                raise TypeError
        if sanity_checking_enabled:
            if not issubclass(differentiable, MolecularProperty):
                # TODO write this error
                raise ValueError
        self._differentiable = differentiable
        if isinstance(base_molecule_or_representation, Molecule):
            # Assume the first InternalRepresentation is the one we want
            self._representation = base_molecule_or_representation.internal_representation
        else:
            self._representation = base_molecule_or_representation
        self.displacements = {}
        self.deltas = deltas
        self.details = details
        if deltas is not None:
            if type_checking_enabled and not isinstance(deltas, Iterable):
                raise TypeError
            if sanity_checking_enabled and len(deltas) != len(self.representation):
                raise ValueError("dimension mismatch.  'deltas' list (length {}) must be the same " \
                                  "as the number of coordinates ({})".format(len(deltas), len(self.representation)))
            self.deltas = []
            for delta, coord in zip(deltas, self.representation):
                if hasunits(delta):
                    if sanity_checking_enabled and not compatible_units(delta.units, coord.units):
                        raise IncompatibleUnitsError("{} is not in valid units for a displacement of a {}".format(
                            delta, coord.__class__.__name__
                        ))
                else:
                    delta = delta * coord.units
                self.deltas.append(delta)
        else:
            self.deltas = Displacement.get_default_deltas(self.representation)

    ##############
    # Properties #
    ##############

    @property
    def base_molecule(self):
        return self.representation.molecule

    @property
    def displaced_molecules(self):
        return [d.displaced_molecule for d in self.displacements.values]

    ###################
    # Special Methods #
    ###################

    def __contains__(self, item):
        if isinstance(item, tuple):
            return item in self.displacements
        elif isinstance(item, Displacement):
            return item.increments(self.deltas) in self.displacements
        else:
            raise TypeError

    ###########
    # Methods #
    ###########

    def displacement_for(self, increments):
        if increments not in self.displacements:
            self.displacements[increments] = Displacement.from_increments(
                increments,
                self.representation,
                self.deltas)
        return self.displacements[increments]

    #----------------------------------------------#
    # Methods abstract in FiniteDifferenceFunction #
    #----------------------------------------------#

    def value_for_displacements(self, pairs):
        # Type checking and sanity checking
        if type_checking_enabled:
            if not all(len(pair) == 2 and isinstance(pair[0], FiniteDifferenceVariable) and isinstance(pair[1], int) for pair in pairs):
                raise TypeError
        if sanity_checking_enabled:
            for pair in pairs:
                if not isinstance(pair[0], Coordinate):
                    raise ValueError(
                        "FiniteDifferenceVariable instances passed to"
                        " DisplacementManager.value_for_displacements"
                        " must be instances of Coordinate subclasses."
                        "  (Got at least one that was a {0})".format(
                            type(pair[0]).__name__
                        ))
        #--------------------------------------------------------------------------------#
        # prepare the increments
        increments = [0]*len(self.representation)
        for coord, ndeltas in pairs:
            i = self.representation.coords.index(coord)
            increments[i] = ndeltas
        increments = tuple(increments)
        #--------------------------------------------------------------------------------#
        # get the displacement, and ask the displaced molecule for the property
        try:
            displacement = self.displacement_for(increments)
            rv = displacement.displaced_molecule.get_property(
                self.differentiable,
                details=self.details
            )
            if rv is None:
                raise ValueError("couldn't get displacement for increments {}".format(increments))
            return rv
        #--------------------------------------------------------------------------------#
        # Old debugging code.  Only applies to LegacyXMLResultGetter
        except RuntimeError as e:
            if all(isinstance(rg, LegacyXMLResultGetter)
                    for rg in self.displacements[increments].displaced_molecule.result_getters):
                legacy_getters = [rg
                                     for rg in self.base_molecule.result_getters
                                         if isinstance(rg, LegacyXMLResultGetter)]
                molecule = self.displacements[increments].displaced_molecule
                smallest_dist = float('inf') * molecule.cartesian_units
                smallest_dist_stub = None
                best_rg = None
                mol = copy(molecule); mol.convert_units(DistanceUnit.default)
                for rg in legacy_getters:
                    for stub in rg.properties.keys():
                        st = copy(stub); st.convert_units(DistanceUnit.default)
                        ldiff = rg.properties.key_representation.value_for_molecule(stub)
                        ldiff -= rg.properties.key_representation.value_for_molecule(mol)
                        ldiff = abs(ldiff)
                        if max(ldiff) < smallest_dist:
                            best_rg = rg
                            smallest_dist = max(ldiff)
                            smallest_dist_stub = st
                if smallest_dist_stub is not None:
                    # Print debugging information for legacy xml parser
                    print("Failed to find match.  Smallest " \
                          "maximum geometric difference was {}".format(smallest_dist),
                        file=sys.stderr)
                    print(indented(
                            "Desired molecule representation " \
                            "was:\n{}".format(indented(str(molecule.internal_representation.copy_with_molecule(mol))))
                        ), file=sys.stderr)
                    print(indented(
                        '-'*40 +
                        "\nClosest stub in the same representation:\n{}".format(
                            indented(str(molecule.internal_representation.copy_with_molecule(smallest_dist_stub))))
                    ), file=sys.stderr)
                    print("Tree leaf: \n{}".format(
                        best_rg.properties._get_possibilities(mol)
                    ), file=sys.stderr)
                    print("Chains: \nMolecule:\n{}\nStub:\n{}".format(
                        indented(str(best_rg.properties._key_chain_str(mol))),
                        indented(str(best_rg.properties._key_chain_str(smallest_dist_stub))),
                    ), file=sys.stderr)
                    sys.stderr.flush()
                raise RuntimeError("couldn't get {} for displacements {}".format(
                    self.differentiable.__name__,
                    increments
                ))
            else:
                raise

    def deltas_for_variables(self, vars):
        return [self.deltas[i] for i in range(len(self.representation)) if self.representation[i] in vars]


#####################
# Dependent Imports #
#####################

from grendel.differentiation.displacement import Displacement

