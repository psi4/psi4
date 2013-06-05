"""
"""
from itertools import combinations_with_replacement as symmetric_product, permutations

from grendel.chemistry.derivative_properties import PropertyDerivative
from grendel.chemistry.molecular_properties import Energy
from grendel.differentiation.derivative_collection import DerivativeCollection
from grendel.differentiation.finite_difference import FiniteDifferenceDerivative
from grendel.util.decorators import typechecked

class ForceField(DerivativeCollection):
    """
    """

    ##############
    # Attributes #
    ##############

    max_order = None
    property_units = None
    molecular_property = None

    ##################
    # Initialization #
    ##################

    def __init__(self, representation, max_order, property=Energy, property_units=None):
        self.property_units = property_units or property.default_units
        self.molecular_property = property
        super(ForceField, self).__init__(representation=representation)
        self.max_order = max_order
        self.filled_for_order = [False] * self.max_order
        ncoord = len(self.representation)
        #TODO units
        for order in xrange(1, self.max_order+1):
            self.tensors[order] = DerivativeTensor(
                representation=self.representation,
                shape=(ncoord,) * order
            )

    ##############
    # Properties #
    ##############

    @property
    def base_molecule(self):
        return self.representation.molecule

    ###########
    # Methods #
    ###########

    @typechecked(order=int)
    def for_order(self, order):
        if order > self.max_order:
            raise ValueError("Requested tensor order ({}) for ForceField is greater"
                             " than max_order ({})".format(
                order, self.max_order
            ))
        else:
            return super(ForceField, self).for_order(order)

    def in_representation(self, representation):
        return self.representation.transform_forcefield(self, representation)

class FiniteDifferenceForceField(ForceField):
    """
    """

    ##############
    # Attributes #
    ##############

    displacement_manager = None
    filled_for_order = None
    computed_property = None
    max_analytic_order = None
    robustness = None
    computation_queue = None

    ######################
    # Private Attributes #
    ######################

    _unit_dict = None

    ##################
    # Initialization #
    ##################

    def __init__(self,
            max_analytic_order=0,
            property=Energy,
            deltas=None,
            details=None,
            run_now=False,
            queue=None,
            robustness=2,
            **kwargs):
        kwargs.update(property=property)
        super(FiniteDifferenceForceField, self).__init__(**kwargs)
        self.robustness = robustness
        self.max_analytic_order = max_analytic_order
        self.filled_for_order = [False] * self.max_order
        if max_analytic_order > 1:
            raise NotImplementedError
        elif max_analytic_order > 0:
            if isinstance(self.representation, CartesianRepresentation):
                self.computed_property = PropertyDerivative(
                    property,
                    representation=self.representation,
                    order=max_analytic_order
                )
                self.computed_property = self.computed_property.in_units(
                    self.property_units / self.representation.units**self.max_analytic_order
                )
            else:
                # TODO figure out why I need to do finite difference of cartesian gradients and then convert to get the correct answer (rather than finite difference of internal coordinate gradients)
                self.computed_property = PropertyDerivative(
                    property,
                    representation=self.representation.molecule.cartesian_representation,
                    order=max_analytic_order
                )
                self.computed_property = self.computed_property.in_units(
                    self.property_units / (
                        self.representation.molecule.cartesian_representation.units**self.max_analytic_order
                    )
                )
        else:
            self.computed_property = property.in_units(self.property_units)
        self.displacement_manager = DisplacementManager(
            self.representation,
            differentiable=self.computed_property,
            deltas=deltas,
            details=details
        )
        if queue is not None:
            self.computation_queue = queue
            queue.enqueue(*self.needed_computations)
            if run_now:
                queue.run()
                self.fill()
        elif run_now:
            self.fill()

    ##############
    # Properties #
    ##############

    @property
    def all_computations(self):
        all_comps = set()
        for order in xrange(self.max_analytic_order+1, self.max_order+1):
            for coords in symmetric_product(self.representation, order-self.max_analytic_order):
                fd = FiniteDifferenceDerivative(
                    self.displacement_manager,
                    *coords,
                    target_robustness=self.robustness
                )
                for incs in fd.needed_increments:
                    increments = [0] * len(self.representation)
                    for coord, inc in zip(fd.variables, incs):
                        increments[coord.index] = inc
                    dmol = self.displacement_manager.displacement_for(tuple(increments)).displaced_molecule
                    comp = dmol.get_computation_for_property(
                        self.computed_property,
                        self.displacement_manager.details
                    )
                    if comp not in all_comps:
                        all_comps.add(comp)
        if self.max_analytic_order == 1:
            comp = self.base_molecule.get_computation_for_property(
                self.computed_property,
                self.displacement_manager.details
            )
            if comp not in all_comps:
                all_comps.add(comp)
        elif self.max_analytic_order > 2:
            # TODO figure out how many analytic computations need to be done
            raise NotImplementedError
        return list(all_comps)

    @property
    def needed_computations(self):
        return [c for c in self.all_computations if not c.completed]


    ###########
    # Methods #
    ###########

    def fill(self):
        if self.computation_queue is not None and not self.computation_queue.is_finished():
            self.computation_queue.run()
        for n in xrange(self.max_order):
            self.fill_order(self.max_order - n)

    def fill_order(self, order):
        if self.filled_for_order[order-1]:
            return
        #----------------------------------------#
        tens = self.tensors[order]
        if order <= self.max_analytic_order:
            #TODO this will need to be significantly modified for hessians and higher, since the transformation to interal coordinates requires all lower derivatives as well...
            tens[...] = self.base_molecule.get_property(
                PropertyDerivative(self.molecular_property, order) if order > 0 else self.molecular_property,
                details=self.displacement_manager.details
            ).value.in_representation(self.representation)
        #----------------------------------------#
        else:
            for f_coords in symmetric_product(self.representation, order-self.max_analytic_order):
                if self.max_analytic_order == 0 or isinstance(self.representation, CartesianRepresentation):
                    spread_val = FiniteDifferenceDerivative(
                        self.displacement_manager,
                        *f_coords,
                        target_robustness=self.robustness
                    ).value
                else:
                    spread_val = RepresentationDependentTensor(
                        FiniteDifferenceDerivative(
                            self.displacement_manager,
                            *f_coords,
                            target_robustness=self.robustness
                        ).value,
                        representation=self.representation.molecule.cartesian_representation
                    )
                    spread_val = spread_val.in_representation(self.representation)
                for perm in permutations(f_coords):
                    tens[perm] = spread_val
        #----------------------------------------#
        self.filled_for_order[order-1] = True

    @typechecked(order=int)
    def for_order(self, order):
        if order <= self.max_order and not self.filled_for_order[order-1]:
            raise ValueError("Order {} tensor of ForceField instance not yet computed.".format(
                order
            ))
        else:
            return super(FiniteDifferenceForceField, self).for_order(order)

#####################
# Dependent Imports #
#####################

from grendel.differentiation.derivative_tensor import DerivativeTensor, RepresentationDependentTensor
from grendel.differentiation.displacement_manager import DisplacementManager
from grendel.representations.internal_representation import InternalRepresentation
from grendel.representations.cartesian_representation import CartesianRepresentation

