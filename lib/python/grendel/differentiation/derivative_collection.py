"""Derivative Collection

"""
import numpy as np

from grendel import sanity_checking_enabled, type_checking_enabled
from grendel.coordinates.cartesian_coordinate import CartesianCoordinate
from grendel.coordinates.coordinate import Coordinate
from grendel.gmath.tensor import Tensor, ComputableTensor
from grendel.gmath.einsum import EinsumTensor
from grendel.representations.representation import Representation
from grendel.util.decorators import typechecked

__author__ = 'dhollman'

# TODO Split this into two classes, a Coordinate-dependent form and a Representation dependent form
class DerivativeCollection(object):
    """A collection of derivatives with indices that are coordinates in a
    given representation.

    """

    ##############
    # Attributes #
    ##############

    tensors = None
    representation = None
    coordinate = None
    tensor_kwargs = None
    first_dimension_different = None
    secondary_representation = None
    einsum_index = None
    uncomputed = None

    ##################
    # Initialization #
    ##################

    @typechecked(
        representation=(Representation, None),
        coordinate=('InternalCoordinate', None),
        first_dimension_different=bool,
        uncomputed=bool,
        secondary_representation=(Representation, None))
    def __init__(self,
            representation=None,
            coordinate=None,
            first_dimension_different=False,
            secondary_representation=None,
            einsum_index=None,
            uncomputed=True,
            **other_tensor_kwargs):
        #--------------------------------------------------------------------------------#
        # miscellanea
        self.tensors = {}
        self.representation = representation
        self.coordinate = coordinate
        self.first_dimension_different = first_dimension_different
        self.secondary_representation = secondary_representation
        self.uncomputed = uncomputed
        #--------------------------------------------------------------------------------#
        # sanity checking
        if sanity_checking_enabled:
            if self.coordinate is None and self.representation is None:
                raise TypeError('DerivativeCollection needs either a Representation or a Coordinate'
                                ' to be based on')
            if self.coordinate is not None and self.representation is not None:
                raise TypeError('DerivativeCollection must be based on either a Representation or a'
                                ' Coordinate, not both.')
            if self.first_dimension_different and self.coordinate is not None:
                raise NotImplementedError('Coordinate-dependent DerivativeCollection cannot have a'
                                ' different first dimension.')
            if self.first_dimension_different and self.secondary_representation is None:
                raise ValueError("a secondary representation must be provided if the first"
                                 " dimension of a DerivativeCollection is to be considered"
                                 " differently.")
            if self.uncomputed and self.coordinate is not None:
                raise NotImplementedError('Coordinate-dependent DerivativeCollection cannot be'
                                          ' uncomputed.')
        #--------------------------------------------------------------------------------#
        # Representation-dependent form
        if self.representation is not None:
            self.tensor_kwargs = other_tensor_kwargs
            self.tensor_kwargs.update(
                representation = self.representation,
                first_dimension_different = self.first_dimension_different,
                secondary_representation = self.secondary_representation,
                uncomputed=self.uncomputed
            )
        #--------------------------------------------------------------------------------#
        # Coordinate-dependent form
        else:
            self.tensor_kwargs = other_tensor_kwargs
            self.einsum_index = einsum_index


    ##############
    # Properties #
    ##############

    @property
    def coordinates(self):
        return self.representation.coords
    coords=coordinates

    ###################
    # Special Methods #
    ###################

    #---------------------#
    # Container Emulation #
    #---------------------#

    def __getitem__(self, item):
        item = (item,) if not isinstance(item, tuple) else item
        #----------------------------------------#
        # Handle einstein summation...
        if all(isinstance(i, basestring) for i in item):
            item = EinsumTensor.split_indices(item)
            #----------------------------------------#
            # Representation dependent form
            if self.representation is not None:
                if self.first_dimension_different:
                    # Act like a Tensor instance
                    return super(ComputableTensor, self.for_order(len(item)-1)).__getitem__(item)
                else:
                    # Act like a Tensor instance
                    return super(ComputableTensor, self.for_order(len(item))).__getitem__(item)
            #----------------------------------------#
            # Coordinate dependent form
            else:
                return self.for_order(len(item)).__getitem__(item)
        #========================================#
        # Representation dependent form
        elif self.representation is not None:
            indices = tuple(i.index if isinstance(i, Coordinate) else i for i in item)
            if self.first_dimension_different:
                # This is a transformation tensor.  The 'order' is one less
                #   than the number of indices.  For instance, the 'first-order'
                #   B tensor has two indices: the internal coordinate and the
                #   cartesian coordinate.
                return self.for_order(len(indices)-1).__getitem__(indices)
            else:
                # This is a different type of derivative tensor, such as a
                #   force tensor.  The order is the same as the number of indices.
                #   For instance, the 'first-order' energy derivative is a vector, with
                #   one index.
                return self.for_order(len(indices)).__getitem__(indices)
        #----------------------------------------#
        # Coordinate dependent form
        else: # self.coordinate is not None
            if sanity_checking_enabled:
                if any(isinstance(c, CartesianCoordinate) for c in item)\
                        and not all(isinstance(c, CartesianCoordinate) for c in item):
                    raise ValueError("mixing of CartesianCoordinates and integers in b tensor"
                                     " element retrieval is confusing and thus no longer allowed."
                                     " indices were ('{}')".format(
                        "', '".join(str(c) for c in item)
                    ))
            if all(isinstance(c, CartesianCoordinate) for c in item):
                item = self.coordinate.internal_indices_for_coordinates(*item)
            return self.for_order(len(item)).__getitem__(item)

    def __setitem__(self, item, value):
        item = (item,) if not isinstance(item, tuple) else item
        #----------------------------------------#
        # Handle odd cases
        if sanity_checking_enabled and item is Ellipsis:
            raise NotImplementedError
        #----------------------------------------#
        # Handle einstein summation...
        if all(isinstance(i, basestring) for i in item):
            item = EinsumTensor.split_indices(item)
            #----------------------------------------#
            # Representation dependent form
            if self.representation is not None:
                if self.first_dimension_different:
                    # Act like a Tensor instance (skip nonsense in Derivative tensor)
                    return super(ComputableTensor, self.for_order(len(item)-1)).__setitem__(item, value)
                else:
                    # Act like a Tensor instance
                    return super(ComputableTensor, self.for_order(len(item))).__setitem__(item, value)
            #----------------------------------------#
            # Coordinate dependent form
            else:
                return self.for_order(len(item)).__setitem__(item, value)
        #----------------------------------------#
        elif self.representation is not None:
            # Representation-dependent form
            # Handle the case where the user is setting things in chunks:
            # Try to figure out what they are trying to do.
            if item == (Ellipsis,) and not self.first_dimension_different:
                try:
                    shp = value.shape
                except AttributeError:
                    raise ValueError(
                        "Setting a DerivativeCollection item by giving an Ellipsis is ambiguous,"
                        " especially if the value doesn't have a 'shape' attribute.  Please set"
                        " items individually or use the 'for_order()' instance method to retrieve"
                        " the individual Tensor to set.")
                tens = self.for_order(len(shp))
                tens[...] = value
                return
            elif self.first_dimension_different and len(item) == 1:
                try:
                    shp = value.shape
                except AttributeError:
                    raise ValueError(
                        "Setting a DerivativeCollection item by giving just one key is ambiguous,"
                        " especially if the value doesn't have a 'shape' attribute.  Please set"
                        " items individually or use the 'for_order()' instance method to retrieve"
                        " the individual Tensor to set.")
                order = len(shp)
                self.for_order(order).__setitem__(item, value)
                return
            elif not np.isscalar(value):
                raise NotImplementedError(
                    "Setting non-scalar chunks of a DerivativeCollection is mostly not implemented;"
                    " the exception is the case of a transformation tensor (a DerivativeCollection"
                    " where 'first_dimension_different' is set to True) where indices for the first"
                    " dimension may be specified along with corresponding chunks.  In other cases,"
                    " it is too difficult to understand to which order Tensor the assignment is being"
                    " made (this may be allowed in the future, but not now).  If you need to do this,"
                    " use the 'for_order()' instance method and explicitly specify the order the"
                    " Tensor to which you wish to make the assignment."
                )
            #----------------------------------------#
            indices = [i.index if isinstance(i, Coordinate) else i for i in item]
            # See notes in __getitem__
            if self.first_dimension_different:
                self.for_order(len(indices)-1).__setitem__(tuple(indices), value)
                return
            else:
                self.for_order(len(indices)).__setitem__(tuple(indices), value)
                return
        else: # self.coordinate is not None
            # Coordinate-dependent form
            if item == (Ellipsis,):
                try:
                    shp = value.shape
                except AttributeError:
                    raise ValueError(
                        "Setting a DerivativeCollection item by giving an Ellipsis is ambiguous,"
                        " especially if the value doesn't have a 'shape' attribute (which is the"
                        " case for the object passed in of type `{ty}`)  Please set items"
                        " individually or use the 'for_order()' instance method to retrieve"
                        " the individual Tensor to set.".format(
                            ty=type(item).__name__
                        ))
                tens = self.for_order(len(shp))
                tens[...] = value
                return
            elif not np.isscalar(value):
                raise NotImplementedError(
                    "Setting non-scalar chunks of a DerivativeCollection is mostly not implemented;"
                    " the exception is the case of a transformation tensor (a DerivativeCollection"
                    " where 'first_dimension_different' is set to True) where indices for the first"
                    " dimension may be specified along with corresponding chunks.  In other cases,"
                    " it is too difficult to understand to which order Tensor the assignment is being"
                    " made (this may be allowed in the future, but not now).  If you need to do this,"
                    " use the 'for_order()' instance method and explicitly specify the order the"
                    " Tensor to which you wish to make the assignment."
                )
            else:
                indices = self.coordinate.internal_indices_for_coordinates(*item)
                self.for_order(len(indices)).__setitem__(indices, value)

    # TODO __str__, __repr__

    ###########
    # Methods #
    ###########

    @typechecked(order=int)
    def for_order(self, order):
        """ Return the Tensor object corresponding to the `order`th order derivative.
        """
        if sanity_checking_enabled:
            if order < 0:
                raise ValueError("don't know about derivatives of order {0}".format(order))
        #--------------------------------------------------------------------------------#
        if order in self.tensors:
            return self.tensors[order]
        else:
            #----------------------------------------#
            # Representation dependent form
            if self.representation is not None:
                if self.first_dimension_different:
                    shape = (len(self.representation),) + (len(self.secondary_representation),) * order
                else:
                    shape = (len(self.representation),) * order
                self.tensors[order] = DerivativeTensor(
                    collection=self,
                    shape=shape,
                    **self.tensor_kwargs)
                return self.tensors[order]
            #----------------------------------------#
            # Coordinate-dependent form
            else: # self.coordinate is not None
                shape = (3*len(self.coordinate.atoms),) * order
                if self.einsum_index is not None:
                    self.tensors[order] = Tensor(
                        indices=','.join((self.einsum_index,) * order),
                        **self.tensor_kwargs)
                else:
                    self.tensors[order] = Tensor(
                        shape=shape,
                        **self.tensor_kwargs)

                return self.tensors[order]

    #-----------------#
    # Inquiry methods #
    #-----------------#

    def has_order(self, order):
        return order in self.tensors


#####################
# Dependent Imports #
#####################

from grendel.differentiation.derivative_tensor import DerivativeTensor
# imports needed for dynamic typechecking, do not delete
if type_checking_enabled:
    from grendel.coordinates.internal_coordinate import InternalCoordinate

