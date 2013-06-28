"""

"""
from collections import Iterable, Sequence
from itertools import product
from grendel.util.iteration import grouper

#TODO better float formatting based on overall magnitude of the tensor
class MatrixFormatter(object):
    """
    """

    ##############
    # Attributes #
    ##############

    float_format_leader = " "
    float_width = 10
    float_type_string = 'g'
    float_digits = 4
    width_overflow_character = '*'
    label_format_leader = "^"
    label_width = 10
    label_overflow_character = 'X'
    col_label_format_leader = None
    col_label_width = None
    row_label_format_leader = None
    row_label_width = None
    line_width = 120
    one_based = True
    name = None

    ##################
    # Initialization #
    ##################

    def __init__(self, **kwargs):
        for key, val in kwargs.items():
            setattr(self, key, val)

    ##############
    # Properties #
    ##############

    @property
    def col_label_format(self):
        return (self.col_label_format_leader or self.label_format_leader)\
                    + str(self.float_width-1)
                   # + str(self.col_label_width or self.label_width)

    @property
    def row_label_format(self):
        return (self.row_label_format_leader or self.label_format_leader)\
                    + str(self.row_label_width or self.label_width)

    @property
    def float_format(self):
        rv = self.float_format_leader
        rv += str(self.float_width-1)
        if self.float_digits is not None:
            rv += "." + str(self.float_digits)
        rv += self.float_type_string
        return rv

    #################
    # Class Methods #
    #################

    @classmethod
    def quick_format(cls, matrix, labels=None, **kwargs):
        formatter = cls(**kwargs)
        return formatter.format(matrix, labels)


    ###########
    # Methods #
    ###########

    def format(self, matrix, labels=None):
        if not isinstance(matrix, Matrix):
            matrix = matrix.view(Matrix)
        #--------------------------------------------------------------------------------#
        ret_val = ''
        #--------------------------------------------------------------------------------#
        if self.name is not None:
            ret_val += '{} x {} Matrix "{}":\n'.format(matrix.shape[0], matrix.shape[1], self.name)
        #--------------------------------------------------------------------------------#
        rwidth = self.row_label_width or self.label_width
        row_fmt = self.row_label_format
        cwidth = self.col_label_width or self.label_width
        col_fmt = self.col_label_format
        # rlf = "row label format"
        def rlf(lbl):
            rv = ("{:" + row_fmt + "}").format(lbl)
            if len(rv) > rwidth:
                return self.label_overflow_character * rwidth
            return rv
        # clf = "column label format"
        def clf(lbl):
            rv = ("{:" + col_fmt + "}").format(lbl)
            if len(rv) > cwidth:
                return self.label_overflow_character * cwidth
            return rv
        # ff = "float format"
        def ff(num):
            rv = ("{:" + self.float_format + "}").format(num)
            if len(rv) > self.float_width:
                return rv[:self.float_width-2] + self.width_overflow_character * 2
            return rv
        #--------------------------------------------------------------------------------#
        shp = matrix.shape
        if labels is None:
            if self.one_based:
                labels = [map(str, range(1, dim+1)) for dim in shp]
            else:
                labels = [map(str, range(dim)) for dim in shp]
        elif isinstance(labels, Sequence) and not isinstance(labels[0], (list, tuple)):
            labels = (tuple(labels),) * len(shp)
        #--------------------------------------------------------------------------------#
        nperline = (self.line_width - rwidth) // max(cwidth+1, self.float_width+1)
        for idx in xrange(matrix.ncols//nperline + 1):
            lbls = labels[1][idx*nperline:(idx+1)*nperline]
            ret_val += rlf("") + " " + " ".join(map(clf, lbls)) + '\n'
            for rownum in xrange(shp[0]):
                data = matrix[rownum][idx*nperline:(idx+1)*nperline]
                ret_val += rlf(labels[0][rownum]) + " " + " ".join(map(ff, data)) + "\n"
            ret_val += "\n"
        #--------------------------------------------------------------------------------#
        return ret_val



class TensorFormatter(MatrixFormatter):
    """
    """

    ##############
    # Attributes #
    ##############

    symmetric = False

    ##################
    # Initialization #
    ##################

    def __init__(self, **kwargs):
        super(TensorFormatter, self).__init__(**kwargs)

    #################
    # Class Methods #
    #################

    @classmethod
    def format(cls, tensor, labels=None, **kwargs):
        formatter = cls(**kwargs)
        return formatter.format(tensor, labels)

    ###########
    # Methods #
    ###########

    #noinspection PyMethodOverriding
    def format(self, tensor, labels=None):
        ret_val = ''
        #--------------------------------------------------------------------------------#
        shp = tensor.shape
        if labels is None:
            if self.one_based:
                labels = [map(str, range(1, dim+1)) for dim in shp]
            else:
                labels = [map(str, range(dim)) for dim in shp]
        elif isinstance(labels, Sequence) and not isinstance(labels[0], (list, tuple)):
            labels = (tuple(labels),) * len(shp)
        #--------------------------------------------------------------------------------#
        for idxs in product(*map(range, shp[:-2])):
            if self.symmetric and idxs != sorted(idxs):
                continue
            mtx = Matrix(tensor[idxs])
            ret_val += ", ".join(labels[i][n] for i, n in enumerate(idxs)) + ":\n"
            ret_val += super(TensorFormatter, self).format(mtx, labels[-2:])
            ret_val += "\n"
        #--------------------------------------------------------------------------------#
        return ret_val



#####################
# Dependent Imports #
#####################

from grendel.gmath.matrix import Matrix








