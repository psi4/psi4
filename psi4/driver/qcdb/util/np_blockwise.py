#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2019 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This file is part of Psi4.
#
# Psi4 is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, version 3.
#
# Psi4 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with Psi4; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

from __future__ import division
from __future__ import print_function

import numpy as np


def blockwise_contract(arr):
    """Undo blockwise_expand, returning `arr` to original 2D array."""

    if len(arr.shape) != 4:
        print('Not appropriate for un_blockwise')
    gr, gc, lr, lc = arr.shape
    ans = np.reshape(arr, (gr * gc, lr, lc))

    def unblockshaped(arr, h, w):
        """
        Return an array of shape (h, w) where
        h * w = arr.size

        If arr is of shape (n, nrows, ncols), n sublocks of shape (nrows, ncols),
        then the returned array preserves the "physical" layout of the sublocks.

        From: https://stackoverflow.com/a/16873755
        """
        n, nrows, ncols = arr.shape
        return (arr.reshape(h // nrows, -1, nrows, ncols).swapaxes(1, 2).reshape(h, w))

    ans = unblockshaped(ans, gr * lr, gc * lc)
    return ans


def blockwise_expand(a, blockshape, aslist=False, require_aligned_blocks=True):
    """
    Return a 2N-D view of the given N-D array, rearranged so each ND block (tile)
    of the original array is indexed by its block address using the first N
    indexes of the output array.

    Note: This function is nearly identical to ``skimage.util.view_as_blocks()``, except:
          - "imperfect" block shapes are permitted (via require_aligned_blocks=False)
          - only contiguous arrays are accepted.  (This function will NOT silently copy your array.)
            As a result, the return value is *always* a view of the input.

    Args:
        a: The ND array
        blockshape: The tile shape

        aslist: If True, return all blocks as a list of ND blocks
                instead of a 2D array indexed by ND block coordinate.
        require_aligned_blocks: If True, check to make sure no data is "left over"
                                in each row/column/etc. of the output view.
                                That is, the blockshape must divide evenly into the full array shape.
                                If False, "leftover" items that cannot be made into complete blocks
                                will be discarded from the output view.

    Sources
    -------
    code: https://github.com/ilastik/lazyflow/blob/master/lazyflow/utility/blockwise_view.py
    so: Inspired by the 2D example shown here: http://stackoverflow.com/a/8070716/162094

    Here's a 2D example (this function also works for ND):

    >>> a = np.arange(1,21).reshape(4,5)
    >>> print a
    [[ 1  2  3  4  5]
     [ 6  7  8  9 10]
     [11 12 13 14 15]
     [16 17 18 19 20]]
    >>> view = blockwise_view(a, (2,2), False)
    >>> print view
    [[[[ 1  2]
       [ 6  7]]
    <BLANKLINE>
      [[ 3  4]
       [ 8  9]]]
    <BLANKLINE>
    <BLANKLINE>
     [[[11 12]
       [16 17]]
    <BLANKLINE>
      [[13 14]
       [18 19]]]]

    """
    assert a.flags['C_CONTIGUOUS'], "This function relies on the memory layout of the array."
    blockshape = tuple(blockshape)
    outershape = tuple(np.array(a.shape) // blockshape)
    view_shape = outershape + blockshape

    if require_aligned_blocks:
        assert (np.mod(a.shape, blockshape) == 0).all(), \
            "blockshape {} must divide evenly into array shape {}"\
            .format( blockshape, a.shape )

    # inner strides: strides within each block (same as original array)
    intra_block_strides = a.strides

    # outer strides: strides from one block to another
    inter_block_strides = tuple(a.strides * np.array(blockshape))

    # This is where the magic happens.
    # Generate a view with our new strides (outer+inner).
    view = np.lib.stride_tricks.as_strided(a, shape=view_shape, strides=(inter_block_strides + intra_block_strides))

    if aslist:
        return list(map(view.__getitem__, np.ndindex(outershape)))
    return view
