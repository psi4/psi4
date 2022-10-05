#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2022 The Psi4 Developers.
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

import numpy as np

def plot_coord(ref, cand=None, orig=None, comment=None):
    """Display target geometry `ref` as black dots in 3D plot. If present, also
    plot candidate geometry `cand` as red dots and starting geometry `orig` as
    pale blue dots. Plot has text `comment`. For assessing alignment, red and
    black should overlap and pale blue shows where red started.
    
    """
    try:
        from matplotlib import pyplot
    except ImportError:
        raise ImportError("""Python module matplotlib not found. Solve by installing it: `conda install matplotlib` or https://matplotlib.org/faq/installing_faq.html""")
    from mpl_toolkits.mplot3d import Axes3D

    fig = pyplot.figure()
    ax = fig.add_subplot(111, projection='3d')

    bound = max(np.amax(ref), -1 * np.amin(ref))
    ax.scatter(ref[:, 0], ref[:, 1], ref[:, 2], c='k', label='goal')
    if cand is not None:
        ax.scatter(cand[:, 0], cand[:, 1], cand[:, 2], c='r', label='post-align')
    if orig is not None:
        ax.scatter(orig[:, 0], orig[:, 1], orig[:, 2], c='lightsteelblue', label='pre-align')

    if comment is not None:
        ax.text2D(0.05, 0.95, comment, transform=ax.transAxes)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_xlim(-bound, bound)
    ax.set_ylim(-bound, bound)
    ax.set_zlim(-bound, bound)
    ax.legend()
    
    pyplot.show()
