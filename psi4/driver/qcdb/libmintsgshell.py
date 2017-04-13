#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2017 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

"""Class to

"""
from __future__ import absolute_import
from __future__ import print_function
import math

#MAX_IOFF = 30000
#extern size_t ioff[MAX_IOFF];
#
#MAX_DF = 500
#extern double df[MAX_DF];
#
#MAX_BC = 20
#extern double bc[MAX_BC][MAX_BC];
#
#MAX_FAC = 100
#extern double fac[MAX_FAC];
#
#
#MAX_DF = 500
#extern double df[MAX_DF];
#
## Globals
#size_t ioff[MAX_IOFF];
#double df[MAX_DF];
#double bc[MAX_BC][MAX_BC];
#double fac[MAX_FAC];
#
#def Wavefunction_initialize_singletons():
#    done = False
#
#    if done:
#        return
#
#    ioff[0] = 0;
#    for (size_t i=1; i<MAX_IOFF; ++i)
#        ioff[i] = ioff[i-1] + i;
#
#    df[0] = 1.0;
#    df[1] = 1.0;
#    df[2] = 1.0;
#    for (int i=3; i<MAX_DF; ++i)
#        df[i] = (i-1)*df[i-2];
#
#    for (int i=0; i<MAX_BC; ++i)
#        for (int j=0; j<=i; ++j)
#            bc[i][j] = combinations(i, j);
#
#    fac[0] = 1.0;
#    for (int i=1; i<MAX_FAC; ++i)
#        fac[i] = i*fac[i-1];
#
#    done = True


def df(n):
    """Gives the double factorial of *n*"""
    return 1.0 if n <= 0 else 1.0 * n * df(n - 2)


def INT_NCART(am):
    """Gives the number of cartesian functions for an angular momentum.
    define INT_NCART(am) ((am>=0) ? ((((am)+2)*((am)+1))>>1) : 0)

    """
    return (((am + 2) * (am + 1)) >> 1) if (am >= 0) else 0


def INT_NPURE(am):
    """Gives the number of spherical functions for an angular momentum.
    #define INT_NPURE(am) (2*(am)+1)

    """
    return 2 * am + 1


def INT_NFUNC(pu, am):
    """Gives the number of functions for an angular momentum based on pu.
    #define INT_NFUNC(pu,am) ((pu)?INT_NPURE(am):INT_NCART(am))

    """
    if pu == 'Cartesian' or pu == False:
        return INT_NCART(am)
    else:
        return INT_NPURE(am)


def INT_CARTINDEX(am, i, j):
    """Computes offset index for cartesian function.
    #define INT_CARTINDEX(am,i,j) (((i) == (am))? 0 : (((((am) - (i) + 1)*((am) - (i)))>>1) + (am) - (i) - (j)))

    """
    return 0 if (i == am) else ((((am - i + 1) * (am - i)) >> 1) + am - i - j)


def INT_ICART(a, b, c):
    """Given a, b, and c compute a cartesian offset.
    #define INT_ICART(a, b, c) (((((((a)+(b)+(c)+1)<<1)-(a))*((a)+1))>>1)-(b)-1)

    """
    return ((((((a + b + c + 1) << 1) - a) * (a + 1)) >> 1) - b - 1)


def INT_IPURE(l, m):
    """Given l and m compute a pure function offset.
    #define INT_IPURE(l, m) ((l)+(m))

    """
    return l + m


# Lookup array that when you index the angular momentum it returns the corresponding letter
PrimitiveType = ['Normalized', 'Unnormalized']
GaussianType = ['Cartesian', 'Pure']  # Cartesian = 0, Pure = 1


class ShellInfo(object):
    """This class has the same behavior as GaussianShell, but implements everything using
    slower data structures, which are easier to construct. These are used to build the
    basis set, which builds more efficient pointer-based GaussianShell objects.
    @param e An array of exponent values.
    @param am Angular momentum.
    @param pure Pure spherical harmonics, or Cartesian.
    @param c An array of contraction coefficients.
    @param nc The atomic center that this shell is located on. Must map
    back to the correct atom in the owning BasisSet molecule. Used
    in integral derivatives for indexing.
    @param center The x, y, z position of the shell. This is passed to
    reduce the number of calls to the molecule.
    @param start The starting index of the first function this shell
    provides. Used to provide starting positions in matrices.
    @param pt Is the shell already normalized?

    """

    def __init__(self, am, c, e, pure, nc, center, start, pt='Normalized', rpowers=None):
        # Angular momentum
        self.l = am
        # Flag for pure angular momentum (Cartesian = 0, Pure = 1)
        self.puream = pure
        # Exponents (of length nprimitives_)
        self.PYexp = e
        # Contraction coefficients (of length nprimitives_)
        self.PYcoef = c
        # ERD normalized contraction coefficients (of length nprimitives_)
        self.PYerd_coef = []
        # Original (un-normalized) contraction coefficients (of length nprimitives)
        self.PYoriginal_coef = [c[n] for n in range(len(c))]
        # Atom number this shell goes to. Needed when indexing integral derivatives.
        self.nc = nc
        # Atomic center number in the Molecule
        self.center = center
        #
        self.start = start
        # How many cartesian functions? (1=s, 3=p, 6=d, ...)
        self.PYncartesian = INT_NCART(self.l)
        # How many functions? (1=s, 3=p, 5/6=d, ...) * Dependent on the value of puream_
        self.PYnfunction = INT_NFUNC(self.puream, self.l)
        # These are the radial factors for ECPs.  They are not defined for regular shells.
        self.rpowers = rpowers

        # Compute the normalization constants
        if pt == 'Unnormalized':
            self.normalize_shell()
            self.erd_normalize_shell()
        else:
            self.PYerd_coef = [0.0 for i in c]

    def primitive_normalization(self, p):
        """Normalizes a single primitive.
        @param p The primitive index to normalize.
        @return Normalization constant to be applied to the primitive.

        """
        tmp1 = self.l + 1.5
        g = 2.0 * self.PYexp[p]
        z = pow(g, tmp1)
        return math.sqrt((pow(2.0, self.l) * z) / (math.pi * math.sqrt(math.pi) * df(2 * self.l)))

    def contraction_normalization(self):
        """Normalizes an entire contraction set. Applies the normalization to the coefficients
        *  @param gs The contraction set to normalize.

        """
        e_sum = 0.0
        for i in range(self.nprimitive()):
            for j in range(self.nprimitive()):
                g = self.PYexp[i] + self.PYexp[j]
                z = pow(g, self.l + 1.5)
                e_sum += self.PYcoef[i] * self.PYcoef[j] / z

        tmp = ((2.0 * math.pi / (2.0 / math.sqrt(math.pi))) * df(2 * self.l)) / pow(2.0, self.l)
        try:
            norm = math.sqrt(1.0 / (tmp * e_sum))
        except ZeroDivisionError:
            self.PYcoef[i] = [1.0 for i in range(self.nprimitive())]
        # Set the normalization
        for i in range(self.nprimitive()):
            self.PYcoef[i] *= norm

    def normalize_shell(self):
        """Handles calling primitive_normalization and
        contraction_normalization for you.

        """
        for i in range(self.nprimitive()):
            normalization = self.primitive_normalization(i)
            self.PYcoef[i] *= normalization
        self.contraction_normalization()

    def erd_normalize_shell(self):
        """

        """
        self.PYerd_coef = []
        tsum = 0.0
        for j in range(self.nprimitive()):
            for k in range(j + 1):
                a1 = self.PYexp[j]
                a2 = self.PYexp[k]
                temp = self.PYoriginal_coef[j] * self.PYoriginal_coef[k]
                temp2 = self.l + 1.5
                temp3 = 2.0 * math.sqrt(a1 * a2) / (a1 + a2)
                temp3 = pow(temp3, temp2)
                temp *= temp3
                tsum += temp
                if j != k:
                    tsum += temp
        prefac = 1.0
        if self.l > 1:
            prefac = pow(2.0, 2 * self.l) / df(2 * self.l)
        norm = math.sqrt(prefac / tsum)
        for j in range(self.nprimitive()):
            self.PYerd_coef.append(self.PYoriginal_coef[j] * norm)

    def copy(self, nc=None, c=None):
        """Make a copy of the ShellInfo"""
        if nc is not None and c is not None:
            return ShellInfo(self.l, self.PYoriginal_coef, self.PYexp,
                self.puream, nc, c,
                self.start, 'Unnormalized', self.rpowers)
        else:
            return ShellInfo(self.l, self.PYoriginal_coef, self.PYexp,
                self.puream, self.nc, self.center,
                self.start, 'Unnormalized', self.rpowers)
        # better to just deepcopy?

    def nprimitive(self):
        """The number of primitive Gaussians"""
        return len(self.PYexp)

    def nfunction(self):
        """Total number of basis functions"""
        return INT_NFUNC(self.puream, self.l)

    def ncartesian(self):
        """Total number of functions if this shell was Cartesian"""
        return self.PYncartesian

    def am(self):
        """The angular momentum of the given contraction"""
        return self.l

    def amchar(self):
        """The character symbol for the angular momentum of the given contraction"""
        return 'spdfghiklmnopqrtuvwxyz'[self.l]

    def AMCHAR(self):
        """The character symbol for the angular momentum of the given contraction (upper case)"""
        return 'SPDFGHIKLMNOPQRTUVWXYZ'[self.l]

    def is_cartesian(self):
        """Returns true if contraction is Cartesian"""
        return True if self.puream == 'Cartesian' else False

    def is_pure(self):
        """Returns true if contraction is pure"""
        return True if self.puream == 'Pure' else False

    def center(self):
        """Returns the center of the Molecule this shell is on"""
        return self.center

    def ncenter(self):
        """Returns the atom number this shell is on. Used by integral derivatives for indexing."""
        return self.nc

    def exp(self, prim):
        """Returns the exponent of the given primitive"""
        return self.PYexp[prim]

    def coef(self, pi):
        """Return coefficient of pi'th primitive"""
        return self.PYcoef[pi]

    def erd_coef(self, pi):
        """Return ERD normalized coefficient of pi'th primitive"""
        return self.PYerd_coef[pi]

    def original_coef(self, pi):
        """Return unnormalized coefficient of pi'th primitive"""
        return self.PYoriginal_coef[pi]

    def rpower(self, pi):
        """Return r exponent (for ECP) of pi'th primitive"""
        return self.rpowers[pi] if self.rpowers else None

    def exps(self):
        """Returns the exponent of the given primitive"""
        return self.PYexp

    def coefs(self):
        """Return coefficient of pi'th primitive and ci'th contraction"""
        return self.PYcoef

    def original_coefs(self):
        """Return unnormalized coefficient of pi'th primitive and ci'th contraction"""
        return self.PYoriginal_coef

    def aslist(self):
        """Return minimal list of shell info"""
        info = [self.l] + [(self.PYexp[K], self.PYoriginal_coef[K], self.rpower(K)) for K in range(self.nprimitive())]
        return info

    def pyprint(self, outfile=None):
        """Print out the shell"""
        text = """    %c %3d 1.00\n""" % (self.AMCHAR(), self.nprimitive())
        for K in range(self.nprimitive()):
            text += """               %20.8f %20.8f\n""" % (self.PYexp[K], self.PYoriginal_coef[K])

        if outfile is None:
            return text
        else:
            with open(outfile, mode='w') as handle:
                handle.write(text)

    def pyprint_gamess(self, outfile=None):
        """Print out the shell in Gamess format"""
        text = """%c %3d\n""" % (self.AMCHAR(), self.nprimitive())
        for K in range(self.nprimitive()):
            text += """%3d %15.8f %15.8f\n""" % (K + 1, self.PYexp[K], self.PYoriginal_coef[K])

        if outfile is None:
            return text
        else:
            with open(outfile, mode='w') as handle:
                handle.write(text)

    def __str__(self):
        """String representation of shell"""
        return self.pyprint(outfile=None)

    def normalize(self, l, m, n):
        """Normalize the angular momentum component"""
        return 1.0

    def function_index(self):
        """Basis function index where this shell starts."""
        return self.start

    def set_function_index(self, i):
        """Set basis function index where this shell starts."""
        self.start = i


class GaussianShell(ShellInfo):
    """Class with same information as :py:class:`ShellInfo`. In C++,
    class uses more efficient data structures, but in Python differences
    minimal.

    """

    def __init__(self, am, nprimitive, oc, c, ec, e, pure, nc, center, start, rpowers=None):
        ShellInfo.__init__(self, am=am, c=c, e=e, pure=pure, nc=nc, center=center, start=start, pt='Normalized', rpowers=rpowers)
        """
        *  @param am Angular momentum.
        *  @param pure Pure spherical harmonics, or Cartesian.
        *  @param oc An array of contraction coefficients.
        *  @param c An array of normalized contraction coefficients.
        *  @param ec An array of ERD normalized contraction coefficients.
        *  @param e An array of exponent values.
        *  @param pure an enum describing whether this shell uses pure or Cartesian functions.
        *  @param nc The atomic center that this shell is located on. Must map back to the correct atom in the owning BasisSet molecule_. Used in integral derivatives for indexing.
        *  @param center The x, y, z position of the shell. This is passed to reduce the number of calls to the molecule.
        *  @param start The starting index of the first function this shell provides. Used to provide starting positions in matrices.
        *  @param pt Is the shell already normalized?

        """
        self.PYnprimitive = nprimitive
        self.PYoriginal_coef = oc
        self.PYerd_coef = ec

    def nprimitive(self):
        """The number of primitive Gaussians"""
        return self.PYnprimitive



#GaussianShell(0, nprimitive_,
#    uoriginal_coefficients_, ucoefficients_, uerd_coefficients_,
#    uexponents_, GaussianType(0), 0, xyz_, 0)
#
#GaussianShell(am, shell_nprim,
#    &uoriginal_coefficients_[ustart+atom_nprim], &ucoefficients_[ustart+atom_nprim], &uerd_coefficients_[ustart+atom_nprim],
#    &uexponents_[ustart+atom_nprim], puream, n, xyz_ptr, bf_count)
#
#GaussianShell(am, shell_nprim,
#    &uoriginal_coefficients_[prim_count], &ucoefficients_[prim_count], &uerd_coefficients_[prim_count],
#    &uexponents_[prim_count], puream, center, xyz_, bf_count)
#
#ShellInfo(am, contractions, exponents, gaussian_type, 0, center, 0, Unnormalized)
