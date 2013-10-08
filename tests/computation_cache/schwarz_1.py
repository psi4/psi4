#!python
from __future__ import print_function, division

from os.path import abspath, dirname, join as path_join

# You probably have to have psi4 installed as a python module for this to work
from computation_cache import ComputationCache

# You will also need a fairly recent version of numpy (which grendel depends on)
import numpy
from grendel import Molecule, Tensor

# Create a (hidden) ComputationCache in a folder named ".comp_cache" in the same
#   directory as this script.
my_dir = dirname(abspath(__file__))
comp_cache = ComputationCache(
    # The path to the root folder that the ComputationCache will store its data in.
    path_join(my_dir, ".comp_cache"),
    # Allow the ComputationCache constructor to create the directory if one does not
    #   already exist (the default for this option is False).
    make_directory=True
)

# Create a CachedComputation object which can retrieve data for a given Molecule
#   in a cached manner.  For larger-scale projects, it is probably better to
#   wrap this as an attribute in a subclass of grendel.Molecule (examples of this
#   will be given elsewhere), but for this simple example, using the bare
#   CachedComputation object will suffice.
comp = comp_cache.get_computation(
    # The first argument can be an xyz string (as shown here) or anything else
    #   for which grendel.Molecule(arg) returns something sensible.  It can
    #   also be a grendel.Molecule instance, which is necessary to access more
    #   complicated features like charge and multiplicity.
    """
    8
    Ethane
    C     0.76500000   0.00000000  -0.00000000
    C    -0.76500000  -0.00000000  -0.00000000
    H     1.12830000  -0.51380000   0.89000000
    H     1.12830000  -0.51380000  -0.89000000
    H     1.12830000   1.02770000   0.00000000
    H    -1.12830000   0.51380000   0.89000000
    H    -1.12830000   0.51380000  -0.89000000
    H    -1.12830000  -1.02770000   0.00000000
    """,
    # The rest of the call is typically just keyword arguments that get passed
    #   in one form or another to the CachedComputation constructor.  Valid keywords
    #   should be one of the values in ComputationCache.required_attributes,
    #   one of the keys in ComputationCache.optional_attributes, or
    #   ComputationCache.available_psi_options.  The details of these three
    #   types of keywords are explained elsewhere.  Currently, the only required
    #   argument is "basis", which should be a string naming a valid basis set.
    basis="6-31G"
)

# Get the two electron AO integrals into a grenel.Tensor object, which is basically
#   a slightly-more-user-friendly subclass of numpy.ndarray.
g = comp.get_datum(
    # There are a number of ways to get a particular datum from a CachedComputation
    #   object.  For common bits of data like two electron integrals, just the name
    #   will suffice.  If you don't want to do anything fancy, you can leave it at
    #   that.  Allowed values in this instance are given as keys in the
    #   CachedComputation.array_getters and CachedComputation.other_getters dict class
    #   attributes.  Other data retrieval methods will described later.
    "ao_eri",
    # In this form of the get_datum call, the string above is the only required
    #   argument.  However, we can give some optional arguments that make it easier
    #   to understand what's going on under the surface.  Here, we give a
    #   pre_computation_callback and a post_computation_callback, both of which
    #   must be callables taking no arguments.  These two optional arguments are
    #   called before and after any actual computation of the datum takes place,
    #   respectively.  If the data is instead retrieved from a cached version of
    #   a previous computation, neither of these callables will be called.
    pre_computation_callback=lambda: print("Computing two electron integrals for ERI kernel"),
    post_computation_callback=lambda: print("Done!"),
    # Side note:  To use print() in a lambda expression like this, you must have
    #     from __future__ import print_function
    #   in the first non-comment line of your script.
)

# The above call to get_datum() returns an instance of the CachedDatum class.
#   Specifically, when the datum is not a single value but a multiply valued (tensorial)
#   datum like two electron integrals, an instance of CachedMemmapArray is returned
#   (both are classes defined in computation_cache.py).  CachedDatum objects
#   have a `value` attribute, which must be retrieved to use them.
g = g.value

# Now `g` is an instance of numpy.ndarray.  Specifically, `g` is an instance of
#   numpy.memmap(), which maps a numpy ndarray-like object in a file to memory
#   (this is used so that CachedComputation can handle tensors that don't fit
#   in memory).  The data is saved in
#       <root of current ComputationCache>/<unique hash>/ao_eri.npy
#   and can be accessed just like any other *.npy file (see the numpy documentation).
#   Here, it is more useful to "view" the memmap object as a grendel.Tensor.
#   We do so by view-casting (see numpy documentation for details), which basically
#   allows us to utilize a new type for the array without copying the data.
g = g.view(Tensor)
# Note (for numpy veterans; if you're not familiar with numpy, ignore this):
#   You should always view-cast from numpy.memmap objects into other
#   np.ndarray subclasses.  The more familiar syntax
#       new = numpy.array(old)
#   copies the data from the memmap `old`, which puts it into memory.  This defeats
#   the purpose of using a memmap in the first place, since the data in the memmap
#   may be larger than memory.  Similar approaches, such as
#       new = numpy.ndarray(shape=old.shape)
#       new[...] = old
#   also don't work.

# Now that we've got `g` as a grendel.Tensor, we want to verify the Schwarz
#   inequality in this example.  To do this, we get the Schwarz-like integrals
#   into a Tensor of their own:
schwarz = Tensor(
    # Make the schwarz matrix the shape of the first two dimensions of g
    shape=g.shape[:2]
)
# We can easily construct the Schwarz matrix from the full `g` tensor.  The
#   tensor `g` defaults to chemists' notation.
for mu in range(g.shape[0]):
    for nu in range(g.shape[1]):
        schwarz[mu, nu] = g[mu, nu, mu, nu]
schwarz = numpy.sqrt(schwarz)

# Now we can construct the Schwarz estimate.  The extended Einstein summation
#   in grendel comes in handy.
g_estimate = Tensor(shape=g.shape)
g_estimate["a,b,c,d"] = schwarz["a,b"] * schwarz["c,d"]
# Note:  If you want to use Greek letters (which usually denote atomic orbital
#   indices), put the following as the first or second line of the script:
#       # -*- coding: utf-8 -*-
#   (See python documentation for more on this).

# Now verify the Schwarz inequality
# See the numpy manual for details of what's going on here, but this should be
#   mostly straightforward.  We set a minimum accuracy because sometimes machine
#   precision causes apparent violations of the Schwarz inequality.
minimum_accuracy = 1e-15
if numpy.all(abs(g) - g_estimate < minimum_accuracy):
    # In other words, |g| < g_estimate to some minimum_accuracy
    print("Schwarz screening works to at least a precision of {0:.1e}!".format(
        minimum_accuracy
    ))
else:
    max_violation = float(numpy.amax( abs(g) - g_estimate))
    print("Schwarz screening failed for minimum accuracy of {0:.1e}."
          "\n  The maximum violation of the Schwarz inequality was {1:.3e}"
          "\n  Perhaps you ask too much of machine precision?".format(
        minimum_accuracy, max_violation
    ))


