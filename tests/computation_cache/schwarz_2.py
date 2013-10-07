#!python
# -*- coding: utf-8 -*-
from __future__ import print_function, division

from os.path import abspath, dirname, join as path_join
from computation_cache import ComputationCache
import numpy
from grendel import Molecule, Tensor, DeclareIndexRange

# The first example in schwarz_1.py shows some basic usage of the
#   computation_cache module.  Some more advanced features will be
#   introduced in this script.

# The parts between these markers is unchanged from schwarz_1.py.
#   For brevity, the documentation has been removed from this part.
#   Read through schwarz_1.py first before reading this script.
#------------- Begin documented in schwarz_1.py -------------#
my_dir = dirname(abspath(__file__))
comp_cache = ComputationCache(
    path_join(my_dir, ".comp_cache"),
    make_directory=True
)

comp = comp_cache.get_computation(
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
    basis="6-31G"
)

g = comp.get_datum(
    "ao_eri",
    pre_computation_callback=lambda: print("Computing two electron integrals for ERI kernel"),
    post_computation_callback=lambda: print("Done!"),
).value.view(Tensor)

schwarz = Tensor(shape=g.shape[:2])
for mu in range(g.shape[0]):
    for nu in range(g.shape[1]):
        schwarz[mu, nu] = g[mu, nu, mu, nu]
schwarz = numpy.sqrt(schwarz)
#-------------- End documented in schwarz_1.py --------------#

# You may have noticed that for larger systems, the computation of
#   the estimates in these lines:
#       g_estimate = Tensor(shape=g.shape)
#       g_estimate["µ,ν,ρ,σ"] = schwarz["µ,ν"] * schwarz["ρ,σ"]
#   is not cached, but could still take quite a bit of time.  It would
#   be better to cache the estimates and only recompute them when
#   something critical to the estimate is different such as the molecular
#   geometry or basis set.  CachedComputation has a mechanism for doing this,
#   in comp.store_datum().  Since what we really want to store is the
#   difference between the estimate and the exact integrals, we can go one
#   step further by storing the difference g - g_estimate.  We'll name
#   the datum "schwarz_estimate_difference" and avoid recomputing it as follows:
if not comp.has_datum("schwarz_estimate_difference"):
    # Only do the computation of g_difference if we don't already have it
    g_difference = Tensor(shape=g.shape)
    # We can use Greek letters now since we set the coding in the second line of
    #   the script!
    g_difference["µ,ν,ρ,σ"] = g["µ,ν,ρ,σ"] - schwarz["µ,ν"] * schwarz["ρ,σ"]
    # Now store the difference so we don't have to compute it next time.
    comp.store_datum(
        # The name can be anything that doesn't conflict with an existing name
        #   in CachedComputation.array_getters or CachedComputation.other_getters
        #   or any previously stored data with a custom getter (described elsewhere).
        "schwarz_estimate_difference",
        # The tensor will be copied into a numpy.memmap wrapped in a
        # CachedMemmapArray instance, just like any other cached tensorial datum
        g_difference
    )
# Retrieve the data, whether it was just computed or computed on a previous run.
#   The process is the same as with the "ao_eri" datum in schwarz_1.py.  Note
#   that you could also do this in an else clause and just keep the g_difference
#   variable when it gets computed, but this way leads to more uniformity
#   and less confusion in more complicated uses.
g_difference = comp.get_datum("schwarz_estimate_difference")

# We can even store simple pieces of data, like the maximum violation, so that
#   we don't have to search through the g_difference for the maximum value every
#   time.
if not comp.has_data("max_schwarz_violation", "max_violation_indices"):
    # Implementation note: The memmap file is not opened and the underlying memmap
    #   array object isn't initialized until the "value" property of the
    #   CachedMemmapArray instance is retrieved.  This can be useful in cases like
    #   this one, where we only want to load g_difference if we need to search
    #   through the array for the maximum value.
    g_difference = g_difference.value.view(numpy.ndarray)
    # This time, let's save the maximum Schwarz violation and the indices of
    #   that violation.  See numpy documentation for more on this next statement
    max_indices = numpy.argmax(g_difference)
    max_indices = numpy.unravel_index(max_indices, g_difference.shape)
    max_violation = g_difference[max_indices]
    # Store the max violation
    comp.store_datum("max_schwarz_violation", max_violation)
    # Any picklable object can be stored in this manner.  For instance, max_indices
    #   is a tuple of integers, with is picklable.  See the Python documentation
    #   for the pickle module for more information.
    comp.store_datum("max_violation_indices", max_indices)
# Retrieve the data as before; remember to get the 'value' attribute
max_indices = comp.get_datum("max_violation_indices").value
max_violation = comp.get_datum("max_schwarz_violation").value

# Now check if it worked, like we did before!
minimum_accuracy = 1e-16
if max_violation < minimum_accuracy:
    print("Schwarz screening works to at least a precision of {0:.1e}!".format(
        minimum_accuracy
    ))
else:
    # This time we can print the indices as well
    print("Schwarz screening failed for minimum accuracy of {0:.1e}."
          "\n  The maximum violation of the Schwarz inequality was {1:.3e}"
          "\n  The indices in g that violated the inequality by this amount were {2}"
          "\n  Perhaps you ask too much of machine precision?".format(
        minimum_accuracy, max_violation, max_indices
    ))
