#!python
# -*- coding: utf-8 -*-
from __future__ import print_function, division

from os.path import abspath, dirname, join as path_join
from computation_cache import ComputationCache, MemmapArrayGetter, ArbitraryBasisDatumGetter
import numpy
import psi4
import itertools
from grendel import Molecule, Tensor, DeclareIndexRange

# As schwarz_2.py did with schwarz_1.py, here we will improve on
#   schwarz_2.py to demonstrate more features of the computation_cache
#   module.  The features introduced here are substantially more
#   advanced than in the first two examples and will mostly only be
#   used by developers looking to expand the computation_cache module
#   itself (though there are plenty of other uses for the features
#   discussed here)

# The parts between these markers is unchanged from schwarz_1.py or
#   schwarz_2.py.  For brevity, the documentation has been removed from this part.
#   Read through schwarz_1.py and schwarz_2.py first before reading this script.
#-------------- Begin documented in previously --------------#
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
    basis="3-21G"
)
#--------------- End documented in previously ---------------#


# You may notice that in the previous scripts, if the two electron integrals are larger
#   than memory, the script will run very slowly the first time through.  In other
#   words, these scripts do not take advantage of the numpy memmap that allows
#   tensors to be stored on disk and accessed in pieces.  To avoid this, we need
#   to write a custom DatumGetter subclass.  DatumGetter instances are what
#   CachedComputation objects use to get data from psi4.  Lots of them have
#   already been written (see computation_cache.py), but in this particular
#   case we have to write our own.
# Note:  This is a pretty contrived case for instructive purposes.  Usually
#   you shouldn't do this unless you're absolutely certain there's no other way.
#   Also, it runs pretty slowly becasue the loop over shells is in python rather
#   than in C++, so this is basically just an example.
# All datum getters should inherit from DatumGetter and implement the __call__
#   method.  In addition, if the datum getter retrieves tensorial data, it
#   should inherit from MemmapArrayGetter (which inherits from DatumGetter) and
#   implement the get_shape() method.  If the getter allows the use of
#   arbitrary basis sets, it should inherit from ArbitraryBasisDatumGetter.
#   CachedComputation looks to see if the datum getter in question is
#   a subclass of MemmapArrayGetter and/or ArbitraryBasisDatumGetter and
#   implements special call behaviors in these cases, as discussed below.
class SchwarzEstimateDifferenceGetter(MemmapArrayGetter):

    # All DatumGetter subclasses must implement __call__.  In the general case,
    #   CachedComputation performs some introspection to look at the argument
    #   names in __call__ and passes the corresponding attributes of the
    #   CachedComputation instance (see the CachedComputation constructor; if
    #   the attribute has not yet been initialized, it is initialized before
    #   it is passed here).  For instance, when this method is called, `basis`
    #   will contain a psi4.BasisSet instance corresponding to the basis name
    #   given in the comp_cache.get_computation() above.  Other attributes
    #   should be similarly self-explanitory.  When a DatumGetter subclass
    #   inherits from MemmapArrayGetter, CachedComputation passes the
    #   output array as the first argument.  The __call__ method should
    #   store the computed data in the out array, which should be a numpy.ndarray
    #   instance (more specifically, a numpy.memmap instance).
    def __call__(self, out, basis):
        # Make a quicky function that gets the slice covered by a given shell
        def shell_slice(ish):
            start = basis.shell_to_basis_function(ish)
            shnbf = basis.shell(ish).nfunction
            return slice(start, start + shnbf)
        #----------------------------------------#
        # Get the number of basis functions and shells
        nbf, nshell = basis.nbf(), basis.nshell()
        # Construct a psi4 integral factory
        factory = psi4.IntegralFactory(basis, basis, basis, basis)
        # Get the TwoBodyInt object
        computer = factory.eri()
        # Enable pybuffer mode (see psi4 source code)
        computer.set_enable_pybuffer(True)
        pybuffer = computer.py_buffer_object
        # Construct the Schwarz array
        schwarz = Tensor(shape=(nbf, nbf))
        #----------------------------------------#
        # Get the Schwarz estimates and the diagonal integrals at once
        # See the itertools module in the python standard library for
        #   a description of the `product` function
        for ish, jsh in itertools.product(range(nshell), repeat=2):
            # Get the shell sizes
            inbf, jnbf = [basis.shell(sh).nfunction for sh in (ish, jsh)]
            # Get the shell slices
            isl, jsl = map(shell_slice, (ish, jsh))
            # Compute the integrals for the shell
            computer.compute_shell(ish, jsh, ish, jsh)
            # Get the buffer and reshape it to the full 4D representation
            #   (see psi4 code for more details of what's going on here)
            buff = numpy.array(pybuffer)
            buff = buff.reshape((inbf, jnbf, inbf, jnbf))
            # Assign the Schwarz integrals to the diagonal (ab|ab) elements
            #   of the buffer
            schwarz[isl, jsl] = numpy.sqrt(
                buff.reshape((inbf*jnbf, inbf*jnbf)).diagonal().reshape(inbf, jnbf)
            )
            # Store the integrals in out for now to avoid recomputing them
            out.value[isl, jsl, isl, jsl] = buff
        #----------------------------------------#
        # Now loop over all the shells, compute the integrals, and store the
        #   difference with the Schwarz estimate.
        for ish, jsh, ksh, lsh in itertools.product(range(nshell), repeat=4):
            # Get the shell sizes
            sizes = [basis.shell(sh).nfunction for sh in (ish, jsh, ksh, lsh)]
            # Get the shell slices
            slices = isl, jsl, ksl, lsl = map(shell_slice, (ish, jsh, ksh, lsh))
            # only compute if we haven't already
            if not (ish == ksh and jsh == lsh):
                computer.compute_shell(ish, jsh, ksh, lsh)
                buff = numpy.array(pybuffer).reshape(sizes)
                out.value[slices] = buff
            # Now compute the Schwarz estimate for the shell
            shell_estimate = Tensor(shape=sizes)
            schwarz_bra = schwarz[isl, jsl]
            schwarz_ket = schwarz[ksl, lsl]
            shell_estimate["a,b,c,d"] = schwarz_bra["a,b"] * schwarz_ket["c,d"]
            # and subtract off the shell estimate to get the difference
            out.value[slices] -= shell_estimate
        #----------------------------------------#
        return out

    # The call sequence for get_shape is similar to that for __call__.  Introspection
    #   is used to determine what attributes should be passed in (only MemmapArrayGetter
    #   subclasses should implement this method).  The shape of the `out` array to
    #   be passed in to the __call__ method should be returned.
    def get_shape(self, basis):
        nbf = basis.nbf()
        return nbf, nbf, nbf, nbf

# Now we can use our custom getter to compute the Schwarz estimate difference
g_difference = comp.get_datum(
    # Again, any unused name can go here
    "schwarz_error",
    # The "custom_getter" keyword is used to pass in an instance of our custom
    #   getter.  Note that in a more general instance, our custom getter might
    #   have initialization arguments, which would obviously be passed to the
    #   constructor, e.g.
    #       SchwarzEstimateDifferenceGetter(arg1, arg2)
    #   One (other) particular use of this keyword is for ArbitraryBasisDatumGetter
    #   instances, which do not have default instances assigned to names in the
    #   CachedComputation.array_getters and CachedComputation.other_getters class
    #   attributes
    custom_getter=SchwarzEstimateDifferenceGetter(),
    # The same callback as before can be used...
    pre_computation_callback=lambda: print("Computing Schwarz difference"),
    post_computation_callback=lambda: print("Done!"),
)

# The rest of this is the same as before
#-------------- Begin documented in previously --------------#
if not comp.has_data("max_schwarz_violation", "max_violation_indices"):
    g_difference = g_difference.value.view(numpy.ndarray)
    max_indices = numpy.argmax(g_difference)
    max_indices = numpy.unravel_index(max_indices, g_difference.shape)
    max_violation = g_difference[max_indices]
    comp.store_datum("max_schwarz_violation", max_violation)
    comp.store_datum("max_violation_indices", max_indices)
max_indices = comp.get_datum("max_violation_indices").value
max_violation = comp.get_datum("max_schwarz_violation").value

minimum_accuracy = 1e-15
if max_violation < minimum_accuracy:
    print("Schwarz screening works to at least a precision of {0:.1e}!".format(
        minimum_accuracy
    ))
else:
    print("Schwarz screening failed for minimum accuracy of {0:.1e}."
          "\n  The maximum violation of the Schwarz inequality was {1:.3e}"
          "\n  The indices in g that violated the inequality by this amount were {2}"
          "\n  Perhaps you ask too much of machine precision?".format(
        minimum_accuracy, max_violation, max_indices
    ))
#--------------- End documented in previously ---------------#
