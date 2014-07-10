from core import asarray
from core import ndarray
from core import me

from mpi4py import MPI

import numpy as np

def random_sample(size=None):
    """random_sample(size=None)

        Return random floats in the half-open interval [0.0, 1.0).

        Results are from the "continuous uniform" distribution over the
        stated interval.  To sample :math:`Unif[a, b), b > a` multiply
        the output of `random_sample` by `(b-a)` and add `a`::

          (b - a) * random_sample() + a

        Parameters
        ----------
        size : int or tuple of ints, optional
            Defines the shape of the returned array of random floats. If None
            (the default), returns a single float.

        Returns
        -------
        out : float or ndarray of floats
            Array of random floats of shape `size` (unless ``size=None``, in which
            case a single float is returned).

        Examples
        --------
        >>> np.random.random_sample()
        0.47108547995356098
        >>> type(np.random.random_sample())
        <type 'float'>
        >>> np.random.random_sample((5,))
        array([ 0.30220482,  0.86820401,  0.1654503 ,  0.11659149,  0.54323428])

        Three-by-two array of random numbers from [-5, 0):

        >>> 5 * np.random.random_sample((3, 2)) - 5
        array([[-3.99149989, -0.52338984],
               [-2.99091858, -0.79479508],
               [-1.23204345, -1.75224494]])

    """
    # this was my first implementation
    # but each process does a lot of unnecessary work
    #a = np.random.random_sample(size)
    #if size is None:
    #    return a
    #else:
    #    return asarray(a)
    a = None
    if size is None:
        a = np.random.random_sample(size)
        a = MPI.COMM_WORLD.bcast(a, 0)
    else:
        a = ndarray(size, dtype=float)
        buf = a.access()
        if buf is not None:
            buf[:] = np.random.random_sample(buf.shape)
        a.release_update()
    return a

def seed(seed=None):
    """seed(seed=None)

        Seed the generator.

        This method is called when `RandomState` is initialized. It can be
        called again to re-seed the generator. For details, see `RandomState`.

        Parameters
        ----------
        seed : int or array_like, optional
            Seed for `RandomState`.

        See Also
        --------
        RandomState

    """
    #np.random.seed(seed)
    if seed is None:
        np.random.seed()
    elif isinstance(seed,int):
        np.random.seed(seed+me())
    else:
        a = np.asarray(seed)
        a += me()
        np.random.seed(a)
