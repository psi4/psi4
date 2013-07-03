""" combinatorics.py

OVERVIEW
========

This module was created to supplement Python's itertools module, filling in gaps
in two important areas of basic combinatorics:

(A) ordered and unordered m-way combinations, and
(B) generalizations of the four basic occupancy problems ('balls in boxes').

Brief descriptions of the included functions and classes follow (more detailed
descriptions and additional examples can be found in the individual doc strings
within the functions):

n_choose_m(n, m): calculate n-choose-m, using a simple algorithm that is less
likely to involve large integers than the direct evaluation of n! / m! / (n-m)!

m_way_ordered_combinations(items, ks): This function returns a generator that
produces all m-way ordered combinations (multinomial combinations) from the
specified collection of items, with with ks[i] items in the ith group, i= 0, 1,
2, ..., m-1, where m= len(ks) is the number of groups. By 'ordered
combinations', we mean that the relative order of equal- size groups is
important; the order of the items within any group is not important. The total
number of combinations generated is given by the multinomial coefficient formula
(see http://en.wikipedia.org/wiki/Multinomial_theorem#Multinomial_coefficients).

m_way_unordered_combinations(items, ks): This function returns a generator that
produces all m-way unordered combinations from the specified collection of
items, with ks[i] items in the ith group, i= 0, 1, 2, ..., m-1, where m= len(ks)
is the number of groups. By 'unordered combinations', we mean that the relative
order of equal-size groups is not important. The order of the items within any
group is also unimportant.

Example of `m_way_unordered_combinations`::

   Issue the following statement from the IPython prompt:

   from combinatorics import *
   list(m_way_unordered_combinations(6,[2,2,2]))

   The output consists of the 15 combinations listed below:

   (0, 1), (2, 3), (4, 5)
   (0, 1), (2, 4), (3, 5)
   (0, 1), (2, 5), (3, 4)
   (0, 2), (1, 3), (4, 5)
   (0, 2), (1, 4), (3, 5)
   (0, 2), (1, 5), (3, 4)
   (0, 3), (1, 2), (4, 5)
   (0, 3), (1, 4), (2, 5)
   (0, 3), (1, 5), (2, 4)
   (0, 4), (1, 2), (3, 5)
   (0, 4), (1, 3), (2, 5)
   (0, 4), (1, 5), (2, 3)
   (0, 5), (1, 2), (3, 4)
   (0, 5), (1, 3), (2, 4)
   (0, 5), (1, 4), (2, 3)

unlabeled_balls_in_labeled_boxes(balls, box_sizes): This function returns a
generator that produces all distinct distributions of indistinguishable balls
among labeled boxes with specified box sizes (capacities). This is a
generalization of the most common formulation of the problem, where each box is
sufficiently large to accommodate all of the balls, and is an important example
of a class of combinatorics problems called 'weak composition' problems.

unlabeled_balls_in_unlabeled_boxes(balls, box_sizes): This function returns a
generator that produces all distinct distributions of indistinguishable balls
among indistinguishable boxes, with specified box sizes (capacities). This is a
generalization of the most common formulation of the problem, where each box is
sufficiently large to accommodate all of the balls. It might be asked, 'In what
sense are the boxes indistinguishable if they have different capacities?' The
answer is that the box capacities must be considered when distributing the
balls, but once the balls have been distributed, the identities of the boxes no
longer matter.

Example of `unlabeled_balls_in_unlabeled_boxes`::

   Issue the following commands from the IPython prompt:

   from combinatorics import *
   list(unlabeled_balls_in_unlabeled_boxes(10,[5,4,3,2,1]))

   The output is as follows:

   [(5, 4, 1, 0, 0),
    (5, 3, 2, 0, 0),
    (5, 3, 1, 1, 0),
    (5, 2, 2, 1, 0),
    (5, 2, 1, 1, 1),
    (4, 4, 2, 0, 0),
    (4, 4, 1, 1, 0),
    (4, 3, 3, 0, 0),
    (4, 3, 2, 1, 0),
    (4, 3, 1, 1, 1),
    (4, 2, 2, 2, 0),
    (4, 2, 2, 1, 1),
    (3, 3, 3, 1, 0),
    (3, 3, 2, 2, 0),
    (3, 3, 2, 1, 1),
    (3, 2, 2, 2, 1)]

labeled_balls_in_unlabeled_boxes(balls, box_sizes): This function returns a
generator that produces all distinct distributions of distinguishable balls
among indistinguishable boxes, with specified box sizes (capacities). This is a
generalization of the most common formulation of the problem, where each box is
sufficiently large to accommodate all of the balls.

labeled_balls_in_labeled_boxes(balls, box_sizes): This function returns a
generator that produces all distinct distributions of distinguishable balls
among distinguishable boxes, with specified box sizes (capacities). This is a
generalization of the most common formulation of the problem, where each box is
sufficiently large to accommodate all of the balls.

Example of `labeled_balls_in_labeled_boxes`::

   Issue the following statements from the IPython prompt:

   from combinatorics import *
   list(labeled_balls_in_labeled_boxes(3,[2,2]))

   The output is as follows:

   [((0, 1), (2,)),
    ((0, 2), (1,)),
    ((1, 2), (0,)),
    ((0,), (1, 2)),
    ((1,), (0, 2)),
    ((2,), (0, 1))]

partitions(n): 'In number theory and combinatorics, a partition of a positive
integer n, also called an integer partition, is a way of writing n as a sum of
positive integers. Two sums that differ only in the order of their summands are
considered to be the same partition.'  We can trivially generate all partitions
of an integer using `unlabeled_balls_in_unlabeled_boxes`.  The quote is from
http://en.wikipedia.org/wiki/Partition_(number_theory) .


AUTHOR
======

Dr. Phillip M. Feldman

Comments and suggestions--especially bug reports--can be communicated to me via
the following e-mail address: Phillip.M.Feldman@gmail.com


REVISION HISTORY
================

04-06-2012, version 1.2.0, Phillip M. Feldman:

I added the function `prod`, which is similar to `numpy.prod` but does all
calculations using large arithmetic when operating on a sequence of integers.

I fixed a bug in `n_choose_m`: We must force the division to be done using
integer arithmetic because otherwise Python attempts to convert the results from
`prod` into floating point numbers, which can fail for n greater than 170.

I added the function `n_choose_m_ln`.  This function calculates the natural
logarithm of choose(n,m), defined as the number of ways in which one can select
m of n distinct objects without regard for order, using SciPy's `gammaln`
function.  For large n, especially for n > 10000, this function is much faster
than `n_choose_m` (computational and memory requirements are both much lower).


10-09-2011, version 1.1.1, Phillip M. Feldman:

I added a function to generate partitions.


10-01-2011, version 1.1.0, Phillip M. Feldman:

I added input error checking to the `labeled_balls_in_unlabeled_boxes` function.

I fixed the function `m_way_ordered_combinations` so that the box order
specified via the input argument `ks` is respected.

I added the function `labeled_balls_in_labeled_boxes`.  This completes the basic
set of functions for solving occupancy functions with capacity limits.

09-24-2011: Initial version.
"""

# NOTE:  MIT LICENSE
# from http://pypi.python.org/pypi/Combinatorics/1.2.0


import collections, itertools, math, operator
from numpy import sort

# Define `fact` as a shorthand name for the `factorial` function:
fact= math.factorial


def prod(seq):
   """
   Because NumPy's `prod` function uses 32-bit integer arithmetic with silent
   handling of overflows, results are wrong if the correct answer would exceed
   the limits of a signed 32-bit integer.  When operating on a sequence of
   integers, the `prod` function that we define here uses large integer
   arithmetic and thus always gives correct results.
   """
   return reduce(operator.mul, seq)


def n_choose_m(n, m):
   """
   OVERVIEW

   This function calculates choose(n,m), defined as the number of ways in which
   one can select m of n distinct objects without regard for order, using only
   integer arithmetic.  The calculation is done as follows:

   1. If m > n-m, we replace m by n-m.

   2. We calculate the answer by evaluating

      prod(range(n-m+1,n+1)) / prod(range(2,m+1)),

   which is equivalent to

      n! / m! / (n-m)!


   NOTE

   Python can handle integers of arbitrary size, but the algorithm tends to bog
   down for very large values of n, partly because of the number of operations
   being performed and partly because of the memory requirements.  For values of
   n above about 10000, use `m_choose_n_ln` instead of `m_choose_n`.
   """

   if not isinstance(n,int) or not isinstance(m,int):
      raise TypeError('The inputs n and m must have type int.')

   if m < 0 or m > n:
      raise ValueError("m (the second argument) must be between 0 and n, "
        "inclusive.")

   if m > n-m: m= n-m

   if m == 0:
      return 1
   elif m == 1:
      return n

   # In the following statement, we force the division to be done using integer
   # arithmetic because otherwise Python attempts to convert the results from
   # `prod` into floating point numbers, which can fail for n greater than 170.
   return prod(range(n-m+1,n+1)) // prod(range(2,m+1))


def n_choose_m_ln(n, m):
   """
   OVERVIEW

   This function calculates the natural logarithm of choose(n,m), defined as the
   number of ways in which one can select m of n distinct objects without regard
   for order, using SciPy's `gammaln` function.  For large n, especially for
   n > 10000, this function is much faster than `n_choose_m` (computational and
   memory requirements are both much lower).


   NOTE

   To obtain a value for choose(n,m), apply the `exp` function to the result
   returned by this function.

   This function works for huge values of `n`, but applying `exp` to the return
   value may produce an overflow.
   """

   return gammaln(n+1) - gammaln(m+1) - gammaln(n-m+1)


def m_way_ordered_combinations(items, ks):
   """
   OVERVIEW

   This function returns a generator that produces all m-way ordered
   combinations (multinomial combinations) from the specified collection of
   items, with ks[i] items in the ith group, i= 0, 1, 2, ..., m-1, where m=
   len(ks) is the number of groups.  By 'ordered combinations', we mean that the
   relative order of equal-size groups is important.  The order of the items
   within any group is not important.  The total number of combinations
   generated is given by the multinomial coefficient formula (see below).


   INPUTS

   `items` must be (A) a list, tuple, or other iterable, or (B) a positive
   integer.  If `items` is an integer, it is replaced by `range(items)`.

   `ks` should be either a list or tuple containing non-negative integers, where
   the sum of these integers does not exceed the length of `items`.


   EXAMPLE

   Let items=[0,1,2,3,4,5] and ks=[2,2,2].  The output includes a total of 90
   combinations.  Two of these are the following:

   ((0, 1), (2, 3), (4, 5))
   ((2, 3), (0, 1), (4, 5))

   These are distinct because the order of the groups, which differs, is
   significant.


   NOTES

   The total number of combinations generated is given by the following
   multinomial coefficient:

                n!
   ----------------------------
   k_0! * k_1! * ... * k_(m-1)!

   where n is the number of items, m is the number of groups, and k_i is the
   number of items in the ith group.
   """
   if isinstance(items,int):
      items= range(items)
   elif not isinstance(items,collections.Iterable):
      raise TypeError("`items` must be a list, tuple, or other iterable.")

   if not isinstance(ks,(list,tuple)):
      raise TypeError("`ks` must be a list or tuple.")

   return _m_way_ordered_combinations(items, ks)

# end def m_way_ordered_combinations

def _m_way_ordered_combinations(items, ks):

   if len(ks) == 1:
      for c in itertools.combinations(items, ks[0]):
         yield (c,)

   else:
      for c_first in itertools.combinations(items, ks[0]):
         items_remaining= set(items) - set(c_first)
         for c_other in _m_way_ordered_combinations(items_remaining, ks[1:]):
            yield (c_first,) + c_other

# end def _m_way_ordered_combinations(items, ns)


def m_way_unordered_combinations(items, ks):
   """
   OVERVIEW

   This function returns a generator that produces all m-way unordered
   combinations from the specified collection of items, with ks[i] items in the
   ith group, i= 0, 1, 2, ..., m-1, where m= len(ks) is the number of groups.
   By 'unordered combinations', we mean that the relative order of equal-size
   groups is not important.  The order of the items within any group is also
   unimportant.


   INPUTS

   `items` must be (A) a list, tuple, or other iterable, or (B) a positive
   integer.  If `items` is an integer, it is replaced by `range(items)`.

   `ks` should be either a list or tuple containing non-negative integers, where
   the sum of these integers does not exceed the length of `items`.


   EXAMPLE

   Issue the following statement issued at the IPython prompt:

   list(m_way_unordered_combinations(6,[2,2,2]))

   The output consists of the 15 combinations listed below:

   (0, 1), (2, 3), (4, 5)
   (0, 1), (2, 4), (3, 5)
   (0, 1), (2, 5), (3, 4)
   (0, 2), (1, 3), (4, 5)
   (0, 2), (1, 4), (3, 5)
   (0, 2), (1, 5), (3, 4)
   (0, 3), (1, 2), (4, 5)
   (0, 3), (1, 4), (2, 5)
   (0, 3), (1, 5), (2, 4)
   (0, 4), (1, 2), (3, 5)
   (0, 4), (1, 3), (2, 5)
   (0, 4), (1, 5), (2, 3)
   (0, 5), (1, 2), (3, 4)
   (0, 5), (1, 3), (2, 4)
   (0, 5), (1, 4), (2, 3)


   NOTES

   When all group sizes are unequal, the total number of combinations generated
   is given by the multinomial coefficient (see above).  When two or more groups
   have equal sizes, the number of combinations is less than the multinomial
   coefficient because combinations that differ only in the relative order of
   equal-size groups are excluded.
   """
   if isinstance(items,int):
      items= range(items)
   elif not isinstance(items,collections.Iterable):
      raise TypeError("`items` must be a list, tuple, or other iterable.")

   if not isinstance(ks,(list,tuple)):
      raise TypeError("`ks` must be a list or tuple.")

   # Sort group sizes from largest to smallest:
   ks= list( sort(ks)[::-1] )

   return _m_way_unordered_combinations(items, ks)

# end def m_way_unordered_combinations

def _m_way_unordered_combinations(items, ks):

   if all(k == 0 for k in ks[1:]):
      for c in itertools.combinations(items, ks[0]):
         yield (c,) + ((),) * (len(ks) - 1)

   else:
      for c_first in itertools.combinations(items, ks[0]):
         items_remaining= set(items) - set(c_first)
         for c_other in \
           _m_way_unordered_combinations(items_remaining, ks[1:]):
            if len(c_first)!=len(c_other[0]) or c_first<c_other[0]:
               yield (c_first,) + c_other

# end def _m_way_unordered_combinations(items, ns)


def unlabeled_balls_in_labeled_boxes(balls, box_sizes):
   """
   OVERVIEW

   This function returns a generator that produces all distinct distributions of
   indistinguishable balls among labeled boxes with specified box sizes
   (capacities).  This is a generalization of the most common formulation of the
   problem, where each box is sufficiently large to accommodate all of the
   balls, and is an important example of a class of combinatorics problems
   called 'weak composition' problems.


   CONSTRUCTOR INPUTS

   n: the number of balls

   box_sizes: This argument is a list of length 1 or greater.  The length of
   the list corresponds to the number of boxes.  `box_sizes[i]` is a positive
   integer that specifies the maximum capacity of the ith box.  If
   `box_sizes[i]` equals `n` (or greater), then the ith box can accommodate all
   `n` balls and thus effectively has unlimited capacity.


   ACKNOWLEDGMENT

   I'd like to thank Chris Rebert for helping me to convert my prototype
   class-based code into a generator function.
   """
   if not isinstance(balls, int):
      raise TypeError("balls must be a non-negative integer.")
   if balls < 0:
      raise ValueError("balls must be a non-negative integer.")

   if not isinstance(box_sizes,list):
      raise ValueError("box_sizes must be a non-empty list.")

   capacity= 0
   for size in box_sizes:
      if not isinstance(size, int):
          raise TypeError("box_sizes must contain only positive integers.")
      if size < 1:
          raise ValueError("box_sizes must contain only positive integers.")
      capacity+= size

   if capacity < balls:
      raise ValueError("The total capacity of the boxes is less than the "
        "number of balls to be distributed.")

   return _unlabeled_balls_in_labeled_boxes(balls, box_sizes)


def _unlabeled_balls_in_labeled_boxes(balls, box_sizes):
   """
   This recursive generator function was designed to be returned by
   `unlabeled_balls_in_labeled_boxes`.
   """

   # If there are no balls, then the only possible distribution is for all boxes
   # to be empty:
   if not balls:
      yield len(box_sizes) * (0,)

   elif len(box_sizes) == 1:

      # If the single available box has sufficient capacity to store the balls,
      # there is only one possible distribution, and we return it to the caller
      # via `yield`.  Otherwise, the flow of control will pass to the end of the
      # function, triggering a `StopIteration` exception.
      if box_sizes[0] >= balls:
          yield (balls,)

   else:

      # Iterate over the number of balls in the first box (from the maximum
      # possible down to zero), recursively invoking the generator to distribute
      # the remaining balls among the remaining boxes.
      for balls_in_first_box in xrange( min(balls, box_sizes[0]), -1, -1 ):
         balls_in_other_boxes= balls - balls_in_first_box

         for distribution_other in _unlabeled_balls_in_labeled_boxes(
           balls_in_other_boxes, box_sizes[1:]):
            yield (balls_in_first_box,) + distribution_other

   # end three alternative blocks

# end def _unlabeled_balls_in_labeled_boxes(balls, box_sizes)


def unlabeled_balls_in_unlabeled_boxes(balls, box_sizes):
   """
   OVERVIEW

   This function returns a generator that produces all distinct distributions of
   indistinguishable balls among indistinguishable boxes, with specified box
   sizes (capacities).  This is a generalization of the most common formulation
   of the problem, where each box is sufficiently large to accommodate all of
   the balls.  It might be asked, 'In what sense are the boxes indistinguishable
   if they have different capacities?' The answer is that the box capacities
   must be considered when distributing the balls, but once the balls have been
   distributed, the identities of the boxes no longer matter.


   CONSTRUCTOR INPUTS

   n: the number of balls

   box_sizes: This argument is a list of length 1 or greater.  The length of
   the list corresponds to the number of boxes.  `box_sizes[i]` is a positive
   integer that specifies the maximum capacity of the ith box.  If
   `box_sizes[i]` equals `n` (or greater), then the ith box can accommodate all
   `n` balls and thus effectively has unlimited capacity.


   NOTE

   For `unlabeled_balls_in_unlabeled_boxes`, the order of the elements of the
   `box_sizes` list is unimportant because the code will sort it into non-
   increasing order before any other processing is done.
   """
   if not isinstance(balls, int):
      raise TypeError("balls must be a non-negative integer.")
   if balls < 0:
      raise ValueError("balls must be a non-negative integer.")

   if not isinstance(box_sizes,list):
      raise ValueError("box_sizes must be a non-empty list.")

   capacity= 0
   for size in box_sizes:
      if not isinstance(size, int):
          raise TypeError("box_sizes must contain only positive integers.")
      if size < 1:
          raise ValueError("box_sizes must contain only positive integers.")
      capacity+= size

   if capacity < balls:
      raise ValueError("The total capacity of the boxes is less than the "
        "number of balls to be distributed.")

   # Sort the box sizes so that the values decrease:
   box_sizes= list( sort(box_sizes)[::-1] )

   return _unlabeled_balls_in_unlabeled_boxes(balls, box_sizes)

# def unlabeled_balls_in_unlabeled_boxes(balls, box_sizes)


def _unlabeled_balls_in_unlabeled_boxes(balls, box_sizes):
   """
   This recursive generator function was designed to be returned by
   `unlabeled_balls_in_unlabeled_boxes`.
   """

   # If there are no balls, then the only possible distribution is for all boxes
   # to be empty:
   if not balls:
      yield len(box_sizes) * (0,)

   elif len(box_sizes) == 1:

      # If the single available box has sufficient capacity to store the balls,
      # there is only one possible distribution, and we return it to the caller
      # via `yield`.  Otherwise, the flow of control will pass to the end of the
      # generator function, triggering a `StopIteration` exception.
      if box_sizes[0] >= balls:
          yield (balls,)

   else:

      # Iterate over the number of balls in the first box (from the maximum
      # possible down to zero), recursively invoking the generator to distribute
      # the remaining balls among the remaining boxes.
      for balls_in_first_box in xrange( min(balls, box_sizes[0]), -1, -1 ):
         balls_in_other_boxes= balls - balls_in_first_box

         for distribution_other in _unlabeled_balls_in_unlabeled_boxes(
           balls_in_other_boxes, box_sizes[1:]):

            # To prevent the possibility of duplicating a distribution that has
            # been obtained previously, we require that the number of balls in
            # the second box (first of the 'other' boxes) not exceed the number
            # of balls in the first box:
            if distribution_other[0] <= balls_in_first_box:
               yield (balls_in_first_box,) + distribution_other

   # end three alternative blocks

# def _unlabeled_balls_in_unlabeled_boxes(balls, box_sizes)


def labeled_balls_in_unlabeled_boxes(balls, box_sizes):
   """
   OVERVIEW

   This function returns a generator that produces all distinct distributions of
   distinguishable balls among indistinguishable boxes, with specified box sizes
   (capacities).  This is a generalization of the most common formulation of the
   problem, where each box is sufficiently large to accommodate all of the
   balls.  It might be asked, 'In what sense are the boxes indistinguishable if
   they have different capacities?'  The answer is that the box capacities must
   be considered when distributing the balls, but once the balls have been
   distributed, the identities of the boxes no longer matter.


   CONSTRUCTOR INPUTS

   n: the number of balls

   box_sizes: This argument is a list of length 1 or greater.  The length of
   the list corresponds to the number of boxes.  `box_sizes[i]` is a positive
   integer that specifies the maximum capacity of the ith box.  If
   `box_sizes[i]` equals `n` (or greater), then the ith box can accommodate all
   `n` balls and thus effectively has unlimited capacity.


   NOTE

   For `labeled_balls_in_unlabeled_boxes`, the order of the elements of the
   `box_sizes` list is unimportant because the code will sort it into non-
   increasing order before any other processing is done.
   """
   if not isinstance(balls, int):
      raise TypeError("balls must be a non-negative integer.")
   if balls < 0:
      raise ValueError("balls must be a non-negative integer.")

   if not isinstance(box_sizes,list):
      raise ValueError("box_sizes must be a non-empty list.")

   capacity= 0
   for size in box_sizes:
      if not isinstance(size, int):
          raise TypeError("box_sizes must contain only positive integers.")
      if size < 1:
          raise ValueError("box_sizes must contain only positive integers.")
      capacity+= size

   if capacity < balls:
      raise ValueError("The total capacity of the boxes is less than the "
        "number of balls to be distributed.")

   for unlabeled_dist in unlabeled_balls_in_unlabeled_boxes(balls, box_sizes):
      for labeled_dist in \
        m_way_unordered_combinations(balls, unlabeled_dist):
         yield labeled_dist

# end def labeled_balls_in_unlabeled_boxes(balls, box_sizes)


def labeled_balls_in_labeled_boxes(balls, box_sizes):
   """
   OVERVIEW

   This function returns a generator that produces all distinct distributions of
   distinguishable balls among distinguishable boxes, with specified box sizes
   (capacities).  This is a generalization of the most common formulation of the
   problem, where each box is sufficiently large to accommodate all of the
   balls.


   CONSTRUCTOR INPUTS

   n: the number of balls

   box_sizes: This argument is a list of length 1 or greater.  The length of
   the list corresponds to the number of boxes.  `box_sizes[i]` is a positive
   integer that specifies the maximum capacity of the ith box.  If
   `box_sizes[i]` equals `n` (or greater), then the ith box can accommodate all
   `n` balls and thus effectively has unlimited capacity.


   EXAMPLE

   Issue the following statement issued at the IPython prompt:

   list(labeled_balls_in_labeled_boxes(3,[2,2]))

   The output is as follows:

   [((0, 1), (2,)),
    ((0, 2), (1,)),
    ((1, 2), (0,)),
    ((0,), (1, 2)),
    ((1,), (0, 2)),
    ((2,), (0, 1))]
   """
   if not isinstance(balls, int):
      raise TypeError("balls must be a non-negative integer.")
   if balls < 0:
      raise ValueError("balls must be a non-negative integer.")

   if not isinstance(box_sizes,list):
      raise ValueError("box_sizes must be a non-empty list.")

   capacity= 0
   for size in box_sizes:
      if not isinstance(size, int):
          raise TypeError("box_sizes must contain only positive integers.")
      if size < 1:
          raise ValueError("box_sizes must contain only positive integers.")
      capacity+= size

   if capacity < balls:
      raise ValueError("The total capacity of the boxes is less than the "
        "number of balls to be distributed.")

   for unlabeled_dist in unlabeled_balls_in_labeled_boxes(balls, box_sizes):
      for labeled_dist in m_way_ordered_combinations(balls, unlabeled_dist):
         yield labeled_dist

# end def labeled_balls_in_labeled_boxes(balls, box_sizes)


def partitions(n):
   """
   'In number theory and combinatorics, a partition of a positive integer n,
   also called an integer partition, is a way of writing n as a sum of positive
   integers.  Two sums that differ only in the order of their summands are
   considered to be the same partition.'  We can trivially generate all
   partitions of an integer using `unlabeled_balls_in_unlabeled_boxes`.  The
   quote is from http://en.wikipedia.org/wiki/Partition_(number_theory) .
   """

   _partitions= unlabeled_balls_in_unlabeled_boxes(n, n*[n])

   for _partition in _partitions:
      yield tuple([p for p in _partition if p])


def partitions2(n):
   """
   This function generates integer partitions; it was written by David Eppstein.
   """
   # base case of recursion: zero is the sum of the empty list
   if n == 0:
      yield []
      return

   # modify partitions of n-1 to form partitions of n
   for p in partitions2(n-1):
      yield [1] + p
      if p and (len(p) < 2 or p[1] > p[0]):
         yield [p[0] + 1] + p[1:]
