.. _Tutorial3:

==========================
 Tutorial 3: Using slices
==========================

This tutorial introduces ``slices``. ``Slices`` are subsets
of data within tracks. To illustrate their use we will count
word sizes of ``.py`` and ``.rst`` files.

***********************
Creating a word counter
***********************

Create the file :file:`Tutorial3.py` in the :file:`python` subdirectory and add 
the following code::

    from SphinxReport.DataTypes import *
    from SphinxReport.Tracker import *

    import os

    class WordCounter(Tracker):
	"""Counting word size."""

	def getTracks( self ):
	    return ( "all", ".py", ".rst" )

	@returnSingleColumnData
	def __call__(self, track, slice = None ):
	    word_sizes = []

	    if track == "all" or track == None:
		tracks = [ ".py", ".rst" ]
	    else:
		tracks = [track]

	    for root, dirs, files in os.walk('.'):
		for f in files:
		    fn, ext = os.path.splitext( f )
		    if ext not in tracks: continue
		    infile = open(os.path.join( root, f),"r")
		    word_sizes.extend( [ len(word) for word in re.split("\s+", "".join(infile.readlines())) ] )
		    infile.close()

	    return word_sizes

This Tracker counts word sizes in ``.py`` and ``.rst`` files in the current directory and returns
a list of word sizes (:class:`DataTypes.SingleColumnData`).

Testing this data source::

   sphinxreport-test -t WordCounter -r histogram-plot -o range=0,100,1"

should produce a line plot. Note the ``-o`` option to :file:`sphinxreport-test` in order to pass 
parameters. In the above example, the histogram is computed in the range from 0 to 100 in steps of size 1.

*************
Adding slices
*************

Let us say we want to examine the length of words starting in with vocals (``AEIOU``) compared to
those starting with consonants. One possibility is to extend the :meth:`getTracks()` method to
return new tracks like .py_vocals, .py_consonants, etc. A better way to do this is to use
slices. Add the following code to :file:`Tutorial3.py`::

    class WordCounterWithSlices(Tracker):
	"""Counting word size."""

	def getTracks( self ):
	    return ( "all", ".py", ".rst" )

	def getSlices( self, subset = None ):
	    return ( "all", "vocals", "consonants")

	@returnSingleColumnData
	def __call__(self, track, slice = None ):
	    word_sizes = []

	    if track == "all" or track == None:
		tracks = [ ".py", ".rst" ]
	    else:
		tracks = [track]

	    if slice == "all" or slice == None:
		test_f = lambda x: True
	    elif slice == "vocals":
		test_f = lambda x: x[0].upper() in "AEIOU"
	    elif slice == "consonants":
		test_f = lambda x: x[0].upper() not in "BCDFGHJKLMNPQRSTVWXYZ"

	    for root, dirs, files in os.walk('.'):
		for f in files:
		    fn, ext = os.path.splitext( f )
		    if ext not in tracks: continue
		    infile = open(os.path.join( root, f),"r")
		    words = [ w for w in re.split("\s+", "".join(infile.readlines())) if len(w) > 0]
		    word_sizes.extend( [ len(w) for w in words if test_f(w)] )
		    infile.close()

	    return word_sizes

This counter again counts word sizes in ``.py`` and ``.rst`` files, but collects counts separately
for words starting with vocals and consonants.

Testing the data source::

   sphinxreport-test -t WordCounterWithSlices -r histogram-plot -o range=0,1,100

will now produce three plots, one for each slice. Per default, plots are grouped by ``slice``, but the grouping
can be changed using the option ``groupby=track``::

   sphinxreport-test -t WordCounterWithSlices -r histogram-plot -o range=0,1,100 -o groubpy=track

Again, three plots are created, but this time there is one plot per ``track``. 

****************************************************
Inserting the graphs in a restructured text document
****************************************************

We can now add these three plots into a restructured text document using
a single report directive block::

    ==========
    Tutorial 3
    ==========

    Using slices

    .. report:: Tutorial3.WordCounterWithSlices
       :render: histogram-plot
       :range: 0,100,1

       Word sizes in .py and .rst files grouped by slice

Additionally you can add the plots grouped by tracks::

    .. report:: Tutorial3.WordCounterWithSlices
       :render: histogram-plot
       :range: 0,100,1
       :groupby: track

       Word sizes in .py and .rst files grouped
       by track.

More fine grained control is possible. The following only shows a single plot::

    .. report:: Tutorial3.WordCounterWithSlices
       :render: histogram-plot
       :range: 0,100,1
       :tracks: .py,.rst
       :slices: vocals

       Word sizes of words starting with vocals in .py and
       .rst files.

See :ref:`Tutorial3Demo` to check how the result should look like.


