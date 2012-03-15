.. _Tutorial3Demo:

==========
Tutorial 3
==========

Using slices

.. report:: Tutorial3.WordCounterWithSlices
   :render: histogram-plot
   :range: 0,100,1

   Word sizes in .py and .rst files. 

.. report:: Tutorial3.WordCounterWithSlices
   :render: histogram-plot
   :range: 0,100,1
   :groupby: track

   Word sizes in .py and .rst files. 

.. report:: Tutorial3.WordCounterWithSlices
   :render: histogram-plot
   :range: 0,100,1
   :tracks: .py,.rst
   :slices: vocals

   Word sizes of words starting with vocals in .py and
   .rst files.
