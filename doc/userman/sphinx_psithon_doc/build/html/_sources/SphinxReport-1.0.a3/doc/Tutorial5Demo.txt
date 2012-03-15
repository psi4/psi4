.. _Tutorial5Demo:

==========
Tutorial 5
==========

The simplest plot:

.. report:: Tutorial4.ExpressionLevelWithSlices
   :render: histogram-plot
   :range: 0,100,4

   Expression level in house-keeping and regulatory genes
   in two experiments.

A customized plot:

.. report:: Tutorial4.ExpressionLevelWithSlices
   :render: histogram-plot
   :range: 0,100,4
   :xtitle: expression level
   :groupby: all
   :as-lines:

   Expression level in house-keeping and regulatory genes
   in two experiments.

The same data in tabular form:

.. report:: Tutorial4.ExpressionLevelWithSlices
   :render: stats

   Expression level in house-keeping and regulatory genes
   in two experiments.

as box plot:

.. report:: Tutorial4.ExpressionLevelWithSlices
   :render: box-plot
   :groupby: all
   :ytitle: expression level

   Expression level in house-keeping and regulatory genes
   in two experiments.

or as literal histogram:

.. report:: Tutorial4.ExpressionLevelWithSlices
   :render: histogram
   :range: 0,100,4

   Expression level in house-keeping and regulatory genes
   in two experiments.
