.. _Introduction:

************
Introduction
************

The :mod:`SphinxReport` module is an extension for :mod:`sphinx`
that provides facilities for data retrieval and data rendering
within reStructured text. 

.. _Installation:

Requirements
************

:mod:`SphinxReport` requires the following software

   * `Python <http://www.python.org>`_ (2.5.2 or higher) 
   * `SQLAlchemy <http://www.sqlalchemy.org/>`_ (0.4.8 or higher)
   * `matplotlib <http://matplotlib.sourceforge.net/>`_ (0.98.1 or higher)
   * `sphinx <http://sphinx.pocoo.org/>`_ (0.5-1 or higher)

Installation
************

In order to install the extension, download the latest sources from *TODO* and unpack:

   tar -xvzf sphinx-report.tar.gz

Alternatively, check out the the latest code from the subversion repository::

   svn checkout http://sphinx-report.googlecode.com/svn/trunk/ sphinx-report-read-only

To install, type::

   python setup.py build
   python setup.py install

First steps
***********

To get running quickly, use the the :file:`sphinxreport-quickstart`` to
create a skeleton project in the directory ``newproject``::

   sphinxreport-quickstart -d newproject

Enter ``newproject`` and build the skeleton report::

   make html

Open :file:`newproject/_build/html/index.html` in your browser 
to view the skeleton documentation. 

At this stage you can review :ref:`Configuration` options
in the file :file:`conf.py` and then start adding content
to your report. See the :ref:`Tutorials` on how to do this.

.. _Configuration:

Configuration
*************

Sphinx and the :mod:`SphinxReport` extension read configuration details
from the file :file:`conf.py` at the top-level of the installation. In order
to use the extension, add the following entries to the variable :data:`extensions`::

   extensions.extend( ['SphinxReport.inheritance_diagram',
              'SphinxReport.only_directives',
              'SphinxReport.render_directive ] )

Further variables can be addded to :file:`conf.py` to customize the extension. The
variables are:

This file containing python code sets mostly :mod:`sphinx` configuration 
options (see the `sphinx documentation <http://sphinx.pocoo.org/config.html>`_
for more information), but some sphinxreport options as well.

































