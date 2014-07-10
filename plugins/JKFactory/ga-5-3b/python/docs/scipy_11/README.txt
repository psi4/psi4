SciPy Proceedings
=================

Paper Format
------------

General Guidelines
``````````````````
- All figures and tables should have captions.
- License conditions on images and figures must be respected
  (Creative Commons, etc.).
- Code snippets should be formatted to fit inside a single column
  without overflow.
- Try to use as little custom LaTeX markup as possible.
- The paper abstract should be a single paragraph.

Authors and affiliations
````````````````````````
Define the fields in the beginning of the paper::

  :author: My Name
  :email: myname@myplace.com
  :institution: Some University

  :author: Author Two
  :email: two@myplace.com
  :institution: Some University

Other markup
------------
Please refer to the example paper in ``papers/00_vanderwalt`` for
examples of how to:

 - Label figures, equations and tables
 - Use math markup
 - Include code snippets

Build Process
-------------
::

  ./make_paper.sh papers/my_paper_dir

Building the entire Proceedings
-------------------------------
::

  ./make_all.sh
  ./build_index.py
  ./make_all.sh
  ./concat_proceedings_pdf.sh

Requirements
------------
 - IEEETran and AMSmath LaTeX classes
 - **Latest** docutils (development version, they haven't released in years)
 - Pygments for code highlighting

