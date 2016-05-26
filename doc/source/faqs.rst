.. _faqs:

FAQS
=======

.. contents::

Installation
-----------------

The :ref:`installation` should guide you to install **GDSCTools** but would you
have any trouble, please let us know and send a ticket 
on `gdsctools github <https://github.com/CancerRxGene/gdsctools/issues>`_.


Usage
--------

Any trouble using **gdsctools_anova** executable, please let us know: send a
on `gdsctools github <https://github.com/CancerRxGene/gdsctools/issues>`_.


I cannot see any plots or figures
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First solution is to import the show function from pylab::

    from pylab import show
    show()

The figures should appear now but you will have to close them to get back to the
shell. To get a better interaction, exit from the shell and start a new one
typing::

    ipython --pylab

Under MAC, the default backend for pylab may cause some trouble with this
message::

    RuntimeError: CGContextRef is NULL

If so, try another backend::

    ipython --pylab=qt


About Python
---------------

Shall I use Python 2 or Python 3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GDSCTools is currently compatible for the version 2.7 but also 3.3, and 3.4 .
We do not have any strong recommendation except to use whatever version is already installed on your system. We will not maintain the code for version below 2.7 although the code may work for those versions. 

How do I get documentation on Python?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The standard documentation for the current stable version of Python is available
at https://docs.python.org/3/. PDF, plain text, and downloadable HTML versions
are also available at https://docs.python.org/3/download.html.

I’ve never used Python before. Is there a Python tutorial?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
There are numerous tutorials and books available. The standard documentation
includes `The Python Tutorial <https://docs.python.org/3/tutorial/index.html#tutorial-index>`_. 


