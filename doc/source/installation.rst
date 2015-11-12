
.. index:: installation, pip, anaconda
.. _installation:

Installation
================

**GDSCTools** is written in Python and depends on a set of established scientific libraries such as `Matplotlib <http://matplotlib.org/>`_ (visualisation) and `Pandas <http://pandas.pydata.org/>`_ (data manipulation) to cite just a few. **GDSCTools** also makes use of the Python ecosystem to provide a robust software. We hope that the installation will work out of the box. Would you have any trouble, please see the :ref:`faqs` or fill an `issue/ticket <https://github.com/CancerRxGene/gdsctools/issues>`_ on github.


Existing Python user and developer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Releases of **GDSCTools** are available on `Pypi <https://pypi.python.org/pypi/gdsctools/0.2.0>`_. Therefore **GDSCTools**
can be installed in a shell using the **pip** command::

    pip install gdsctools

Dependencies (e.g., Pandas, Matplotlib) should be taken care of automatically.

If you are new to Python and IPython notebook
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are not familiar with Python, or have issues with compilation, or do not have root access, we would recommend to use the `Anaconda <https://www.continuum.io/downloads>`_ solution. 

Anaconda is a free Python distribution. It includes most popular Python packages for science and data analysis. Anaconda will install the software required by  **GDSCTools**. Since it does not require root access, it should not interfere with your system.

Please, visit the `Anaconda <https://www.continuum.io/downloads>`_ webiste and follow the instructions. You may need to choose between 2 versions of Python (2.X or 3.X). We do not recommend any specific version but note that **GDSCTools** is tested under Python 2.7, 3.3 and 3.4 (not 3.5 for now) and is mainly developed under Python 2.7.



You have now two options. You can either (as an **end-user**) type the following command that will install the latest release of **GDSCTools**::

    pip install gdsctools

or you can (as **developer**) get the latest source code from github and install **GDSCTools** from source as follows (in a shell) ::

    # go in a working directory and type:
    git clone https://github.com/CancerRxGene/gdsctools
    cd gdsctools
    python setup.py install


Testing your installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You should now be ready to use **GDSCTools**. A good test is to check
that the following executable is available. In a shell, type::

    gdsctools_anova --test

or ::

    gdsctools_anova --help

or for developers, start an IPython shell, and type e.g.::

    from gdsctools import *
    an = ANOVA(ic50_test)

Please, go to the next section for a :ref:`quickstart` session.



