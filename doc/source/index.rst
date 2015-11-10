############################################################################
GDSCtools
############################################################################

.. contents::

Overview
---------

**GDSCTools** contains utilities used to find significant associations between
drug response and genomic features.



Installation
---------------

You know Python and/or **pip** utility ?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Releases of **GDSCTools** are available on `Pypi <https://pypi.python.org/pypi/gdsctools/0.2.0>`_. Therefore **GDSCTools**
can be installed in a shell using **pip**::

    pip install gdsctools

Dependencies (e.g., Pandas, matplotlib) should be taken care of automatically.

You don't know Python  or have compilation issues ?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you are not familiar with Python, or have issues with compilation, or do not have root access, we would recommend to use the `Anaconda <https://www.continuum.io/downloads>`_ solution. It will install binary packages required to install **GDSCTools** itself; since it does not require root access, it should not interfere with your system.

Please, visit the link above and install Anaconda following the
instructions.

Once you have installed Anaconda, open a new shell.

As a **developer**, you can get the latest source code from github and install GDSCTools from source as follows (in a shell) ::

    # go in a working directory and type:
    git clone https://github.com/CancerRxGene/gdsctools
    cd gdsctools
    python setup.py install

As an **end-user**, open a shell and type. This will install the latest official release of **GDSCTools**::

    pip install gdsctools

Testing your installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You should now be ready to use **GDSCTools**. A good test is to check
that the following executable is available. In a shell, type::

    gdsctools_anova --test

or ::

    gdsctools_anova --help

or for developers, starts an IPython shell, and type e.g.::

    from gdsctools import *
    an = ANOVA(ic50_test)

Please, see the quickstart session for the usage.


User Guide
------------

.. toctree::
    :maxdepth: 1

    userguide.rst

Developers Guide
------------------

.. toctree::
    :maxdepth: 1

    references.rst
    developers.rst
    ChangeLog.rst

Issues
-----------

Please fill bug report in https://github.com/CancerRxGene/gdsctools/issues

Contributions
----------------

Please join https://github.com/CancerRxGene/gdsctools

.. toctree::
    :hidden::

    glossary
