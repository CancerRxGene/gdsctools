############################################################################
GDSCtools
############################################################################

.. contents::

Overview
---------

**gdsctools** contains utilities used to find significant associations between
drug response and genomic features.



Installation
---------------

Releases of **gdsctools** are available on 
`Pypi <https://pypi.python.org/pypi/gdsctools/0.2.0>`_. Therefore **gdsctools**
can be installed in a shell using **pip**::

    pip install gdsctools

There are a few dependencies (e.g., Pandas, matplotlib) that will need to have
compilers. You may also need root access. If this is an issue or you are not a developer and just wish to install and try gdsctools, we would recommend to 
use the `Anaconda <https://www.continuum.io/downloads>`_. This page will
guide you through a simple installation procedure. It will install the binary
packages required to use gdsctools. Once done, type the command above (pip
install gdsctools), or get the latest version from github. In a shell::

    git clone https://github.com/CancerRxGene/gdsctools
    cd gdsctools
    python setup.py install
    
You should be ready now to use gdsctools. A good test is to test that the
executable is available. In a shell, type::

    gdsctools_anova --help 




Issues
-----------

Please fill bug report in https://github.com/CancerRxGene/gdsctools/issues


Contributions
----------------

Please join https://github.com/CancerRxGene/gdsctools

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
    ChangeLog.rst


.. toctree::
    :hidden::

    glossary
