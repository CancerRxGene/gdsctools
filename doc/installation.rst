
.. index:: installation, pip, anaconda
.. _installation:

Installation
================

You know **bioconda**, then you can install **GDSCTools** as follows::


    conda install gdsctools

Otherwise, please read this section in details.

**GDSCTools** is written in Python and depends on a set of established scientific libraries such as `Matplotlib <http://matplotlib.org/>`_ (visualisation) and `Pandas <http://pandas.pydata.org/>`_ (data manipulation) to cite just a few. We post releases on the `Python repository  <https://pypi.python.org/pypi/gdsctools>`_ and make use of the Python ecosystem to provide a robust software. Would you have any troubles, please see the :ref:`faqs` or fill an `issue/ticket <https://github.com/CancerRxGene/gdsctools/issues>`_ on github.


Installation using Conda
----------------------------

Install conda executable
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In practice, we do use `Anaconda <https://conda.readthedocs.io/>`_ . We
recommend to
install **conda** executable via the manual installer (`download
<https//continuum.io/downloads>_`). 
You may have the choice
between Python 2 and 3. We recommend to choose a Python version 3.

Add conda channels
~~~~~~~~~~~~~~~~~~~~~~~~~~

When you want to install a new package, you have to use this type of syntax::

    conda install ipython

where **ipython** is the package you wish to install. Note that by default,
**conda** looks on the official Anaconda website (channel). However, there are
many channels available. We will use the **bioconda** channel. To use it, type
these commands (once for all)::

    conda config --add channels r
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda

.. warning:: **it is important to add them in this order**, as mentionned on
   bioconda webpage (https://bioconda.github.io/).


Create an environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once **conda** is installed, open a new shell.
Although this is not required strictly speaking, we would
recomment to create an environment dedicated to Sequana. This environment can
later be removed without affecting your system or conda installation. A
**conda** environment is nothing else than a directory and can be created as
follows::

    conda create --name gdsctools_env python=3.5

Then, since you may have several environments, you must activate the **gdsctools**
environment itself::

    source activate gdsctools_env


Install gdsctools via conda (bioconda)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finally, just type::

    conda install gdsctools

This should install most of the required dependencies. However, you may need to
install more packages depending on the pipeline used.


Installation using **pip**
---------------------------------------

Python users and developers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Releases of **GDSCTools** are available on `Pypi <https://pypi.python.org/pypi/gdsctools/>`_ so **GDSCTools** can be installed in a :term:`Terminal` using the **pip** command::

    pip install gdsctools

Dependencies (e.g., Pandas, Matplotlib) should be taken care of automatically.
However, compilation are needed making the installation process much longer. 


Installation from source
--------------------------------

:: 


    git clone https://github.com/CancerRxGene/gdsctools
    cd github_gdsctools
    git pull
    python setup.py install


Testing your installation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You should now be ready to use **GDSCTools**. A good test is to check
that the following executable is available. In a shell, type::

    gdsctools_anova --test

or ::

    gdsctools_anova --help

or for developers, start an IPython shell, and type for example::

    from gdsctools import *
    an = ANOVA(ic50_test)

Please, go to the next section for a :ref:`quickstart` session.

Open an IPython shell
~~~~~~~~~~~~~~~~~~~~~~~~~

Under Windows, got to All Programs-->Anaconda -->Anaconda Prompt.

A shell will be opened where you can type **ipython** command.

Or alternatively, under Windows, got to All Programs-->Anaconda -->IPython

Notes for windows/mac/linux
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Anaconda method was tested successfully on the following systems: MAC,
Windows 7 Pack1, Fedora 19 (Nov 2015) with version 0.16.5 of GDSCTools.

Under Windows, an error was raised due to scipy. This was fixed by typing::

    conda remove scipy scikit-learn -y
    conda install scipy scikit-learn -y

https://github.com/scikit-learn/scikit-learn/issues/4830
