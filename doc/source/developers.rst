For developers
=================

The GDSCTools test suite
--------------------------

The GDSCTools package has a large set of tests in the ./test/gdsctools directory.  In version 0.10, about 80% of the functionalities are covered. If you add a 
new module or function, please add a corresponding test file in
  ./test/gdsctools directory.

The test suite uses the **nosetests** utility and is integrated within the
  Travis continuous integration (see below).

To run the tests locally, you will need nose and coverage packages::

    pip install nose coverage

Then, in the directory where is the **setup.py** file, type ::

    python setup.py nosetests 

You may also go directly in the ./tests/gdsctools directory  but some tests may
required extra files or permission.

Documentation
----------------

The documentation is based on Sphinx. This means that all code is documentated
using the REST syntax. Docstring are added in classes and
functions as much as possible with code examples. 

See for example the file in gdsctools/anova.py and in particular the class
ANOVA. 

In addition, in the ./doc directory there are a set of files encoded using 
the REST syntax. If you go to that directory and type::

    make html

then the entire documentation included Tutorial and Developer guide 
will be parsed and interpreted. The resulting html documentation can then be found in doc/build/html. For the Sphinx documentation to be generated, the **gdsctools** package must be installed.

.. note:: the command above will work only if gdsctools source code is 
    available and installed in the shell where the command is executed.


The documentation can then be uploaded on pypi. Go to the directory were is the
file **setup.py** and type::

    python setup.py upload_docs


Continuous Integration
---------------------------

This is based on Travis. The reports should be available here: https://travis-ci.org/CancerRxGene/gdsctools and the status is also reported in the main github page (https://github.com/CancerRxGene/gdsctools) as an icon (**build**)  that will be green or red depending  on the status of the build within Travis. 

