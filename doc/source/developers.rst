For developers
=================

The GDSCTools test suite
--------------------------

* The entire GDSCTools packages is tested. Tests are in the ./test/gdsctools 
  directory of the source code. 

* As of version 0.3, about 80% of the code is covered. 

* If you add a new module, please add a corresponding test file in
  ./test/gdsctools

Documentation
----------------

The documentation is based on Sphinx. This means that all code is documentated
using the REST syntax. Docstring are added in classes and
functions as much as possible with code examples. 

See for example the file in gdsctools/anova.py and in particular the class
ANOVA. 

In addition, in the ./doc directory there are a set of files usng REST syntax.
If you go to that directory and type::

    make html

The entire documentation included Tutorial and Developer guide will be parsed
and interpreted. The resulting html documentation can then be found in
doc/build/html. 

.. note:: the command above will work only if gdsctools source code is 
    available and installed in the shell where the command is executed.


The documentation can then be uploaded on pypi. Go to the directory were is the
file **setup.py** and type::

    python setup.py upload_docs


    

