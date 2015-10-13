# -*- coding: utf-8 -*-
__revision__ = "$Id$"
import sys
import os
from setuptools import setup, find_packages
import glob



_MAJOR               = 0
_MINOR               = 0
_MICRO               = 2
version              = '%d.%d.%d' % (_MAJOR, _MINOR, _MICRO)
release              = '%d.%d' % (_MAJOR, _MINOR)

metainfo = {
    'authors': {
        'Cokelaer':('Thomas Cokelaer','cokelaer@gmail.com'),
        },
    'version': version,
    'license' : 'BSD',
    'download_url' : ['http://pypi.python.org/pypi/gdsctools'],
    'url' : ['http://pypi.python.org/pypi/gdsctools'],
    'description':'Access to Biological Web Services from Python' ,
    'platforms' : ['Linux', 'Unix', 'MacOsX', 'Windows'],
    'keywords' : ['corrplot', 'heatmap'],
    'classifiers' : [
          'Development Status :: 1 - Planning',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Topic :: Software Development :: Libraries :: Python Modules',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Information Analysis',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Physics']
    }



with open('README.rst') as f:
    readme = f.read()

from distutils.core import setup, Extension


setup(
    name             = 'gdsctools',
    version          = version,
    maintainer       = metainfo['authors']['Cokelaer'][0],
    maintainer_email = metainfo['authors']['Cokelaer'][1],
    author           = metainfo['authors']['Cokelaer'][0],
    author_email     = metainfo['authors']['Cokelaer'][1],
    long_description = readme,
    keywords         = metainfo['keywords'],
    description = metainfo['description'],
    license          = metainfo['license'],
    platforms        = metainfo['platforms'],
    url              = metainfo['url'],      
    download_url     = metainfo['download_url'],
    classifiers      = metainfo['classifiers'],

    zip_safe=False,
    packages = find_packages(),
    install_requires = [ 'numpy', 'pandas', 'easydev', 'scipy', 'biokit'],

    )


