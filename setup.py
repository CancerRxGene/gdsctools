# -*- coding: utf-8 -*-
__revision__ = "$Id$"
import sys
import os
from setuptools import setup, find_packages
import glob


_MAJOR               = 0
_MINOR               = 18
_MICRO               = 0
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
    'description':'Set of tools and pipelines to analyse GDSC data (cancerrxgene.org)' ,
    'platforms' : ['Linux', 'Unix', 'MacOsX', 'Windows'],
    'keywords' : ['gdsc', 'drug response', 'anova', 'genomic features'],
    'classifiers' : [
          'Development Status :: 1 - Planning',
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: BSD License',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
          'Topic :: Software Development :: Libraries :: Python Modules',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Information Analysis',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Physics']
    }


with open('README.rst') as f:
    readme = f.read()

from distutils.core import setup, Extension


on_rtd = os.environ.get('READTHEDOCS', None) == 'True'

# on rtd, pandas is 1.3.1 cannot use something better for now (march 2016)

# jinja <= 2.9.1 because there is a py36 dependencies otherwise
if on_rtd is True:  # only import and set the theme if we're building docs
    install_requires = ['colormap', 'easydev', "sphinx-gallery", "numpydoc",
                        "colorlog"]
else:
    install_requires = ['numpy', "numexpr", 'matplotlib>=1.4.3',
        'pandas>=0.19', 'easydev>=0.9.34', 'scipy', "colormap>=1.0",
        'jinja2', 'statsmodels', "scikit-learn", "reports>=0.3.0",
        "biokit>=0.4", "colorlog", "xlrd"],


packages = find_packages()
packages = [this for this in packages if this not in ['test']]


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
    packages = packages,

    # here below '': pattern means include that pattern in all packages
    # so '' :['README.rst'] will include all README.rst recursively
    package_data = {
        '': ["*.js", "*txt", "*.csv", "*tsv", "*.gz"],
        'gdsctools.pipelines': ['*.rules', '*.yaml'],
        'gdsctools.data.css': ['*.css'],
        'gdsctools.data.images' : ['*.png', '*.ico'],
        'gdsctools.data.templates' : ['*.html'],
        'gdsctools.data.javascript' : ['*.js'],
        },

    # comment the requirements otherwise RTD fails
    # but we then need a requirements.txt file !
    install_requires = install_requires,

    zip_safe=False,
    entry_points = {
        'console_scripts': [
        'gdsctools_anova=gdsctools.scripts.anova:main',
        'gdsctools_regression=gdsctools.scripts.regression:main']
    },
)


