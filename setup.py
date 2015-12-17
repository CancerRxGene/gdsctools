# -*- coding: utf-8 -*-
__revision__ = "$Id$"
import sys
import os
from setuptools import setup, find_packages
import glob


_MAJOR               = 0
_MINOR               = 10
_MICRO               = 1 
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


# on rtd, pandas is 1.3.1 cannot use something better for now
if on_rtd is True:  # only import and set the theme if we're building docs
    install_requires = ['colormap', 'easydev']
else:
    install_requires = ['numpy', 'matplotlib>=1.4.3',
        'pandas>=0.16.2', 'easydev>=0.9.5', 'scipy', "colormap>=0.9.7",
        'mpld3', 'jinja2', 'statsmodels'],



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
    include_package_data = True,

    # tells discutils extra packages are under share/data
    package_dir={'share.data': 'share/data'},

    # here below '': pattern means include that pattern in all packages
    # so '' :['README.rst'] will include all README.rst recursively
    package_data = {
        'share.data' : ['*.css', '*.js', '*.txt', '*.csv', '*.tsv', 
            '*.gz', 'README.rst'],
        '' : ['README.rst'],
        'share.data.images' : ['*.png'],
        'share.data.templates' : ['*.html'],
        'share.data.javascript' : ['*.js'],
        },

    # comment the requirements otherwise RTD fails
    # but we then need a requirements.txt file !
    install_requires = install_requires,
    entry_points = {
        'console_scripts': [
        'gdsctools_anova=gdsctools.pipelines:anova_pipeline',]
        },


    )


