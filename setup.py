#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pyGTF
from os.path import join, abspath, dirname
from setuptools import setup, find_packages

with open(join(abspath(dirname(__file__)), 'README.md')) as fi:
    long_description = fi.read()

setup(
    name = pyGTF.__package__,
    version = pyGTF.__version__,
    license = pyGTF.__licence__,
    url = pyGTF.__url__,
    author = pyGTF.__author__,
    description = pyGTF.__description__,
    long_description=long_description,
    long_description_content_type='text/markdown',

    keywords='Fastx, GTF, GFF, RefSeq',
    py_modules=['pyGTF'],
    # packages = find_packages(),
    include_package_data = False,
    # platforms = 'any',
    # install_requires = None,
)