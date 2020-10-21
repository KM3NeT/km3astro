#!/usr/bin/env python
# vim:set ts=4 sts=4 sw=4 et:

from setuptools import setup
from km3astro import __version__

with open('requirements.txt') as fobj:
    requirements = [l.strip() for l in fobj.readlines()]


setup(
    name='km3astro',
    version=__version__,
    description='Astro Utils',
    url='http://git.km3net.de/km3py/km3astro',
    author='Tamas Gal and Moritz Lotze',
    author_email='tgal@km3net.de',
    license='BSD-3',
    packages=['km3astro', ],
    install_requires=requirements,
    python_requires='>=3.5',
    include_package_data=True,
)
