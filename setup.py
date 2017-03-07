#!/usr/bin/env python
# vim:set ts=4 sts=4 sw=4 et:

from setuptools import setup


setup(
    name='km3astro',
    version='0.2',
    description='Astro Utils',
    url='http://git.km3net.de/moritz/km3astro',
    author='Moritz Lotze',
    author_email='mlotze@km3net.de',
    license='BSD-3',
    packages=['km3astro', ],
    install_requires=[
        'astropy',
        'numpy',
        'pandas',
    ]
)
