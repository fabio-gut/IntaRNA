#!/usr/bin/env python


from setuptools import setup

setup(
    name='IntaRNApvalue',
    version='0.1',
    description='Calculates pvalues to IntaRNA scores',
    author='Fabio Gutmann',
    author_email='fabio.gutmann@jupiter.uni-freiburg.de',
    python_requires='>=3.6.0',
    url='https://github.com/fabio-gut/IntaRNA/tree/pvalue/',
    packages=['IntaRNApvalue'],
    # py_modules=['IntaRNApvalue'],
    install_requires=['scipy', 'numpy'],
    include_package_data=True,
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy'
    ]
)
