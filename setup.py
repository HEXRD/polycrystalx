import os
import sys

from setuptools import setup, find_packages

entry_points = {
    'console_scripts': [
        "pxx_job = polycrystalx.scripts.run_job:main",
        "pxx_suite = polycrystalx.scripts.run_suite:main"
    ]
}

setup(
    name = 'polycrystalx',
    author = 'Donald E. Boyce',
    author_email = 'donald.e.boyce@gmail.com',
    description = 'polycrystal material modeling',
    classifiers = [
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering',
        ],
    packages = find_packages(),
    entry_points = entry_points,
    )
install_reqs = [
    'dolfinx',
    'ufl',
    'numpy',
    'scipy',
    'h5py',
    'pytest',
]
