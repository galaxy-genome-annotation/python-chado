# http://bugs.python.org/issue15881#msg170215
from setuptools import setup, find_packages

from os import listdir

setup(
    name="chado",
    version='2.1',
    description="Chado library",
    author="Anthony Bretaudeau",
    author_email="anthony.bretaudeau@inra.fr",
    url="https://github.com/galaxy-genome-annotation/python-chado",
    install_requires=['sqlalchemy', 'psycopg2', 'biopython==1.67', 'bcbio-gff', 'wrapt', 'click', 'pyyaml', 'future'],
    packages=find_packages(),
    license='MIT',
    platforms="Posix; MacOS X; Windows",
    entry_points='''
        [console_scripts]
        chakin=chakin.cli:chakin
    ''',
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.6",
        "Programming Language :: Python :: 2.7"
    ])
