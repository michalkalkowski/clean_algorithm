import os
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="clean",
    version="0.0.1",
    description=("An implementation of the CLEAN algorithm for extracting"
                 + "individual components from a multi-component signal"),
    author="Michal K Kalkowski",
    author_email="m.kalkowski@imperial.ac.uk",
    packages=find_packages(exclude=['data', 'references', 'output', 'notebooks']),
    long_description=read('README.md'),
    license='MIT',
)
