from setuptools import setup, find_packages
import sys

install_requires = ["numpy>=1.17",
                    "matplotlib>=3.1.1",
                    "geopandas>=0.6.1",
                    "netcdf4>=1.4.2",
                    "ipython>=7.9.0",
                    "scipy>=1.3.1",
                    "rasterio>=1.1.0",
                    "jupyterlab>=3.0.9"]

setup(name='synthetic_offsets',
      version='0.0.3',
      description='Simulation of offsets',
      author='Andy Howell',
      author_email='andrew.howell@canterbury.ac.nz',
      packages=find_packages()
      )
