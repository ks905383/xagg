from setuptools import setup
import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
	name="xagg",
    version="0.2.6",
    author="Kevin Schwarzwald",
    author_email="kschwarzwald@iri.columbia.edu",
    description="Aggregating raster data over polygons",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ks905383/xagg",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy',
        'scipy',
        'xarray',
        'pandas',
        'netcdf4',
        'geopandas',
        'shapely',
        'xesmf>=0.5.2',
        'cf_xarray>=0.5.1',
        'esmpy>=8.1.0' 
    ],
)
