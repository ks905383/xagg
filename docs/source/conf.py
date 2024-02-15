# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

import os
import sys
#sys.path.insert(0,os.path.abspath('../../xagg/'))
sys.path.insert(0,os.path.abspath('../..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'xagg'
copyright = '2021-2024, Kevin Schwarzwald'
author = 'Kevin Schwarzwald'
release = '0.3.2.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc','sphinx.ext.viewcode','numpydoc','nbsphinx','sphinx_rtd_theme',]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
# (allegedly, adding '_build' (or 'build'?) and "**.ipynb_checkpoints' gets
# rid of a bunch of "Pygments lexer name 'ipython3' not known" warnings...)
exclude_patterns = ['build', '_build','html','_html', 'Thumbs.db', '.DS_Store','**.ipynb_checkpoints']

# Don't actually build environment for docs (this was the reason why
# the API in the docs was coming up empty)
autodoc_mock_imports = ['pytest','xarray','numpy','scipy','shapely',
'netcdf4','pandas','geopandas','cartopy','xesmf','matplotlib','cmocean','cf_xarray','pytables']


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
