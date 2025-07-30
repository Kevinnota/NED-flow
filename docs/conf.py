# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'NED-flow'
copyright = '2025, Kevin Nota, Benjamin Vernot'
author = 'Kevin Nota, Benjamin Vernot'
release = 'stable-beta_1'

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx_rtd_dark_mode',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'

# -- Dark mode toggle --------------------------------------------------------

default_dark_mode = True

# -- Static files (custom CSS, images, etc.) ---------------------------------

html_static_path = ['_static']
