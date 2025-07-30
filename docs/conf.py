# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------

project = 'NED-flow'
copyright = '2025, Kevin Nota, Benjamin Vernot'
author = 'Kevin Nota, Benjamin Vernot'
release = 'stable-beta_1'
language = 'en'

# -- General configuration ---------------------------------------------------


extensions = [
    'sphinx.ext.intersphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']
source_suffix = '.rst'
exclude_patterns = []
locale_dirs = ['locale/']
gettext_compact = False

master_doc = 'index'
suppress_warnings = ['image.nonlocal_uri']
pygments_style = 'default'

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': True,
    'navigation_depth': 5,
}
html_context = {}

# -- Dark mode toggle --------------------------------------------------------

default_dark_mode = True

# -- Static files (custom CSS, images, etc.) ---------------------------------

html_static_path = ['_static']
