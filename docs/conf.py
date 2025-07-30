# Configuration file for the Sphinx documentation builder.

import os
import sys
import re
sys.path.insert(0, os.path.abspath('..'))

import sphinx_rtd_theme

# -- Project information -----------------------------------------------------

project = 'NED-flow'
slug = re.sub(r'\W+', '-', project.lower())
copyright = '2025, Kevin Nota, Benjamin Vernot'
author = 'Kevin Nota, Benjamin Vernot'
release = 'stable-beta_1'
language = 'en'

# -- General configuration ---------------------------------------------------

templates_path = ['_templates']
exclude_patterns = []

extensions = [
    'sphinx_rtd_theme',
]

html_theme = 'sphinx_rtd_theme'

html_static_path = ['_static']

    
