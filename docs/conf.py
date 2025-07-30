# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'NED-flow'
copyright = '2025, Kevin Nota, Benjamin Vernot'
author = 'Kevin Nota, Benjamin Vernot'
release = '0.1.0_beta_1'

#import sphinx_rtd_theme

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

templates_path = ['_templates']
exclude_patterns = []
#exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_css_files = [
    'custom.css',
]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-out

html_theme = "pydata_sphinx_theme"

extensions = ["sphinx.ext.autodoc"] 

html_static_path = ['_static']

html_theme_options = {
    "navbar_start": ["navbar-logo"],
    "navbar_end": ["navbar-icon-links"],
    "navbar_persistent": ["search-button"],
    "logo": {
        "text": "My Project",
    },
    "use_edit_page_button": False,
    "github_url": "https://github.com/yourname/yourrepo",  # optional
}
