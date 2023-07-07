# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'GRSS'
copyright = '2023, Rahil Makadia'
author = 'Rahil Makadia'
# get release and version from version.txt
with open('../../grss/version.txt', 'r', encoding='utf-8') as f:
    release = version = f.read().strip()

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx_copybutton', # for adding copy button to code blocks
    'sphinx.ext.autosummary', # for generating documentation from docstrings
    'sphinx.ext.duration', # for printing duration when building docs
    'sphinx.ext.napoleon', # for parsing numpy style docstrings
    'nbsphinx', # for parsing jupyter notebooks
    "sphinx_gallery.load_style", # for displaying jupyter notebook thumbnails
]
autosummary_generate = True

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
exclude_patterns = ['**.ipynb_checkpoints']
