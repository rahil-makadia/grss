# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
# find directory of current file
current_dir = os.path.dirname(os.path.abspath(__file__))
# add one level up to path
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

nbsphinx_thumbnails = {
    'tests/prop/apophis': '_images/tests_prop_apophis_9_0.png',
    'tests/prop/didymos': '_images/tests_prop_didymos_9_0.png',
    'tests/prop/eggl': '_images/tests_prop_eggl_9_0.png',
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
