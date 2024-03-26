"""Configuration file for the Sphinx documentation builder."""
# pylint: disable=invalid-name
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "GRSS"
project_copyright = "2024, GRSS Development Team"
author = "Rahil Makadia"
# get release and version from version.txt
with open("../../grss/version.txt", "r", encoding="utf-8") as f:
    release = version = f.read().strip()

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx_copybutton",  # for adding copy button to code blocks
    "sphinx.ext.autosummary",  # for generating documentation from docstrings
    "sphinx.ext.duration",  # for printing duration when building docs
    "sphinx.ext.napoleon",  # for parsing numpy style docstrings
    "nbsphinx",  # for parsing jupyter notebooks
    "sphinx_favicon",  # for adding full favicon support
    "IPython.sphinxext.ipython_console_highlighting",  # for syntax highlighting
    "sphinx_gallery.load_style",  # for displaying jupyter notebook thumbnails
    "breathe",  # for linking to C++ documentation
]
autosummary_generate = True

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_context = {
    "default_mode": "light",
}
html_sidebars = {
    "*": [
        "navbar-logo.html",
        "search-field.html",
        "sbt-sidebar-nav.html",
    ]
}
html_theme_options = {
    "repository_url": "https://github.com/rahil-makadia/grss",
    "repository_branch": "main",
    "path_to_docs": "docs",
    "use_download_button": True,
    "use_edit_page_button": False,
    "use_fullscreen_button": True,
    "use_issues_button": True,
    "use_repository_button": True,
    "use_source_button": False,
    "use_sidenotes": False,
    "home_page_in_toc": False,
    "show_navbar_depth": 1,
    "show_toc_level": 1,
    "logo": {
        # Because the logo is also a homepage link,
        # including "home" in the alt text is good practice
        "alt_text": "GRSS - Home",
        "text": f"{project} v{version} documentation",
    },
    "navigation_with_keys": True,
    "extra_footer": ("Created using "
                        "<a href=https://sphinx-book-theme.readthedocs.io/>"
                        "The Sphinx Book Theme</a>."),
}
html_last_updated_fmt = "%b %d, %Y"
html_logo = "_static/grss-cropped.svg"
html_static_path = ["_static"]
favicons = [
    # generic icons compatible with most browsers
    "favicon.ico",
    "favicon-32x32.png",
    "favicon-16x16.png",
    {"rel": "shortcut icon", "sizes": "any", "href": "favicon.ico"},
    # chrome specific
    "android-chrome-192x192.png",
    # apple icons
    {"rel": "mask-icon", "color": "#67b934", "href": "safari-pinned-tab.svg"},
    {"rel": "apple-touch-icon", "href": "apple-touch-icon.png"},
    # msapplications
    {"name": "msapplication-TileColor", "content": "#da532c"},
    {"name": "theme-color", "content": "#ffffff"},
    {"name": "msapplication-TileImage", "content": "mstile-310x310.png"},
]
exclude_patterns = ["**.ipynb_checkpoints"]

breathe_projects = {
	"GRSS": "../doxygen/xml"
}
breathe_default_project = "GRSS"
breathe_default_members = ('members', 'undoc-members')
