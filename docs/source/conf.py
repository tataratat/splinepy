# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os
import sys

#sys.path.insert(0, os.path.abspath("../../"))

import splinepy

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "splinepy"
copyright = "2022, Jaewook Lee"
author = "Jaewook Lee"
release = splinepy.__version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "markdown",
    ".md": "markdown",
}

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.inheritance_diagram",
    "sphinx_mdinclude",
]


templates_path = ["_templates"]
exclude_patterns = []

pygments_style = "sphinx"


# -- Options for HTML output -------------------------------------------------

html_theme = "piccolo_theme"
# html_static_path = ['_static']

# def skip(app, what, name, obj, would_skip, options):
#    if name == "__init__":
#        return False
#
#    return would_skip
#
# def setup(app):
#    app.connect("autodoc-skip-member", skip)
