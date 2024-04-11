# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import glob
import os

import splinepy

# sys.path.insert(0, os.path.abspath("../../"))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "splinepy"
copyright = "2022, Jaewook Lee"  # noqa A001
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
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_mdinclude",
    "breathe",
]


templates_path = ["_templates"]
exclude_patterns = []

pygments_style = "sphinx"


include_path = os.path.join(os.path.dirname(__file__), "../../include")
files = glob.glob(
    os.path.join(include_path, "**/**/**/**/*.hpp"), recursive=True
)
l_inc = len(include_path)
files_path = []
for f in files:
    files_path.append(os.path.relpath(f, include_path))

doxygenxml_path = os.path.join(
    os.path.dirname(__file__), "../build/doxy/doxygenxml"
)
breathe_projects = {"splinepy": doxygenxml_path}
breathe_default_project = "splinepy"
breathe_projects_source = {"splinepy": (include_path, files_path)}

# -- Options for HTML output -------------------------------------------------

html_theme = "pydata_sphinx_theme"
html_theme_options = {
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/tataratat/splinepy",
            "icon": "fa-brands fa-square-github",
        },
        {
            "name": "PyPI",
            "url": "https://pypi.org/project/splinepy/",
            "icon": "fa-solid fa-box",
        },
    ],
    "navigation_with_keys": False,
}
# html_favicon = "_static/favicon.ico"
html_static_path = ["_static"]
html_css_files = ["style.css"]

autodoc_default_options = {
    "autosummary": True,
}

autosummary_context = {
    "skipmethods": ["__init__"],
}
