from cgitb import html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'NekRS'
copyright = '2024, ANL et al'
author = 'Jez Swann & ??????'

extensions = [
    'breathe'
]

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'venv']

html_theme = "sphinx_rtd_theme"

# Breathe configuration
breathe_projects = {"nekrs": "doxygen/xml"}
breathe_default_project = "nekrs"
breathe_default_members = ('members', 'undoc-members', 'private-members')