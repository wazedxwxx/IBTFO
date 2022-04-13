import sys

extensions = [ 'sphinx.ext.mathjax']
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'
project = u'IBTFO'
copyright = u'IBTFO Copyright (c) 2022, The Regents of the National University of Defense technology. All rights reserved.'
author = u'Xinxin Wang, Jianhan Liang'
version = u'0.1'
release = u'0.1'
language = None
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
todo_include_todos = False
numfig = True
numfig_format = {'figure': '%s', 'table': '%s', 'code-block': '%s'}
html_theme = 'sphinx_rtd_theme'
#html_static_path = ['_static']
htmlhelp_basename = 'IBTFOdoc'
latex_elements = { }
latex_documents = [
    (master_doc, 'IBTFO.tex', u'IBTFO Documentation',
     author, 'manual'),
]
texinfo_documents = [
    (master_doc, 'IBTFO', u'IBTFO Documentation',
     author, 'IBTFO', 'One line description of project.',
     'Miscellaneous'),
]
