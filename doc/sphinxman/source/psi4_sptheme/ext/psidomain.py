"""Extension to format and index PSI variables."""

#Sphinx.add_object_type(psivar, rolename, indextemplate='', parse_node=None, ref_nodeclass=None, objname='', doc_field_types=[])

def setup(app):

    app.add_object_type('psivar', 'psivar', indextemplate='single: %s')

