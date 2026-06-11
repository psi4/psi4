from docutils import nodes
from sphinx.ext.intersphinx import InventoryAdapter

INTERSPHINX_ALIASES = {
    ("py", "class", "qcelemental.models.v1.AtomicResult"): (
        "qcelemental_v1",
        "py:pydantic_model",
        "qcelemental.models.results.AtomicResult",
    ),
    ("py", "class", "qcelemental.models.v1.results.AtomicResult"): (
        "qcelemental_v1",
        "py:pydantic_model",
        "qcelemental.models.results.AtomicResult",
    ),
    ("py", "class", "qcelemental.models.v1.basemodels.ProtoModel"): (
        "qcelemental_v1",
        "py:pydantic_model",
        "qcelemental.models.basemodels.ProtoModel",
    ),
    ("py", "class", "pydantic.v1.main.BaseModel"): (
        "pydantic",
        "py:class",
        "pydantic.main.BaseModel",
    ),
    ("py", "class", "pathlib._local.Path"): (
        "python",
        "py:class",
        "pathlib.Path",
    ),
    ("py", "class", "qcelemental.models.v1.results.AtomicInput"): (
        "qcelemental_v1",
        "py:pydantic_model",
        "qcelemental.models.results.AtomicInput",
    ),
    ("py", "class", "qcelemental.models.v1.common_models.FailedOperation"): (
        "qcelemental_v1",
        "py:pydantic_model",
        "qcelemental.models.common_models.FailedOperation",
    ),
}


def resolve_intersphinx_aliases(app, env, node, contnode):
    alias = INTERSPHINX_ALIASES.get((
        node.get("refdomain"),
        node.get("reftype"),
        node.get("reftarget"),
    ))
    if alias is None:
        return None

    inv_name, inv_type, target = alias
    inventory = InventoryAdapter(env).named_inventory
    item = inventory.get(inv_name, {}).get(inv_type, {}).get(target)
    if item is None:
        return None

    project, version, uri, _dispname = item
    reftitle = f"(in {project} v{version})" if version else f"(in {project})"
    newnode = nodes.reference("", "", internal=False, refuri=uri, reftitle=reftitle)
    newnode.append(contnode)
    return newnode


def setup(app):
    app.connect("missing-reference", resolve_intersphinx_aliases)
