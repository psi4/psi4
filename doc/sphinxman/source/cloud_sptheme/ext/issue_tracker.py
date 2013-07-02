"""cloud_sptheme.ext.issue_tracker - added ":issue:" text role to link to issue tracker"""
#===========================================================================
# imports
#===========================================================================
#core
import re
#site
from docutils import nodes, utils
from docutils.parsers.rst.roles import set_classes
#local

#===========================================================================
# issue role
#===========================================================================

def make_error(inliner, rawtext, line, value):
    "generate error node and msg"
    msg = inliner.reporter.error(value, line=line)
    node = inliner.problematic(rawtext, rawtext, msg)
    return [node], [msg]

def get_issue_tracker_title(config):
    "retreive issue_tracker_title template"
    return getattr(config, "issue_tracker_title", None) or "issue {issue}"

def get_issue_tracker_url(config):
    "retrieve issue_tracker_url template, replacing aliases"
    template = getattr(config, "issue_tracker_url", None)
    if not template:
        # causes :issue:`xx` to be replaced with label instead of url.
        return None

    elif template.startswith("bb:"):
        # parse "bb:<project>", and replace with bitbucket url
        project = template[3:].strip("/")
        return "http://bitbucket.org/" + project + "/issue/{issue}/"

    elif template.startswith("gc:"):
        # parse "gc:<project>", and replace with google code url
        project = template[3:].strip("/")
        return "http://code.google.com/p/" + project + \
                    "/issues/detail?id={issue}"

    else:
        # assume it contains {issue} and possibly {title}
        return template

# pattern allows inside :issue: text roles
issue_re = re.compile("""
    ^
    (?:
        (?P<title>[^<]+)
        \s*
        <
        (?P<issue1>\d+)
        >
    |
        (?P<issue2>\d+)
    )
    $
    """, re.X)

def issue_role(name, rawtext, text, line, inliner, options={}, content=[]):
    "generate link to an issue"
    #NOTE:
    #   name - role name in doc, should be 'issue'
    #   rawtext - text of entire node
    #   text - contents of role
    #   lineno
    #   inliner - ???
    #   options - ???
    #   content - ???
    # returns (nodes, messages)

    # extract title & issue number from text
    m = issue_re.match(text)
    if m:
        issue = int(m.group("issue1") or m.group("issue2"))
        title = m.group("title")
    else:
        return make_error(inliner, rawtext, line,
                          "Invalid issue identifier: %r" % (text,))

    # get url template from config, resolve aliases
    config = inliner.document.settings.env.app.config
    url_template = get_issue_tracker_url(config)
    title_template = get_issue_tracker_title(config)

    # generate replacement node
    if not title:
        title = title_template.format(issue=issue)
    set_classes(options)
    clist = options.setdefault('classes',[])
    clist.append("issue")
    if url_template:
        url = url_template.format(issue=issue, title=title)
        node = nodes.reference(rawtext, title, refuri=url, **options)
    else:
        node = nodes.label(rawtext, title, **options)
    return [node], []

#===========================================================================
# init
#===========================================================================

def setup(app):
    "install plugin"

    app.add_config_value('issue_tracker_url', None, 'env')
    app.add_config_value('issue_tracker_title', None, 'env')
    app.add_role('issue', issue_role)

#===========================================================================
# eof
#===========================================================================
