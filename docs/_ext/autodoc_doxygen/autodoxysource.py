"""autodoxysource directive: render doxygen <programlisting> as
syntax-highlighted source with clickable cross-references."""

from docutils import nodes
from docutils.parsers.rst import Directive
from sphinx import addnodes
from sphinx.util import logging

from . import get_doxygen_root

logger = logging.getLogger(__name__)

# Map doxygen highlight classes to CSS classes for the source listing.
_HL_CSS = {
    'comment': 'f-hl-comment',
    'normal': 'f-hl-normal',
    'keyword': 'f-hl-keyword',
    'keywordtype': 'f-hl-keywordtype',
    'keywordflow': 'f-hl-keywordflow',
    'stringliteral': 'f-hl-stringliteral',
    'preprocessor': 'f-hl-preprocessor',
}

# Lazily-built index: doxygen refid -> (refdomain, reftype, reftarget)
_REFID_INDEX = None

# Index: absolute file path -> doxygen file-ID (for [source] links)
_FILE_PATH_INDEX = None


def _build_refid_index():
    """Walk the merged doxygen tree once, building refid -> xref info."""
    global _REFID_INDEX
    _REFID_INDEX = {}
    root = get_doxygen_root()

    for cd in root.findall('./compounddef'):
        cid = cd.get('id')
        kind = cd.get('kind')
        cn = cd.find('compoundname')
        if cn is None or cn.text is None:
            continue
        name = cn.text

        if kind == 'namespace':
            _REFID_INDEX[cid] = ('f', 'mod', name)
        elif kind in ('type', 'struct'):
            # compoundname is "parent::typename" -> "parent/typename"
            _REFID_INDEX[cid] = ('f', 'type', name.replace('::', '/'))
        elif kind == 'interface':
            _REFID_INDEX[cid] = ('f', 'func', name.replace('::', '/'))
        # kind='file' / 'page' / 'dir' -> no xref

        # Index memberdefs inside this compound
        for md in cd.iter('memberdef'):
            mid = md.get('id')
            mk = md.get('kind')
            mn_el = md.find('name')
            if mn_el is None or mn_el.text is None:
                continue
            qualified = name + '/' + mn_el.text
            if mk == 'function':
                _REFID_INDEX[mid] = ('f', 'func', qualified)
            # variables / enums aren't useful xref targets here


def _resolve_ref(refid):
    """O(1) lookup of a doxygen refid to (refdomain, reftype, reftarget)."""
    global _REFID_INDEX
    if _REFID_INDEX is None:
        _build_refid_index()
    return _REFID_INDEX.get(refid)


def _build_file_path_index():
    """Build filepath -> file_id index from compounddef[@kind='file']."""
    global _FILE_PATH_INDEX
    _FILE_PATH_INDEX = {}
    root = get_doxygen_root()
    for cd in root.findall('./compounddef[@kind="file"]'):
        fid = cd.get('id')
        loc = cd.find('location')
        if loc is not None and loc.get('file'):
            _FILE_PATH_INDEX[loc.get('file')] = fid


def get_source_link(xml_node):
    """Given a doxygen XML node (compounddef or memberdef), return
    (docname, line) for a [source] link, or None if unavailable.

    *docname* is the Sphinx document name (e.g. 'api/generated/source/MOM_8F90').
    *line* is the source line number string.
    """
    global _FILE_PATH_INDEX
    if _FILE_PATH_INDEX is None:
        _build_file_path_index()

    loc = xml_node.find('location')
    if loc is None:
        return None

    filepath = loc.get('file')
    line = loc.get('line', '1')
    if not filepath:
        return None

    file_id = _FILE_PATH_INDEX.get(filepath)
    if file_id is None:
        return None

    docname = 'api/generated/source/' + file_id
    return (docname, line)


class AutoDoxySourceDirective(Directive):
    """.. autodoxysource:: <file-id>

    Render the source listing for a doxygen file compound, with
    per-line anchors, syntax highlighting, and clickable identifiers
    that link to the Sphinx API documentation.
    """
    required_arguments = 1
    optional_arguments = 0
    has_content = False

    def run(self):
        file_id = self.arguments[0]
        root = get_doxygen_root()

        compounddef = root.find('./compounddef[@id="%s"]' % file_id)
        if compounddef is None:
            logger.warning('autodoxysource: compounddef not found for %s',
                           file_id)
            return [nodes.paragraph('', 'Source listing not available.')]

        programlisting = compounddef.find('.//programlisting')
        if programlisting is None:
            logger.warning('autodoxysource: no programlisting in %s', file_id)
            return [nodes.paragraph('', 'Source listing not available.')]

        # Outer wrapper
        table = nodes.container(classes=['autodoxysource'])

        for codeline in programlisting.findall('codeline'):
            lineno = codeline.get('lineno', '')

            line_node = nodes.container(classes=['source-line'])

            # Anchor target for #L<N> links
            target = nodes.target('', '', ids=['L' + lineno])
            line_node += target

            # Line number display
            ln = nodes.inline(lineno, lineno, classes=['source-lineno'])
            line_node += ln

            # Code content container
            code_node = nodes.inline(classes=['source-code'])

            for hl in codeline:
                if hl.tag != 'highlight':
                    continue
                css = _HL_CSS.get(hl.get('class', 'normal'), 'f-hl-normal')
                _walk_highlight(hl, css, code_node)

            line_node += code_node
            table += line_node

        return [table]


def _walk_highlight(hl, css_class, parent):
    """Walk mixed content of a <highlight> element, appending nodes to *parent*."""
    # Leading text of the <highlight> element itself
    if hl.text:
        parent += nodes.inline(hl.text, hl.text, classes=[css_class])

    for child in hl:
        if child.tag == 'sp':
            parent += nodes.Text(' ')
        elif child.tag == 'ref':
            _emit_ref(child, css_class, parent)
        else:
            # Unknown child — render as plain text
            if child.text:
                parent += nodes.inline(child.text, child.text,
                                       classes=[css_class])

        # Tail text after the closing tag of this child
        if child.tail:
            parent += nodes.inline(child.tail, child.tail,
                                   classes=[css_class])


def _emit_ref(ref_el, css_class, parent):
    """Emit a pending_xref (or plain text fallback) for a <ref> element."""
    ref_text = ref_el.text or ''
    refid = ref_el.get('refid', '')

    resolved = _resolve_ref(refid)
    if resolved:
        refdomain, reftype, reftarget = resolved
        inner = nodes.inline(ref_text, ref_text, classes=[css_class])
        xref = addnodes.pending_xref(
            '', inner,
            refdomain=refdomain,
            reftype=reftype,
            reftarget=reftarget,
        )
        parent += xref
    else:
        parent += nodes.inline(ref_text, ref_text, classes=[css_class])
