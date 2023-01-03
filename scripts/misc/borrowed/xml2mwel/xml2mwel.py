from __future__ import division, print_function, unicode_literals
from collections import defaultdict, OrderedDict
import sys
from xml.etree import ElementTree
from xml.parsers import expat

sys.path.insert(0, '/Library/Application Support/MWorks/MWEL')
import mwel
import mwel.parser


def parse_xml(filename):
    tb = ElementTree.TreeBuilder()

    def xml_decl(version, encoding, standalone):
        pass

    def start_element(name, attrs):
        attrs = OrderedDict(attrs[i:i+2] for i in range(0, len(attrs), 2))
        tb.start(name, attrs)

    def end_element(name):
        tb.end(name)

    def comment(data):
        tag = ElementTree.Comment
        tb.start(tag, {})
        tb.data(data)
        tb.end(tag)

    def default(data):
        if data.strip():
            raise RuntimeError('Unhandled XML data: %r', data)

    p = expat.ParserCreate()

    p.XmlDeclHandler = xml_decl
    p.StartElementHandler = start_element
    p.EndElementHandler = end_element
    p.CommentHandler = comment
    p.DefaultHandlerExpand = default

    p.buffer_text = True
    p.ordered_attributes = True

    with open(filename, 'rb') as fp:
        p.ParseFile(fp)

    return ElementTree.ElementTree(tb.close())


_str = type('')


def _is_true(value):
    return (_str(value).lower() in ('true', 'yes', '1'))


def _is_false(value):
    return (_str(value).lower() in ('false', 'no', '0'))


class Converter(object):

    _groups = {
        'io_devices': 'I/O Devices',
        'variables': 'Variables',
        'sounds': 'Sounds',
        'stimuli': 'Stimuli',
        'filters': 'Filters',
        'optimizers': 'Optimizers',
        'calibrators': 'Calibrators',
        'resources': 'Resources',
        'experiment': 'Experiment',
        }

    _ignored = (
        'action_marker',
        'transition_marker',
        )

    _preferred_aliases = {
        'task_system': 'task',
        'task_system_state': 'state',
        'variable': 'var',
        }

    _omit_tag = (
        'action',
        'list_replicator',
        'range_replicator',
        'transition',
        )

    _ignored_attrib = (
        '_id',
        '_line_number',
        '_location',
        'full_name',
        )

    _tab = '    '

    def __init__(self):
        self._expr_parser = mwel.parser.ExpressionParser(self._log_error)

        self._params = {}
        for info in mwel.get_component_info().values():
            if not info.get('abstract', False):
                signature = info['signature']
                params = dict((p['name'], p)
                              for p in info.get('parameters', []))
                self._params[signature] = params
                for alias in info.get('alias', []) + info.get('mwel_alias', []):
                    self._params[alias] = params

        short_type_counts = defaultdict(lambda: 0)
        for signature in self._params:
            short_type_counts[signature.split('/')[-1]] += 1
        self._unique_short_types = set(k for k, v in short_type_counts.items()
                                       if v == 1)

        self._single_param_names = {}
        for signature, params in self._params.items():
            params = list(params.values())
            if len(params) > 1:
                params = [p for p in params if p.get('required', False)]
            if len(params) == 1:
                self._single_param_names[signature] = params[0]['name']

    def _log_error(self, msg, token=None, lineno=None, colno=None, filename=''):
        self._error_count += 1

    def _is_valid_expr(self, expr):
        self._error_count = 0
        self._expr_parser.parse(expr)
        return (self._error_count == 0)

    def _quote_tag(self, tag):
        if self._is_valid_expr(tag):
            return tag
        return "'%s'" % tag

    def _quote_param_value(self, expr):
        if self._is_valid_expr('[%s]' % expr):
            return expr
        return "'%s'" % expr

    def _print_line(self, text, indent=0):
        sys.stdout.write(self._tab * indent)
        sys.stdout.write(text)
        sys.stdout.write('\n')

    def _convert_group(self, elem, indent):
        self._print_line('')
        self._print_line('//', indent)
        self._print_line('// ' + self._groups[elem.tag], indent)
        self._print_line('//', indent)
        self._print_line('')
        for child in elem:
            self._convert(child, indent)

    def _convert_component(self, elem, indent):
        signature = elem.tag
        attrib = elem.attrib
        children = list(elem)

        type = attrib.get('type', '').lower()
        if type and ((signature != 'variable') or (type == 'selection')):
            signature += '/' + type
            del attrib['type']
        signature = self._preferred_aliases.get(signature, signature)

        if signature == 'action/assignment':
            var = attrib['variable']
            if self._is_valid_expr(var):
                value = attrib['value']
                self._print_line('%s = %s' % (var, value), indent)
                return

        short_type = signature.split('/')[-1]
        if short_type in self._unique_short_types:
            line = short_type
        else:
            line = signature

        tag = attrib.pop('tag', None)
        if tag and (elem.tag not in self._omit_tag):
            line += ' ' + self._quote_tag(tag)

        if attrib:
            for name in self._ignored_attrib:
                attrib.pop(name, None)
            for name, value in list(attrib.items()):
                if value == '':
                    del attrib[name]
            for name, info in self._params.get(signature, {}).items():
                if (name in attrib) and ('default' in info):
                    value = attrib[name]
                    default = info['default']
                    if ((value == default) or
                        (_is_true(value) and _is_true(default)) or
                        (_is_false(value) and _is_false(default))):
                        attrib.pop(name)

        if not (attrib or children):
            line += ' ()'
            self._print_line(line, indent)
            return

        if ((len(attrib) == 1) and
            (self._single_param_names.get(signature) in attrib)):
            line += ' (%s)' % self._quote_param_value(attrib.popitem()[1])
        elif attrib:
            line += ' ('
            self._print_line(line, indent)
            for param, value in attrib.items():
                self._print_line('%s = %s' % (param,
                                              self._quote_param_value(value)),
                                 indent+1)
            line = ')'

        if children:
            line += ' {'
            self._print_line(line, indent)
            for child in children:
                self._convert(child, indent+1)
            line = '}'

        self._print_line(line, indent)

    def _convert(self, elem, indent):
        tag = elem.tag
        if tag is ElementTree.Comment:
            # TODO
            pass
        elif tag in self._ignored:
            pass
        elif tag in self._groups:
            self._convert_group(elem, indent)
        else:
            self._convert_component(elem, indent)

    def convert(self, filename):
        root = parse_xml(filename).getroot()
        assert root.tag == 'monkeyml'
        for elem in root:
            self._convert(elem, 0)


if __name__ == '__main__':
    Converter().convert(sys.argv[1])
