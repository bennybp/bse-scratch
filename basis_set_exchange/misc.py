'''
Miscellaneous helper functions
'''

import re
from . import lut


def _Z_from_str(s):
    if s.isdecimal():
        return int(s)
    else:
        return lut.element_Z_from_sym(s)


def transpose_matrix(mat):
    '''Transposes a matrix (list of lists) commonly done do coefficients'''

    return list(map(list, zip(*mat)))


def max_am(shells):
    '''Determine the maximum angular momentum of a list of shells or potentials'''
    all_am = [max(x['angular_momentum']) for x in shells]
    return max(all_am)


def contraction_string(element, compact=False):
    """
    Forms a string specifying the contractions for an element

    ie, (16s,10p) -> [4s,3p]
    or, if compact=True, 16s10p.4s3p
    """

    # Does not have electron shells (ECP only?)
    if 'electron_shells' not in element:
        return ""

    cont_map = dict()
    for sh in element['electron_shells']:
        nprim = len(sh['exponents'])
        ngeneral = len(sh['coefficients'])

        # is a combined general contraction (sp, spd, etc)
        is_spdf = len(sh['angular_momentum']) > 1

        for am in sh['angular_momentum']:
            # If this a general contraction (and not combined am), then use that
            ncont = ngeneral if not is_spdf else 1

            if am not in cont_map:
                cont_map[am] = (nprim, ncont)
            else:
                cont_map[am] = (cont_map[am][0] + nprim, cont_map[am][1] + ncont)

    primstr = ""
    contstr = ""
    for am in sorted(cont_map.keys()):
        nprim, ncont = cont_map[am]

        if am != 0 and not compact:
            primstr += ','
            contstr += ','
        primstr += str(nprim) + lut.amint_to_char([am])
        contstr += str(ncont) + lut.amint_to_char([am])

    if compact:
        return "{}.{}".format(primstr, contstr)

    return "({}) -> [{}]".format(primstr, contstr)


def compact_elements(elements):
    """
    Create a string (with ranges) given a list of element numbers

    For example, [1, 2, 3, 6, 7, 8, 10] will return "H-Li,C-O,Ne"
   """

    if not elements:
        return

    # We have to convert to integers for this function
    elements = [int(el) for el in elements]

    # Just to be safe, sort the list
    el = sorted(set(elements))

    ranges = []
    i = 0
    while i < len(el):
        start_el = el[i]
        end_el = start_el

        i += 1
        while i < len(el):
            if el[i] != end_el + 1:
                break

            end_el += 1
            i += 1

        if start_el == end_el:
            ranges.append([start_el])
        else:
            ranges.append([start_el, end_el])

    # Convert to elemental symbols
    range_strs = []
    for r in ranges:
        sym = lut.element_sym_from_Z(r[0], True)

        if len(r) == 1:
            range_strs.append(sym)
        elif len(r) == 2 and r[1] == r[0] + 1:
            sym2 = lut.element_sym_from_Z(r[1], True)
            range_strs.append(sym + "," + sym2)
        else:
            sym2 = lut.element_sym_from_Z(r[1], True)
            range_strs.append(sym + "-" + sym2)

    return ",".join(range_strs)


def expand_elements(compact_el, as_str=False):
    """
    Create a list of integers given a string or list of compacted elements

    This is partly the opposite of compact_elements, but is more flexible.

    compact_el can be a list or a string. If compact_el is a list, each element is processed individually
    as a string (meaning list elements can contain commas, ranges, etc)
    If compact_el is a string, it is split by commas and then each section is processed.

    In all cases, element symbols (case insensitive) and Z numbers (as integers or strings)
    can be used interchangeably. Ranges are also allowed in both lists and strings.

    Some examples:
        "H-Li,C-O,Ne" will return [1, 2, 3, 6, 7, 8, 10]
        "H-N,8,Na-12" will return [1, 2, 3, 4, 5, 6, 7, 8, 11, 12]
        ['C', 'Al-15,S', 17, '18'] will return [6, 13, 14, 15, 16, 17, 18]

    If as_str is True, the list will contain strings of the integers
    (ie, the first example above will return ['1', '2', '3', '6', '7', '8', '10']
    """

    # If an integer, just return it
    if isinstance(compact_el, int):
        if as_str:
            return [str(compact_el)]
        else:
            return [compact_el]

    # If compact_el is a list, make it a comma-separated string
    if isinstance(compact_el, list):
        compact_el = [str(x) for x in compact_el]
        compact_el = [x for x in compact_el if x]
        compact_el = ','.join(compact_el)

    # Find multiple - or ,
    # Also replace all whitespace with spaces
    compact_el = re.sub(r',+', ',', compact_el)
    compact_el = re.sub(r'-+', '-', compact_el)
    compact_el = re.sub(r'\s+', '', compact_el)

    # Find starting with or ending with comma and strip them
    compact_el = compact_el.strip(',')

    # Check if I was passed an empty string or list
    if not compact_el:
        return []

    # Find some erroneous patterns
    # -, and ,-
    if '-,' in compact_el:
        raise RuntimeError("Malformed element string")
    if ',-' in compact_el:
        raise RuntimeError("Malformed element string")

    # Strings ends or begins with -
    if compact_el.startswith('-') or compact_el.endswith('-'):
        raise RuntimeError("Malformed element string")

    # x-y-z
    if re.search(r'\w+-\w+-\w+', compact_el):
        raise RuntimeError("Malformed element string")

    # Split on commas
    tmp_list = compact_el.split(',')

    # Now go over each one and replace elements with ints
    el_list = []
    for el in tmp_list:
        if '-' not in el:
            el_list.append(_Z_from_str(el))
        else:
            begin, end = el.split('-')
            begin = _Z_from_str(begin)
            end = _Z_from_str(end)
            el_list.extend(list(range(begin, end + 1)))

    if as_str:
        return [str(x) for x in el_list]
    else:
        return el_list


def transform_basis_name(name):
    """
    Transforms the name of a basis set to an internal representation

    This makes comparison of basis set names easier by, for example,
    converting the name to all lower case.
    """

    name = name.lower()
    name = name.replace('/', '_sl_')
    name = name.replace('*', '_st_')
    return name


def basis_name_to_filename(name):
    '''
    Given a basis set name, transform it into a valid filename

    This makes sure filenames don't contain invalid characters
    '''

    return transform_basis_name(name)


def basis_name_from_filename(filename):
    '''
    Given a basis set name that was part of a filename, determine the basis set name

    This is opposite of :func:`transform_basis_name`

    Pass only the part of the filename that contains the basis set name
    '''

    name = filename.lower()
    name = name.replace('_sl_', '/')
    name = name.replace('_st_', '*')
    return name
