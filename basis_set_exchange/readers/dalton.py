import re
from .. import lut, manip
from . import helpers
from .nwchem import _parse_ecp_lines

all_element_names = [x.lower() for x in lut.all_element_names()]
all_element_names.append('wolffram')  # For tungsten

# Lines beginning a shell can be:
#   H {nprim} {ngen}
# or
#   {nprim} {ngen} 0
# This regex handles both cases
shell_begin_re = re.compile(r'^(?:[hH]\s+)?(\d+)\s+(\d+)(?: +0)?$')

# Typical format for starting an element with a comment
# Not sure this is absolutely required though
element_begin_re = re.compile(r'^!\s+([a-z]+)\s+\(.*\)\s*->\s*\[.*\]$')


#############################################################################
# There seems to be several dalton formats. One splits out by element, with
# each element block starting with "a {element_Z}". Then each shell
# starts with "{nprim} {ngen} 0"
# The second is also split by elements, but the element name
# is in a comment...
#############################################################################
def _line_begins_element(line):
    if not line:
        return False

    line = line.lower()

    if line.startswith('a '):
        return True

    if element_begin_re.match(line):
        s = line[1:].split()
        if s[0] in all_element_names:
            return True
        else:
            raise RuntimeError("Line looks to start an element, but element name is unknown to me. Line: " + line)

    return False


def _parse_electron_lines(basis_lines, bs_data):
    '''Parses lines representing all the electron shells for all elements

    Resulting information is stored in bs_data
    '''
    # fix common spelling mistakes
    basis_lines = [re.sub('PHOSPHOROUS', 'PHOSPHORUS', line) for line in basis_lines]
    # A little bit of a hack here
    # If we find the start of an element, remove all the following comment lines
    new_basis_lines = []
    i = 0
    while i < len(basis_lines):
        if _line_begins_element(basis_lines[i]):
            new_basis_lines.append(basis_lines[i])
            i += 1
            while basis_lines[i].startswith('!'):
                i += 1
        else:
            new_basis_lines.append(basis_lines[i])
            i += 1

    basis_lines = helpers.prune_lines(new_basis_lines, '$')

    # Now split out all the element blocks
    element_blocks = helpers.partition_lines(basis_lines, _line_begins_element, min_size=3)

    # For each block, split out all the shells
    for el_lines in element_blocks:
        # Figure out which type of block this is (does it start with 'a ' or a comment
        header = el_lines[0].lower()
        if header.startswith('a '):
            element_Z = helpers.parse_line_regex(r'a +(\d+) *$', header, "a {element_z}", convert_int=False)
        elif header.startswith('!'):
            element_name = helpers.parse_line_regex(element_begin_re, header, '! {element_name}')
            element_Z = lut.element_Z_from_name(element_name, as_str=True)
        else:
            raise RuntimeError("Unable to parse block in dalton: header line is \"{}\"".format(header))

        element_data = manip.create_element_data(bs_data, element_Z, 'electron_shells')
        el_lines.pop(0)

        # Remove all the rest of the comment lines
        el_lines = helpers.prune_lines(el_lines, '!')

        # Now partition again into blocks of shells for this element
        shell_blocks = helpers.partition_lines(el_lines, lambda x: shell_begin_re.match(x))

        # Shells are written in increasing angular momentum
        shell_am = 0

        for sh_lines in shell_blocks:
            nprim, ngen = helpers.parse_line_regex(shell_begin_re, sh_lines[0], 'nprim, ngen')
            bas_lines = sh_lines[1:]
            # fix for split over newline
            if nprim > 0 and bas_lines:
                num_line_splits = len(sh_lines[1:]) // nprim
                if num_line_splits * nprim == len(sh_lines[1:]):
                    bas_lines = [
                        ' '.join([sh_lines[1 + num_line_splits * i + offset] for offset in range(num_line_splits)])
                        for i in range(nprim)
                    ]
            exponents, coefficients = helpers.parse_primitive_matrix(bas_lines, nprim=nprim, ngen=ngen)

            func_type = lut.function_type_from_am([shell_am], 'gto', 'spherical')

            shell = {
                'function_type': func_type,
                'region': '',
                'angular_momentum': [shell_am],
                'exponents': exponents,
                'coefficients': coefficients
            }

            element_data['electron_shells'].append(shell)
            shell_am += 1


def read_dalton(basis_lines):
    '''Reads Dalton-formatted file data and converts it to a dictionary with the
       usual BSE fields
    '''

    ###########################################################################
    # We need to leave in comments until later, since they can be significant
    # (one format allows "! {ELEMENT}" to start an element block)
    ###########################################################################

    # But we still prune blank lines
    basis_lines = helpers.prune_lines(basis_lines)

    bs_data = {}

    # Skip forward until either:
    # 1. Line begins with 'a'
    # 2. Line begins with 'ecp'
    # 2. Lines begins with '!', with an element name following
    while basis_lines and not _line_begins_element(basis_lines[0]) and basis_lines[0].lower() != 'ecp':
        basis_lines.pop(0)

    # Empty file?
    if not basis_lines:
        return bs_data

    # Partition into ECP and electron blocks
    # I don't think Dalton supports ECPs, but the original BSE
    # Used the NWChem output format for the ECP part
    basis_sections = helpers.partition_lines(basis_lines, lambda x: x.lower() == 'ecp', min_blocks=1, max_blocks=2)

    for s in basis_sections:
        if s[0].lower() == 'ecp':
            _parse_ecp_lines(s, bs_data)
        else:
            _parse_electron_lines(s, bs_data)

    return bs_data
