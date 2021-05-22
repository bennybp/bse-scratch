'''
Parser for OpenMolcas' RICDlib format

Written by Susi Lehtola, 2021
'''

import re
from .. import lut, misc
from . import helpers

# Start for new basis entry
basis_head_re = re.compile(r'^/([a-zA-Z]+).({}|)....(aCD|acCD)-aux-basis.\s*$'.format(helpers.basis_name_re_str))
# Elemental charge, lmax, number of basis set blocks
charge_line_re = re.compile(r'^\s*({})\s+(\d+)\s+(\d+)\s*$'.format(helpers.floating_re_str))
dummy_line_re = re.compile(r'^\s*Dummy reference line.\s*$')
# nprim, ncontr, functype
shell_start_re = re.compile('^\s*(\d+)\s+(\d+)\s+(\d+)\s*$')

# Floating point data
array_data_re = re.compile(r'^\s*(?:\s({}))+\s*$'.format(helpers.floating_re_str))


def _parse_basis(basis_lines, bs_data):
    # Parse symbol
    element_symbol, basis_name, basis_type = helpers.parse_line_regex(basis_head_re, basis_lines[0],
                                                                      'Symbol.Basis....a(c)CD-aux-basis.')

    # Initialize BSE basis
    element_Z = str(lut.element_Z_from_sym(element_symbol))
    if element_Z not in bs_data or 'electron_shells' not in bs_data[element_Z]:
        element_data = helpers.create_element_data(bs_data, element_Z, 'electron_shells')
    else:
        element_data = bs_data[element_Z]

    # Read charge line (we don't care about the number of basis sets
    # since we parse everything anyway)
    charge, lmax, _ = helpers.parse_line_regex(charge_line_re, basis_lines[1], 'charge, lmax, foo')

    # We should have two dummy lines
    if not dummy_line_re.match(basis_lines[2]):
        raise RuntimeError('Expected dummy line, got: "{}"!'.format(basis_lines[2]))
    if not dummy_line_re.match(basis_lines[3]):
        raise RuntimeError('Expected dummy line, got: "{}"!'.format(basis_lines[3]))

    # Shell data
    shell_data = basis_lines[4:]

    # Parse the shells
    for l in range(lmax + 1):
        nprim, ncontr, amtype = helpers.parse_line_regex(shell_start_re, shell_data[0], 'charge, lmax, foo')
        shell_data = shell_data[1:]
        if nprim == 0 or ncontr == 0:
            # Skip over dummy entries
            continue

        # Read the exponents
        exponents, shell_data = helpers.read_n_floats(shell_data, nprim)

        # Read the contraction coefficients
        coefficients, shell_data = helpers.read_n_floats(shell_data, nprim * ncontr)
        coefficients = helpers.chunk_list(coefficients, nprim, ncontr)
        coefficients = misc.transpose_matrix(coefficients)

        # Store the shell. It is cartesian, unless the amtype flag is
        # 3 which means spherical transform and removal of the
        # lower-angular-momentum contaminants.
        func_type = helpers.function_type_from_am([l], 'gto', 'spherical' if amtype == 3 else 'cartesian')
        shell = {
            'function_type': func_type,
            'region': '',
            'angular_momentum': [l],
            'exponents': exponents,
            'coefficients': coefficients
        }
        element_data['electron_shells'].append(shell)


def read_ricdlib(basis_lines):
    '''Reads ricdlib-formatted file data and converts it to a dictionary with the
       usual BSE fields

       Note that the molcas format does not store all the fields we
       have, so some fields are left blank
    '''

    bs_data = {}

    # Split into elements
    element_blocks = helpers.partition_lines(basis_lines, basis_head_re.match)

    # Inside the loop, check that all blocks refer to the same basis set
    basis_names_found = set()

    for element_lines in element_blocks:
        element_symbol, basis_name, basis_type = helpers.parse_line_regex(basis_head_re, element_lines[0],
                                                                          'Symbol.Basis....a(c)CD-aux-basis.')
        basis_names_found.add(basis_name.lower())
        _parse_basis(element_lines, bs_data)

    # Check for multiple basis sets
    if len(basis_names_found) > 1:
        raise RuntimeError("Multiple basis sets found in file: " + ','.join(basis_names_found))

    return bs_data
