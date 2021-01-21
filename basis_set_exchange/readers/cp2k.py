'''
Reader for CP2K basis sets

Written by Susi Lehtola, 2021
'''

import regex
from .. import lut
from . import helpers

# Shell entry: 'element names'
element_shell_re = regex.compile(r'^\s*(?P<sym>\w+)(?:\s+(?P<name>\w+))+\s*$')
# Number of blocks
nblocks_re = regex.compile(r'^\s*(?P<nblocks>\d+)\s*$')
# Block entry: n lmin lmax nexp nshell(lmin) nshell(lmin+1) ... nshell(lmax-1) nshell(lmax)
block_re = regex.compile(r'^\s*(?P<n>\d+)\s+(?P<lmin>\d+)\s+(?P<lmax>\d+)\s+(?P<nexp>\d+)\s+(?P<nshell>\d+)+\s*$')

# Function type: spherical by default
_func_type = 'gto_spherical'


def _read_shell(basis_lines, bs_data):
    '''Read in a shell from the input'''
    # Read the shell entry
    iline = 0
    print('Block:\n{}'.format(basis_lines))
    element_shell = helpers.parse_line_regex_dict(element_shell_re, basis_lines[iline], 'element (aliases)')
    iline += 1

    # Create entry
    element_sym = element_shell['sym'][0]
    element_Z = lut.element_Z_from_sym(element_sym, as_str=True)
    if element_Z not in bs_data or 'electron_shells' not in bs_data[element_Z]:
        element_data = helpers.create_element_data(bs_data, element_Z, 'electron_shells')
    else:
        element_data = bs_data[element_Z]

    # Read the number of blocks
    nblocks = helpers.parse_line_regex_dict(nblocks_re, basis_lines[iline], 'nblocks')['nblocks'][0]
    iline += 1

    # Read blocks
    for iblock in range(nblocks):
        block = helpers.parse_line_regex_dict(
            block_re, basis_lines[iline],
            'n lmin lmax nexp nshell(lmin) nshell(lmin+1) ... nshell(lmax-1) nshell(lmax)')
        iline += 1

        lmin = block['lmin'][0]
        lmax = block['lmax'][0]
        nprim = block['nprim'][0]
        nshell = block['nshell']
        assert (len(nshell) == lmax - lmin + 1)

        # Total number of contractions is
        ncontr = sum(nshell)

        # Parse contractions
        exps, coeffs = helpers.parse_primitive_matrix(basis_lines[iline:iline + nprim], nprim=nprim, ngen=ncontr)
        iline += nprim

        # Add the shells
        for l in range(lmin, lmax + 1):
            col_offset = sum(nshell[:l - lmin])
            ncontr_l = nshell[l - lmin]
            l_coeff = coeffs[:, col_offset:col_offset + ncontr_l]
            shell = {
                'function_type': func_type,
                'region': '',
                'angular_momentum': shell_am,
                'exponents': exps,
                'coefficients': l_coeff
            }
            element_data['electron_shells'].append(shell)


def read_cp2k(basis_lines):
    '''Reads basis set in CP2K format and converts it to a dictionary with
       the usual BSE fields

       Note that the CP2K format does not store all the fields we
       have, so some fields are left blank

    '''

    # Removes comments
    basis_lines = helpers.prune_lines(basis_lines, '!')

    bs_data = {}

    # Empty file?
    if not basis_lines:
        return bs_data

    # split into element sections
    element_sections = helpers.partition_lines(basis_lines, element_shell_re.match)
    for es in element_sections:
        print('section\n{}'.format(es))
    for es in element_sections:
        _read_shell(es, bs_data)

    return bs_data
