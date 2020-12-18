'''
Reader for Molpro system library basis sets

Written by Susi Lehtola, 2020
'''

import re
import regex
from .. import lut
from . import helpers

# Shell entry: 'element am (aliases) : nprim ncontr start1.end1 start2.end2 ... startn.endn' allowing whitespace
element_shell_re = regex.compile(r'^\s*(?P<sym>\w+)\s+(?P<am>[spdfghikSPDFGHIK])\s*(?:\s*(?P<alias>\S+)\s*)+\s*:\s*(?P<nprim>\d+)\s*(?P<ncontr>\d+)\s*(?:\s*(?P<range>\d+.\d+)\s*)+\s*$')
# Exponent / coefficient entry: val1 val2 ... valn, allowing whitespace
entry_re = regex.compile(r'^\s*(?:\s*(?P<val>({}|{}))\s*)+\s*$'.format(helpers.floating_re_str, helpers.integer_re_str))

# Function type: spherical by default
_func_type = 'gto_spherical'

def _parse_electron_lines(basis_lines, bs_data):
    '''Parses lines representing all the electron shells for a single element

    Resulting information is stored in bs_data
    '''

    # Read the data one line at a time
    iline=0

    assert(element_shell_re.match(basis_lines[iline]))

    while element_shell_re.match(basis_lines[iline]):
        # Read the shell entry
        shell = helpers.parse_line_regex_dict(element_shell_re, basis_lines[iline], 'element am (aliases) : nprim ncontr start.end')
        # Read the comment
        iline += 1
        comment = basis_lines[iline]

        # Angular momentum
        assert(len(shell['am']) == 1)
        shell_am = lut.amchar_to_int(shell['am'][0])
        # Element
        assert(len(shell['sym']) == 1)
        element_sym = shell['sym'][0]
        # Number of primitives
        assert(len(shell['nprim']) == 1)
        nprim = shell['nprim'][0]
        assert(nprim>0)
        # Number of contractions
        assert(len(shell['ncontr']) == 1)
        ncontr = shell['ncontr'][0]
        assert(ncontr>0)
        # Contraction ranges
        cranges = shell['range']
        assert(len(cranges) == ncontr)
        # Parse contraction ranges
        cranges = [ r.split('.') for r in cranges ]
        cranges = [ [int(v) for v in r] for r in cranges ]

        # Create entry
        element_Z = lut.element_Z_from_sym(element_sym, as_str=True)
        if element_Z not in bs_data:
            element_data = helpers.create_element_data(bs_data, element_Z, 'electron_shells')
        else:
            element_data = bs_data[element_Z]

        # Count the number of entries we're supposed to read
        nread = nprim
        for r in cranges:
            nread += r[1]-r[0]+1

        # Read in data
        rawdata = []
        nreadin = 0
        while nreadin < nread:
            iline += 1
            entry = helpers.parse_line_regex_dict(entry_re, basis_lines[iline], 'exponent / contraction', convert_int = False)
            # Ensure all parsed values have a decimal point
            readval = [k if k.find('.')!=-1 else k+'.0' for k in entry['val']]
            rawdata = rawdata + readval
            nreadin += len(readval)

        # Collect the exponents
        exponents = rawdata[:nprim]
        # Collect the contraction coefficients
        coefficients = []
        offset = nprim
        for i in range(ncontr):
            # Start and end
            start = cranges[i][0]
            end = cranges[i][1]
            # Number of entries
            nentries = end-start+1
            # Contraction coefficients
            cc = rawdata[offset:offset+nentries]
            offset += nentries

            # Pad coefficients with zeros
            if start > 1:
                cc = ['0.0' for _ in range(1,start)] + cc
            if end < nprim:
                cc = cc + ['0.0' for _ in range(end,nprim)]

            # Add to contraction
            coefficients.append(cc)

        # Function type
        func_type = 'gto' if shell_am[0] < 2 else _func_type
        # Store the data
        shell = {
            'function_type': func_type,
            'region': '',
            'angular_momentum': shell_am,
            'exponents': exponents,
            'coefficients': coefficients
        }
        element_data['electron_shells'].append(shell)


def _parse_ecp_lines(basis_lines, bs_data):
    '''Parses lines representing all the ECP potentials for a single element

    Resulting information is stored in bs_data
    '''

    # Nothing to do
    x = 1
    return


def read_libmol(basis_lines):
    '''Reads basis set from Molpro system library data and converts it to
       a dictionary with the usual BSE fields

       Note that the Molpro system library format does not store all
       the fields we have, so some fields are left blank

    '''

    # Removes comments
    basis_lines = helpers.prune_lines(basis_lines, '!')

    # Go through input and check basis type
    for line in basis_lines:
        if line.strip().lower() == 'spherical':
            _func_type = 'gto_spherical'
        if line.strip().lower() == 'cartesian':
            _func_type = 'gto_cartesian'

    bs_data = {}

    # Empty file?
    if not basis_lines:
        return bs_data

    _parse_electron_lines(basis_lines, bs_data)

    return bs_data
