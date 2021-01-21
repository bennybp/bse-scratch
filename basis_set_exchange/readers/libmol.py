'''
Reader for Molpro system library basis sets

Written by Susi Lehtola, 2020
'''

import regex
from .. import lut
from . import helpers

# Shell entry: 'element am (aliases) : nprim ncontr start1.end1 start2.end2 ... startn.endn' allowing whitespace
element_shell_re = regex.compile(
    r'^\s*(?P<sym>\w+)\s+(?P<am>[spdfghikSPDFGHIK])\s*(?:\s*(?P<alias>\w+)\s*)+\s*:\s*(?P<nprim>\d+)\s*(?P<ncontr>\d+)\s*(?:\s*(?P<range>\d+.\d+)\s*)+\s*$'
)
# Exponent / coefficient entry: val1 val2 ... valn, allowing whitespace
entry_re = regex.compile(r'^\s*(?:\s*(?P<val>({}|{}))\s*)+\s*$'.format(helpers.floating_re_str,
                                                                       helpers.integer_re_str))

# ECP entry: symbol ECP names : ncore lmax lmaxso ndata
ecp_re = regex.compile(
    r'\s*(?P<sym>\w+)\s+ECP\s+(?:\s*(?P<alias>\w+)\s*)\s*:\s*(?P<ncore>\d+)\s+(?P<lmax>\d+)\s+(?P<lmaxso>\d+)\s+(?P<ndata>\d+)\s*$'
)
# ECP data block beginning: nterms (rexp1 expn1 coeff1) (rexp2 expn2 coeff2) ...
ecp_block_start_re = regex.compile(
    r'\s*(?P<nterms>\d+)(?:\s+(?P<rexp>\d+)\s+(?P<expn>{0})\s+(?P<coeff>{0}))+\s*$'.format(helpers.floating_re_str))
# ECP data block continuation lines: (rexp3 expn3 coeff3) (rexp4 expn4 coeff4) ...
ecp_block_cont_re = regex.compile(r'\s*(?:(?P<rexp>\d+)\s+(?P<expn>{0})\s+(?P<coeff>{0})\s*)+\s*$'.format(
    helpers.floating_re_str))

# Function type: spherical by default
_func_type = 'gto_spherical'


def _read_shell(basis_lines, bs_data, iline):
    '''Read in a shell from the input'''
    # Read the shell entry
    shell = helpers.parse_line_regex_dict(element_shell_re, basis_lines[iline],
                                          'element am (aliases) : nprim ncontr start.end')
    # Skip the comment line
    iline += 1

    # Angular momentum
    assert (len(shell['am']) == 1)
    shell_am = lut.amchar_to_int(shell['am'][0])
    # Element
    assert (len(shell['sym']) == 1)
    element_sym = shell['sym'][0]
    # Number of primitives
    assert (len(shell['nprim']) == 1)
    nprim = shell['nprim'][0]
    assert (nprim > 0)
    # Number of contractions
    assert (len(shell['ncontr']) == 1)
    ncontr = shell['ncontr'][0]
    assert (ncontr > 0)
    # Contraction ranges
    cranges = shell['range']
    assert (len(cranges) == ncontr)
    # Parse contraction ranges
    cranges = [r.split('.') for r in cranges]
    cranges = [[int(v) for v in r] for r in cranges]

    # Create entry
    element_Z = lut.element_Z_from_sym(element_sym, as_str=True)
    if element_Z not in bs_data or 'electron_shells' not in bs_data[element_Z]:
        element_data = helpers.create_element_data(bs_data, element_Z, 'electron_shells')
    else:
        element_data = bs_data[element_Z]

    # Count the number of entries we're supposed to read
    nread = nprim
    for r in cranges:
        nread += r[1] - r[0] + 1

    # Read in data
    rawdata = []
    nreadin = 0
    while nreadin < nread:
        iline += 1
        entry = helpers.parse_line_regex_dict(entry_re,
                                              basis_lines[iline],
                                              'exponent / contraction',
                                              convert_int=False)
        # Ensure all parsed values have a decimal point
        readval = [k if k.find('.') != -1 else k + '.0' for k in entry['val']]
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
        nentries = end - start + 1
        # Contraction coefficients
        cc = rawdata[offset:offset + nentries]
        offset += nentries

        # Pad coefficients with zeros
        if start > 1:
            cc = ['0.0' for _ in range(1, start)] + cc
        if end < nprim:
            cc = cc + ['0.0' for _ in range(end, nprim)]

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

    return iline


def _read_ecp(basis_lines, bs_data, iline):
    '''Reads an ECP from the input'''

    # Read the ECP entry
    shell = helpers.parse_line_regex_dict(ecp_re, basis_lines[iline], 'symbol ECP names : ncore lmax lmaxso ndata')
    element_sym = shell['sym'][0]
    ncore = shell['ncore'][0]
    lmax = shell['lmax'][0]
    lmaxso = shell['lmaxso'][0]
    ndata = shell['ndata'][0]

    # Skip the comment line
    iline += 1

    # Data entries read in
    nread = 0

    def parse_ecp(basis_lines, iline, lmax, so=False):
        '''Reads in a block of ECP data'''

        if lmax == 0:
            # Nothing to do if lmax is 0
            return [], iline

        # Initialize storage
        ecpblocks = [[] for _ in range(lmax + 1)]
        # SO ECP runs from l=1 to lmaxso
        lmin = 0 if not so else 1

        for l in range(lmin, lmax + 1):
            # ECP data is in the order lmax, 0, 1, ..., lmax-1
            if not so:
                target_l = lmax if (l == 0) else l - 1
            else:
                # SO runs from l=1 to l=lmaxso
                target_l = l

            # Get the entry
            iline += 1
            entry = helpers.parse_line_regex_dict(ecp_block_start_re, basis_lines[iline],
                                                  'nterms (rexp1 expn1 coeff1) (rexp2 expn2 coeff2) ...')
            # Number of terms
            nterms = entry['nterms'][0]
            # Parse the data
            coeffs = entry['coeff']
            gexps = entry['expn']
            rexps = entry['rexp']

            # Read any continuation lines
            while len(rexps) < nterms:
                iline += 1
                entry = helpers.parse_line_regex_dict(ecp_block_cont_re, basis_lines[iline],
                                                      '(rexp3 expn3 coeff3) (rexp4 expn4 coeff4) ...')
                coeffs += entry['coeff']
                gexps += entry['expn']
                rexps += entry['rexp']

            ecpblocks[target_l] = [rexps, gexps, coeffs]

        # Return the parsed data
        return ecpblocks, iline

    # Read the ECP data
    ecpdata, iline = parse_ecp(basis_lines, iline, lmax, so=False)
    # Read the SO data
    ecpsodata, iline = parse_ecp(basis_lines, iline, lmaxso, so=True)

    # Count the number of values read in
    nread = 0
    for b in ecpdata:
        if not len(b):
            continue
        nread += 1 + 3 * len(b[0])
    for b in ecpsodata:
        if not len(b):
            continue
        nread += 1 + 3 * len(b[0])
    assert (nread == ndata)

    # Create entry
    element_Z = lut.element_Z_from_sym(element_sym, as_str=True)
    if element_Z not in bs_data or 'ecp_potentials' not in bs_data[element_Z]:
        element_data = helpers.create_element_data(bs_data, element_Z, 'ecp_potentials')
    else:
        element_data = bs_data[element_Z]

    # Store data
    element_data['ecp_electrons'] = ncore
    for il in range(len(ecpdata)):
        ecp_l = il
        r_exp = ecpdata[il][0]
        g_exp = ecpdata[il][1]
        coeff = ecpdata[il][2]
        ecp_pot = {
            'angular_momentum': [ecp_l],
            'ecp_type': 'scalar_ecp',
            'r_exponents': r_exp,
            'gaussian_exponents': g_exp,
            'coefficients': [coeff]
        }
        element_data['ecp_potentials'].append(ecp_pot)

    if lmaxso != 0:
        print('Warning: spin-orbit ECPs not yet supported in reader')

    return iline


def _parse_lines(basis_lines, bs_data):
    '''Parses lines representing all the electron shells for a single element

    Resulting information is stored in bs_data
    '''

    # Read the data one line at a time
    iline = 0
    while iline < len(basis_lines):
        if element_shell_re.match(basis_lines[iline]):
            iline = _read_shell(basis_lines, bs_data, iline)
        elif ecp_re.match(basis_lines[iline]):
            iline = _read_ecp(basis_lines, bs_data, iline)
        else:
            iline += 1


def read_libmol(basis_lines):
    '''Reads basis set from Molpro system library data and converts it to
       a dictionary with the usual BSE fields

       Note that the Molpro system library format does not store all
       the fields we have, so some fields are left blank

    '''

    # Removes comments
    basis_lines = helpers.prune_lines(basis_lines, '!')

    bs_data = {}

    # Empty file?
    if not basis_lines:
        return bs_data

    _parse_lines(basis_lines, bs_data)

    return bs_data
