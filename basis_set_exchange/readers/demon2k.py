# Written by Susi Lehtola, 2021

import re
from .. import lut, manip
from . import helpers

# Orbital entry: O-ELEMENT [possible repetition(s) of element symbol] (name)
orbital_re = re.compile(r'^O-([A-Za-z]+)(?: ([A-Za-z]+))* \(({})\)\s*$'.format(helpers.basis_name_re_str))
# Shell entry: formal_quantum_number am nprim
shell_re = re.compile(r'^\s*(\d+)\s+(\d+)\s+(\d+)\s*$')

# ECP block start: ECP
ecp_start_re = re.compile(r'^\s*ECP\s*$')
# ECP entry: symbol nelec N
ecp_entry_re = re.compile(r'^([A-Za-z]+)\s+nelec\s+(\d+)\s*$')
# ECP shell entry: symbol am
ecp_shell_re = re.compile(r'^([A-Za-z]+)\s+([A-Za-z]+)\s*$')
# ECP data entry: rexp gexp gcoeff
ecp_data_re = re.compile(r'^\s*(\d+)\s+({0})\s+({0})\s*$'.format(helpers.floating_re_str))

basis_end_re = re.compile(r'^\s*END\s*$')


def _parse_electron_lines(basis_lines, bs_data):
    '''Parses lines representing all the electron shells for a single element

    Resulting information is stored in bs_data
    '''

    # First line should be element and name of basis
    assert orbital_re.match(basis_lines[0])
    element_name, repetitions, basis_name = helpers.parse_line_regex(orbital_re, basis_lines[0],
                                                                     "O-element (el?) name")

    # Create basis
    element_Z = lut.element_Z_from_name(element_name, as_str=True)
    element_data = manip.create_element_data(bs_data, element_Z, 'electron_shells')

    # Second line is number of shells
    n_shells = int(basis_lines[1].split()[0])

    # After that come the shells.
    shell_blocks = helpers.partition_lines(basis_lines[2:], shell_re.match)
    if len(shell_blocks) != n_shells:
        raise RuntimeError('Expected {} shells, got {}'.format(n_shells, len(shell_blocks)))

    for sh_lines in shell_blocks:
        # Check that this is a shell block
        if not shell_re.match(sh_lines[0]):
            continue
        formal_n, shell_am, nprim = helpers.parse_line_regex(shell_re, sh_lines[0], "formal_n am nprim")

        # BSE expects this as an array
        shell_am = [shell_am]

        # Read the exponents and coefficients
        ngen = 1
        exponents, coefficients = helpers.parse_primitive_matrix(sh_lines[1:1 + nprim], nprim, ngen)

        # Function type (assuming always spherical)
        func_type = lut.function_type_from_am(shell_am, 'gto', 'spherical')
        shell = {
            'function_type': func_type,
            'region': '',
            'angular_momentum': shell_am,
            'exponents': exponents,
            'coefficients': coefficients
        }

        element_data['electron_shells'].append(shell)


def _parse_ecp_lines(ecp_lines, bs_data):
    '''Parses lines representing all the ECP potentials for a single element

    Resulting information is stored in bs_data
    '''

    # Split by symbol
    element_ecps = helpers.partition_lines(ecp_lines, ecp_entry_re.match)
    for element_lines in element_ecps:
        # Check that this is an ECP block
        if not ecp_entry_re.match(element_lines[0]):
            continue
        # Block must have at least three lines
        if len(element_lines) < 3:
            continue

        # First line is "{element} nelec N"
        element_sym, ecp_electrons = helpers.parse_line_regex(ecp_entry_re, element_lines[0], "symbol nelec n")

        # Initialize basis
        element_Z = lut.element_Z_from_sym(element_sym, as_str=True)
        element_data = manip.create_element_data(bs_data, element_Z, 'ecp_potentials')
        element_data['ecp_electrons'] = ecp_electrons

        # Read angular momentum blocks
        amblocks = []
        current_block = None
        iline = 1
        while iline < len(element_lines):
            if ecp_shell_re.match(element_lines[iline]):
                # Start new shell
                if current_block is not None:
                    amblocks.append(current_block)

                # Parse entry
                ecp_element_sym, ecp_am = helpers.parse_line_regex(ecp_shell_re, element_lines[iline], "symbol ecp-am")
                if ecp_element_sym != element_sym:
                    raise RuntimeError('Expected ECP for {} not {}\n'.format(element_sym, ecp_element_sym))

                current_block = [ecp_am, []]

            elif ecp_data_re.match(element_lines[iline]):
                # Append to current block
                rexp, gexp, gcoeff = helpers.parse_line_regex(ecp_data_re, element_lines[iline], "rexp gexp gcoeff")
                current_block[1].append([rexp, gexp, gcoeff])

            elif basis_end_re.match(element_lines[iline]):
                break

            else:
                raise RuntimeError('Unexpected format of ECP block!\n')

            iline += 1

        # The last entry has not been added yet
        if current_block is not None:
            amblocks.append(current_block)

        # Loop over the angular momentum blocks
        Nblocks = len(amblocks)
        for iblock in range(Nblocks):
            block = amblocks[iblock]

            # First entry is highest projector, then S, P, ...
            if iblock == 0:
                assert block[0] == 'ul'
                current_am = Nblocks - 1
            else:
                current_am = iblock - 1
                assert lut.amint_to_char([current_am])[0].lower() == block[0].lower()

            # Collect the r exponents
            r_exp = [x[0] for x in block[1]]
            g_exp = [x[1] for x in block[1]]
            g_coeff = [[x[2] for x in block[1]]]

            # Store the ECP
            ecp_pot = {
                'angular_momentum': [current_am],
                'ecp_type': 'scalar_ecp',
                'r_exponents': r_exp,
                'gaussian_exponents': g_exp,
                'coefficients': g_coeff
            }
            element_data['ecp_potentials'].append(ecp_pot)


def read_demon2k(basis_lines):
    '''Reads deMon2k formatted file data and converts it to a dictionary
       with the usual BSE fields

       Note that the deMon2k format does not store all the fields we
       have, so some fields are left blank

    '''

    # Removes comments
    basis_lines = helpers.prune_lines(basis_lines, '#')

    bs_data = {}

    # Empty file?
    if not basis_lines:
        return bs_data

    # split orbital basis into element sections
    orbital_sections = helpers.partition_lines(basis_lines, orbital_re.match, min_size=3)
    for es in orbital_sections:
        _parse_electron_lines(es, bs_data)

    # split orbital basis into element sections
    ecp_sections = helpers.partition_lines(basis_lines, ecp_start_re.match, min_size=3)
    for ecp in ecp_sections:
        _parse_ecp_lines(ecp, bs_data)

    # Last line should be "END"
    last = basis_lines.pop()
    if last != 'END':
        raise RuntimeError("Basis set is missing terminating END")

    return bs_data
