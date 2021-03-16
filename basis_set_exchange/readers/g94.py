import re
from .. import lut, manip
from . import helpers

element_re = re.compile(r'^-?([A-Za-z]{1,3})(?:\s+0)?$')
ecp_shell_start_re = re.compile(r'^([A-Za-z]{1,3})\s+(\d+)\s+(\d+)$')
ecp_am_nelec_re = re.compile(r'^\S+\s+(\d+)\s+(\d+)$')

# This beast of a regex captures all the scaling factors as a single group:
#    '(?:{floating_re_str})+ ' is a non-capturing group of a number of floating point numbers.
#    Then, we capture all of them.
am_line_re = re.compile(r'^([A-Za-z]+)\s+(\d+)((?:\s+{})+)$'.format(helpers.floating_re_str))


def _parse_electron_lines(basis_lines, bs_data):
    '''Parses lines representing all the electron shells for a single element

    Resulting information is stored in bs_data
    '''

    # Last line should be "****"
    last = basis_lines.pop()
    if last != '****':
        raise RuntimeError("Electron shell is missing terminating ****")

    # First line is "{element} 0"
    element_sym = basis_lines[0].split()[0]

    # In the format, if the element symbol starts with a dash, then Gaussian does not crash
    # with an error in case this element does not exist in the molecule input
    # (this is what you need for system basis set libraries).
    element_sym = element_sym.lstrip('-')

    element_Z = lut.element_Z_from_sym(element_sym, as_str=True)
    element_data = manip.create_element_data(bs_data, element_Z, 'electron_shells')

    # After that come shells. We determine the start of a shell
    # by if the line starts with an angular momentum (a non-numeric character
    shell_blocks = helpers.partition_lines(basis_lines[1:], lambda x: x[0].isalpha())

    for sh_lines in shell_blocks:
        # Shells start with AM nprim scaling
        shell_am, nprim, scaling_factors = helpers.parse_line_regex(am_line_re, sh_lines[0],
                                                                    "Shell AM, nprim, scaling")
        shell_am = lut.amchar_to_int(shell_am, hij=True)
        func_type = lut.function_type_from_am(shell_am, 'gto', 'spherical')

        # Handle gaussian scaling factors
        # The square of the scaling factor is applied to exponents.
        # Typically they are 1.0, but not always
        scaling_factors = helpers.replace_d(scaling_factors)
        scaling_factors = [float(x) for x in scaling_factors.split()]

        # Remove any scaling factors that are 0.0
        scaling_factors = [x for x in scaling_factors if x != 0.0]

        # We should always have at least one scaling factor
        if len(scaling_factors) == 0:
            raise RuntimeError("No scaling factors given for element {}: Line: {}".format(element_sym, sh_lines[0]))

        # There can be multiple scaling factors, but we don't handle that. It seems to be very rare
        if len(scaling_factors) > 1:
            raise NotImplementedError("Number of scaling factors > 1")

        scaling_factor = float(scaling_factors[0])**2
        has_scaling = scaling_factor != 1.0

        # How many columns of coefficients do we have?
        # Gaussian doesn't support general contractions, so only >1 if
        # you have a fused shell
        ngen = len(shell_am)

        # Now read the exponents and coefficients
        exponents, coefficients = helpers.parse_primitive_matrix(sh_lines[1:], nprim, ngen)

        # If there is a scaling factor, apply it
        # But we keep track of some significant-figure type stuff (as best we can)
        if has_scaling:
            new_exponents = []
            for ex in exponents:
                ex = float(ex) * scaling_factor
                ex = '{:.16E}'.format(ex)

                # Trim useless zeroes
                ex_splt = ex.split('E')
                ex = ex_splt[0].rstrip('0')
                if ex[-1] == '.':  # Stripped all the zeroes...
                    ex += '0'
                ex += 'E' + ex_splt[1]
                new_exponents.append(ex)

            exponents = new_exponents

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

    # First line is "{element} 0", with the zero being optional
    element_sym = basis_lines[0].split()[0]
    element_Z = lut.element_Z_from_sym(element_sym, as_str=True)
    element_data = manip.create_element_data(bs_data, element_Z, 'ecp_potentials')

    # Second line is information about the ECP
    max_am, ecp_electrons = helpers.parse_line_regex(ecp_am_nelec_re, basis_lines[1], 'ECP max_am, nelec')

    element_data['ecp_electrons'] = ecp_electrons

    # Partition all the potentials
    # Look for lines containing only an integer, but include the line
    # before it (is a comment line)
    ecp_blocks = helpers.partition_lines(basis_lines[2:], helpers.is_integer, before=1)

    for pot_lines in ecp_blocks:
        # first line is comment
        # second line is number of lines

        # Check that the number of lines is consistent
        if not helpers.is_integer(pot_lines[1]):
            raise RuntimeError("Number of lines for potential is not an integer: " + pot_lines[1])

        nlines = int(pot_lines[1])
        if nlines <= 0:
            raise RuntimeError("Number of lines for potential is <= 0")

        if len(pot_lines) != (nlines + 2):
            raise RuntimeError("Number of lines is incorrect. Expected {}, got {}".format(nlines, len(pot_lines) - 2))

        ecp_data = helpers.parse_ecp_table(pot_lines[2:])

        ecp_pot = {
            'angular_momentum': None,
            'ecp_type': 'scalar_ecp',
            'r_exponents': ecp_data['r_exp'],
            'gaussian_exponents': ecp_data['g_exp'],
            'coefficients': ecp_data['coeff']
        }

        element_data['ecp_potentials'].append(ecp_pot)

    # Determine the AM of the potentials
    # Highest AM first, then the rest in order
    all_pot_am = helpers.potential_am_list(max_am)

    # Were there as many potentials as we thought there should be?
    if len(all_pot_am) != len(element_data['ecp_potentials']):
        raise RuntimeError("Found incorrect number of potentials for {}: Expected {}, got {}".format(
            element_sym, len(all_pot_am), len(element_data['ecp_potentials'])))

    for idx, pot in enumerate(element_data['ecp_potentials']):
        pot['angular_momentum'] = [all_pot_am[idx]]


def read_g94(basis_lines):
    '''Reads G94-formatted file data and converts it to a dictionary with the
       usual BSE fields

       Note that the gaussian format does not store all the fields we
       have, so some fields are left blank
    '''

    # Removes comments
    basis_lines = helpers.prune_lines(basis_lines, '!')

    bs_data = {}

    # Empty file?
    if not basis_lines:
        return bs_data

    # split into element sections (may be electronic or ecp)
    element_sections = helpers.partition_lines(basis_lines, element_re.match, min_size=3)

    for es in element_sections:
        # Try to guess if this is an ecp
        # Each element block starts with the element symbol
        # If the number of lines > 3, and the 4th line is just an integer, then it is an ECP
        if len(es) > 3 and helpers.is_integer(es[3]):
            _parse_ecp_lines(es, bs_data)
        else:
            _parse_electron_lines(es, bs_data)

    return bs_data
