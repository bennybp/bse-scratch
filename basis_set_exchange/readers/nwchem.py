import re
from .. import lut, manip
from . import helpers

am_line_re = re.compile(r'^([A-Za-z]+)\s+([A-Za-z]+)$')
nelec_re = re.compile(r'^([a-z]+)\s+nelec\s+(\d+)$', flags=re.IGNORECASE)


def _parse_electron_lines(basis_lines, bs_data):
    ''' Parses lines representing all the electron shells for all elements

    Resulting information is stored in bs_data
    '''

    # Remove 'end' from the lines (if they exist)
    # They may exist when this is called from other readers
    basis_lines = [x for x in basis_lines if x.lower() != 'end']

    # Basis entry needs to start with 'basis'
    if not basis_lines[0].lower().startswith('basis'):
        raise RuntimeError("Basis entry must start with 'basis'")

    # Is the basis set spherical or cartesian?
    am_type = 'cartesian' if basis_lines[0].lower().find('spherical') == -1 else 'spherical'

    # Start at index 1 in order to strip of the first line ('BASIS AO PRINT' or something)
    shell_blocks = helpers.partition_lines(basis_lines[1:], lambda x: x[0].isalpha(), min_size=2)

    for sh_lines in shell_blocks:
        element_sym, shell_am = helpers.parse_line_regex(am_line_re, sh_lines[0], "Element sym, shell am")
        shell_am = lut.amchar_to_int(shell_am)

        element_Z = lut.element_Z_from_sym(element_sym, as_str=True)
        element_data = manip.create_element_data(bs_data, element_Z, 'electron_shells', key_exist_ok=True)

        func_type = lut.function_type_from_am(shell_am, 'gto', am_type)

        # How many columns of coefficients do we have?
        # Only if this is a fused shell do we know
        ngen = len(shell_am) if len(shell_am) > 1 else None

        exponents, coefficients = helpers.parse_primitive_matrix(sh_lines[1:], ngen=ngen)

        shell = {
            'function_type': func_type,
            'region': '',
            'angular_momentum': shell_am,
            'exponents': exponents,
            'coefficients': coefficients
        }

        element_data['electron_shells'].append(shell)


def _parse_ecp_lines(basis_lines, bs_data):
    ''' Parses lines representing all the ECP potentials for all elements

    Resulting information is stored in bs_data
    '''

    # Remove 'end' from the lines (if they exist)
    # They may exist when this is called from other readers
    basis_lines = [x for x in basis_lines if x.lower() != 'end']

    # Start at index 1 in order to strip of the first line ('ECP' or something)
    # This splits out based on starting with an alpha character. A block can be either
    #   a potential or the nelec line
    ecp_blocks = helpers.partition_lines(basis_lines[1:], lambda x: x[0].isalpha())

    for pot_lines in ecp_blocks:
        # Check if this is the nelec line
        if len(pot_lines) == 1:
            # Check for this and give a better error message
            if not nelec_re.match(pot_lines[0]):
                raise RuntimeError("Unknown block with single line. Perhaps it's an empty block? Line: " +
                                   pot_lines[0])

            element_sym, n_elec = helpers.parse_line_regex(nelec_re, pot_lines[0].lower(), "ECP: Element sym, nelec")

            element_Z = lut.element_Z_from_sym(element_sym, as_str=True)
            element_data = manip.create_element_data(bs_data, element_Z, 'ecp_electrons', create=int)
            element_data['ecp_electrons'] = n_elec
        else:
            element_sym, pot_am = helpers.parse_line_regex(am_line_re, pot_lines[0], "ECP: Element sym, pot AM")

            element_Z = lut.element_Z_from_sym(element_sym, as_str=True)
            element_data = manip.create_element_data(bs_data, element_Z, 'ecp_potentials', key_exist_ok=True)

            # See if this is the 'ul' angular momentum
            if pot_am.lower() == 'ul':
                pot_am = []  # Placeholder - leave for later
            else:
                pot_am = lut.amchar_to_int(pot_am)

            ecp_data = helpers.parse_ecp_table(pot_lines[1:])
            ecp_pot = {
                'angular_momentum': pot_am,
                'ecp_type': 'scalar_ecp',
                'r_exponents': ecp_data['r_exp'],
                'gaussian_exponents': ecp_data['g_exp'],
                'coefficients': ecp_data['coeff']
            }

            element_data['ecp_potentials'].append(ecp_pot)

    # Fix ecp angular momentum now that everything has been read
    # Specifically, we can set the 'ul' potential to be the max am + 1
    for el, v in bs_data.items():
        if 'ecp_potentials' not in v:
            continue

        all_ecp_am = []
        for x in v['ecp_potentials']:
            all_ecp_am.extend(x['angular_momentum'])

        max_ecp_am = max(all_ecp_am)

        for s in v['ecp_potentials']:
            if not s['angular_momentum']:
                s['angular_momentum'] = [max_ecp_am + 1]

    # Make sure the number of electrons replaced by the ECP was specified for all elements
    for el, v in bs_data.items():
        if 'ecp_potentials' in v and 'ecp_electrons' not in v:
            raise RuntimeError("Number of ECP electrons not specified for element {}".format(el))


def read_nwchem(basis_lines):
    '''Reads NWChem-formatted file data and converts it to a dictionary with the
       usual BSE fields

       Note that the nwchem format does not store all the fields we
       have, so some fields are left blank
    '''

    basis_lines = helpers.prune_lines(basis_lines, '#')

    bs_data = {}

    # split into basis and ecp
    basis_sections = helpers.partition_lines(basis_lines,
                                             lambda x: x.lower() == 'end',
                                             min_blocks=1,
                                             max_blocks=2,
                                             include_match=False)

    for s in basis_sections:
        if s[0].lower().startswith('basis'):
            _parse_electron_lines(s, bs_data)
        elif s[0].lower().startswith('ecp'):
            _parse_ecp_lines(s, bs_data)
        else:
            raise RuntimeError("Unknown section: " + s[0])

    return bs_data
