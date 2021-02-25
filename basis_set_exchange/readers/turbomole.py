import re
from .. import lut
from . import helpers

section_re = re.compile(r'^\$(basis|ecp|cbas|jbas|jkbas)$')
element_re = re.compile(r'^([a-zA-Z]{1,3})\s+(.*)$')
shell_re = re.compile(r'^(\d+) +([a-zA-Z])$')
ecp_info_re = re.compile(r'^ncore\s*=\s*(\d+)\s+lmax\s*=\s*(\d+)$', flags=re.IGNORECASE)
ecp_pot_am_re = re.compile(r'^([a-z])(-[a-z])?$')


def _parse_electron_lines(basis_lines, bs_data):
    # Strip all lines beginning with $
    basis_lines = helpers.prune_lines(basis_lines, '$')

    # Last line should be *
    # We don't need it
    if basis_lines[-1] != '*':
        raise RuntimeError("Missing terminating * line")
    basis_lines.pop()

    # Partition based on lines beginning with a character
    element_blocks = helpers.partition_lines(basis_lines, element_re.match, before=1, min_size=4)

    # Element lines should be surrounded by *
    # Check all first. the partition_lines above will eat part of the previous
    # element if the * is missing
    for element_lines in element_blocks:
        if element_lines[0] != '*':
            raise RuntimeError("Element line not preceded by *")
        if element_lines[2] != '*':
            raise RuntimeError("Element line not followed by *")

        # Check for any other lines starting with *
        for line in element_lines[3:]:
            if line.startswith('*'):
                raise RuntimeError("Found line starting with * that probably doesn't belong: " + line)

    # Now process them all
    for element_lines in element_blocks:
        element_sym, _ = helpers.parse_line_regex(element_re, element_lines[1], 'Element line')

        element_Z = lut.element_Z_from_sym(element_sym, as_str=True)
        element_data = helpers.create_element_data(bs_data, element_Z, 'electron_shells')

        # Partition into shells
        shell_blocks = helpers.partition_lines(element_lines[3:], shell_re.match, min_size=2)

        for sh_lines in shell_blocks:
            nprim, shell_am = helpers.parse_line_regex(shell_re, sh_lines[0], 'shell nprim, am')
            shell_am = lut.amchar_to_int(shell_am)

            func_type = helpers.function_type_from_am(shell_am, 'gto', 'spherical')

            # Primitive matrix
            prim_mat = sh_lines[1:1+nprim]
            # The turbomole format appears to allow an enumeration of the exponents; if there are three fields, we drop the first one
            for iline in range(len(prim_mat)):
                entries = prim_mat[iline].split()
                if len(entries) == 3:
                    if int(entries[0]) != 1+iline:
                        raise RuntimeError('Enumerated format detected, but numbering is off on line: ', prim_mat[iline])
                    prim_mat[iline] = ' '.join(entries[1:])

            exponents, coefficients = helpers.parse_primitive_matrix(prim_mat, nprim=nprim, ngen=1)

            shell = {
                'function_type': func_type,
                'region': '',
                'angular_momentum': shell_am,
                'exponents': exponents,
                'coefficients': coefficients
            }

            element_data['electron_shells'].append(shell)


def _parse_ecp_potential_lines(element_lines, bs_data):
    #########################################################
    # This is split out because the turbomole ECP format is
    # almost identical to the genbas ECP format
    #########################################################
    element_sym, _ = helpers.parse_line_regex(element_re, element_lines[0], 'Element line')

    element_Z = lut.element_Z_from_sym(element_sym, as_str=True)

    # We don't need the return value - we will use the one from creating ecp_electrons
    helpers.create_element_data(bs_data, element_Z, 'ecp_potentials')

    # 4th line should be ncore and lmax
    n_elec, max_am = helpers.parse_line_regex(ecp_info_re, element_lines[1], 'ECP ncore, lmax')

    element_data = helpers.create_element_data(bs_data, element_Z, 'ecp_electrons', key_exist_ok=False, create=int)
    element_data['ecp_electrons'] = n_elec

    # split the remaining lines by lines starting with a character
    ecp_potentials = helpers.partition_lines(element_lines[2:], lambda x: x[0].isalpha(), min_size=2)

    # Keep track of what the max AM we actually found is
    found_max = False
    for pot_lines in ecp_potentials:
        pot_am, pot_base_am = helpers.parse_line_regex(ecp_pot_am_re, pot_lines[0], 'ECP potential am')

        pot_am = lut.amchar_to_int(pot_am)

        if pot_base_am:
            pot_base_am = lut.amchar_to_int(pot_base_am[1:])  # Strip the - from the beginning

            if pot_base_am[0] != max_am:
                raise RuntimeError("Potential does not use max_am of {}. Uses {}".format(max_am, pot_base_am[0]))
        else:
            if found_max:
                raise RuntimeError("Found multiple potentials with single AM")

            if pot_am[0] != max_am:
                raise RuntimeError("Potential with single AM {} is not the same as lmax = {}".format(
                    pot_am[0], max_am))

            found_max = True

        ecp_data = helpers.parse_ecp_table(pot_lines[1:], order=['coeff', 'r_exp', 'g_exp'])
        ecp_pot = {
            'angular_momentum': pot_am,
            'ecp_type': 'scalar_ecp',
            'r_exponents': ecp_data['r_exp'],
            'gaussian_exponents': ecp_data['g_exp'],
            'coefficients': ecp_data['coeff']
        }

        element_data['ecp_potentials'].append(ecp_pot)


def _parse_ecp_lines(basis_lines, bs_data):
    # Strip all lines beginning with $
    basis_lines = helpers.prune_lines(basis_lines, '$')

    # Last line should be *
    # We don't need it
    if basis_lines[-1] != '*':
        raise RuntimeError("Missing terminating * line")
    basis_lines.pop()

    # Partition based on lines beginning with a character
    element_blocks = helpers.partition_lines(basis_lines, element_re.match, before=1)

    # Element lines should be surrounded by *
    # Check all first. the partition_lines above will eat part of the previous
    # element if the * is missing
    for element_lines in element_blocks:
        if element_lines[0] != '*':
            raise RuntimeError("Element line not preceded by *")
        if element_lines[2] != '*':
            raise RuntimeError("Element line not followed by *")

        # Check for any other lines starting with *
        for line in element_lines[3:]:
            if line.startswith('*'):
                raise RuntimeError("Found line starting with * that probably doesn't belong: " + line)

    # Now process all elements
    for element_lines in element_blocks:
        # Remove the two * lines and parse using the separate function
        element_lines = element_lines[1:2] + element_lines[3:]
        _parse_ecp_potential_lines(element_lines, bs_data)


def read_turbomole(basis_lines):
    '''Reads turbomole-formatted file data and converts it to a dictionary with the
       usual BSE fields

       Note that the turbomole format does not store all the fields we
       have, so some fields are left blank
    '''

    basis_lines = helpers.prune_lines(basis_lines, '#')

    # first line must begin with $, last line must be $end
    if basis_lines and not basis_lines[0][0].startswith('$'):
        raise RuntimeError("First line does not begin with $. Line: " + basis_lines[0])
    if basis_lines and basis_lines[-1] != '$end':
        raise RuntimeError("Last line of basis is not $end. Line: " + basis_lines[-1])

    bs_data = {}

    # Split into basis and ecp
    # Just split based on lines beginning with $
    basis_sections = helpers.partition_lines(basis_lines,
                                             lambda x: x.startswith('$') and x != '$end',
                                             min_blocks=1,
                                             max_blocks=2)

    for s in basis_sections:
        # Check if section is empty. If so, all lines start with $
        if all(x.startswith('$') for x in s):
            continue

        if len(s) == 0:  # Empty section
            continue
        elif s[0].lower() == '$ecp':
            _parse_ecp_lines(s, bs_data)
        elif section_re.match(s[0]):
            _parse_electron_lines(s, bs_data)
        else:
            raise RuntimeError("Unknown section " + s[0])

    return bs_data
