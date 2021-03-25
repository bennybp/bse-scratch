'''
Conversion of basis sets to Gaussian format
'''

from .. import lut, manip, sort, printing


def _write_g94_common(basis, add_harm_type, psi4_am, system_library):
    '''Converts a basis set to Gaussian format
    '''

    s = ''

    basis = manip.uncontract_general(basis, True)
    basis = manip.uncontract_spdf(basis, 1, False)
    basis = sort.sort_basis(basis, False)

    # Elements for which we have electron basis
    electron_elements = [k for k, v in basis['elements'].items() if 'electron_shells' in v]

    # Elements for which we have ECP
    ecp_elements = [k for k, v in basis['elements'].items() if 'ecp_potentials' in v]

    # Electron Basis
    if electron_elements:
        for z in electron_elements:
            data = basis['elements'][z]

            sym = lut.element_sym_from_Z(z, True)
            s += '{}{}     0\n'.format('-' if system_library else '', sym)

            for shell in data['electron_shells']:
                exponents = shell['exponents']
                coefficients = shell['coefficients']
                ncol = len(coefficients) + 1
                nprim = len(exponents)

                am = shell['angular_momentum']
                amchar = lut.amint_to_char(am, hij=True).upper()

                if psi4_am and len(am) == 1 and am[0] >= 7:
                    # For am=7 and above, use explicit L={am} notation
                    amchar = "L=" + str(am[0])

                harm = ''
                if add_harm_type and shell['function_type'] == 'gto_cartesian':
                    harm = ' c'
                s += '{:4} {}   1.00{}\n'.format(amchar, nprim, harm)

                point_places = [8 * i + 15 * (i - 1) for i in range(1, ncol + 1)]
                s += printing.write_matrix([exponents, *coefficients], point_places, convert_exp=True)

            s += '****\n'

    # Write out ECP
    if ecp_elements:
        s += '\n'
        for z in ecp_elements:
            data = basis['elements'][z]
            sym = lut.element_sym_from_Z(z).upper()
            max_ecp_am = max([x['angular_momentum'][0] for x in data['ecp_potentials']])
            max_ecp_amchar = lut.amint_to_char([max_ecp_am], hij=True)

            # Sort lowest->highest, then put the highest at the beginning
            ecp_list = sorted(data['ecp_potentials'], key=lambda x: x['angular_momentum'])
            ecp_list.insert(0, ecp_list.pop())

            s += '{}     0\n'.format(sym)
            s += '{}-ECP     {}     {}\n'.format(sym, max_ecp_am, data['ecp_electrons'])

            for pot in ecp_list:
                rexponents = pot['r_exponents']
                gexponents = pot['gaussian_exponents']
                coefficients = pot['coefficients']
                nprim = len(rexponents)

                am = pot['angular_momentum']
                amchar = lut.amint_to_char(am, hij=True)

                if am[0] == max_ecp_am:
                    s += '{} potential\n'.format(amchar)
                else:
                    s += '{}-{} potential\n'.format(amchar, max_ecp_amchar)

                s += '  ' + str(nprim) + '\n'

                point_places = [0, 9, 32]
                s += printing.write_matrix([rexponents, gexponents, *coefficients], point_places, convert_exp=True)

    return s


def write_g94(basis):
    '''Converts a basis set to Gaussian format
    '''
    return _write_g94_common(basis, False, False, False)


def write_g94lib(basis):
    '''Converts a basis set to Gaussian system library format
    '''
    return _write_g94_common(basis, False, False, True)


def write_xtron(basis):
    '''Converts a basis set to xTron format

    xTron uses a modified gaussian format that puts 'c' on the same
    line as the angular momentum if the shell is cartesian.
    '''
    return _write_g94_common(basis, True, False, False)


def write_psi4(basis):
    '''Converts a basis set to Psi4 format

    Psi4 uses the same output as gaussian94, except
    that the first line must be cartesian/spherical,
    and it prefers to have a starting asterisks

    The cartesian/spherical line is added later, since it must
    be the first non-blank line.
    '''

    s = '****\n'
    s += _write_g94_common(basis, False, True, False)
    return s
