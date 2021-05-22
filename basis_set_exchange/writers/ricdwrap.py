'''
This is a wrapper for generating acCD basis sets with OpenMolcas
'''

from .. import lut, manip, printing, misc, sort


def write_ricdwrap(basis):
    '''Generates an input file for OpenMolcas that generates the acCD basis
    '''

    basis = manip.make_general(basis, False, True)
    basis = sort.sort_basis(basis, False)

    s = '''
&GATEWAY
  ricd
  accd
  cdthreshold=1.0d-4
'''

    for z, data in basis['elements'].items():
        s += 'Basis set\n'
        has_electron = 'electron_shells' in data

        el_name = lut.element_name_from_Z(z).upper()
        el_sym = lut.element_sym_from_Z(z, normalize=True)
        s += '* {}  {}\n'.format(el_name, misc.contraction_string(data))

        s += ' {}    / inline\n'.format(el_sym)

        if has_electron:
            max_am = misc.max_am(data['electron_shells'])

            # number of electrons
            # should be z - number of ecp electrons
            nelectrons = int(z)

            s += '{:>7}.00   {}\n'.format(nelectrons, max_am)

            for shell in data['electron_shells']:
                exponents = shell['exponents']
                coefficients = shell['coefficients']
                nprim = len(exponents)
                ngen = len(coefficients)

                amchar = lut.amint_to_char(shell['angular_momentum']).upper()
                s += '* {}-type functions\n'.format(amchar)
                s += '{:>6}    {}\n'.format(nprim, ngen)

                s += printing.write_matrix([exponents], [17])

                point_places = [8 * i + 15 * (i - 1) for i in range(1, ngen + 1)]
                s += printing.write_matrix(coefficients, point_places)

        # Make a nucleus; use 10 angstrom distance
        s += '{} 0.0 0.0 {:.1f}\n'.format(el_sym, 10.0 * (int(z) - 1))

        # Are there cartesian shells?
        if has_electron:
            cartesian_shells = set()
            for shell in data['electron_shells']:
                if shell['function_type'] == 'gto_cartesian':
                    for am in shell['angular_momentum']:
                        cartesian_shells.add(lut.amint_to_char([am]))
            if len(cartesian_shells):
                s += 'cartesian {}\n'.format(' '.join(cartesian_shells))

        s += 'End of basis set\n\n'

    return s
