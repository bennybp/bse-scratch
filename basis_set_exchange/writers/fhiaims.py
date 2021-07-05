'''
Conversion of basis sets to FHI-aims format
'''

from .. import lut, manip, sort, printing

def write_fhiaims(basis):
    '''Converts a basis set to FHI-aims format
    '''

    # Angular momentum type
    types = basis['function_types']
    pure = False if 'gto_cartesian' in types else True

    # Set up
    s = ''
    basis = manip.uncontract_general(basis, True)
    basis = manip.uncontract_spdf(basis, 1, False)
    basis = sort.sort_basis(basis, False)

    # Elements for which we have electron basis
    electron_elements = [k for k, v in basis['elements'].items() if 'electron_shells' in v]

    # Electron Basis
    if electron_elements:
        for z in electron_elements:
            data = basis['elements'][z]

            # FHI-aims defines elements in species default files.
            # The options need to be specified for every element basis.
            s += '''
            # The default minimal basis should not be included
            include_min_basis .false.
            # Use spherical functions?
            pure_gauss {}
            '''.format('.true.' if pure else '.false.')
            
            sym = lut.element_sym_from_Z(z, True)
            s += '# {} {}\n'.format(sym, basis['name'])

            for shell in data['electron_shells']:
                exponents = shell['exponents']
                coefficients = shell['coefficients']
                ncol = len(coefficients) + 1
                nprim = len(exponents)
                am = shell['angular_momentum']
                assert len(am) == 1
                
                if nprim == 1:
                    s += 'gaussian {} {} {}\n'.format(am[0], nprim, exponents[0])
                else:
                    s += 'gaussian {} {}\n'.format(am[0], nprim)
                    point_places = [8 * i + 15 * (i - 1) for i in range(1, ncol + 1)]
                    s += printing.write_matrix([exponents, *coefficients], point_places, convert_exp=False)

    return s
