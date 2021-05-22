'''
Conversion of basis sets to Molcas basis_library format
'''

from .. import lut, manip, printing, misc, sort, api
import unidecode


def write_molcas_library(basis):
    '''Converts a basis set to Molcas basis_library format
    '''

    basis = manip.make_general(basis, False, True)
    basis = sort.sort_basis(basis, False)

    s = ''

    ref_data = api.get_reference_data(None)

    for z, data in basis['elements'].items():
        has_electron = 'electron_shells' in data
        has_ecp = 'ecp_potentials' in data

        el_name = lut.element_name_from_Z(z).upper()
        el_sym = lut.element_sym_from_Z(z, normalize=True)
        if 'names' in basis:
            bs_name = basis['names'][0]
        else:
            bs_name = basis['name']
        bs_name = bs_name.replace(' ', '_')
        cont = misc.contraction_string(data, compact=True)
        try:
            ref = data['references'][-1]['reference_keys'][-1]
        except (IndexError, KeyError):
            ref = None
        author = first_author(ref, ref_data)
        ecp = ''

        # number of electrons
        # should be z - number of ecp electrons
        nelectrons = int(z)
        if has_ecp:
            nelectrons -= data.get('ecp_electrons', 0)
            ecp = 'ECP.{}el.'.format(nelectrons)
        s += '/{}.{}.{}.{}.{}\n'.format(el_sym, bs_name, author, cont, ecp)
        s += '{}\n'.format(format_reference(ref, ref_data))
        s += '{} {}\n'.format(el_name, misc.contraction_string(data))

        # Are there cartesian shells?
        cartesian_shells = set()
        for shell in data['electron_shells']:
            if shell['function_type'] == 'gto_cartesian':
                for am in shell['angular_momentum']:
                    cartesian_shells.add(lut.amint_to_char([am]))
        if len(cartesian_shells):
            s += 'Options\n'
            s += 'Cartesian {}\n'.format(' '.join(cartesian_shells))
            s += 'EndOptions\n'

        if has_electron:
            max_am = misc.max_am(data['electron_shells'])

            s += '{:>7}.0   {}\n'.format(nelectrons, max_am)

            for shell in data['electron_shells']:
                exponents = shell['exponents']
                coefficients = shell['coefficients']
                nprim = len(exponents)
                ngen = len(coefficients)

                amchar = lut.amint_to_char(shell['angular_momentum']).lower()
                s += '* {}-type functions\n'.format(amchar)
                s += '{:>6}    {}\n'.format(nprim, ngen)

                s += printing.write_matrix([exponents], [17])

                point_places = [8 * i + 15 * (i - 1) for i in range(1, ngen + 1)]
                s += printing.write_matrix(coefficients, point_places)

        if has_ecp:
            max_ecp_am = max([x['angular_momentum'][0] for x in data['ecp_potentials']])

            # Sort lowest->highest, then put the highest at the beginning
            ecp_list = sorted(data['ecp_potentials'], key=lambda x: x['angular_momentum'])
            ecp_list.insert(0, ecp_list.pop())

            s += 'PP, {}, {}, {} ;\n'.format(el_sym, data['ecp_electrons'], max_ecp_am)

            for pot in ecp_list:
                rexponents = pot['r_exponents']
                gexponents = pot['gaussian_exponents']
                coefficients = pot['coefficients']

                am = pot['angular_momentum']
                amchar = lut.amint_to_char(am)
                s += '{};'.format(len(rexponents))

                if am[0] == max_ecp_am:
                    s += ' !  ul potential\n'
                else:
                    s += ' !  {}-ul potential\n'.format(amchar)

                for p in range(len(rexponents)):
                    s += '{},{},{};\n'.format(rexponents[p], gexponents[p], coefficients[0][p])

            s += 'Spectral Representation Operator\n'
            s += 'End of Spectral Representation Operator\n'
        s += '\n'

    return s


def first_author(ref, ref_data):
    if ref is None:
        return ''
    author = ''
    if ref == 'gaussian09e01':
        author = 'Gaussian'
    else:
        author = ref_data[ref].get('authors', [''])[0]
    author = author.split(',')[0]
    author = unidecode.unidecode(author)
    author = author.replace(' ', '_')
    return author


def format_author(author):
    text = ''
    parts = [i.strip() for i in author.split(',')]
    sep = ''
    first_name = parts[-1].split()
    if len(first_name) == 1:
        sep = '-'
        first_name = parts[-1].split(sep)
    for i,n in enumerate(first_name):
        if i > 0:
            text += sep
        if n.endswith('.'):
            text += n
        else:
            text += n[0] + '.'
    text += ' ' + ', '.join(parts[0:-1])
    return text


def format_reference(ref, ref_data):
    if ref is None:
        return 'Unknown reference'
    text = ''
    authors = []
    if len(ref_data[ref]['authors']) > 9:
        author = ref_data[ref]['authors'][0]
        authors.append(format_author(author))
        authors.append('et al.')
    else:
        for author in ref_data[ref]['authors']:
            authors.append(format_author(author))
    text += ', '.join(authors)
    if not text.endswith('.'):
        text += '.'
    if 'journal' in ref_data[ref]:
        text += ' ' + ref_data[ref]['journal']
        text += ' ' + ref_data[ref]['volume']
        text += ' (' + ref_data[ref]['year'] + ')'
        text += ' ' + ref_data[ref]['pages']
    elif 'booktitle' in ref_data[ref]:
        text += ' In "' + ref_data[ref]['booktitle'] + '"'
        text += ' (' + ref_data[ref]['year'] + ')'
        text += ' ' + ref_data[ref]['pages']
    else:
        text += ' ' + ref_data[ref]['title']
    if not text.endswith('.'):
        text += '.'
    if 'doi' in ref_data[ref]:
        text += ' doi:' + ref_data[ref]['doi'].lower()
    text = unidecode.unidecode(text)
    return text
