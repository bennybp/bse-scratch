'''
Functions for reading and writing the standard JSON-based
basis set format
'''

import json
import os
import collections

# Determine the path to the data directory
my_path = os.path.dirname(os.path.abspath(__file__))
data_path = os.path.join(my_path, 'data')


def sort_basis_dict(bs):
    '''Sorts a basis set dictionary into a standard order

    This allows the written file to be more easily read by humans by,
    for example, putting the name and description before more detailed fields
    '''

    keyorder = [
        'molssi_bse_magic',
        'basisSetName', 'basisSetDescription',
        'basisSetRole', 'basisSetElements', 'elementReferences', 'elementECPElectrons',
        'elementElectronShells', 'elementECP', 'elementComponents', 'elementEntry',
        'shellFunctionType', 'shellHarmonicType', 'shellRegion', 'shellAngularMomentum',
        'shellExponents', 'shellCoefficients', 'potentialECPType', 'potentialAngularMomentum',
        'potentialRExponents', 'potentialGaussianExponents', 'potentialCoefficients'
    ]

    # Add integers for the elements
    keyorder.extend(list(range(150)))

    bs_sorted = sorted(bs.items(), key=lambda x: keyorder.index(x[0]))
    bs_sorted = collections.OrderedDict(bs_sorted)

    for k, v in bs_sorted.items():
        if isinstance(v, dict):
            bs_sorted[k] = sort_basis_dict(v)
        elif k == 'elementElectronShells':
            bs_sorted[k] = [sort_basis_dict(x) for x in v]

    return bs_sorted


def check_compatible_merge(dest, source):
    '''TODO - check for any incompatibilities between the two elements
    '''
    pass


def merge_element_data(dest, sources):

    # return a shallow copy
    ret = {k: v for k, v in dest.items()}

    for s in sources:
        check_compatible_merge(dest, s)

        if 'elementElectronShells' in s:
            if not 'elementElectronShells' in ret:
                ret['elementElectronShells'] = []
            ret['elementElectronShells'].extend(s['elementElectronShells'])
        if 'elementECP' in s:
            if 'elementECP' in ret:
                raise RuntimeError('Cannot overwrite existing ECP')
            ret['elementECP'] = s['elementECP']
            ret['elementECPElectrons'] = s['elementECPElectrons']
        if 'elementReferences' in s:
            if not 'elementReferences' in ret:
                ret['elementReferences'] = []
            ret['elementReferences'].extend(s['elementReferences'])

    # Sort the shells
    ret['elementElectronShells'].sort(key=lambda x: x['shellAngularMomentum'])

    return ret


def read_json_by_path(file_path):

    if not os.path.isfile(file_path):
        raise RuntimeError('Basis set file \'{}\' does not exist, is not '
                           'readable, or is not a file'.format(file_path))

    with open(file_path, 'r') as f:
        js = json.loads(f.read())

    # Check for magic key/number
    if not 'molssi_bse_magic' in js:
        raise RuntimeError('This file does not appear to be a BSE JSON file')

    # change the element keys to integers
    js['basisSetElements'] = {int(k): v for k, v in js['basisSetElements'].items()}

    return js


def read_json_by_name(name, filetype=None):
    if filetype:
        name += '.' + filetype
    name += ".json"

    file_path = os.path.join(data_path, name)

    return read_json_by_path(file_path)


def dump_basis(bs):
    '''Returns a string with all the basis information (pretty-printed)
    '''
    return json.dumps(sort_basis_dict(bs), indent=4)


def write_basis_file(filepath, bs):
    '''Read a JSON basis set file to a given path

       The keys are first sorted into a standard order
    '''
    with open(filepath, 'w') as f:
        json.dump(sort_basis_dict(bs), f, indent=4)


def read_component_by_name(name):
    return read_json_by_name(name)


def read_elemental_basis_by_name(name):
    js = read_json_by_name(name, 'element')

    # construct a list of all files to read
    # TODO - can likely be replaced by memoization
    component_names = set()
    for k, v in js['basisSetElements'].items():
        component_names.update(set(v['elementComponents']))

    component_map = {k: read_json_by_name(k) for k in component_names}

    for k, v in js['basisSetElements'].items():

        components = v['elementComponents']

        v['elementElectronShells'] = []

        # all of the component data for this element
        el_data = [component_map[c]['basisSetElements'][k] for c in components]

        v = merge_element_data(v, el_data)

        # Remove the 'elementComponents' now that the data has been inserted
        v.pop('elementComponents')

        # Set it in the actual dict (v was a reference before)
        js['basisSetElements'][k] = v

    return js


def read_table_basis_by_name(name):
    js = read_json_by_name(name, 'table')

    # construct a list of all files to read
    # TODO - can likely be replaced by memoization
    component_names = set()
    for k, v in js['basisSetElements'].items():
        component_names.add(v['elementEntry'])

    component_map = {k: read_elemental_basis_by_name(k) for k in component_names}

    for k, v in js['basisSetElements'].items():
        entry = v['elementEntry']
        data = component_map[entry]

        # Replace the basis set for this element with the one
        # from the elemental basis
        js['basisSetElements'][k] = data['basisSetElements'][k]

    return js


def get_available_names():
    all_files = [x for x in os.listdir(data_path) if x.endswith('.table.json')]
    all_names = [x.replace('.table.json', '') for x in all_files]
    return sorted(all_names)
