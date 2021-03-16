'''
Main interface to Basis Set Exchange internal basis sets

This module contains the interface for getting basis set data
and references from the internal data store of basis sets. As much
as possible, this is being kept separate from the typical reading/writing
functionality of the library.
'''

import os
import textwrap

from . import compose
from . import writers
from . import fileio
from . import manip
from . import memo
from . import notes
from . import refconverters
from . import references
from . import sort
from . import misc
from . import lut
from ._version import get_versions

__version__ = get_versions()['version']


def version():
    '''Obtain the version of the basis set exchange library (as a string)'''
    return __version__


# Determine the path to the data directory that is part
# of this installation
_my_dir = os.path.dirname(os.path.abspath(__file__))
_default_data_dir = os.path.join(_my_dir, 'data')

# Main URL of the project
_main_url = 'https://www.basissetexchange.org'


def _get_basis_metadata(name, data_dir):
    '''Get metadata for a single basis set

    If the basis doesn't exist, an exception is raised
    '''

    # Transform the name into an internal representation
    tr_name = misc.transform_basis_name(name)

    # Get the metadata for all basis sets
    metadata = get_metadata(data_dir)

    if tr_name not in metadata:
        raise KeyError("Basis set {} does not exist".format(name))

    return metadata[tr_name]


def _header_string(basis_dict):
    '''Creates a header with information about a basis set

    Information includes description, revision, etc, but not references
    '''

    tw = textwrap.TextWrapper(initial_indent='', subsequent_indent=' ' * 20)

    header = '-' * 70 + '\n'
    header += ' Basis Set Exchange\n'
    header += ' Version ' + version() + '\n'
    header += ' ' + _main_url + '\n'
    header += '-' * 70 + '\n'
    header += '   Basis set: ' + basis_dict['name'] + '\n'
    header += tw.fill(' Description: ' + basis_dict['description']) + '\n'
    header += '        Role: ' + basis_dict['role'] + '\n'
    header += tw.fill('     Version: {}  ({})'.format(basis_dict['version'],
                                                      basis_dict['revision_description'])) + '\n'
    header += '-' * 70 + '\n'

    return header


def fix_data_dir(data_dir):
    '''
    If data_dir is None, returns the default data_dir. Otherwise,
    returns the data_dir parameter unmodified
    '''

    return _default_data_dir if data_dir is None else data_dir


def get_basis(name,
              elements=None,
              version=None,
              fmt=None,
              uncontract_general=False,
              uncontract_spdf=False,
              uncontract_segmented=False,
              remove_free_primitives=False,
              make_general=False,
              optimize_general=False,
              augment_diffuse=0,
              augment_steep=0,
              get_aux=0,
              data_dir=None,
              header=True):
    '''Obtain a basis set

    This is the main function for getting basis set information.
    This function reads in all the basis data and returns it either
    as a string or as a python dictionary.

    Parameters
    ----------
    name : str
        Name of the basis set. This is not case sensitive.
    elements : str or list
        List of elements that you want the basis set for.
        Elements can be specified by Z-number (int or str) or by symbol (str).
        If this argument is a str (ie, '1-3,7-10'), it is expanded into a list.
        Z numbers and symbols (case insensitive) can be used interchangeably
        (see :func:`basis_set_exchange.misc.expand_elements`)

        If an empty string or list is passed, or if None is passed (the default),
        all elements for which the basis set is defined are included.
    version : int or str
        Obtain a specific version of this basis set. By default,
        the latest version is returned.
    fmt: str
        The desired output format of the basis set. By default,
        basis set information is returned as a python dictionary. Otherwise,
        if a format is specified, a string is returned.
        Use :func:`basis_set_exchange.api.get_formats` to programmatically obtain the available
        formats.  The `fmt` argument is not case sensitive.

        Available formats are

            * bdf
            * gamess_us
            * gaussian94
            * json
            * nwchem
            * psi4
            * turbomole

    uncontract_general : bool
        If True, remove general contractions by duplicating the set
        of primitive exponents with each vector of coefficients.
        Primitives with zero coefficient are removed, as are duplicate shells.
    uncontract_spdf : bool
        If True, remove general contractions with combined angular momentum (sp, spd, etc)
        by duplicating the set of primitive exponents with each vector of coefficients.
        Primitives with zero coefficient are removed, as are duplicate shells.
    uncontract_segmented : bool
        If True, remove segmented contractions by duplicating each primitive into new shells.
        Each coefficient is set to 1.0
    make_general : bool
        If True, make the basis set as generally-contracted as possible. There will be one
        shell per angular momentum (for each element)
    optimize_general : bool
        Optimize by removing general contractions that contain uncontracted
        functions (see :func:`basis_set_exchange.manip.optimize_general`)
    augment_diffuse : int
        Add n diffuse functions by even-tempered extrapolation
    augment_steep : int
        Add n steep functions by even-tempered extrapolation
    get_aux : int
        Instead of the orbital basis, get an auxiliary basis
        set. Options 0 (return orbital basis), 1 (return AutoAux
        basis), 2 (return Auto-ABS Coulomb fitting basis)
    data_dir : str
        Data directory with all the basis set information. By default,
        it is in the 'data' subdirectory of this project.

    Returns
    -------
    str or dict
        The basis set in the desired format. If `fmt` is **None**, this will be a python
        dictionary. Otherwise, it will be a string.

    '''

    data_dir = fix_data_dir(data_dir)
    bs_data = _get_basis_metadata(name, data_dir)

    # If version is not specified, use the latest
    if version is None:
        version = bs_data['latest_version']
    else:
        version = str(version)  # Version may be an int

    if version not in bs_data['versions']:
        raise KeyError("Version {} does not exist for basis {}".format(version, name))

    # Compose the entire basis set (all elements)
    file_relpath = bs_data['versions'][version]['file_relpath']
    basis_dict = compose.compose_table_basis(file_relpath, data_dir)

    # Set the name (from the global metadata)
    # Only the list of all names will be returned from compose_table_basis
    basis_dict['name'] = bs_data['display_name']

    # Handle optional arguments
    if elements is not None:
        # Convert to purely a list of strings that represent integers
        elements = misc.expand_elements(elements, True)

        # Did the user pass an empty string or empty list? If so, include
        # all elements
        if elements:
            bs_elements = basis_dict['elements']

            # Are elements part of this basis set?
            for el in elements:
                if el not in bs_elements:
                    elsym = lut.element_sym_from_Z(el)
                    raise KeyError("Element {} (Z={}) not found in basis {} version {}".format(
                        elsym, el, name, version))

            # Set to only the elements we want
            basis_dict['elements'] = {k: v for k, v in bs_elements.items() if k in elements}

            # Since we only grab some of the elements, we need to
            # update the function types used, too
            basis_dict['function_types'] = compose._whole_basis_types(basis_dict)

    # Note that from now on, the pipleline is going to modify basis_dict. That is ok,
    # since we are returned a unique instance from compose_table_basis

    needs_pruning = False

    if remove_free_primitives:
        basis_dict = manip.remove_free_primitives(basis_dict, False)
        needs_pruning = True

    if optimize_general:
        basis_dict = manip.optimize_general(basis_dict, False)
        needs_pruning = True

    # uncontract_segmented implies uncontract_general
    if uncontract_segmented:
        basis_dict = manip.uncontract_segmented(basis_dict, False)
        needs_pruning = True

    elif uncontract_general:
        basis_dict = manip.uncontract_general(basis_dict, False)
        needs_pruning = True

    if uncontract_spdf:
        basis_dict = manip.uncontract_spdf(basis_dict, 0, False)
        needs_pruning = True

    if make_general:
        basis_dict = manip.make_general(basis_dict, False, False)
        needs_pruning = True

    # Remove dead and duplicate shells
    if needs_pruning:
        basis_dict = manip.prune_basis(basis_dict, False)

    # Augment
    if augment_diffuse > 0:
        basis_dict = manip.geometric_augmentation(basis_dict,
                                                  augment_diffuse,
                                                  use_copy=False,
                                                  as_component=False,
                                                  steep=False)
    if augment_steep > 0:
        basis_dict = manip.geometric_augmentation(basis_dict,
                                                  augment_steep,
                                                  use_copy=False,
                                                  as_component=False,
                                                  steep=True)
        # Need to sort to get added steep functions first
        basis_dict = sort.sort_basis(basis_dict)
    # Re-make general
    if (augment_diffuse > 0 or augment_steep > 0) and make_general:
        basis_dict = manip.make_general(basis_dict, False, False)

    # Did we actually want an auxiliary basis set?
    if get_aux == 1:
        basis_dict = manip.autoaux_basis(basis_dict)
    elif get_aux == 2:
        basis_dict = manip.autoabs_basis(basis_dict)

    # If fmt is not specified, return as a python dict
    if fmt is None:
        return basis_dict

    if header:
        header_str = _header_string(basis_dict)
    else:
        header_str = None

    return writers.write_formatted_basis_str(basis_dict, fmt, header_str)


def lookup_basis_by_role(primary_basis, role, data_dir=None):
    '''Lookup the name of an auxiliary basis set given a primary basis set and role

    Parameters
    ----------
    primary_basis : str
        The primary (orbital) basis set that we want the auxiliary
        basis set for. This is not case sensitive.
    role: str
        Desired role/type of auxiliary basis set.
        Use :func:`get_roles` to programmatically obtain the available
        roles.  The `role` argument is not case sensitive.

        Available roles are

            * jfit
            * jkfit
            * rifit
            * optri
            * admmfit
            * dftxfit
            * dftjfit
            * guess

    data_dir : str
        Data directory with all the basis set information. By default,
        it is in the 'data' subdirectory of this project.

    Returns
    -------
    str
        The name of the auxiliary basis set for the given primary basis
        and role.
    '''

    data_dir = fix_data_dir(data_dir)

    role = role.lower()

    if role not in get_roles():
        raise RuntimeError("Role {} is not a valid role".format(role))

    bs_data = _get_basis_metadata(primary_basis, data_dir)
    auxdata = bs_data['auxiliaries']

    if role not in auxdata:
        raise RuntimeError("Role {} doesn't exist for {}".format(role, primary_basis))

    return auxdata[role]


@memo.BSEMemoize
def get_metadata(data_dir=None):
    '''Obtain the metadata for all basis sets

    The metadata includes information such as the display name of the basis set,
    its versions, and what elements are included in the basis set

    The data is read from the METADATA.json file in the `data_dir` directory.

    Parameters
    ----------
    data_dir : str
        Data directory with all the basis set information. By default,
        it is in the 'data' subdirectory of this project.
    '''

    data_dir = fix_data_dir(data_dir)
    metadata_file = os.path.join(data_dir, "METADATA.json")
    return fileio.read_metadata(metadata_file)


@memo.BSEMemoize
def get_reference_data(data_dir=None):
    '''Obtain information for all stored references

    This is a nested dictionary with all the data for all the references

    The reference data is read from the REFERENCES.json file in the given
    `data_dir` directory.
    '''

    data_dir = fix_data_dir(data_dir)
    reffile_path = os.path.join(data_dir, 'REFERENCES.json')

    return fileio.read_references(reffile_path)


def get_all_basis_names(data_dir=None):
    '''Obtain a list of all basis set names

    The returned list is the internal representation of the basis set name.

    The data is read from the METADATA.json file in the data directory.

    Parameters
    ----------
    data_dir : str
        Data directory with all the basis set information. By default,
        it is in the 'data' subdirectory of this project.
    '''

    md = get_metadata(data_dir)
    names = list(v['display_name'] for v in md.values())
    return sorted(names)


def get_references(basis_name, elements=None, version=None, fmt=None, data_dir=None):
    '''Get the references/citations for a basis set

    Parameters
    ----------
    basis_name : str
        Name of the basis set. This is not case sensitive.
    elements : list
        List of element numbers that you want the basis set for. By default,
        all elements for which the basis set is defined are included.
    version : int
        Obtain a specific version of this basis set. By default,
        the latest version is returned.
    fmt: str
        The desired output format of the basis set references. By default,
        basis set information is returned as a list of dictionaries. Use
        get_reference_formats() to programmatically obtain the available formats.
        The `fmt` argument is not case sensitive.

        Available reference formats are

            * bib
            * txt
            * json

    data_dir : str
        Data directory with all the basis set information. By default,
        it is in the 'data' subdirectory of this project.

    Returns
    -------
    str or dict
        The references for the given basis set in the desired format. If `fmt` is **None**, this will be a python
        dictionary. Otherwise, it will be a string.
    '''

    data_dir = fix_data_dir(data_dir)
    basis_dict = get_basis(basis_name, elements=elements, version=version, data_dir=data_dir)

    all_ref_data = get_reference_data(data_dir)
    ref_data = references.compact_references(basis_dict, all_ref_data)

    if fmt is None:
        return ref_data

    return refconverters.convert_references(ref_data, fmt)


def get_basis_family(basis_name, data_dir=None):
    '''Lookup a family by a basis set name
    '''

    data_dir = fix_data_dir(data_dir)
    bs_data = _get_basis_metadata(basis_name, data_dir)
    return bs_data['family']


@memo.BSEMemoize
def get_families(data_dir=None):
    '''Return a list of all basis set families'''
    data_dir = fix_data_dir(data_dir)
    metadata = get_metadata(data_dir)

    families = set()
    for v in metadata.values():
        families.add(v['family'])

    return sorted(families)


def filter_basis_sets(substr=None, family=None, role=None, elements=None, data_dir=None):
    '''Filter basis sets by some criteria

    All parameters are ANDed together and are not case sensitive.

    Parameters
    ----------
    substr : str
        Substring to search for in the basis set name
    family : str
        Family the basis set belongs to
    role : str
        Role of the basis set
    elements : str or list
        List of elements that the basis set must include.
        Elements can be specified by Z-number (int or str) or by symbol (str).
        If this argument is a str (ie, '1-3,7-10'), it is expanded into a list.
        Z numbers and symbols (case insensitive) can be used interchangeably
        (see :func:`basis_set_exchange.misc.expand_elements`)
    data_dir : str
        Data directory with all the basis set information. By default,
        it is in the 'data' subdirectory of this project.

    Returns
    -------
    dict
        Basis set metadata that matches the search criteria
    '''

    data_dir = fix_data_dir(data_dir)
    metadata = get_metadata(data_dir)

    # family and role are required to be lowercase (via schema and validation functions)

    if family is not None:
        family = family.lower()
        if family not in get_families(data_dir):
            raise RuntimeError("Family '{}' is not a valid family".format(family))
        metadata = {k: v for k, v in metadata.items() if v['family'] == family}
    if role is not None:
        role = role.lower()
        if role not in get_roles():
            raise RuntimeError("Role '{}' is not a valid role".format(role))
        metadata = {k: v for k, v in metadata.items() if v['role'] == role}
    if elements is not None:
        elements = misc.expand_elements(elements, True)
        elements = set(elements)

        for basis_name, basis_data in metadata.items():
            ver_data = basis_data['versions']
            basis_data['versions'] = {k: v for k, v in ver_data.items() if elements <= set(v['elements'])}

        # There will be basis sets with no versions. So clean that up
        metadata = {k: v for k, v in metadata.items() if v['versions']}
    if substr:
        substr = substr.lower()
        metadata = {k: v for k, v in metadata.items() if substr in k or substr in v['display_name']}

    return metadata


@memo.BSEMemoize
def _family_notes_path(family, data_dir):
    '''Form a path to the notes for a family'''

    data_dir = fix_data_dir(data_dir)

    family = family.lower()
    if family not in get_families(data_dir):
        raise RuntimeError("Family '{}' does not exist".format(family))

    file_name = 'NOTES.' + family.lower()
    file_path = os.path.join(data_dir, file_name)
    return file_path


@memo.BSEMemoize
def _basis_notes_path(name, data_dir):
    '''Form a path to the notes for a basis set'''

    data_dir = fix_data_dir(data_dir)
    bs_data = _get_basis_metadata(name, data_dir)

    # the notes file is the same as the base file name, with a .notes extension
    filebase = bs_data['basename']
    file_path = os.path.join(data_dir, filebase + '.notes')
    return file_path


@memo.BSEMemoize
def get_family_notes(family, data_dir=None):
    '''Return a string representing the notes about a basis set family

    If the notes are not found, an empty string is returned
    '''

    file_path = _family_notes_path(family, data_dir)
    notes_str = fileio.read_notes_file(file_path)

    if notes_str is None:
        notes_str = ""

    ref_data = get_reference_data(data_dir)
    return notes.process_notes(notes_str, ref_data)


@memo.BSEMemoize
def has_family_notes(family, data_dir=None):
    '''Check if notes exist for a given family

    Returns True if they exist, false otherwise
    '''

    file_path = _family_notes_path(family, data_dir)
    return os.path.isfile(file_path)


@memo.BSEMemoize
def get_basis_notes(name, data_dir=None):
    '''Return a string representing the notes about a specific basis set

    If the notes are not found, an empty string is returned
    '''

    file_path = _basis_notes_path(name, data_dir)
    notes_str = fileio.read_notes_file(file_path)

    if notes_str is None:
        return ""

    ref_data = get_reference_data(data_dir)
    return notes.process_notes(notes_str, ref_data)


@memo.BSEMemoize
def has_basis_notes(family, data_dir=None):
    '''Check if notes exist for a given basis set

    Returns True if they exist, false otherwise
    '''

    file_path = _basis_notes_path(family, data_dir)
    return os.path.isfile(file_path)


def get_formats(function_types=None):
    '''Return information about the basis set formats available

    The returned data is a map of format to display name. The format
    can be passed as the fmt argument to :func:`get_basis()`

    If a list is specified for function_types, only those formats
    supporting the given function types will be returned.
    '''

    #####################################################################################
    # This function is being kept for two reasons. One, for backwards compatibility.
    # Two, the idea is that this module deals exclusively with obtaining data
    # from the internal basis set files, and we want to know what format we can
    # get them in. So it is semantically clear what get_formats means in that context.
    #####################################################################################
    return writers.get_writer_formats(function_types)


def get_reference_formats():
    '''Return information about the reference/citation formats available

    The returned data is a map of format to display name. The format
    can be passed as the fmt argument to :func:`get_references`
    '''
    return refconverters.get_reference_formats()


def get_roles():
    '''Return information about the available basis set roles available

    The returned data is a map of role to display name. The format
    can be passed as the role argument to fmt argument to :func:`lookup_basis_by_role`
    '''

    # yapf: disable
    return { 'orbital': 'Orbital basis',
             'jfit': 'J-fitting',
             'jkfit': 'JK-fitting',
             'rifit': 'RI-fitting',
             'optri': 'Optimized RI-fitting',
             'admmfit': 'Auxiliary-Density Matrix Method Fitting',
             'dftxfit': 'DFT Exchange Fitting',
             'dftjfit': 'DFT Correlation Fitting',
             'guess': 'Initial guess'
            }
    # yapf: enable


def get_data_dir():
    '''Get the default data directory of this installation'''
    return _default_data_dir
