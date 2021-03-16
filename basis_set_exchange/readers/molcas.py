import re
from .. import lut, manip, misc
from . import helpers

element_head_re = re.compile(r'^/([a-zA-Z]{1,3})\.(?:ECP\.)?([^.]+)\..*$')
electron_head_re = re.compile(r'^/([a-zA-Z]{1,3})\.([^.]+)\..*$')
ecp_head_re = re.compile(r'^/([a-zA-Z]{1,3})\.ECP\.([^.]+)\..*$')

electron_z_max_am_re = re.compile(r'^(\d+|{})(?:\s+(\d+))?$'.format(helpers.floating_re_str))
shell_nprim_ngen_re = re.compile(r'^(\d+)(?:\s+(\d+))?$')

ecp_info_re = re.compile(r'^[Pp]{2}\s*,\s*([a-zA-Z]+)\s*,\s*(\d+)\s*,\s*(\d+)\s*;$')
ecp_pot_begin_re = re.compile(r'^(\d+)\s*;.*$')  # Sometime comments are after the semicolon


def _parse_electron_lines(basis_lines, bs_data, element_Z):
    element_data = manip.create_element_data(bs_data, element_Z, 'electron_shells')

    # Handle the options block
    # This specifies the kind of data that might be found at the end of the element block
    # We ignore that data, but we need to know it is there
    options_lines, basis_lines = helpers.remove_block(basis_lines, 'Options', 'EndOptions')

    # Every option adds another block after each shell.
    # I don't care about the contents, just how many options there are
    n_option_blocks = len(options_lines)

    # Next is the nuclear charge and the max_am
    # max_am may be missing
    nuc_charge, max_am = helpers.parse_line_regex(electron_z_max_am_re, basis_lines[0], 'Electron: Z, max_am')

    # If the nuclear charge is not equal to the element Z, then this must be an ECP
    # If the number of ECP electrons already exists, check it.
    # If not, put it into the element data, and the ecp reader below will double check it
    nuc_charge = float(nuc_charge)

    # Is this actually an integer?
    if float(int(nuc_charge)) != nuc_charge:
        raise RuntimeError("Non-integer specified for nuclear charge: " + nuc_charge)

    ecp_electrons = int(element_Z) - int(nuc_charge)
    if 'ecp_electrons' in element_data and element_data['ecp_electrons'] != ecp_electrons:
        raise RuntimeError("Element Z: {} with charge {} does not match ECP electrons {}".format(
            element_Z, nuc_charge, ecp_electrons))
    elif ecp_electrons > 0:
        element_data['ecp_electrons'] = ecp_electrons

    # Partition the remaining lines. Blocks start with one or two integers.
    # This includes typical electron shell blocks, and blocks corresponding to options
    shell_blocks = helpers.partition_lines(basis_lines[1:], lambda x: x.split()[0].isdecimal())

    # There should be max_am + 1 shell blocks
    # Each shell will also have n_option_blocks additional blocks
    if max_am is not None and len(shell_blocks) != (max_am + 1) * (n_option_blocks + 1):
        raise RuntimeError("Found {} shell blocks. Expected {}".format(len(shell_blocks),
                                                                       (max_am + 1) * (n_option_blocks + 1)))

    # Shells are simply increasing AM
    shell_am = 0

    # Loop over all shell lines before the option block data
    # Option blocks are after each shell, so just slice and skip
    # every n_option_blocks
    if n_option_blocks > 0:
        shell_blocks = shell_blocks[::n_option_blocks + 1]

    for shell_lines in shell_blocks:
        nprim, ngen = helpers.parse_line_regex(shell_nprim_ngen_re, shell_lines[0], 'Shell nprim, ngen')

        if nprim <= 0:
            raise RuntimeError("Cannot have {} primitives in a shell".format(nprim))
        if ngen is not None and ngen <= 0:
            raise RuntimeError("Cannot have {} general contractions in a shell".format(ngen))

        exponents, shell_lines = helpers.read_n_floats(shell_lines[1:], nprim)

        # Read the coefficient matrix from the next block
        # We may or may not know the dimensions of the matrix since &@*!# is optional
        # in this format.
        coefficients = helpers.read_all_floats(shell_lines)

        # Now we can check if this makes sense
        n_coefs = len(coefficients)
        if n_coefs == 0:
            raise RuntimeError("Have zero coefficients?")
        if n_coefs % nprim != 0:
            raise RuntimeError("Number of coefficients is not a multiple of nprim: {} % {} = {}".format(
                n_coefs, nprim, n_coefs % nprim))

        # If we do actually have the number of general contractions, does it match?
        if ngen is not None and ngen != n_coefs // nprim:
            raise RuntimeError("Expected {} general contractions, but found {}".format(ngen, n_coefs // nprim))
        else:
            ngen = n_coefs // nprim

        # Turn the coefficients into a matrix
        coefficients = helpers.chunk_list(coefficients, nprim, ngen)
        coefficients = misc.transpose_matrix(coefficients)

        # Now add to the bs_data
        func_type = lut.function_type_from_am([shell_am], 'gto', 'spherical')

        shell = {
            'function_type': func_type,
            'region': '',
            'angular_momentum': [shell_am],
            'exponents': exponents,
            'coefficients': coefficients
        }

        element_data['electron_shells'].append(shell)
        shell_am += 1


def _parse_ecp_lines(basis_lines, bs_data, element_Z):
    # Remove "Spectral" Stuff
    _, basis_lines = helpers.remove_block(basis_lines, r'^Spectral.*', r'^End\s*Of\s*Spectral.*')

    # Parse the ecp info line
    element_sym, ecp_electrons, max_am = helpers.parse_line_regex(ecp_info_re, basis_lines[0],
                                                                  "ECP Info: pp,sym,nelec,max_am")
    element_Z_ecp = lut.element_Z_from_sym(element_sym, as_str=True)

    # Does this block match the element symbol from the main element header?
    if element_Z_ecp != element_Z:
        raise RuntimeError("ECP element Z={} found in block for element Z={}".format(element_Z, element_Z_ecp))

    element_data = manip.create_element_data(bs_data, element_Z, 'ecp_potentials')

    # Does the ecp_electrons key exist? This may have been determined when reading the
    # electron shells above
    if 'ecp_electrons' in element_data and element_data['ecp_electrons'] != ecp_electrons:
        raise RuntimeError(
            "No. of electrons specified in ECP block do not match already-determined number of electrons: {} vs {}".
            format(ecp_electrons, element_data['ecp_electrons']))
    else:
        element_data['ecp_electrons'] = ecp_electrons

    # Now split into potentials
    # The beginning of each potential is a number followed by a semicolon
    pot_blocks = helpers.partition_lines(basis_lines[1:], ecp_pot_begin_re.match, min_size=2)

    if len(pot_blocks) != max_am + 1:
        raise RuntimeError("Expected {} potentials, but got {}".format(max_am + 1, len(pot_blocks)))

    # Set up the AM for the potentials
    all_pot_am = helpers.potential_am_list(max_am)

    for pot_lines in pot_blocks:
        pot_am = all_pot_am.pop(0)

        nlines = helpers.parse_line_regex(ecp_pot_begin_re, pot_lines[0], "ECP Potential: # of lines")
        if nlines != len(pot_lines) - 1:
            raise RuntimeError("Expected {} lines in potential, but got {}".format(nlines, len(pot_lines) - 1))

        # Strip trailing semicolon
        pot_lines = [x.rstrip(';') for x in pot_lines[1:]]
        ecp_data = helpers.parse_ecp_table(pot_lines, split=r'\s*,\s*')
        ecp_pot = {
            'angular_momentum': [pot_am],
            'ecp_type': 'scalar_ecp',
            'r_exponents': ecp_data['r_exp'],
            'gaussian_exponents': ecp_data['g_exp'],
            'coefficients': ecp_data['coeff']
        }

        element_data['ecp_potentials'].append(ecp_pot)


def read_molcas(basis_lines):
    '''Reads molcas-formatted file data and converts it to a dictionary with the
       usual BSE fields

       Note that the molcas format does not store all the fields we
       have, so some fields are left blank
    '''

    basis_lines = helpers.prune_lines(basis_lines, '*#$')

    bs_data = {}

    # Split into elements. Every start of an element is /
    element_blocks = helpers.partition_lines(basis_lines, lambda x: x.startswith('/'), min_size=4)

    # Inside this loop, check that all blocks refer to the same basis set
    basis_names_found = set()

    for element_lines in element_blocks:
        # The start of the element block is:
        #     /{element_sym}.{bs_name}.stuff
        #     Reference line 1
        #     Reference line 2
        # The next line should either be "options" or the element Z number
        element_sym, basis_name = helpers.parse_line_regex(element_head_re, element_lines[0], 'Start of element line')
        element_Z = lut.element_Z_from_sym(element_sym, as_str=True)
        basis_names_found.add(basis_name.lower())

        # remove the header and comments
        element_lines = element_lines[3:]

        # Split based on PP (pseudopotential)
        element_split = helpers.partition_lines(element_lines,
                                                lambda x: x.lower().startswith('pp,'),
                                                min_blocks=1,
                                                max_blocks=2)

        for block_lines in element_split:
            if block_lines[0].lower().startswith('pp'):
                _parse_ecp_lines(block_lines, bs_data, element_Z)
            else:
                _parse_electron_lines(block_lines, bs_data, element_Z)

    # Check for multiple basis sets
    if len(basis_names_found) > 1:
        raise RuntimeError("Multiple basis sets found in file: " + ','.join(basis_names_found))

    return bs_data
