import re
from .. import lut, manip
from . import helpers

element_entry_re = re.compile(r'^([a-zA-Z]{1,3}):(.*):(.*)$')
shell_info_re = re.compile(r'^([a-zA-Z])\s+(\d+)\s+(\d+)$')


def read_gbasis(basis_lines):
    '''Reads gbasis-formatted file data and converts it to a dictionary with the
       usual BSE fields

       GBASIS only supports electronic shells (no ecp)

       Note that the gbasis format does not store all the fields we
       have, so some fields are left blank
    '''

    basis_lines = helpers.prune_lines(basis_lines, '!#')

    bs_data = {}

    # The file just contains sections separated headed by
    # a line looking like "Al:aug-cc-pV5+dZ:(21s13p6d4f3g2h) -> [8s7p6d4f3g2h]"
    element_sections = helpers.partition_lines(basis_lines, element_entry_re.match, min_size=4)

    found_basis = set()
    for element_lines in element_sections:
        # First line is element + basis name (already checked)
        # Second is highest angular momentum
        element_sym, basis_name, _ = helpers.parse_line_regex(element_entry_re, element_lines[0],
                                                              "Element entry: sym:basis:pattern")

        element_Z = lut.element_Z_from_sym(element_sym, as_str=True)

        # Add what basis we found
        # We only support one basis per file
        found_basis.add(basis_name)

        element_data = manip.create_element_data(bs_data, element_Z, 'electron_shells')
        max_am = helpers.parse_line_regex(r'^(\d+)$', element_lines[1], 'Highest AM')

        # Split all the shells based on lines beginning with alpha character
        shell_blocks = helpers.partition_lines(element_lines[2:], lambda x: x[0].isalpha())

        # We know how many blocks there should be
        if (max_am + 1) != len(shell_blocks):
            raise RuntimeError("Different number of blocks for element {}. Expected {}, found {}".format(
                element_sym, max_am + 1, len(shell_blocks)))

        # Now loop over the blocks
        found_am = []

        for shell_lines in shell_blocks:
            shell_am, nprim, ngen = helpers.parse_line_regex(shell_info_re, shell_lines[0], 'Shell: AM, nprim, ngen')
            shell_am = lut.amchar_to_int(shell_am)

            func_type = lut.function_type_from_am(shell_am, 'gto', 'spherical')

            if len(shell_am) > 1:
                raise RuntimeError("Fused AM not supported by gbasis")

            if shell_am in found_am:
                raise RuntimeError("Duplicate AM for element {}: AM {} already found".format(element_sym, shell_am[0]))

            found_am.append(shell_am)

            if shell_am[0] > max_am:
                raise RuntimeError("Found AM greater than max AM: {}".format(shell_am[0]))

            # Now parse the exponents & coefficients
            exponents, coefficients = helpers.parse_primitive_matrix(shell_lines[1:], nprim=nprim, ngen=ngen)

            shell = {
                'function_type': func_type,
                'region': '',
                'angular_momentum': shell_am,
                'exponents': exponents,
                'coefficients': coefficients
            }

            element_data['electron_shells'].append(shell)

    if len(found_basis) > 1:
        raise RuntimeError("Multiple basis sets in a single file: " + str(found_basis))

    return bs_data
