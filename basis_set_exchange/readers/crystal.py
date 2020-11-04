import re
from .. import lut
from . import helpers

# Element entry starts out with two integers
element_re = re.compile(r'^([1-9]{1,3})\s+([1-9]+)?$')
# The shell definition has three integers and two (potentially floating point, maybe integer) numbers
shell_re = re.compile(r'^([1-9]+)\s+([1-9]+)\s+([1-9]+)\s+({0})\s+({0})$'.format(helpers.floating_re_str))

def _parse_electron_lines(basis_lines, bs_data):
    '''Parses lines representing all the electron shells for a single element

    Resulting information is stored in bs_data
    '''

    # First line is "{element} {number of shells}"
    element_Z = int(basis_lines[0].split()[0])
    element_sym = lut.element_name_from_Z(element_Z)
    num_shells = int(basis_lines[0].split()[1])

    element_data = helpers.create_element_data(bs_data, element_Z, 'electron_shells')

    # After that come the shells.
    shell_blocks = helpers.partition_lines(basis_lines[1:], shell_re.match)
    for sh_lines in shell_blocks:
        # Shell starts with five integers
        ityb, raw_shell_am, ngen, formal_charge, scale = helpers.parse_line_regex(shell_re, sh_lines[0], "ityb, lat, ng, che, scal")
        assert(ityb == 0) # other choices 1 or 2 are Pople STO-nG and 3(6)-21G

        # Parse angular momentum
        if raw_shell_am == 0:
            shell_am = [0]
        elif raw_shell_am == 1:
            # SP shell
            shell_am = [0, 1]
        else:
            # offset by 1
            shell_am = [raw_shell_am-1]

        # CRYSTAL uses a spherical basis
        func_type = helpers.function_type_from_am(shell_am, 'gto', 'spherical')

        # Handle gaussian scaling factors
        # The square of the scaling factor is applied to exponents.
        # Typically they are 1.0, but not always
        scaling_factors = helpers.replace_d(scaling_factors)
        scaling_factors = [float(x) for x in scaling_factors.split()]

        # Remove any scaling factors that are 0.0
        scaling_factors = [x for x in scaling_factors if x != 0.0]

        # We should always have at least one scaling factor
        if len(scaling_factors) == 0:
            raise RuntimeError("No scaling factors given for element {}: Line: {}".format(element_sym, sh_lines[0]))

        # There can be multiple scaling factors, but we don't handle that. It seems to be very rare
        if len(scaling_factors) > 1:
            raise NotImplementedError("Number of scaling factors > 1")

        scaling_factor = float(scaling_factors[0])**2
        has_scaling = scaling_factor != 1.0

        # How many columns of coefficients do we have?
        # Gaussian doesn't support general contractions, so only >1 if
        # you have a fused shell
        ngen = len(shell_am)

        # Now read the exponents and coefficients
        exponents, coefficients = helpers.parse_primitive_matrix(sh_lines[1:], nprim, ngen)

        # If there is a scaling factor, apply it
        # But we keep track of some significant-figure type stuff (as best we can)
        if has_scaling:
            new_exponents = []
            for ex in exponents:
                ex = float(ex) * scaling_factor
                ex = '{:.16E}'.format(ex)

                # Trim useless zeroes
                ex_splt = ex.split('E')
                ex = ex_splt[0].rstrip('0')
                if ex[-1] == '.':  # Stripped all the zeroes...
                    ex += '0'
                ex += 'E' + ex_splt[1]
                new_exponents.append(ex)

            exponents = new_exponents

        shell = {
            'function_type': func_type,
            'region': '',
            'angular_momentum': shell_am,
            'exponents': exponents,
            'coefficients': coefficients
        }

        element_data['electron_shells'].append(shell)


def read_crystal(basis_lines):
    '''Reads Crystal-formatted file data and converts it to a dictionary
       with the usual BSE fields

       Note that the Crystal format does not store all the fields we
       have, so some fields are left blank
    '''

    # Removes comments
    basis_lines = helpers.prune_lines(basis_lines, '!')

    bs_data = {}

    # Empty file?
    if not basis_lines:
        return bs_data

    # split into element sections (may be electronic or ecp)
    element_sections = helpers.partition_lines(basis_lines, element_re.match, min_size=3)
    for es in element_sections:
        _parse_electron_lines(es, bs_data)

    return bs_data
