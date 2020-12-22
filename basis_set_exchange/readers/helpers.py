'''
Some helper functions for parsing basis set files
'''

import re
import regex
from ..misc import transpose_matrix

floating_re_str = r'[-+]?\d*\.\d*(?:[dDeE][-+]?\d+)?'
floating_re = re.compile(floating_re_str)
floating_only_re = re.compile('^' + floating_re_str + '$')
integer_re_str = r'[-+]?\d+'
integer_re = re.compile(integer_re_str)
integer_only_re = re.compile('^' + integer_re_str + '$')


def _convert_str_int(s):
    '''Optionally convert a string to an integer

    If string s represents an integer, returns an int. Otherwise, returns s unchanged
    '''

    # May throw ValueError if string is not an int
    # May throw TypeError if it is NoneType (some captures may be optional)
    try:
        return int(s)
    except ValueError:
        return s
    except TypeError:
        return s


def is_floating(s):
    '''Tests if a string is a floating point number'''
    return floating_only_re.match(s)


def is_integer(s):
    '''Tests if a string is an integer'''
    return integer_only_re.match(s)


def replace_d(s):
    '''Replaces fortran-style 'D' with 'E'''
    s = s.replace('D', 'E')
    s = s.replace('d', 'e')
    return s


def function_type_from_am(shell_am, base_type, spherical_type):
    if max(shell_am) <= 1:
        return base_type
    else:
        return base_type + '_' + spherical_type


def potential_am_list(max_am):
    '''Creates a canonical list of AM for use with ECP potentials

    The list is [max_am, 0, 1, ..., max_am-1]
    '''

    am_list = list(range(max_am + 1))
    am_list.insert(0, am_list.pop())
    return am_list


def chunk_list(lst, rows, cols):
    '''Turns a list into a matrix of the given dimensions'''

    n_elements = len(lst)
    if n_elements != rows * cols:
        raise RuntimeError("Cannot partition {} elements into a {}x{} matrix".format(n_elements, rows, cols))

    mat = [lst[i:i + cols] for i in range(0, n_elements, cols)]
    assert len(mat) == rows
    return mat


def remove_expected_line(lines, expected='', position=0):
    '''Tests the first element of the list to see if it is an expected string, and removes it

    If line does not match, or lines is empty, an exception is raised
    '''
    if not lines:
        raise RuntimeError("No lines to test for expected line")
    if position >= 0 and len(lines) <= position:
        raise RuntimeError("Not enough lines. Can't test line {} when there are {} lines".format(position, len(lines)))
    if position < 0 and len(lines) < -position:
        raise RuntimeError("Not enough lines. Can't test line {} when there are {} lines".format(position, len(lines)))
    if lines[position] != expected:
        raise RuntimeError("Expected line '{}', but found '{}'".format(expected, lines[0]))

    new_lines = lines[:]
    new_lines.pop(position)
    return new_lines


def create_element_data(bs_data, element_Z, key, key_exist_ok=False, element_exist_ok=True, create=list):
    '''Creates an element and a subkey of the element in bs_data

    Note that bs_data is modified!
    '''

    if element_Z not in bs_data:
        bs_data[element_Z] = {}
    elif not element_exist_ok:
        raise RuntimeError("Element {} already exists in basis data".format(element_Z))

    if key not in bs_data[element_Z]:
        bs_data[element_Z][key] = create()
    elif not key_exist_ok:
        raise RuntimeError("Key {} already exists in basis data for element {}".format(key, element_Z))

    return bs_data[element_Z]


def parse_line_regex(rex, line, description=None, convert_int=True):
    if isinstance(rex, str):
        rex = re.compile(rex)

    r = rex.match(line)
    if not r:
        if description:
            raise RuntimeError("Regex '{}' does not match line: '{}'. Regex is '{}'".format(
                description, line, rex.pattern))
        else:
            raise RuntimeError("Regex '{}' does not match line: '{}'".format(rex.pattern, line))

    g = r.groups()
    if convert_int:
        g = [_convert_str_int(x) for x in g]

    if len(g) == 1:
        return g[0]
    else:
        return g


def parse_line_regex_dict(rex, line, description=None, convert_int=True):
    if isinstance(rex, str):
        rex = regex.compile(rex)

    r = rex.match(line)
    if not r:
        if description:
            raise RuntimeError("Regex '{}' does not match line: '{}'. Regex is '{}'".format(
                description, line, rex.pattern))
        else:
            raise RuntimeError("Regex '{}' does not match line: '{}'".format(rex.pattern, line))

    g = r.capturesdict()
    if convert_int:
        g = {k: [_convert_str_int(x) for x in g[k]] for k in g}

    return g


def partition_lines(lines,
                    condition,
                    before=0,
                    min_after=None,
                    min_blocks=None,
                    max_blocks=None,
                    min_size=1,
                    include_match=True):
    '''Partition a list of lines based on some condition

    Parameters
    ----------
    lines : list
        List of strings representing the lines in the file
    condition : function
        Function or lambda that takes a line as an argument and returns True if
        that line is the start of a section
    before : int
        Number of lines prior to the splitting line (where condition(line) == True)
        to include
`   min_after : int
        Minimum number of lines to include after the match. This number of lines is forced
        even if a match is found again with that set of lines.
    min_blocks : int
        Minimum number of blocks to find. If fewer are found, an exception is thrown
    max_blocks : int
        Maximum number of blocks to find. If more are found, an exception is thrown
    min_size : int
        Minimum size/length of each block. If one is found that is smaller, an exception is thrown

    Returns
    -------
    list of list
        The original list of strings partitioned into blocks
    '''

    # First, just partition into blocks
    all_blocks = []
    cur_block = []

    i = 0
    while i < len(lines):
        line = lines[i]
        if condition(line):
            # This is beginning of a new block. If we have anything,
            # append it to block list
            if cur_block:
                all_blocks.append(cur_block)
                cur_block = []

            if include_match:
                cur_block.append(line)
            if min_after:
                cur_block.extend(lines[i + 1:i + 1 + min_after])
                i += min_after
        else:
            cur_block.append(line)

        i += 1

    # Handle the last block, if needed
    if cur_block:
        all_blocks.append(cur_block)

    # Now handle if the caller specified the 'before' parameter
    if before > 0:
        # Each block after the first steals from the previous
        # If we are going to do this, then:
        #    1.) There must be more than one block
        #    2.) The first block must be only the part to be moved to the second block.
        if len(all_blocks) <= 1:
            raise RuntimeError("Cannot partition lines with before = {}: have {} blocks".format(
                before, len(all_blocks)))

        if len(all_blocks[0]) != before:
            raise RuntimeError("Cannot partition lines with before = {}: first block has {} lines".format(
                before, len(all_blocks[0])))

        # Starting with the second block (index 1), move 'before' lines from the end of the previous block
        # to the beginning of this block
        for idx in range(1, len(all_blocks)):
            all_blocks[idx] = all_blocks[idx - 1][-before:] + all_blocks[idx]
            all_blocks[idx - 1] = all_blocks[idx - 1][:-before]

        # Double check that the first block is now empty (we moved it all)
        first_block = all_blocks.pop(0)
        assert len(first_block) == 0

    # Check to make sure all blocks for min_size
    if min_size > 0:
        for idx, block in enumerate(all_blocks):
            if len(block) < min_size:
                raise RuntimeError("Block {} does not have minimum number of lines ({})".format(idx, min_size))

    if min_blocks and len(all_blocks) < min_blocks:
        raise RuntimeError("Found {} blocks, but need at least {}".format(len(all_blocks), min_blocks))
    if max_blocks and len(all_blocks) > max_blocks:
        raise RuntimeError("Found {} blocks, but need at most {}".format(len(all_blocks), max_blocks))

    return all_blocks


def read_n_floats(lines, n_numbers, convert=False, split=r'\s+'):
    '''Reads in a number of space-separated floating-point numbers

    These numbers may span multiple lines.

    An exception will be thrown if there are remaining numbers on the last line.

    If a non-floating point entry is found, an exception is also thrown.

    If convert is True, then float objects are created from the strings

    Returns the found floating point numbers (as str), and the remaining lines
    '''

    found_numbers = []
    while len(found_numbers) < n_numbers:
        if not lines:
            raise RuntimeError("Wanted {} numbers but ran out of lines after {}".format(n_numbers, len(found_numbers)))
        if not lines[0]:
            raise RuntimeError("Wanted {} numbers but found empty line after {}".format(n_numbers, len(found_numbers)))

        l = replace_d(lines[0])
        s = re.split(split, l.strip())
        found_numbers.extend(s)
        lines = lines[1:]

    if len(found_numbers) != n_numbers:
        raise RuntimeError("Wanted {} numbers, but found {}".format(n_numbers, found_numbers))

    # Make sure these are floating point
    if not all(is_floating(x) for x in found_numbers):
        raise RuntimeError("Non-floating-point value found in numbers: " + ' '.join(found_numbers))

    if convert:
        found_numbers = [float(x) for x in found_numbers]

    return found_numbers, lines


def read_all_floats(lines, convert=False, split=r'\s+'):
    '''Reads in all floats on all lines

    This function takes a block of numbers and splits them all, for all lines in the block.

    If a non-floating point entry is found, an exception is also thrown.

    If convert is True, then float objects are created from the strings
    '''

    found_numbers = []
    for l in lines:
        l = replace_d(l)
        s = re.split(split, l.strip())
        found_numbers.extend(s)

    # Make sure these are floating point
    if not all(is_floating(x) for x in found_numbers):
        raise RuntimeError("Non-floating-point value found in numbers: " + ' '.join(found_numbers))

    if convert:
        found_numbers = [float(x) for x in found_numbers]

    return found_numbers


def read_n_integers(lines, n_ints, convert=False, split=r'\s+'):
    '''Reads in a number of space-separated integers

    These numbers may span multiple lines.

    An exception will be thrown if there are remaining numbers on the last line.

    If a non-integer entry is found, an exception is also thrown.

    If convert is True, then int objects are created from the strings

    Returns the found integers point numbers (as str), and the remaining lines
    '''

    found_numbers = []
    while len(found_numbers) < n_ints:
        s = re.split(split, lines[0].strip())
        found_numbers.extend(s)
        lines = lines[1:]

    if len(found_numbers) != n_ints:
        raise RuntimeError("Wanted {} numbers, but found {}".format(n_ints, found_numbers))

    # Make sure these are integers
    if not all(is_integer(x) for x in found_numbers):
        raise RuntimeError("Non-integer value found in numbers: " + ' '.join(found_numbers))

    if convert:
        found_numbers = [int(x) for x in found_numbers]

    return found_numbers, lines


def parse_fixed_matrix(lines, rows, cols, split=r'\s+'):
    '''Parses a simple matrix of numbers with a predefined number of rows/columns

    This will read in a matrix of the given number of rows and columns, even if the
    rows span multiple lines. There must be a newline at the very end of a row.

    Returns the matrix and the remaining lines
    '''
    mat = []

    for i in range(rows):
        row_data, lines = read_n_floats(lines, cols)
        mat.append(row_data)

    return mat, lines


def parse_matrix(lines, rows=None, cols=None, split=r'\s+'):
    '''Parses a simple matrix of numbers

    The lines parameter must specify a list of strings containing the entire matrix.

    If rows and/or cols is specified, and the found number of rows/cols does not
    match, an exception is raised.
    '''
    mat = []

    for l in lines:
        l = replace_d(l)
        s = re.split(split, l.strip())

        if not all(is_floating(x) for x in s):
            raise RuntimeError("Non-floating-point value found in matrix: " + ' '.join(s))

        mat.append(s)

    # Make sure all rows have the same number of columns
    for x in mat:
        if len(x) == 0:
            raise RuntimeError("Matrix row has zero values")
        if len(x) != len(mat[0]):
            raise RuntimeError("Inconsistent number of columns: {} vs {}".format(len(x), len(mat[0])))

    if len(mat) == 0:
        raise RuntimeError("Empty matrix?")

    # if rows/cols are given, sanity check
    if rows is not None:
        rows = int(rows)

        if len(mat) != rows:
            raise RuntimeError("Inconsistent number of rows: {} vs {}".format(rows, len(mat)))

    # We already checked that all rows have the same number of columns
    # so just check the first row
    if cols is not None and len(mat[0]) != cols:
        raise RuntimeError("Inconsistent number of columns: {} vs {}".format(cols, len(mat[0])))

    return mat


def parse_primitive_matrix(lines, nprim=None, ngen=None, split=r'\s+'):
    '''Parses a matrix/table of exponents and coefficients

    The first column of the matrix contains exponents, and the remaining
    columns contain the coefficients for all general contractions.

    The lines parameter must specify a list of strings containing the entire matrix.

    If nprim and/or ngen are specified, and the found number of primitives/contractions 
    match, an exception is raised.
    '''
    exponents = []
    coefficients = []

    for l in lines:
        l = replace_d(l)
        s = re.split(split, l.strip())

        e = s[0]
        c = s[1:]

        # Make sure all exponents and coefficients are floating point
        if not is_floating(e):
            raise RuntimeError("Non-floating-point value found in exponents: " + e)

        if not all(is_floating(x) for x in c):
            raise RuntimeError("Non-floating-point value found in coefficients: " + ' '.join(c))

        exponents.append(e)
        coefficients.append(c)

    # Make sure all primitives have the same number of coefficients
    for c in coefficients:
        if len(c) == 0:
            raise RuntimeError("Missing contraction coefficients")
        if len(c) != len(coefficients[0]):
            raise RuntimeError("Inconsistent number of coefficients: {} vs {}".format(len(c), len(coefficients[0])))

    coefficients = transpose_matrix(coefficients)

    if len(exponents) == 0:
        raise RuntimeError("No exponents found")
    if len(coefficients) == 0:
        raise RuntimeError("No coefficients found")

    # if nprim and/or ngen are given, sanity check
    if nprim is not None:
        nprim = int(nprim)

        if len(exponents) != nprim:
            raise RuntimeError("Inconsistent number of primitives in exponents: {} vs {}".format(
                nprim, len(exponents)))

        # We already checked that all coefficients have the same number of primitives
        # so just check the first general contraction
        if len(coefficients[0]) != nprim:
            raise RuntimeError("Inconsistent number of primitives in coefficients: {} vs {}".format(
                nprim, len(coefficients[0])))

    if ngen is not None:
        ngen = int(ngen)
        if len(coefficients) != ngen:
            raise RuntimeError("Inconsistent number of general contractions: {} vs {}".format(ngen, len(coefficients)))

    return exponents, coefficients


def parse_ecp_table(lines, order=['r_exp', 'g_exp', 'coeff'], split=r'\s+'):
    ecp_data = {'r_exp': [], 'g_exp': [], 'coeff': []}

    for l in lines:
        l = replace_d(l)
        s = re.split(split, l.strip())

        if len(s) != 3:
            raise RuntimeError("Expected 3 values in ecp table")

        for i, k in enumerate(order):
            ecp_data[k].append(s[i])

    # Note: This function does not handle multiple coefficients
    # So we only have to add another layer to the coefficient list
    ecp_data['coeff'] = [ecp_data['coeff']]

    # Change r exponents to integers
    # But do a check first
    if not all(is_integer(x) for x in ecp_data['r_exp']):
        raise RuntimeError("Non-integer value found in r exponents")

    # Make sure g exponents and coefficients are floating point
    if not all(is_floating(x) for x in ecp_data['g_exp']):
        raise RuntimeError("Non-floating-point value found in g exponents")

    for c in ecp_data['coeff']:
        if not all(is_floating(x) for x in c):
            raise RuntimeError("Non-floating-point value found in coefficients")

    # Now convert r exponents to integers
    ecp_data['r_exp'] = [int(x) for x in ecp_data['r_exp']]

    return ecp_data


def prune_lines(lines, skipchars='', prune_blank=True, strip_end_blanks=True):
    '''Remove comment and blank lines

    Also strips all lines of beginning/ending whitespace.

    Parameters
    ----------
    lines : list of str
        List of lines to prune
    skipchars : str
        Comment characters designating lines to be skipped. Lines starting
        with any character in this string will be removed.
    prune_blank: bool
        Remove blank lines
    strip_end_blank_lines: bool
        Remove starting/ending blank lines (even if prune_blank=False)

    Returns
    -------
    list of str
        Pruned lines
    '''

    lines = [l.strip() for l in lines]

    if skipchars:
        lines = [l for l in lines if len(l) == 0 or l[0] not in skipchars]
    if prune_blank:
        lines = [l for l in lines if l]

    if not lines:
        return lines

    # No need to strip end blank lines if we got rid of all blank lines
    # Note that we can use pop. The lines variable has already been copied
    if strip_end_blanks and not prune_blank:
        while lines and len(lines[0]) == 0:
            lines.pop(0)
        while lines and len(lines[-1]) == 0:
            lines.pop()

    return lines


def remove_block(lines, start_re, end_re):
    '''Removes a block of data from the lines of text
    
       For example, there may be an optional block of options (like in molcas)

       This will only remove a single block

       Parameters
       ----------
       lines : list of str
           Line of text to parse
       start_re : str
           Regex string representing the start of the block (case insensitive)
       end_re : str
           Regex string representing the end of the block (case insensitive)

       Returns
       -------
       list of str, list of str
           The block found (may be empty), the input lines without the block
    '''

    start = re.compile(start_re, flags=re.IGNORECASE)
    end = re.compile(end_re, flags=re.IGNORECASE)

    start_idx = -1
    for idx, l in enumerate(lines):
        if start.match(l):
            if start_idx != -1:
                raise RuntimeError("Multiple blocks starting with '{}' found".format(start))
            start_idx = idx
            break
    else:
        return [], lines  # Block not found

    i = start_idx + 1
    block_lines = []
    while i < len(lines) and not end.match(lines[i]):
        block_lines.append(lines[i])
        i += 1

    # Did we find the end of the block
    if i == len(lines):
        raise RuntimeError("Cannot find end of block. Looking for '{}' to close '{}'".format(end_re, start_re))

    lines = lines[:start_idx] + lines[i + 1:]
    return block_lines, lines
