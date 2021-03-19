"""
Common basis set manipulations

This module contains functions for uncontracting and merging basis set
data, as well as some other small functions.
"""

import copy
from . import skel, misc


def merge_element_data(dest, sources, use_copy=True):
    """
    Merges the basis set data for an element from multiple sources
    into dest.

    The destination is not modified, and a (shallow) copy of dest is returned
    with the data from sources added.

    If use_copy is True, then the data merged into dest will be a (deep)
    copy of that found in sources. Otherwise, data may be shared between dest
    and sources
    """

    if dest is not None:
        ret = dest.copy()
    else:
        ret = {}

    if use_copy:
        sources = copy.deepcopy(sources)

    # Note that we are not copying notes/data_sources
    for s in sources:
        if 'electron_shells' in s:
            if 'electron_shells' not in ret:
                ret['electron_shells'] = []
            ret['electron_shells'].extend(s['electron_shells'])
        if 'ecp_potentials' in s:
            if 'ecp_potentials' in ret:
                raise RuntimeError('Cannot overwrite existing ECP')
            ret['ecp_potentials'] = s['ecp_potentials']
            ret['ecp_electrons'] = s['ecp_electrons']
        if 'references' in s:
            if 'references' not in ret:
                ret['references'] = []
            for ref in s['references']:
                ret['references'].append(ref)

    return ret


def prune_shell(shell, use_copy=True):
    """
    Removes exact duplicates of primitives, and condenses duplicate exponents
    into general contractions

    Also removes primitives if all coefficients are zero
    """

    new_exponents = []
    new_coefficients = []

    exponents = shell['exponents']
    nprim = len(exponents)

    # transpose of the coefficient matrix
    coeff_t = list(map(list, zip(*shell['coefficients'])))

    # Group by exponents
    ex_groups = []
    for i in range(nprim):
        for ex in ex_groups:
            if float(exponents[i]) == float(ex[0]):
                ex[1].append(coeff_t[i])
                break
        else:
            ex_groups.append((exponents[i], [coeff_t[i]]))

    # Now collapse within groups
    for ex in ex_groups:
        if len(ex[1]) == 1:
            # only add if there is a nonzero contraction coefficient
            if not all([float(x) == 0.0 for x in ex[1][0]]):
                new_exponents.append(ex[0])
                new_coefficients.append(ex[1][0])
            continue

        # ex[1] contains rows of coefficients. The length of ex[1]
        # is the number of times the exponent is duplicated. Columns represent general contractions.
        # We want to find the non-zero coefficient in each column, if it exists
        # The result is a single row with a length representing the number
        # of general contractions

        new_coeff_row = []

        # so take yet another transpose.
        ex_coeff = list(map(list, zip(*ex[1])))
        for g in ex_coeff:
            nonzero = [x for x in g if float(x) != 0.0]
            if len(nonzero) > 1:
                raise RuntimeError("Exponent {} is duplicated within a contraction".format(ex[0]))

            if not nonzero:
                new_coeff_row.append(g[0])
            else:
                new_coeff_row.append(nonzero[0])

        # only add if there is a nonzero contraction coefficient anywhere for this exponent
        if not all([float(x) == 0.0 for x in new_coeff_row]):
            new_exponents.append(ex[0])
            new_coefficients.append(new_coeff_row)

    # take the transpose again, putting the general contraction
    # as the slowest index
    new_coefficients = list(map(list, zip(*new_coefficients)))

    shell['exponents'] = new_exponents
    shell['coefficients'] = new_coefficients

    return shell


def prune_basis(basis, use_copy=True):
    """
    Removes primitives that have a zero coefficient, and
    removes duplicate primitives and shells

    This only finds EXACT duplicates, and is meant to be used
    after other manipulations

    If use_copy is True, the input basis set is not modified.
    """

    if use_copy:
        basis = copy.deepcopy(basis)

    for k, el in basis['elements'].items():
        if 'electron_shells' not in el:
            continue

        shells = el.pop('electron_shells')
        shells = [prune_shell(sh, False) for sh in shells]

        # Remove any duplicates
        el['electron_shells'] = []

        for sh in shells:
            if sh not in el['electron_shells']:
                el['electron_shells'].append(sh)

    return basis


def uncontract_spdf(basis, max_am=0, use_copy=True):
    """
    Removes sp, spd, spdf, etc, contractions from a basis set

    The general contractions are replaced by uncontracted versions

    Contractions up to max_am will be left in place. For example,
    if max_am = 1, spd will be split into sp and d

    The input basis set is not modified. The returned basis
    may have functions with coefficients of zero and may have duplicate
    shells.

    If use_copy is True, the input basis set is not modified.
    """

    if use_copy:
        basis = copy.deepcopy(basis)

    for k, el in basis['elements'].items():

        if 'electron_shells' not in el:
            continue
        newshells = []

        for sh in el['electron_shells']:

            # am will be a list
            am = sh['angular_momentum']
            coeff = sh['coefficients']

            # if this is an sp, spd,...  orbital
            if len(am) > 1:
                newsh = sh.copy()
                newsh['angular_momentum'] = []
                newsh['coefficients'] = []

                ngen = len(sh['coefficients'])
                for g in range(ngen):
                    if am[g] > max_am:
                        newsh2 = sh.copy()
                        newsh2['angular_momentum'] = [am[g]]
                        newsh2['coefficients'] = [coeff[g]]
                        newshells.append(newsh2)
                    else:
                        newsh['angular_momentum'].append(am[g])
                        newsh['coefficients'].append(coeff[g])

                newshells.insert(0, newsh)

            else:
                newshells.append(sh)

        el['electron_shells'] = newshells

    return basis


def uncontract_general(basis, use_copy=True):
    """
    Removes the general contractions from a basis set

    The input basis set is not modified. The returned basis
    may have functions with coefficients of zero and may have duplicate
    shells.

    If use_copy is True, the input basis set is not modified.
    """

    if use_copy:
        basis = copy.deepcopy(basis)

    for k, el in basis['elements'].items():

        if 'electron_shells' not in el:
            continue

        newshells = []

        for sh in el['electron_shells']:
            # See if we actually have to uncontract
            # Also, don't uncontract sp, spd,.... orbitals
            #      (leave that to uncontract_spdf)
            if len(sh['coefficients']) == 1 or len(sh['angular_momentum']) > 1:
                newshells.append(sh)
            else:
                if len(sh['angular_momentum']) == 1:
                    for c in sh['coefficients']:
                        # copy, them replace 'coefficients'
                        newsh = sh.copy()
                        newsh['coefficients'] = [c]
                        newshells.append(newsh)

        el['electron_shells'] = newshells

    # If use_copy is True, we already made our deep copy
    return prune_basis(basis, False)


def uncontract_segmented(basis, use_copy=True):
    """
    Removes the segmented contractions from a basis set

    This implicitly removes general contractions as well,
    but will leave sp, spd, ... orbitals alone

    The input basis set is not modified. The returned basis
    may have functions with coefficients of zero and may have duplicate
    shells.

    If use_copy is True, the input basis set is not modified.
    """

    if use_copy:
        basis = copy.deepcopy(basis)

    for k, el in basis['elements'].items():

        if 'electron_shells' not in el:
            continue

        newshells = []

        for sh in el['electron_shells']:
            exponents = sh['exponents']
            nam = len(sh['angular_momentum'])

            for i in range(len(exponents)):
                newsh = sh.copy()
                newsh['exponents'] = [exponents[i]]
                newsh['coefficients'] = [["1.00000000E+00"] * nam]

                # Remember to transpose the coefficients
                newsh['coefficients'] = list(map(list, zip(*newsh['coefficients'])))

                newshells.append(newsh)

        el['electron_shells'] = newshells

    return basis


def make_general(basis, skip_spdf=False, use_copy=True):
    """
    Makes one large general contraction for each angular momentum

    If use_copy is True, the input basis set is not modified.

    The output of this function is not pretty. If you want to make it nicer,
    use sort_basis afterwards.
    """

    zero = '0.00000000'

    if use_copy:
        basis = copy.deepcopy(basis)

    if not skip_spdf:
        basis = uncontract_spdf(basis, 0, False)

    for k, el in basis['elements'].items():
        if 'electron_shells' not in el:
            continue

        newshells = []

        # See what we have
        all_am = []
        for sh in el['electron_shells']:
            am = sh['angular_momentum']

            # Skip sp shells
            if len(am) > 1:
                newshells.append(sh)
                continue

            if am not in all_am:
                all_am.append(am)

        all_am = sorted(all_am)

        for am in all_am:
            newsh = {
                'angular_momentum': am,
                'exponents': [],
                'coefficients': [],
                'region': '',
                'function_type': None,
            }

            # Do exponents first
            for sh in el['electron_shells']:
                if sh['angular_momentum'] == am:
                    newsh['exponents'].extend(sh['exponents'])

            # Number of primitives in the new shell
            nprim = len(newsh['exponents'])

            cur_prim = 0
            for sh in el['electron_shells']:
                if sh['angular_momentum'] != am:
                    continue

                if newsh['function_type'] is None:
                    newsh['function_type'] = sh['function_type']

                # Make sure the shells we are merging have the same function types
                ft1 = newsh['function_type']
                ft2 = sh['function_type']

                # Check if one function type is the subset of another
                # (should handle gto/gto_spherical, etc)
                if ft1 not in ft2 and ft2 not in ft1:
                    raise RuntimeError("Cannot make general contraction of different function types")

                ngen = len(sh['coefficients'])

                for g in range(ngen):
                    coef = [zero] * cur_prim
                    coef.extend(sh['coefficients'][g])
                    coef.extend([zero] * (nprim - len(coef)))
                    newsh['coefficients'].append(coef)

                cur_prim += len(sh['exponents'])

            newshells.append(newsh)

        el['electron_shells'] = newshells

    # If the basis was read in from a segmented format, it will have
    # duplicate primitives, and so a pruning is necessary
    prune_basis(basis, False)

    return basis


def _is_single_column(col):
    return sum(float(x) != 0.0 for x in col) == 1


def _free_primitives(coeffs):
    # Find which columns represent free primitives
    single_columns = [c for c in coeffs if _is_single_column(c)]
    if len(single_columns) == 0:
        return []
    # Now dig out the functions on those columns
    csum = [0.0 for k in range(len(single_columns[0]))]
    for c in single_columns:
        for k in range(len(csum)):
            csum[k] += abs(float(c[k]))
    # Since we're only looking at columns that represent free
    # primitives, the rows that have non-zero sums correspond to free
    # exponents.
    free_prims = [k for k, s in enumerate(csum) if s != 0.0]
    return free_prims


def remove_free_primitives(basis, use_copy=True):
    """
    Removes any free primitives from a basis set as a way to generate a minimal basis

    The input basis set is not modified. The returned basis
    may have functions with coefficients of zero and may have duplicate
    shells.

    If use_copy is True, the input basis set is not modified.
    """

    if use_copy:
        basis = copy.deepcopy(basis)

    for k, el in basis['elements'].items():

        if 'electron_shells' not in el:
            continue

        newshells = []
        for sh in el['electron_shells']:
            # Find contractions
            coefficients = sh['coefficients']
            contracted_columns = [idx for idx, c in enumerate(coefficients) if not _is_single_column(c)]
            coefficients = [coefficients[c] for c in contracted_columns]
            if len(coefficients):
                sh['coefficients'] = coefficients
                newshells.append(sh)
        el['electron_shells'] = newshells

    # We can now have exponents that aren't contracted so we need to
    # prune. If use_copy is True, we already made our deep copy
    return prune_basis(basis, False)


def optimize_general(basis, use_copy=True):
    """
    Optimizes the general contraction using the method of Hashimoto et al

    .. seealso :: | T. Hashimoto, K. Hirao, H. Tatewaki
                  | 'Comment on Dunning's correlation-consistent basis set'
                  | Chemical Physics Letters v243, Issues 1-2, pp, 190-192 (1995)
                  | https://doi.org/10.1016/0009-2614(95)00807-G

    """

    if use_copy:
        basis = copy.deepcopy(basis)

    # Make as generally-contracted as possible first
    basis = make_general(basis, skip_spdf=True, use_copy=False)

    for eldata in basis['elements'].values():

        if 'electron_shells' not in eldata:
            continue

        elshells = eldata['electron_shells']
        for sh in elshells:
            exponents = sh['exponents']
            coefficients = sh['coefficients']
            nprim = len(exponents)
            nam = len(sh['angular_momentum'])

            # Skip sp shells and shells with only one general contraction
            if nam > 1 or len(coefficients) < 2:
                continue

            # First, find columns (general contractions) with a single non-zero value
            single_columns = [idx for idx, c in enumerate(coefficients) if _is_single_column(c)]

            # Find the corresponding rows that have a value in one of these columns
            # Note that at this stage, the row may have coefficients in more than one
            # column. That is what we are looking for

            # Also, test to see that each row is only represented once. That is, there should be
            # no rows that are part of single columns (this would represent duplicate shells).
            # This can happen in poorly-formatted basis sets and is an error
            row_col_pairs = []
            all_row_idx = []
            for col_idx in single_columns:
                col = coefficients[col_idx]
                for row_idx in range(nprim):
                    if float(col[row_idx]) != 0.0:
                        if row_idx in all_row_idx:
                            raise RuntimeError("Badly-formatted basis. Row {} makes duplicate shells".format(row_idx))

                        # Store the index of the nonzero value in single_columns
                        row_col_pairs.append((row_idx, col_idx))
                        all_row_idx.append(row_idx)

            # Now for each row/col pair, zero out the entire row
            # EXCEPT for the column that has the single value
            for row_idx, col_idx in row_col_pairs:
                for idx, col in enumerate(coefficients):
                    if float(col[row_idx]) != 0.0 and col_idx != idx:
                        col[row_idx] = '0.0000000E+00'

    return basis


def geometric_augmentation(basis, nadd, use_copy=True, as_component=False, steep=False):
    '''Extends a basis set by adding extrapolated diffuse or steep functions.

    For augmented Dunning sets (aug), the diffuse augmentation
    corresponds to multiple augmentation (aug -> daug, taug, ...).

    In order for the augmentation to make sense, the two outermost
    primitives have to be free i.e. uncontracted.

    Parameters
    ----------
    basis : dict
        Basis set dictionary to work with
    nadd: int
        Number of functions to add (must be >=1). For diffuse augmentation on an augmented set: 1 -> daug, 2 -> taug, etc
    use_copy: bool
        If True, the input basis set is not modified.
    steep: bool
        If True, the augmentation is done for steep functions instead of diffuse functions.

    '''

    if nadd < 1:
        raise RuntimeError("Adding {} functions makes no sense for geometric_augmentation".format(nadd))

    # We need to combine shells by AM
    # I guess we don't NEED to, but it makes things a lot easier
    basis_copy = make_general(basis, use_copy=True)

    # We will now store the new shells in 'basis'
    # If as_component is specified, then create a new empty component basis
    if as_component:
        basis = skel.create_skel('component')
        basis['elements'] = {k: {} for k, v in basis_copy['elements'].items() if 'electron_shells' in v}
    elif use_copy:
        basis = copy.deepcopy(basis)

    # From Woon & Dunning, Jr
    # J. Chem. Phys. v100, No. 4, p. 2975 (1994)
    # DOI: 10.1063/1.466439
    #
    # The exponent for d-aug-cc-pVXZ: alpha*beta
    #                  t-aug-cc-pVXZ: alpha*(beta**2)
    # and so on.
    #
    # alpha = smallest exponent in aug-cc-pVXZ
    # beta = ratio of two most diffuse functions (beta < 1)
    #
    # This applies to all angular momentum (which have been combined into
    # shells already)

    for el_z, eldata in basis_copy['elements'].items():
        if 'electron_shells' not in eldata:
            continue

        new_shells = []

        for shell in eldata['electron_shells']:
            # Find the two smallest exponents. The smallest is alpha
            # beta is the ratio alpha/(next smallest)
            # Keep track of the indices as well
            exponents = [(float(x), idx) for idx, x in enumerate(shell['exponents'])]
            sorted_exponents = sorted(exponents)

            if len(sorted_exponents) < 2:
                # Need at least two exponents to perform augmentation
                continue

            if steep:
                # If we're augmenting by steep functions, the
                # references are the steepest and second-steepest
                # function.
                ref_idx = -1
                next_idx = -2
            else:
                # If we're augmenting by diffuse functions, the
                # references are the diffusemost and second-most
                # diffuse function.
                ref_idx = 0
                next_idx = 1

            ref_exp, ref_idx = sorted_exponents[ref_idx]
            next_exp, next_idx = sorted_exponents[next_idx]
            # Even-tempered spacing parameter
            beta = ref_exp / next_exp

            if ref_exp == next_exp:
                raise RuntimeError(
                    "The two outermost exponents are the same. Duplicate exponents are not a good thing here. Exponent: {}"
                    .format(ref_exp))

            # Test that the primitives for the references are free.
            free_prims = _free_primitives(shell['coefficients'])
            if (ref_idx not in free_prims) or (next_idx not in free_prims):
                # The shell does not have enough free primitives so
                # skip the extrapolation. (The alternative would be to
                # add in free primitives anyway, but then you'd have a
                # problem with contraction errors.)
                continue

            # Form new exponents
            new_exponents = []
            for i in range(1, nadd + 1):
                new_exponents.append(ref_exp * (beta**i))

            new_exponents = ['{:.6e}'.format(x) for x in new_exponents]

            # add the new exponents as new uncontracted shells
            for ex in new_exponents:
                new_shells.append({
                    'function_type': shell['function_type'],
                    'region': shell['region'],
                    'angular_momentum': shell['angular_momentum'],
                    'exponents': [ex],
                    'coefficients': [['1.00000000']]
                })

        # add the shells to the original basis set
        # (if this is a component basis, then 'electron_shells' may not exist)
        if 'electron_shells' not in basis['elements'][el_z]:
            basis['elements'][el_z]['electron_shells'] = []

        basis['elements'][el_z]['electron_shells'].extend(new_shells)

    return basis


def remove_primitive(electron_shell, idx_to_remove):
    '''Removes a primitive (by index) of a electron_shell

    The electron_shell passed in is modified
    '''

    electron_shell['exponents'].pop(idx_to_remove)
    for c in electron_shell['coefficients']:
        c.pop(idx_to_remove)

    # We may have removed the only primitive of a general contraction
    # So remove general contractions if all coefficients are zero
    electron_shell['coefficients'] = [
        gen for gen in electron_shell['coefficients'] if any([float(c) != 0.0 for c in gen])
    ]


def _element_remove_diffuse(eldata, nremove):
    if 'electron_shells' not in eldata:
        pass

    shells = eldata['electron_shells']
    max_am = misc.max_am(shells)
    if nremove == 'all':
        nremove = max_am + 1

    am_toremove = list(range(max_am, max_am - nremove, -1))

    # We may be asked to remove more than we have. That is ok
    am_toremove = [x for x in am_toremove if x >= 0]

    for shell in shells:
        shell_am = shell['angular_momentum']

        if len(shell_am) > 1:
            raise RuntimeError("Cannot not remove diffuse functions from fused shell: {}".format(shell_am))

        shell_am = shell_am[0]
        if shell_am not in am_toremove:
            continue

        # Find the most diffuse function and remove it
        exponents = [(float(x), idx) for idx, x in enumerate(shell['exponents'])]
        sorted_exponents = sorted(exponents)
        diffuse_exponent, diffuse_idx = sorted_exponents[0]

        remove_primitive(shell, diffuse_idx)


def truhlar_calendarize(basis, month, use_copy=True):
    '''Create the truhlar "month" basis sets from the corresponding aug basis set

    In Papajak 2011, removal of diffuse function stopped before removal of s and p functions.
    This, conceivably, extends to d functions for transition metals.
    However, you can keep extending to further in the year by removing these functions,
    although it may affect stability. This is implemented in Gaussian with the t(month)
    basis sets - tjul, tjun, etc. The 'tmonth' and regular 'month' basis sets are equivalent
    until the 'maug' basis set is reached (containg no diffuse functions on H,He, s,p,d on transition
    metals, and s, p on other elements).

    Since the regular 'month' basis sets are equivalent until maug, we do not adopt the
    t(month) nomenclature. Instead, you can just go further backwards through the months
    until you run out of diffuse functions or run out of months.

    Parameters
    ----------
    basis : dict
        Basis set dictionary to work with
    month: str
        Month to create ('apr', 'jul', etc). Not case sensitive
    use_copy: bool
        If True, the input basis set is not modified.
    '''

    valid_months = ['jul', 'jun', 'may', 'apr', 'mar', 'feb', 'jan']
    month = month.lower()
    if month not in valid_months:
        raise RuntimeError("Month {} is not valid for truhlar calendarization")

    # We need to combine shells by AM
    # I guess we don't NEED to, but it makes things a lot easier
    basis = make_general(basis, use_copy=use_copy)

    # Find out he max am for the entire basis
    all_am = [misc.max_am(x['electron_shells']) for x in basis['elements'].values()]
    basis_max_am = max(all_am)

    # Figure out what to remove
    # We do this by finding the offset from jul
    # The index represents the offset, since index('jul') = 0
    # This offset represents how many diffuse functions to remove from each element
    month_offset = valid_months.index(month)

    # Are we removing all diffuse functions from all elements?
    # If so, we will end up with all elements corresponding to the cc-pVXZ basis.
    # While it is ok for some elements to end up that way, we don't want the entire
    # basis set to be like that
    # Note that for a given basis_max_am, there will be (basis_max_am+1) diffuse functions
    if month_offset > basis_max_am:  # or month_offset >= (basis_max_am+1)
        raise RuntimeError("Will not remove {} diffuse functions for basis with max am = {}".format(
            month_offset, basis_max_am))

    # jul = same as aug, except remove all diffuse from H,He
    # First, remove diffuse functions from H,He
    # This will always be true - we don't handle 'aug'
    for el_z in ['1', '2']:
        if el_z not in basis['elements']:
            continue

        eldata = basis['elements'][el_z]
        _element_remove_diffuse(eldata, 'all')

    # short circuit if we are doing jul
    if month_offset == 0:
        basis = prune_basis(basis, False)
        return basis

    # Now the rest
    for el_z, eldata in basis['elements'].items():
        if el_z == '1' or el_z == '2':
            continue

        _element_remove_diffuse(eldata, month_offset)

    basis = prune_basis(basis, False)
    return basis


def split_blocked_contractions(basis, use_copy=True):
    '''Checks if the contraction coefficients in a general contraction can be made block diagonal and thereby split into two or more distinct shells


    Parameters
    ----------
    basis : dict
        Basis set dictionary to work with
    use_copy: bool
        If True, the input basis set is not modified.
    '''

    if use_copy:
        basis = copy.deepcopy(basis)

    for eldata in basis['elements'].values():

        if 'electron_shells' not in eldata:
            continue

        orig_shells = eldata['electron_shells']
        new_shells = []
        for sh in orig_shells:
            coefficients = sh['coefficients']
            ncontr = len(coefficients)
            nam = len(sh['angular_momentum'])
            # Skip sp shells and shells with only one general contraction
            if nam > 1 or ncontr == 1:
                new_shells.append(sh)
                continue

            exponents = sh['exponents']
            nprim = len(exponents)

            # Figure out which contractions share primitives between them
            shared_primitives = [[False for _ in range(ncontr)] for _ in range(ncontr)]
            for icontr in range(ncontr):
                for jcontr in range(icontr + 1):
                    for iprim in range(nprim):
                        if float(coefficients[icontr][iprim]) != 0.0 and float(coefficients[jcontr][iprim]) != 0.0:
                            shared_primitives[icontr][jcontr] = True
                            shared_primitives[jcontr][icontr] = True
                            break

            # Which contractions have been processed
            contraction_processed = [False for _ in range(ncontr)]

            # Indices of the contractions that are coupled
            blocks = []
            for icontr in range(ncontr):
                if not contraction_processed[icontr]:
                    # List of contracted functions in the block
                    block = [icontr]
                    contraction_processed[icontr] = True

                    # Form the list
                    for jcontr in range(icontr + 1, ncontr):
                        # No need to analyze functions that have already been processed
                        if contraction_processed[jcontr]:
                            continue
                        if shared_primitives[icontr][jcontr]:
                            block.append(jcontr)
                            contraction_processed[jcontr] = True
                    blocks.append(block)

            # Do we need to do anything?
            if len(blocks) == 1:
                # All functions are in a single block; we keep the shell as it is
                new_shells.append(sh)
                continue

            # Create new shells
            for block in blocks:
                # Identify the used primitives
                used_primitives = []
                for iprim in range(nprim):
                    for icontr in block:
                        if float(coefficients[icontr][iprim]) != 0.0:
                            if iprim not in used_primitives:
                                used_primitives.append(iprim)
                            continue

                # Form submatrices
                reduced_exponents = [exponents[p] for p in used_primitives]
                reduced_coefficients = [[coefficients[b][p] for p in used_primitives] for b in block]
                redsh = sh.copy()
                redsh['exponents'] = reduced_exponents
                redsh['coefficients'] = reduced_coefficients
                new_shells.append(redsh)

        # Replace the shells in the basis
        eldata['electron_shells'] = new_shells

    return basis
