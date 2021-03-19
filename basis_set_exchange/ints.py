'''
One-center integrals for Gaussian-type and Slater-type orbitals

Written by Susi Lehtola, 2020
'''

from math import gamma, sqrt

try:
    import numpy
    _use_numpy = True
except ImportError:
    _use_numpy = False


def _transform_numpy(C0, P0):
    '''Transforms the primitive integrals P into the contracted basis C using NumPy'''
    C = numpy.asarray(C0)
    P = numpy.asarray(P0)
    np_result = numpy.dot(numpy.dot(C, P), C.T)
    return [[np_result[i][j] for i in range(np_result.shape[0])] for j in range(np_result.shape[1])]


def _inner_product_numpy(Cleft0, P0, Cright0):
    '''Transforms the primitive integrals P into the contracted basis C using NumPy'''
    Cleft = numpy.asarray(Cleft0)
    Cright = numpy.asarray(Cright0)
    P = numpy.asarray(P0)
    ip = numpy.dot(Cleft, numpy.dot(P, Cright))
    return ip


def _matmul(A, B):
    '''Matrix multiply'''
    assert (len(A[0]) == len(B))
    result = [[0.0 for i in range(len(B[0]))] for j in range(len(A))]
    for i in range(len(A)):
        for j in range(len(B[0])):
            for k in range(len(B)):
                result[i][j] += A[i][k] * B[k][j]
    return result


def _transform_python(C, P):
    '''Transforms the primitive integrals P into the contracted basis C in pure Python'''

    # C transpose
    nrows = len(C)
    ncols = len(C[0])
    Ct = [[C[j][i] for j in range(nrows)] for i in range(ncols)]
    # Transform
    return _matmul(C, _matmul(P, Ct))


def _inner_product_python(Cleft, P, Cright):
    '''Calculates the inner product of Cleft and Cright against metric P in pure Python'''

    assert len(Cleft) == len(P)
    assert len(P[0]) == len(Cright)
    assert len(Cleft) == len(Cright)
    nprim = len(Cleft)

    CP = [sum([Cleft[i] * P[i][j] for i in range(nprim)]) for j in range(nprim)]
    return sum([CP[i] * Cright[i] for i in range(nprim)])


def _to_float(M):
    '''Transforms a vector or matrix from string representation to floating point representation'''

    nrows = len(M)
    if isinstance(M[0], list):
        # Matrix
        ncols = len(M[0])
        if isinstance(M[0][0], str):
            return [[float(M[i][j]) for j in range(ncols)] for i in range(nrows)]
        else:
            return M.copy()
    else:
        # Vector
        if isinstance(M[0], str):
            return [float(M[i]) for i in range(nrows)]
        else:
            return M.copy()


def _to_int(M):
    '''Transforms a vector from string representation to int representation'''
    nrows = len(M)
    if isinstance(M[0], str):
        return [int(M[i]) for i in range(nrows)]
    else:
        return M.copy()


def _transform(C, P):
    '''Transforms the primitive integrals P into the contracted basis C in numpy if it's available'''

    if _use_numpy:
        return _transform_numpy(C, P)
    else:
        return _transform_python(C, P)


def inner_product(Cleft, P, Cright):
    '''Calculates the inner product of Cleft and Cright in metric P'''

    if _use_numpy:
        return _inner_product_numpy(Cleft, P, Cright)
    else:
        return _inner_product_python(Cleft, P, Cright)


def _zero_matrix(N):
    '''Allocates a zero matrix of size N x N'''
    return [[0.0 for _ in range(N)] for _ in range(N)]


def _normalize_contraction(contr, ovl):
    '''Returns a normalized contraction matrix given the overlap matrix'''

    # Number of contractions
    nprim = len(contr[0])
    # Number of primitives
    ncontr = len(contr)

    # Size check: we want the primitive overlap
    assert (len(ovl) == nprim)
    assert (len(ovl[0]) == nprim)

    # Transform the overlap to the contracted basis
    ovl = _transform(contr, ovl)
    normfac = [1.0 / sqrt(ovl[i][i]) for i in range(ncontr)]

    # Form the normalized contraction matrix
    norm_contr = contr.copy()
    for icontr in range(ncontr):
        for iprim in range(nprim):
            norm_contr[icontr][iprim] *= normfac[icontr]
    return norm_contr


def _gto_overlap(exps, l):
    '''Computes the primitive overlap matrix for the given exponents,
    assuming the basis functions are of the spherical form r^l exp(-z
    r^2).

    '''
    assert (isinstance(l, int) and l >= 0)

    def ovl(a, b, l):
        ab = 0.5 * (a + b)
        return ((a * b) / (ab * ab))**(l / 2 + 3 / 4)

    # Initialize memory
    overlaps = _zero_matrix(len(exps))
    # Compute overlaps
    for i in range(len(exps)):
        for j in range(i + 1):
            overlaps[i][j] = ovl(exps[i], exps[j], l)
            overlaps[j][i] = overlaps[i][j]
    return overlaps


def gto_overlap_contr(exps0, contr0, l):
    '''Computes the overlap matrix in the contracted basis, assuming the
    basis functions are of the spherical form r^l \sum_i c_i exp(-z_i
    r^2).

    '''

    # Convert exponents and contractions to floating point
    exps = _to_float(exps0)
    contr = _to_float(contr0)
    # Get primitive integrals
    ovl = _gto_overlap(exps, l)
    return _transform(contr, ovl)


def _gto_Rsq(exps, l):
    '''Computes the r^2 matrix for the given exponents, assuming the basis
    functions are of the normalized spherical form r^l exp(-z r^2).

    '''
    assert (isinstance(l, int) and l >= 0)

    def rsq(a, b, l):
        ab = 0.5 * (a + b)
        return ((a * b) / (ab * ab))**(l / 2 + 3 / 4) / ab

    # Initialize memory
    rsqs = _zero_matrix(len(exps))
    # Compute overlaps
    prefactor = (3 / 4 + l / 2)
    for i in range(len(exps)):
        for j in range(i + 1):
            rsqs[i][j] = prefactor * rsq(exps[i], exps[j], l)
            rsqs[j][i] = rsqs[i][j]
    return rsqs


def gto_Rsq_contr(exps0, contr0, l):
    '''Computes the r^2 matrix in the contracted basis, assuming the basis
    functions are of the spherical form r^l \sum_i c_i exp(-z_i r^2).
    The function also takes care of proper normalization.

    '''

    # Convert exponents and contractions to floating point
    exps = _to_float(exps0)
    contr = _to_float(contr0)

    # Get primitive integrals
    rsqs = _gto_Rsq(exps, l)
    ovl = _gto_overlap(exps, l)

    # Normalize the contraction
    contr = _normalize_contraction(contr, ovl)
    # Transform to normalized contracted form
    return _transform(contr, rsqs)


def _sto_overlap(exps, ns):
    '''Computes the primitive overlap matrix for the given exponents and
    primary quantum numbers, assuming the basis functions are of the
    spherical form r^(n-1) exp(-z r).

    '''
    assert (len(exps) == len(ns))

    def ovl(za, zb, na, nb):
        return gamma(na + nb + 1) / sqrt(
            gamma(2 * na + 1) * gamma(2 * nb + 1)) * za**(na + 0.5) * zb**(nb + 0.5) / (0.5 * (za + zb))**(na + nb + 1)

    # Initialize memory
    overlaps = _zero_matrix(len(exps))
    # Compute overlaps
    for i in range(len(exps)):
        for j in range(i + 1):
            overlaps[i][j] = ovl(exps[i], exps[j], ns[i], ns[j])
            overlaps[j][i] = overlaps[i][j]
    return overlaps


def sto_overlap_contr(exps0, contr0, ns0):
    '''Computes the overlap matrix in the contracted basis, assuming the
    basis functions are of the spherical form r^(n-1) exp(-z r).

    '''
    # Convert exponents and contractions to floating point
    exps = _to_float(exps0)
    contr = _to_float(contr0)
    ns = _to_float(ns0)

    # Get primitive integrals
    ovl = _sto_overlap(exps, ns)
    return _transform(contr, ovl)


def _sto_Rsq(exps, ns):
    '''Computes the primitive r^2 matrix for the given exponents and
    primary quantum numbers, assuming the basis functions are of the
    spherical form r^(n-1) exp(-z r).

    '''
    assert (len(exps) == len(ns))

    def rsq(za, zb, na, nb):
        return gamma(na + nb + 3) / sqrt(
            gamma(2 * na + 1) * gamma(2 * nb + 1)) * za**(na + 0.5) * zb**(nb + 0.5) / (0.5 * (za + zb))**(na + nb + 3)

    # Initialize memory
    Rsqs = _zero_matrix(len(exps))
    # Compute Rsqs
    for i in range(len(exps)):
        for j in range(i + 1):
            Rsqs[i][j] = rsq(exps[i], exps[j], ns[i], ns[j])
            Rsqs[j][i] = Rsqs[i][j]
    return Rsqs


def sto_Rsq_contr(exps0, contr0, ns0):
    '''Computes the r^2 matrix in the contracted basis, assuming the basis
    functions are of the spherical form r^(n-1) exp(-z r).

    '''
    # Convert exponents and contractions to floating point
    exps = _to_float(exps0)
    contr = _to_float(contr0)
    ns = _to_int(ns0)

    # Get primitive integrals
    rsqs = _sto_Rsq(exps, ns)
    ovl = _sto_overlap(exps, ns)

    # Normalize the contraction
    contr = _normalize_contraction(contr, ovl)
    # Transform to normalized contracted form
    return _transform(contr, rsqs)
