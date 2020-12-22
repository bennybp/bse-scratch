'''
Helper functions for writing out basis set in various formats
'''

from math import ceil


def find_range(coeffs):
    '''
    Find the range in a list of coefficients where the coefficient is nonzero
    '''

    coeffs = [float(x) != 0 for x in coeffs]
    first = coeffs.index(True)
    coeffs.reverse()
    last = len(coeffs) - coeffs.index(True) - 1
    return first, last


def reshape(data, block_size):
    '''
    Reshape the input array as a matrix
    '''
    output_data = []
    for iblock in range(ceil(len(data) / block_size)):
        start = iblock * block_size
        end = min(len(data), (iblock + 1) * block_size)
        output_data.append(data[start:end])
    return output_data
