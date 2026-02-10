import os
import gzip

rORForise_VERSION = "v0.0.4"


def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))


def fix_path(path):
    fixed_path = os.path.normpath(path)
    fixed_path = os.path.realpath(fixed_path)
    return fixed_path


def open_file(file):
    if file.endswith('.gz'):
        return gzip.open(file, 'rb')
    else:
        return open(file, 'rb')