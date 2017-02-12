#! /usr/bin/env python
# -*- coding: utf-8 -*-

from itertools import tee, izip
import math

import sysutils

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

def merge_interval_list(ivs, dist=1):
    """ Merge intervals

    Args:
        ivs (list): List of intervals. Each interval is represented by a tuple of
            integers (start, end) where end > start.
        dist (int): Distance between intervals to be merged. Setting dist=1 (default) will
            merge adjacent ("bookended") intervals

    Returns:
        list: Merged list of intervals

    Examples:
        >>> merge_interval_list([])
        []
        >>> merge_interval_list([(1,10)])
        [(1, 10)]
        >>> merge_interval_list([(4, 9), (10, 14), (1, 3)])
        [(1, 3), (4, 9), (10, 14)]
        >>> merge_interval_list([(4, 9), (10, 14), (1, 3)], dist=1)
        [(1, 14)]
    """
    if len(ivs)<= 1: return ivs
    ivs.sort(key=lambda x:x[0])
    ret = [ivs[0]]
    for iv in ivs[1:]:
        if iv[0] - ret[-1][1] > dist:
            ret.append(iv)
        else:
           ret[-1] = (ret[-1][0], max(iv[1],ret[-1][1]))
    return ret

def overlaps(iv1, iv2):
    """ Returns True if iv1 overlaps iv2
    """
    return (min(iv1[1], iv2[1]) - max(iv1[0], iv2[0])) > 0

def guess_encoding(fh, nsamp=100):
    """ Guess the encoding of a fastq file

    S - Sanger        Phred+33,  raw reads typically (0, 40).  ASCII 33-73
    X - Solexa        Solexa+64, raw reads typically (-5, 40). ASCII 59-104
    I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40).  ASCII 64-104
    J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40).  ASCII 66-104
        with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator (bold)
    L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)   ASCII 33-74
    
    Args:
        fh (str or file): Fastq file to guess encoding
        
    """
    fh = sysutils.get_filehandle(fh)

    # Initialize min and max
    minq, maxq = ('z', '!')

    for i,l in enumerate(fh):
        if i % 4 == 3:
            """ Solexa+64 can be as low as -5, so if there are any ASCII characters below
                59 (64-5), it is definitely phred-33.
                chr(64-5) == ';'
            """
            if any(c < ';' for c in l.strip('\n')):
                return 'Phred+33'
            """ Phred+33 maxes out at 74 (33+41) for Illumina data, but PacBio QV can be
                above 60. We'll assume that we will not see any QVs above 64. If there are
                ASCII characters above 97 (33+64) we will assume this is Phred+64.
                chr(33+64) == 'a'
            """
            if any(c >= 'a' for c in l.strip('\n')):
                return 'Phred+64'
            """ Otherwise just set the overall min and max values """
            minq = min(minq+l)
            maxq = max(maxq+l)
        if i >= nsamp*4:
            break
    
    """ All sample lines have been processed and no characters below chr(64-5) == ';'
        (indicating Phred+33) or above chr(33+64) == 'a' (indicating Phred+64) were found.
        We will assume that the file is Phred+64 if qualities greater than chr(75) == 'K'
        were observed, since 1.8+ Illumina reads max out at Q41 ('J').
        
    """
    return 'Phred+64' if maxq > 'K' else 'Phred+33'

def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)
    
def percentile(d, p):
    """ Returns percentile from sorted list """
    k = (len(d)-1) * p
    return (float(d[int(math.floor(k))] + d[int(math.ceil(k))]) / 2)
