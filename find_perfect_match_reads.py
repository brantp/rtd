#!/usr/bin/env python
'''
given a uniqued file and set of new reads, passes only those reads that perfectly match a uniqued entry present at or above some copy number
'''

import os,sys
import numpy
from editdist import distance
from preprocess_radtag_lane import next_read_from_fh, smartopen, get_read_count

idx_bp = 5
cut_bp = 5

lnum = 4
min_seqs = 7

uniqued, fastq = sys.argv[1:]

readlen = len(next_read_from_fh(smartopen(fastq),4)[1])

print >> sys.stderr, 'readlen: %s' % readlen

num_reads = get_read_count(fastq,4)
tickon = num_reads/200

useqs = []
for l in open(uniqued):
    s,cntstr = l.strip().split()[0], l.strip().split()[4]
    cnt = numpy.mean([int(i) for i in cntstr.split(',')])
    if cnt >= min_seqs:
        useqs.append(s[cut_bp:readlen-idx_bp])

useqs = list(set(useqs))
print >> sys.stderr, '%s unique %sbp sequences in uniqued file' % (len(useqs),len(s[cut_bp:readlen-idx_bp]))

fh = smartopen(fastq)

for i in range(num_reads):
    if i%tickon==0: print >> sys.stderr, '%s / %s' % (i,num_reads)
    line = next_read_from_fh(fh)
    found = False
    for seq in useqs:
        if line[1][idx_bp+cut_bp:] == seq:
            found = True
            break

    if found:
        print '\n'.join(['@'+line[0],line[1],'+',line[2]])
