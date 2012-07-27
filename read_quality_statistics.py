#!/usr/bin/env python

import preprocess_radtag_lane
smartopen = preprocess_radtag_lane.smartopen


import os,sys,numpy

def get_fastq_properties(fq):
    if smartopen(fq).read(1) == '@':
        lnum = 4
    else:
        lnum = 1
    print >> sys.stderr, 'fastq format lnum: %s' % lnum

    baseQ = None
    qfh = smartopen(fq)
    while baseQ is None:
        t,r,q = preprocess_radtag_lane.next_read_from_fh(qfh,lnum)
        baseQ = preprocess_radtag_lane.get_baseQ(q)
    qfh.close()
    print >> sys.stderr, 'fastq format baseQ: %s' % baseQ

    readlen = len(r)
    print >> sys.stderr, 'fastq format readlen: %s' % readlen

    return lnum,baseQ,readlen

if __name__ == "__main__":

    if len(sys.argv) == 2:
        fq = sys.argv[1]
        boundstr = "0:"
    else:
        fq, boundstr = sys.argv[1:]

    start,end = boundstr.split(':')
    start = int(start)

    lnum,baseQ,readlen = get_fastq_properties(fq)

    if end == '':
        end = readlen

    readcount = preprocess_radtag_lane.get_read_count(fq)

    qsc_n = 0
    qsc_tot = numpy.zeros(readlen)
    qsc_by_read = []

    fh = smartopen(fq)

    tickon = readcount/1000
    for i in range(readcount):
        if i % tickon == 0:
            print >> sys.stderr, '\r%0.1f' % ((i/float(readcount)) * 100),
        t,r,q = preprocess_radtag_lane.next_read_from_fh(fh,lnum)
        qsc = [ord(c)-baseQ for c in q]
        qsc_n += 1
        qsc_tot += qsc
        qsc_by_read.append(numpy.mean(qsc[start:end]))

    qsc_by_base = list(qsc_tot/qsc_n)

    print >> sys.stderr, 'write per-base mean qual ...',
    open(fq+'-per_base_qual.list','w').write(qsc_by_base.__repr__())
    print >> sys.stderr, 'done'
    print >> sys.stderr, 'write per-read qual ..',
    open(fq+'-per_read_qual.list','w').write(qsc_by_read.__repr__())
    print >> sys.stderr, 'done'
