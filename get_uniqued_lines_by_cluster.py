#!/usr/bin/env python

import numpy,os,sys,re
from preprocess_radtag_lane import get_baseQ


def load_cluster_data(gr,tab,mID=None):
    '''given an mcl graph output file, an mcl tab file of sequence labels, and an mID file of individuals per label,

    returns (clusters,labels,mID_by_label) dicts
    '''

    #extract clusters from mcl output
    print >>sys.stderr, 'load graph...',
    fc = open(gr).read()
    print >> sys.stderr, 'done\nparse graph...',
    body = re.search('begin(.+?)\)',fc.replace('\n',' ')).groups()[0]
    clusters = dict([(s.strip().split()[0],s.strip().split()[1:]) for s in body.strip().split('$') if len(s.strip().split()) > 1])
    print >> sys.stderr, 'done\nload labels...',

    #cluster labels
    labels = dict([l.strip().split() for l in open(tab).readlines()])

    print >> sys.stderr, 'done\nload mIDs...',
    #individuals by labels
    if mID is None:
        mID_by_label = None
    else:
        mID_by_label = dict([(l.split()[0],l.split()[1:]) for l in open(mID).readlines()])
    print >> sys.stderr, 'done'

    return clusters,labels,mID_by_label

def load_lines_from_uniqued(source_uniques,rv_sort = True, sort_key = lambda x: (len(x[0]),int(x[1])), keep_source_id = False):
    '''
    if keep_source_id is True
        returns list of 2-tuples uniqued_id (eg 100617_lane6_PE for "data/100617/100617_lane6_PE.uniqued")
        tuples are (parsed_lines,uniqued_id)

    else list of lines.
    '''
    uniquedlines = []
    for f in source_uniques:
        lines = []

        print >> sys.stderr, 'load %s ...' % f,
        lines = tuple([l.strip().split() for l in open(f).readlines()])
        print >> sys.stderr, '%s lines' % len(lines)

        #get qual base
        baseQ = None
        for l in lines:
        	baseQ = get_baseQ(l[2])
        	if baseQ is not None:
        		break
        print >> sys.stderr, 'qual base: %s' % baseQ
        
        if baseQ == 64: 
        	print >> sys.stderr, 'Translate quality encoding to base 33 ...',
        	for l in lines:
        		l[2] = ''.join([chr(ord(c)-64+33) for c in l[2]])
        	print >> sys.stderr, 'done'
        
        if keep_source_id:
            uniqued_id = os.path.basename(os.path.splitext(f)[0])
            uniquedlines.extend( zip( lines,[uniqued_id]*len(lines) ) )
        else:
            uniquedlines.extend(lines)

    print >> sys.stderr, 'sort',
    if keep_source_id:
        uniquedlines.sort(reverse = rv_sort,key = lambda x: sort_key(x[0]))
    else:
        uniquedlines.sort(reverse = rv_sort,key = sort_key)
    print >> sys.stderr, 'done'
    return uniquedlines


if __name__ == '__main__':

    gr,tab = sys.argv[1:3]
    source_uniques = sys.argv[3:]

    nticks = 100

    clusters,labels,mID_by_label = load_cluster_data(gr,tab)

    clid_by_lid = {}
    for clid in clusters.keys():
        for lblid in clusters[clid]:
            clid_by_lid[lblid] = clid

    lid_by_seq = {}
    for k in labels.keys():
        lid_by_seq[labels[k].split('.')[2]] = k

    try_lens = sorted(list(set([len(k) for k in lid_by_seq])),reverse=True)

    seqlen = try_lens[0]

    uniquedlines = load_lines_from_uniqued(source_uniques, rv_sort = False, sort_key = lambda x: x[0][:seqlen],keep_source_id=True)

    tickon = len(uniquedlines) / nticks

    print >> sys.stderr, 'process clusters'

    for i,(l,uniqued_id) in enumerate(uniquedlines):
        if i % tickon == 0: print >> sys.stderr, '%s / %s lines' % (i,len(uniquedlines))
        lid = None
        try:
            lid = lid_by_seq[l[0][:seqlen]]
        except:
            for sl in try_lens:
                try:
                    lid = lid_by_seq[l[0][:sl]]
                    break
                except:
                    pass
        if lid is None: 
            pass
        else:
            print '\t'.join([clid_by_lid[lid],'%s.%s' % (lid,uniqued_id)]+l)

