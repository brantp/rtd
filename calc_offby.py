#!/usr/bin/env python

import os,sys,re
import preprocess_radtag_lane,rtd_run
from collections import defaultdict

def reject_pair(s1,s2,dist):
    mm = 0
    for c1,c2 in zip(s1,s2):
        if c1 != c2:
            mm += 1
            if mm>dist:
                return True
    return False

if __name__ == "__main__":
    seqlen,dist,uniqueds,outfile = sys.argv[1:]

    all_quality = defaultdict(dict)
    
    for u in uniqueds.split(','):
        rtd_run.load_uniqued(all_quality,u,count_by_ind=True)
    
    seq_by_len = defaultdict(list)

    for k in all_quality.keys():
        seq_by_len[len(k)].append(k)

    seqs = seq_by_len[int(seqlen)]
    seqs.sort()

    offby = {}
    for i,s in enumerate(seqs):
        offby[s] = [si for si in seqs if si != s and not si in offby.keys() and not reject_pair(s,si,int(dist))]
        print >> sys.stderr, '\r%s / %s' % (i,len(seqs)),

    open(outfile,'w').write(offby.__repr__())
    
