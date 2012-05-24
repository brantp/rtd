#!/usr/bin/env python

import os,sys,numpy,Util
from collections import defaultdict
import preprocess_radtag_lane

def get_pool(db,ind,fc,lane,index=None):

    hits = [d['pool'] for d in db if d['sampleid'] == ind and d['flowcell'] == fc and d['lane'] == lane and (index is None or d['index'] == index)]

    return hits

def get_pool_lookup(db,fc,lane,index=None):
    all_ind = [d['sampleid'] for d in db if d.get('flowcell','') == fc and d['lane'] == lane and (index is None or d['index'] == index)]
    pool_lookup = {}
    for ind in all_ind:
        pool_lookup[ind] = get_pool(db,ind,fc,lane,index)
        if len(pool_lookup[ind]) != 1:
            raise ValueError, 'individual %s non-unique: %s' % (ind,pool_lookup[ind])
        else:
            pool_lookup[ind] = pool_lookup[ind][0]

    return pool_lookup

def get_counts_by_pool(uniqued,db):
    ufields = get_uniqued_info(uniqued)
    pool_lookup = get_pool_lookup(db,ufields[0],ufields[1],ufields[3])
    counts_by_pool = {}
    fh = preprocess_radtag_lane.smartopen(uniqued)
    for l in fh:
        f = l.split()
        for ind,ct in zip(f[3].split(','),[int(i) for i in f[4].split(',')]):
            pool = pool_lookup[ind]
            
            try:
                counts_by_pool[pool][ind] += ct
            except:
                counts_by_pool[pool] = defaultdict(int)
                counts_by_pool[pool][ind] += ct

    return counts_by_pool
            
def get_uniqued_info(uniqued):
    if 'index' in uniqued:
        ufields = os.path.splitext(os.path.basename(uniqued))[0].rsplit('_',3)
        ufields[3] = ufields[3][5:]
    else:
        ufields = os.path.splitext(os.path.basename(uniqued))[0].rsplit('_',2)
        ufields.append(None)

    ufields[1] = ufields[1][4:]

    return ufields


if __name__ == "__main__":

    db = preprocess_radtag_lane.get_table_as_dict('DB_library_data',suppress_fc_check=True)
    uniqued = sys.argv[1]

    ufields = get_uniqued_info(uniqued)

    counts_by_pool = get_counts_by_pool(uniqued,db)

    for k,v in counts_by_pool.items():
        print '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%0.1f\t%d' % (ufields[0],ufields[1],ufields[2],ufields[3],k,sum(v.values()),len(v),numpy.mean(v.values()),numpy.median(v.values()))
    
    
