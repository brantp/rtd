#!/usr/bin/env python

import os, sys
from config import RTDROOT
from preprocess_radtag_lane import smartopen
from rtd_run import readlen_from_uniqued
from glob import glob

if __name__ == "__main__":

    outfile, cutsite, engine, cores, parts, threads, radius = sys.argv[1:]

    readlen = readlen_from_uniqued(outfile)
    
    err_clust_root = outfile.rstrip('.gz') + '-rtd'
    cmd = os.path.join(RTDROOT,'rtd_run.py --cleanup -pe %s -np %s -nc %s -I %s -te %s -s %s -cs %s %s' % (engine,parts, cores, radius, threads, cutsite, err_clust_root, outfile))
    print >> sys.stderr, cmd
    ret = os.system(cmd)
    if ret != 0:
        raise OSError, 'rtd_run failed'

    cdest_file = glob(os.path.join(err_clust_root,'*I%s*.clstats.cdest' % radius))[0]
    cdest = float(open(cdest_file).read())
    errest = cdest/readlen
    print >> sys.stderr, 'sequencing error estimates: cluster dirt %0.4f, per-base error: %0.4f' % (cdest,errest)
    open(os.path.splitext(cdest_file)[0]+'.errest','w').write(errest.__repr__()) 
