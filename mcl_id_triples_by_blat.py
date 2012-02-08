#!/usr/bin/env python

'''
given a list of query files and a single subject, generates a sparse id matrix in the label input format specified by mcl

ALL SEQUENCE HEADERS MUST BEGIN WITH A UNIQUE INTEGER ID FOLLOWED BY A PERIOD

'''

import os, sys, re
from config import SCRATCH as scratch

keep_blat = False
pct_id_cut = 0.9

if __name__ == '__main__':
    s,q,args,outbase = sys.argv[1:]

    if keep_blat:
        outf = outbase+'.blat'
    else:
        try:
            os.makedirs(scratch)
        except:
            pass
        outf = os.path.join(scratch,'%s-%s.blat' % (os.path.basename(s),os.path.basename(q)))
        ret = os.system('touch %s' % outf)
        if not os.path.exists(outf):
            outf = outbase+'.blat'
        else:
            os.unlink(outf)

    matf = outbase+'.label'

    if os.path.exists(matf):
        print >> sys.stderr, 'output %s already present' % matf
        sys.exit(0)
    
    exitcode = os.system('blat %s %s %s %s' % (s,q,args,outf))
    if not exitcode == 0:
        if not keep_blat:
            os.unlink(outf)
        raise OSError, 'blat failed; exit'

    print >> sys.stderr, 'process blat output %s (%s bytes)' % (outf,os.path.getsize(outf))
    matfh = open(matf,'w')
    for l in open(outf):
        fields = l.strip().split()
        try:
            m = int(fields[0])
            topstrand = (fields[8] == '+')
        except:
            m = None
        if m and topstrand:
            s1 = fields[9]
            s2 = fields[13]
            l1 = len(s1.split('.')[2])
            l2 = len(s2.split('.')[2])
            if s1 != s2 and m/float(min((l1,l2))) >= pct_id_cut:
                matfh.write('%s\t%s\t%s\n' % (fields[9],fields[13],m))
    matfh.close()

    if not keep_blat:
        os.unlink(outf)
        
