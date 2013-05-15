#!/usr/bin/env python

'''
given a list of query files and a single subject, generates a sparse id matrix in the label input format specified by mcl

ALL SEQUENCE HEADERS MUST BEGIN WITH A UNIQUE INTEGER ID FOLLOWED BY A PERIOD

'''

<<<<<<< HEAD
import os, sys, re, gzip
=======
import os, sys, re, numpy
>>>>>>> 9746a9d04d8803a838c0229517a1dacbdb12cb1b
from config import SCRATCH as scratch

fail_to_local = False #must set to "True" to allow use of working directory for blat output
                      #note that this can result in large .blat files on shared volumes!
keep_blat = False
pct_id_cut = 0.8
min_space_on_scratch = 100000 #in mb; fail if less than this amount of free space

def splitpath_rec(path, maxdepth=20):
     ( head, tail ) = os.path.split(path)
     return splitpath_rec(head, maxdepth - 1) + [ tail ] \
         if maxdepth and head and head != path \
         else [ head or tail ]

def splitpath(path):
    pathli = splitpath_rec(path)
    if pathli[0] == '/':
        pathli = pathli[1:]
        pathli[0] = '/'+pathli[0]
    return pathli

def space_free_on_volume(vol,unit='M',verbose=False):
    '''returns free space on the specified volume in units <unit>
    '''
    from subprocess import Popen,PIPE
    if verbose:
        print >> sys.stderr, 'checking free space on volume %s ...' % vol,
    try:
        free = int(Popen('df -P --sync -B %s %s' % (unit,vol), shell=True, stdout=PIPE).stdout.readlines()[-1].strip().split()[3].rstrip(unit))
    except:
        print >> sys.stderr, 'free space check failed; proceeding.  MONITOR AVAILABLE SPACE ON %s' % vol
        free = numpy.inf
    if verbose:
        print >> sys.stderr, '%s %sB' % (free,unit)
    return free

if __name__ == '__main__':
    s,q,args,outbase = sys.argv[1:]

    scratch_space = space_free_on_volume(splitpath(scratch)[0],verbose=True)
    if keep_blat:
        if fail_to_local:
             outf = outbase+'.blat'
        else:
             print >> sys.stderr, 'working directory failover not allowed'
             raise OSError, 'no space for blat output'
    elif scratch_space < min_space_on_scratch:
        print >> sys.stderr, 'available space on scratch (%s MB) is less than specified minimum (%s MB); use working directory' % ( scratch_space,min_space_on_scratch)
        if fail_to_local:
             outf = outbase+'.blat'
        else:
             print >> sys.stderr, 'working directory failover not allowed'
             raise OSError, 'no space for blat output'
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

    matf = outbase+'.label.gz'

    if os.path.exists(matf):
        print >> sys.stderr, 'output %s already present' % matf
        sys.exit(0)
    
    exitcode = os.system('blat %s %s %s %s' % (s,q,args,outf))
    if not exitcode == 0:
        if not keep_blat:
            os.unlink(outf)
        raise OSError, 'blat failed; exit'

    print >> sys.stderr, 'process blat output %s (%s bytes)' % (outf,os.path.getsize(outf))
    matfh = gzip.open(matf,'w')
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
        
