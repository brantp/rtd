#!/usr/bin/env python
'''
provides simulation for ddRAD recovery modeled on pippin prep automated size
selection (assuming normally distributed size selection).
In our experience, the pippin prep lower and upper set bounds correspond to
a 4SD range, i.e. pippin settings of 260 - 340 are modeled by mean of 300,
standard deviation of 20 (300+/-40bp = 300+/- 2SD; 1SD = 20)
'''
import matplotlib
matplotlib.use('Agg')

import pylab,numpy, os, sys, re
from preprocess_radtag_lane import smartopen as open

def fraglen_from_seq(seq,e1,e2):
    sites_by_pos = sorted(reduce(lambda x,y:x+y,[[(m.start(),e) for m in re.finditer(e,seq)] for e in [e1,e2]]))
    fraglens = []
    lastpos = None
    lastsite = None
    for pos,site in sites_by_pos:
        if lastpos and lastsite and lastsite != site:
            fraglens.append(pos-lastpos)
        lastpos = pos
        lastsite = site
    return fraglens


def fa_cuts(fname,e1,e2):
    listfile = '%s_%s-%s.list' % (fname,e1,e2)
    if os.path.exists(listfile):
        print >> sys.stderr, '%s found, load fragment lengths ...' % listfile,
        fraglens = eval(open(listfile).read())
        print >> sys.stderr, '%s fragments loaded' % len(fraglens)
        return fraglens
    
    fraglens = []
    this = []
    seqnum = 0
    for l in open(fname):
        if l[0] == '>':
            if this:
                seq = ''.join(this)
                fraglens.extend(fraglen_from_seq(seq,e1,e2))
                this = []
            sname = l[1:].strip()
            seqnum += 1
            outn = ' (sequence %s)' % (seqnum)
            om = '\r%s %s' % (sname[:80-len(outn)], outn)
            print >> sys.stderr, om, ' '*(80-len(om)+1),
        else:
            this.append(l.strip().upper())
        
    seq = ''.join(this)
    fraglens.extend(fraglen_from_seq(seq,e1,e2))
    print >> sys.stderr, 'store %s fragment lengths to %s' % (len(fraglens),listfile)
    open(listfile,'w').write(fraglens.__repr__())
    return fraglens

def countdict(li):
    '''return a dictionary where keys are unique list values and values are counts of items of that value'''
    d = {}.fromkeys(set(li),0)
    for l in li:
        d[l] += 1
    return d


def simulate_coverage(reads,size_mean,size_sd,genome_frag_counts,max_frag_size=10000):
    #simulate read distribution from sizes
    rc = countdict(numpy.random.normal(size_mean,size_sd,size=reads).astype('int'))
    #simulate breakdown of reads within a size
    coverage = []
    for fragsize in xrange(min(max(genome_frag_counts.keys()),max_frag_size)):
        if genome_frag_counts.get(fragsize,0):
            counts_this_size = countdict(numpy.random.randint(0,genome_frag_counts.get(fragsize,0),size=rc.get(fragsize,0)))
        coverage.extend([counts_this_size.get(s,0) for s in xrange(genome_frag_counts.get(fragsize,0))])
    return coverage

def coverage_curve(X,size_mean,size_sd,genome_frag_counts,cov_cut,max_frag_size=10000):
    '''
    given:
    vector of read counts X
    size mean and SD
    dictionary of fragment counts by size in basepairs
    cutoff for cluster coverage to consider "sampled"

    returns number of sampled clusters
    '''
    Y = []
    print >> sys.stderr, 'simulate coverage for %s readcount conditions (%s - %s)' % (len(X),min(X),max(X))
    for i,x in enumerate(X):
        print >> sys.stderr, '\r%s (%s)' % (x,i),
        Y.append(len(numpy.nonzero(numpy.array(simulate_coverage(x,size_mean,size_sd,genome_frag_counts,max_frag_size=max_frag_size))>=cov_cut)[0]))
    return Y

def mean_sd_from_bounds(lo,hi,adapt_len=76,num_sd=2):
    inslo = lo-adapt_len
    inshi = hi-adapt_len
    m = inslo + ((inshi-inslo)/2)
    s = (m-inslo)/2
    print >> sys.stderr, '%s to %s size selection; inserts are %s to %s; %s +/- %s; SD %s' % (lo,hi,inslo,inshi,m,m-inslo,s)
    return m,s

if __name__ == "__main__":

    import argparse

    ds =  ' [%(default)s]'
    #create command line parser
    parser = argparse.ArgumentParser(description='simulate ddRAD coverage from genome sequence.  Supply either pippin bounds (-pl, -ph) or mean and SD (-m, -s). Required arguments: two enzyme recognition sequences and genome fasta file')

    parser.add_argument('-t','--threshold',default=10,type=int,help='read coverage condition for cluster inclusion'+ds)

    parser.add_argument('-m','--mean',type=int,default=0,help='mean of size selection'+ds)
    parser.add_argument('-s','--sd',type=float,default=0,help='standard deviation of size selection'+ds)

    parser.add_argument('-pl','--pippin_low',type=int,default=0,help='lower bound for pippin collection range'+ds)
    parser.add_argument('-ph','--pippin_high',type=int,default=0,help='upper bound for pippin collection range'+ds)

    parser.add_argument('-xm','--xmax',type=int,default=2000000,help='max readcount to simulate'+ds)
    parser.add_argument('-xs','--xstep',type=int,default=10000,help='step size of readcounts to simulate'+ds)

    parser.add_argument('-mfs','--max_frag_size',type=int,default=10000,help='maximum fragment size to consider in simulation'+ds)

    parser.add_argument('e1',help='recognition sequence for enzyme 1')
    parser.add_argument('e2',help='recognition sequence for enzyme 2')
    
    parser.add_argument('fa',help='genome fasta file')

    opts = parser.parse_args()

    X = range(0,opts.xmax,opts.xstep)

    frags_by_len = countdict(fa_cuts(opts.fa,opts.e1,opts.e2))

    if opts.mean and opts.sd:
        m = opts.mean
        s = opts.sd
    elif opts.pippin_low and opts.pippin_high:
        m,s = mean_sd_from_bounds(opts.pippin_low,opts.pippin_high)
    else:
        raise ValueError, 'either -m and -s OR -pl and -ph must be set'

    Y = numpy.array(coverage_curve(X,m,s,frags_by_len,opts.threshold,opts.max_frag_size))

    plotname = '%s_%s-%s_mean%s_sd%s_threshold%s.pdf' % (opts.fa,opts.e1,opts.e2,m,s,opts.threshold)
    

    pylab.plot(X,Y,'k')
    pylab.plot(X,Y*0.75,'--k')
    pylab.plot(X,Y*0.5,':k')
    pylab.savefig(plotname)

    print >> sys.stderr, '\nplot stored as %s' % plotname
