#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')

import preprocess_radtag_lane, sam_from_clust_uniqued
from collections import defaultdict
import os,sys,numpy,pylab
from glob import glob

open = preprocess_radtag_lane.smartopen

#taken from py_util/iplot 20120708

histcolors = {'b':'blue',
              'g':'green',
              'r':'red',
              'c':'cyan',
              'm':'magenta',
              'y':'yellow',
              'k':'black',
              'w':'white',}

dig = range(9)+['A','B','C','D','E','F']

spectrum = ['#FF%s%s00' % (i,j) for i in dig for j in dig] + ['#%s%sFF00' % (i,j) for i in dig[::-1] for j in dig[::-1]] + ['#00FF%s%s' % (i,j) for i in dig for j in dig] +['#00%s%sFF' % (i,j) for i in dig[::-1] for j in dig[::-1]] + ['#%s%s00FF' % (i,j) for i in dig for j in dig]

def subspectrum(numcol):
    '''given a number of colors desired (<= len spectrum)
    returns a list of hex codes for colors'''
    #pad spectrum if undersized
    spec = spectrum
    while len(spec) < numcol:
        spec = list(reduce(lambda x,y:x+y,zip(spec,spec)))
    try:
        step = len(spec)/numcol
        return spec[::step][:numcol]
    except:
        return []

def subspec_enum(iterable):
    return zip(subspectrum(len(iterable)),iterable)

def readcounts_from_uniqueds(uniqueds):
    '''given a list of .uniqued files
    returns dict of readcounts for individuals'''

    readcounts = defaultdict(int)

    for uni in uniqueds:
        fh = open(uni)
        print >> sys.stderr, 'get counts for %s ...' % uni,
        lines = 0
        for l in fh:
            lines += 1
            fields = l.split()
            inds = fields[3].split(',')
            counts = map(int,fields[4].split(','))
            for ind,count in zip(inds,counts):
                readcounts[ind] += count
        print >> sys.stderr, '%s lines read' % lines

    return readcounts

def rc_from_bam(bam,rc_bam=None):
    from subprocess import Popen, PIPE
    rg_lookup = {}
    for l in Popen('samtools view -H %s  | grep "@RG"' % bam, shell=True,stdout=PIPE).stdout.readlines():
        f = l.strip().split()
        rg_lookup[f[1][3:]] = f[4][3:]
    if rc_bam is None:
        rc_bam = defaultdict(int)
    for l in Popen('samtools view -F 0x0080 %s' % bam, shell=True,stdout=PIPE).stdout:
        rc_bam[rg_lookup[[f for f in l.strip().split() if f.startswith('RG:Z')][0][5:]]] += 1
    return rc_bam

def readcounts_from_bams(bams):
    readcounts = defaultdict(int)

    for b in bams:
        print >> sys.stderr, 'get counts for %s' % b
        readcounts = rc_from_bam(b,readcounts)

    return readcounts

def readcounts_from_cl_lines(cl_lines):
    '''given a list of lines corresponding to a single cluster from a .cluni file
    returns a dictionary of readcounts for each individual in that cluster'''

    readcounts = defaultdict(int)
    for l in cl_lines:
        fields = l.split()
        inds = fields[5].split(',')
        counts = map(int,fields[6].split(','))
        for ind,count in zip(inds,counts):
            readcounts[ind] += count

    return readcounts

def readcounts_from_cluni(cluni):

    fh = open(cluni)
    readcounts = []
    cl_dirt = []
    all_indiv = set([])
    print >> sys.stderr, 'load individual readcounts from %s ...' % cluni,
    cl_lines = sam_from_clust_uniqued.next_cluster_lines(fh)
    while cl_lines:
        rc_this = readcounts_from_cl_lines(cl_lines)
        readcounts.append(rc_this)
        all_indiv = all_indiv.union(set(rc_this.keys()))
        cl_dirt.append(sam_from_clust_uniqued.calc_cluster_dirt(cl_lines))
        cl_lines = sam_from_clust_uniqued.next_cluster_lines(fh)
    print >> sys.stderr, 'data for %s individuals loaded' % len(all_indiv)
    
    return readcounts,cl_dirt,all_indiv

def load_sorted_data(uniqueds,cluni,sortfn=sum):
    rc,cd,inds = readcounts_from_cluni(cluni)
    rc_orig = readcounts_from_uniqueds(uniqueds)
    rc_sort,cd_sort = sort_lists(rc,cd,sortfn)

    cd_cut = None
    print >> sys.stderr, 'attempt to load cluster dirt threshold:'
    for f in glob(os.path.splitext(cluni.rstrip('.gz'))[0]+'*.clstats.cdest'):
        print >> sys.stderr, '\t%s ...' % f,
        if os.path.getsize(f) > 0:
            cd_cut = float(open(f).read())
            print >> sys.stderr, '%s' % cd_cut
            break
        print >> sys.stderr, 'failed'
    
    return rc_sort,cd_sort,list(inds),rc_orig,cd_cut

def lol_by_segment(li,seglen):
    lol = []
    for i in xrange(0,len(li),seglen):
        lol.append(li[i:i+seglen])
    return lol

def mat_from_list_of_dict(list_of_dict,labels=None):
    if labels is None:
        labels = reduce(lambda x,y: x.union(y), [set(d.keys()) for d in list_of_dict])

    return numpy.array([numpy.array([d.get(lab,0) for lab in labels]) for d in list_of_dict])

def sort_lists(readcounts,cl_dirt,fn):
    '''sort readcounts and dirt by function fn (descending) of readcounts list-of-dict'''
    sort_order = [fn(d.values()) for d in readcounts]
    rc_sort = zip(*sorted(zip(sort_order,readcounts),reverse=True))[1]
    cd_sort = zip(*sorted(zip(sort_order,cl_dirt),reverse=True))[1]
    return rc_sort,cd_sort

def draw_efficiency_plots(rc_sort,cd_sort,inds,rc_orig,cd_cut=None,fignum=1,figsize=(8,10),filename=None):
    if cd_cut is None:
        cd_cut = 0.1
        
    pylab.figure(fignum,figsize)

    print >> sys.stderr, 'draw plots'
    
    pylab.subplot(4,1,1)
    for c,ind in subspec_enum(sorted(inds,key=lambda x: rc_orig[x])):
        pylab.plot(numpy.cumsum([d.get(ind,0) for d in rc_sort]),c=c)

    pylab.subplot(4,1,2)
    for c,ind in subspec_enum(sorted(inds,key=lambda x: rc_orig[x])):
        pylab.plot(numpy.cumsum([d.get(ind,0) for d in rc_sort])/float(rc_orig[ind]),c=c)
    pylab.ylim(0,1)
    
    pylab.subplot(4,1,3)
    for c,ind in subspec_enum(sorted(inds,key=lambda x: rc_orig[x])):
        pylab.plot(numpy.cumsum([d.get(ind,0) for this_cd,d in zip(cd_sort,rc_sort) if this_cd <= cd_cut]),c=c)

    pylab.subplot(4,1,4)
    for c,ind in subspec_enum(sorted(inds,key=lambda x: rc_orig[x])):
        pylab.plot(numpy.cumsum([d.get(ind,0) for this_cd,d in zip(cd_sort,rc_sort) if this_cd <= cd_cut])/float(rc_orig[ind]),c=c)
    pylab.ylim(0,1)

    if filename is not None:
        print >> sys.stderr, 'store plots'
        try:
            pylab.savefig(filename)
        except IOError:
            print >> sys.stderr, 'unable to write %s, output not stored' % filename
    
def draw_ind_by_clust_plots(rc_sort,cd_sort,inds,rc_orig,cd_cut=None,win=1000,rc_low=4,rc_hi=10,fignum=1,figsize=(8,10),filename=None):
    if cd_cut is None:
        cd_cut = 0.1
        
    pylab.figure(fignum,figsize)

    lol = lol_by_segment(rc_sort,win)
    cdlol = lol_by_segment(cd_sort,win)
    ncat = len(lol)
    step = ncat/20
    print >> sys.stderr, 'draw boxplots'

    pylab.subplot(4,1,1)
    pylab.boxplot([[len([v for v in d.values() if v>=rc_low]) for d in li] for li in lol])
    pylab.xticks(numpy.arange(0,ncat,step),(numpy.arange(0,ncat,step)*win)/1000,rotation=90)
    
    pylab.subplot(4,1,2)
    pylab.boxplot([[len([v for v in d.values() if v>=rc_hi]) for d in li] for li in lol])
    pylab.xticks(numpy.arange(0,ncat,step),(numpy.arange(0,ncat,step)*win)/1000,rotation=90)
    
    pylab.subplot(4,1,3)
    pylab.boxplot([[len([v for v in d.values() if v>=rc_low]) for this_cd,d in zip(cdli,li) if this_cd <= cd_cut] for cdli,li in zip(cdlol,lol)])
    pylab.xticks(numpy.arange(0,ncat,step),(numpy.arange(0,ncat,step)*win)/1000,rotation=90)
    
    pylab.subplot(4,1,4)
    pylab.boxplot([[len([v for v in d.values() if v>=rc_hi]) for this_cd,d in zip(cdli,li) if this_cd <= cd_cut] for cdli,li in zip(cdlol,lol)])
    pylab.xticks(numpy.arange(0,ncat,step),(numpy.arange(0,ncat,step)*win)/1000,rotation=90)

    if filename is not None:
        print >> sys.stderr, 'store boxplots'
        try:
            pylab.savefig(filename)
        except IOError:
            print >> sys.stderr, 'unable to write %s, output not stored' % filename

def draw_clust_by_reads_scatter(rc_sort,cd_sort,inds,rc_orig,cd_cut=None,rc_low=4,rc_hi=10,simstep=10000,xmax=None,fignum=1,figsize=(8,10),filename=None):
    if cd_cut is None:
        cd_cut = 0.1

    if xmax is None:
        xmax = min(max(rc_orig.values())*2,sum(rc_orig.values()),numpy.mean(rc_orig.values())*5)
    else:
        xmax = int(xmax)

    pylab.figure(fignum,figsize)

    allreads = sum(rc_orig.values())
    sumrc = [sum(d.values()) for d in rc_sort]
    sumrc1 = numpy.array(sumrc)/float(sum(sumrc))
    sumy = len([d for d in rc_sort if sum(d.values()) >= rc_low])

    try:
        ymax = [y for x,y in zip(range(0,allreads,simstep),[len(numpy.nonzero((sumrc1*c)>=rc_low)[0]) for c in range(0,allreads,simstep)]) if x>xmax][0]
    except:
        ymax = sumy

    print >> sys.stderr, 'draw scatters'

    ax = pylab.subplot(2,1,1)
    ax.set_xticks(range(0,int(xmax),100000),minor=True)
    ax.set_yticks(range(0,int(ymax*1.2),1000),minor=True)
    ax.grid(which='minor')
    pylab.scatter([rc_orig[ind] for ind in inds],[len([d for d in rc_sort if d.get(ind,0) >= rc_low]) for ind in inds],facecolor='none',edgecolor='y')
    pylab.scatter([rc_orig[ind] for ind in inds],[len([d for this_cd,d in zip(cd_sort,rc_sort) if d.get(ind,0) >= rc_low and this_cd <= cd_cut]) for ind in inds],facecolor='none',edgecolor='b')
    pylab.scatter([sum(rc_orig.values())],[sumy],facecolor='k')
    pylab.plot(range(0,allreads,simstep), [len(numpy.nonzero((sumrc1*c)>=rc_low)[0]) for c in range(0,allreads,simstep)],'--k')
    pylab.plot(range(0,allreads,simstep), [len(numpy.nonzero((sumrc1*c)>=rc_low)[0])*0.75 for c in range(0,allreads,simstep)],'--k')
    pylab.plot(range(0,allreads,simstep), [len(numpy.nonzero((sumrc1*c)>=rc_low)[0])*0.5 for c in range(0,allreads,simstep)],'--k')
    pylab.text(20000,ymax,'mean: %0.2f, cv: %0.2f\n%s indiv, cd: %0.2f' % (numpy.mean([rc_orig[ind] for ind in inds]),numpy.std([rc_orig[ind] for ind in inds])/numpy.mean([rc_orig[ind] for ind in inds]),len(inds),cd_cut))
    pylab.ylim(0,ymax*1.2)
    pylab.xlim(0,xmax)

    
    sumy = len([d for d in rc_sort if sum(d.values()) >= rc_hi])
    try:
        ymax = [y for x,y in zip(range(0,allreads,simstep),[len(numpy.nonzero((sumrc1*c)>=rc_hi)[0]) for c in range(0,allreads,simstep)]) if x>xmax][0]
    except:
        ymax = sumy

    ax = pylab.subplot(2,1,2)
    ax.set_xticks(range(0,int(xmax),100000),minor=True)
    ax.set_yticks(range(0,int(ymax*1.2),1000),minor=True)
    ax.grid(which='minor')
    pylab.scatter([rc_orig[ind] for ind in inds],[len([d for d in rc_sort if d.get(ind,0) >= rc_hi]) for ind in inds],facecolor='none',edgecolor='y')
    pylab.scatter([rc_orig[ind] for ind in inds],[len([d for this_cd,d in zip(cd_sort,rc_sort) if d.get(ind,0) >= rc_hi and this_cd <= cd_cut]) for ind in inds],facecolor='none',edgecolor='b')
    pylab.scatter([sum(rc_orig.values())],[sumy],facecolor='k')
    pylab.plot(range(0,allreads,simstep), [len(numpy.nonzero((sumrc1*c)>=rc_hi)[0]) for c in range(0,allreads,simstep)],'--k')
    pylab.plot(range(0,allreads,simstep), [len(numpy.nonzero((sumrc1*c)>=rc_hi)[0])*0.75 for c in range(0,allreads,simstep)],'--k')
    pylab.plot(range(0,allreads,simstep), [len(numpy.nonzero((sumrc1*c)>=rc_hi)[0])*0.5 for c in range(0,allreads,simstep)],'--k')
    pylab.ylim(0,ymax*1.2)
    pylab.xlim(0,xmax)

    if filename is not None:
        print >> sys.stderr, 'store scatters'
        try:
            pylab.savefig(filename)
        except IOError:
            print >> sys.stderr, 'unable to write %s, output not stored' % filename


def main(uniqueds, cluni, set_cd_cut=None, rc_low=4, rc_hi=10, simstep=10000, win=1000, xmax=None, startfig=0, figext='.pdf'):
    '''generate three plots (see module source)
    set_cd_cut sets cluster dirt if .clstats.cdest not found (see preprocess_radtag_lane "-e" option)
    '''

    rc_sort,cd_sort,inds,rc_orig,cd_cut = load_sorted_data(uniqueds,cluni)

    if cd_cut is None and set_cd_cut is not None:
        cd_cut = set_cd_cut

    try:
        draw_efficiency_plots(rc_sort,cd_sort,inds,rc_orig,cd_cut,fignum=startfig+1,filename='%s.efficiency%s' % (cluni,figext))
    except:
        print >> sys.stderr, 'efficiency plot draw failed'
    try:
        draw_ind_by_clust_plots(rc_sort,cd_sort,inds,rc_orig,cd_cut,rc_low=rc_low,rc_hi=rc_hi,win=win,fignum=startfig+2,filename='%s.indByClust%s' % (cluni,figext))
    except:
        print >> sys.stderr, 'ind_by_clust plot draw failed'
    try:
        draw_clust_by_reads_scatter(rc_sort,cd_sort,inds,rc_orig,cd_cut,rc_low=rc_low,rc_hi=rc_hi,simstep=simstep,xmax=xmax,fignum=startfig+3,filename='%s.clustByReads%s' % (cluni,figext))
    except:
        print >> sys.stderr, 'clust_by_reads plot draw failed'
    

if __name__ == '__main__':

    import argparse

    ds =  ' [%(default)s]'
    #create command line parser
    parser = argparse.ArgumentParser(description='plot run metrics for ddRAD de novo analysis')

    parser.add_argument('-lt','--low_threshold',default=4,type=int,help='"low" read coverage condition for cluster inclusion'+ds)
    parser.add_argument('-ht','--high_threshold',default=10,type=int,help='"high" read coverage condition for cluster inclusion'+ds)

    parser.add_argument('-s','--simulation_step',default=10000,type=int,help='x-axis values to interpolate predicted coverage'+ds)
    parser.add_argument('-x','--scatter_xmax',default=None,help='x-axis maximum value fof scatter'+ds)

    parser.add_argument('-w','--window',default=1000,type=int,help='cluster bin size for boxplot analyses'+ds)

    parser.add_argument('cluni',help='.cluni clustering result')
    parser.add_argument('uniqueds',nargs='+',help='.uniqued files incorporated in clustering analysis')

    opts = parser.parse_args()

    main(opts.uniqueds, opts.cluni, set_cd_cut=None, rc_low=opts.low_threshold, rc_hi=opts.high_threshold, simstep=opts.simulation_step, win=opts.window, xmax=opts.scatter_xmax, startfig=0, figext='.pdf')
    

    
