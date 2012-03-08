#!/usr/bin/env python

if __name__ == "__main__":
    import matplotlib
    matplotlib.use('Agg')

import os, sys, numpy, pylab, math
from preprocess_radtag_lane import dezip
from collections import defaultdict

def plot_cd(cd_dict,fignum=1,xtickspace = 5.0,fig_title = None):
    xmax = (math.ceil(max(cd_dict.keys())/float(xtickspace)) * xtickspace) + 1
    fig = pylab.figure(fignum)
    ax = pylab.subplot(2,1,1)
    pylab.boxplot([v for k,v in sorted(cd_dict.items())])
    x,y = dezip([(k,numpy.percentile(v,99)) for k,v in sorted(cd_dict.items())])
    pylab.plot(x,y,'r')
    ax.set_xticks(range(0,int(xmax),int(xtickspace)))
    ax.set_yticks(numpy.arange(0,1,0.1))
    pylab.xlabel('number of individuals in cluster')
    pylab.ylabel('cluster dirt')
    if fig_title:
        pylab.title(fig_title)

    ax = pylab.subplot(2,1,2)
    x,y = dezip([(k,len(v)) for k,v in sorted(cd_dict.items())])
    pylab.plot(x,y,'g')
    ax.set_xticks(range(0,int(xmax),int(xtickspace)))
    pylab.xlabel('number of individuals in cluster')
    pylab.ylabel('number of clusters')
    
    return fig

def cd_dict_from_clstats(clstats):
    fh = open(clstats)
    cdd = defaultdict(list)
    for l in fh:
        f = l.strip().split()
        cdd[int(f[2])].append(float(f[3]))

    return cdd

def cd_cut_by_pctile_range(cd_dict,pct_lo=50,pct_hi=75,ret_pct=99):
    ctv = []
    for k,v in cd_dict.items():
        ctv.extend([k]*len(v))

    iqc = []
    for k,v in cd_dict.items():
        if numpy.percentile(ctv,pct_lo) < k < numpy.percentile(ctv,pct_hi):
            iqc.extend(v)

    return numpy.percentile(iqc,ret_pct)

if __name__ == "__main__":
    clstats = sys.argv[1]

    cd_dict = cd_dict_from_clstats(clstats)
    fig = plot_cd(cd_dict, fig_title = os.path.split(os.path.split(clstats)[0])[-1])
    fig.savefig(clstats+'.pdf')
    print >> sys.stderr, 'plot saved as',(clstats+'.pdf')
    print cd_cut_by_pctile_range(cd_dict)



