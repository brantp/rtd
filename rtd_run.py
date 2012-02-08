#!/usr/bin/env python

'''
pipeline script generates reference-sorted, indexed BAM from uniqued reads from radtag sequencing lanes.

To generate uniqued reads, see preprocess_radtag_lane.py

four accessory programs and three python libraries are used, listed below.
for parallel execution, GNU parallel is also HIGHLY recommended

REQUIREMENTS:
-        PATH must contain: blat mcl mcxload muscle samtools [parallel]
-  PYTHONPATH must contain: numpy gdata editdist

see (at the time of this writing, March 09 2011)
  blat         http://hgdownload.cse.ucsc.edu/downloads.html
  mcl/mcxload  http://www.micans.org/mcl/
  muscle       http://www.drive5.com/muscle/
  samtools     http://samtools.sourceforge.net/
  GNU parallel http://savannah.gnu.org/projects/parallel/

  numpy *      http://sourceforge.net/projects/numpy/files/
  gdata        http://code.google.com/p/gdata-python-client/downloads/list
  editdist     http://www.mindrot.org/projects/py-editdist/

* N.B. numpy is also available as part of the excellent Enthought Python Distribution,
available free for academic/non-profit use at http://www.enthought.com/products/epd.php
'''

from config import RTDROOT as radtag_denovo
from preprocess_radtag_lane import get_baseQ

from collections import defaultdict
from itertools import groupby
from operator import itemgetter

from subprocess import Popen, PIPE
from editdist import distance
from glob import glob
from copy import deepcopy

import os, sys, re, numpy
import gdata.spreadsheet.service

def load_uniqued(all_quality,uniqued,readlen=None,nticks=20,baseQ=None):
    '''given a .uniqued file produced by preprocess_radtag_lane.py

    loads data into all_quality, ensuring sequences remain unique

    all_quality per 20101114 - UPDATE below    
    '''

    print >> sys.stderr, '%s readcount: ' % (uniqued),
    #number of sequences
    nreads = int(Popen('wc -l %s' % uniqued,shell=True,stdout=PIPE).stdout.read().split()[0])
    print >> sys.stderr, nreads
    
    qfh = open(uniqued)
    while baseQ is None:
		line = qfh.next()
		qstr = line.strip().split()[2]
		baseQ = get_baseQ(qstr)
    qfh.close()
    
    print >> sys.stderr, 'uniqued qualities base %s' % (baseQ)

    
    tickon = nreads/nticks
    if tickon < 1:
    	tickon = 1
    print >> sys.stderr, '\tloading'


    for i,line in enumerate(open(uniqued)):
        if i % tickon == 0: print >> sys.stderr, '\t\t%s / %s (%d%%)' % (i,nreads,(float(i)/nreads)*100)

        try:
            s,c,qstr,indiv,indcnt,r2,r2cnt = line.strip().split()
        except ValueError:
            print >> sys.stderr, 'line %s split: incorrect element number (%s) line:\n%ssplit:\n%s\n' % (i,len(line.strip().split()),line,line.strip().split())
        q = numpy.array([ord(ch)-baseQ for ch in qstr])
        c = int(c)
        indiv = set(indiv.split(','))

        if readlen is not None:
            s = s[:readlen]
            q = q[:readlen]

        if all_quality.has_key(s):
            all_quality[s]['mIDs'] = list(set(all_quality[s]['mIDs']).union(indiv))
            all_quality[s]['sum_quality'] += q*c
            all_quality[s]['tot'] += c
        else:
            all_quality[s]['mIDs'] = list(indiv)
            all_quality[s]['sum_quality'] = q*c
            all_quality[s]['tot'] = c

def preprocess_sequence_for_match(all_quality, cutsite, mIDfile, subject, queries, minlen=20):
    '''given a quality dictionary
    {
    20101114 - UPDATE:
    modified dict structure:
    
    
    <sequence> : {
                 "tot" : int
                 "mIDs" : [  <sampleID> ,<sampleID>,  ]
                 "sum_quality" : array([int,int,int ...])
                 }
    }

    generates three types of files:
    1x mIDlookup file containing header\tmID\tmID ... for each sequence
    1x "subject" contains all sequences that start with <cutsite>
    Nx "query" each contain a partition (<nparts> total) of fasta formatted sequence.  All seqs greater than <minlen> included

    '''

    mID_fh = open(mIDfile,'w')
    subj_fh = open(subject,'w')

    print >> sys.stderr, 'write query sequences'

    gen_queries = []
    
    if len(queries) == 1: #write all queries to single file
        this_outfile = queries[0]
        this_query_fh = open(this_outfile,'w')
        print >> sys.stderr, this_outfile
    else: #write to multiple query files for parallel execution
        this_outfile = None
        break_at = int(len(all_quality)/(len(queries)))
        qcopy = deepcopy(queries)
    
    for i,s in enumerate(sorted(all_quality.keys())):
        c = all_quality[s]['tot']
        qsum = all_quality[s]['sum_quality']
        q = qsum / c

        if len(queries) > 1 and i%break_at==0 and len(qcopy) > 0: #move to the next query chunk
            if this_outfile:
                gen_queries.append(this_outfile)
                this_query_fh.close()
            this_outfile = qcopy.pop(0)
            print >> sys.stderr, i,this_outfile
            this_query_fh = open(this_outfile,'w')

        if 2 in q:
            first2 = numpy.arange(len(q))[q==2][0]
        else:
            first2 = len(q)

        if first2 > minlen:
            header = '%s.%s.%s.%s' % (i,c,s[:first2],''.join([chr(int(n)+64) for n in q[:first2]]))
            this_query_fh.write('>%s\n%s\n' % (header,s[:first2]))
            mID_fh.write(header+'\t'+('\t'.join(all_quality[s]['mIDs']))+'\n')
            if s.startswith(cutsite) and c > 1:
                subj_fh.write('>%s\n%s\n' % (header,s[:first2]))

    gen_queries.append(this_outfile)
    this_query_fh.close()
    gen_subject = subject
    subj_fh.close()
    return gen_subject, gen_queries

def get_shortest_readlen(unifiles):
    readlen = numpy.inf
    for uniqued in unifiles:
        rl = len(open(uniqued).readline().strip().split()[0])
        if rl < readlen:
            readlen = rl
    return readlen

def make_similarity_calc_inputs(unifiles,set_mincycles,nticks,cutsite,mIDfile,subject,queries):
    if set_mincycles:
        readlen = set_mincycles
        print >> sys.stderr, 'read length set manually:',readlen
    else:
        readlen = get_shortest_readlen(unifiles)
        print >> sys.stderr, 'read length selected:',readlen

    all_quality = defaultdict(dict)

    for i, uniqued in enumerate(unifiles):
        print >> sys.stderr, '%s of %s, load' % (i+1,len(unifiles))
        load_uniqued(all_quality,uniqued,readlen=readlen,nticks=nticks)
        print >> sys.stderr, 'Total number of unique sequences loaded:',len(all_quality)

    gen_subject,gen_queries = preprocess_sequence_for_match(all_quality,cutsite,mIDfile,subject,queries)
    if gen_subject == subject and set(gen_queries) == set(queries):
        print >> sys.stderr, 'similarity calculation inputs successfully created (subject %s and %s queries)' % (subject, len(queries))
    else:
        raise ValueError, 'failure to match subjects: given %s generated %s, queries: (given - generated) %s, (generated - given) %s' % (subject,gen_subject,set(queries)-set(gen_queries),set(gen_queries)-set(queries))

    #drop sequences to save on ram for graph processing
    del(all_quality)


def run_lsf_blat(subject,queries,blattile,blatargstr='',num_batches=100):
    '''submits mcl_id_triples_by_blat.py jobs to LSF

    intended as an example of parallelization over a compute grid;
    uses a module LSF.py for interaction with scheduler

    '''
    import LSF
    
    blatargstr += ' -tileSize=%s' % blattile
    blatargstr += ' -stepSize=%s' % (int(blattile)/2)

    cmds = []
    labf = []
    for q in queries:
        outbase = q.rstrip('.fa').rstrip('_query')+'_blat'+blatargstr.replace('=','').replace(' ','')
        labf.append(outbase+'.label')
        cmds.append('%smcl_id_triples_by_blat.py %s %s \\"%s\\" %s' % (radtag_denovo,subject,q,blatargstr,outbase))

    logfile = os.path.join(os.path.dirname(subject),'blat-log')
    try:
        os.unlink(logfile)
    except:
        pass
    #print >> sys.stderr, 'LSF %s\nlog: %s' % (cmds,logfile)
    import time
    while len(cmds) > 0:
        jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'normal_serial',jobname_base='blat2mat',num_batches=num_batches)
        time.sleep(20)
        LSF.lsf_wait_for_jobs(jobids,logfile,namedict=namedict)

        cmds = LSF.lsf_no_success_from_log(logfile)

    return labf

def run_local_blat(subject,queries,blattile,blatargstr='',num_cores=1):
    '''
    runs blat commands using os.system()
    runs all jobs as a single batch, to run on multiple cores/computers, consider run_parallel_blat()
    '''

    blatargstr += ' -tileSize=%s' % blattile
    blatargstr += ' -stepSize=%s' % (int(blattile)/2)

    cmds = []
    labf = []
    for q in queries:
        outbase = q.rstrip('.fa').rstrip('_query')+'_blat'+blatargstr.replace('=','').replace(' ','')
        labf.append(outbase+'.label')
        cmds.append('%smcl_id_triples_by_blat.py %s %s "%s" %s' % (radtag_denovo,subject,q,blatargstr,outbase))

    shscr = os.path.join(os.path.dirname(subject) , 'runblat.sh')
    open(shscr, 'w').writelines([cmd+';\n' for cmd in cmds])
    os.system('chmod +x '+shscr)
    os.system(shscr)

    return labf

def run_parallel_blat(subject,queries,blattile,blatargstr='',num_cores='+0'):
    '''
    runs blat commands using GUN parallel.

    '''

    blatargstr += ' -tileSize=%s' % blattile
    blatargstr += ' -stepSize=%s' % (int(blattile)/2)

    cmds = []
    labf = []
    for q in queries:
        outbase = q.rstrip('.fa').rstrip('_query')+'_blat'+blatargstr.replace('=','').replace(' ','')
        labf.append(outbase+'.label')
        cmds.append('%smcl_id_triples_by_blat.py %s %s "%s" %s' % (radtag_denovo,subject,q,blatargstr,outbase))

    shscr = os.path.join(os.path.dirname(subject) , 'runblat.sh')
    open(shscr, 'w').writelines([cmd+';\n' for cmd in cmds])
    os.system('chmod +x '+shscr)
    os.system('parallel --progress -j %s < %s' % (num_cores,shscr))

    return labf


def run_match(match_engine,parallel_engine,ncores,subject,queries,minlen,other_argstr):
    if match_engine == 'blat':
        if parallel_engine == 'local':
            return run_local_blat(subject,queries,minlen,other_argstr,ncores)
        elif parallel_engine == 'lsf':
            return run_lsf_blat(subject,queries,minlen,other_argstr,ncores)
        elif parallel_engine == 'parallel':
            return run_parallel_blat(subject,queries,minlen,other_argstr,ncores)
    
    
if __name__ == '__main__':

    import argparse

    ds =  ' [%(default)s]'
    #create command line parser
    parser = argparse.ArgumentParser(description='performs distance matrix calculation, graph clustering and SAM generation')

    parser.add_argument('-me','--match_engine',default='blat',choices=['blat',],help='specifies match engine to use for distance matrix calculation'+ds)
    parser.add_argument('-pe','--parallel_engine',default='local',choices=['local','parallel','lsf'],help='specifies parallelization engine to use for distance matrix calculation'+ds)
    parser.add_argument('-l','--minlen',default=8,type=int,help='minimum tile/string match to calculate alignment in distance matrix calculation'+ds)
    parser.add_argument('--blatargstr',default='-minScore=20 -minMatch=2 -repMatch=1000000',help='additional arguments passed to blat if match_engine == "blat"'+ds)
    parser.add_argument('-np','--nparts',default=40,type=int,help='number of match_engine runs to split distance matrix calculation into'+ds)
    parser.add_argument('-nc','--ncores',default=4,type=int,help='number of cores to run distance matrix sub-parts on'+ds)
     
    parser.add_argument('-I','--mclradius',default=2,type=float,help='radius term for mcl clustering (given as -I term to mcl; see mcl documention)'+ds)

    parser.add_argument('-c','--set_mincycles',default=0,type=int,help='arbitrarily truncate reads at set length, rather than shortest read in infiles'+ds)
    parser.add_argument('-ts','--truncate_seqs',action='store_true',help='only output final sequences to the length of cluster seed (as set by --set_mincycles or automatically determined from input data)'+ds)
    parser.add_argument('-s','--cutsite',default='AATTC',help='sequence left behind by restriction enzyme at read1 end of library NOT NECESSARILY FULL R.E. SITE'+ds)
    
    parser.add_argument('-cd','--clustdirt',default='0.05',help='fraction of invalid reads in a putatively orthologous sequence set to permit before removing sequence set. See sam_from_clust_uniqued.py'+ds)
    parser.add_argument('-mi','--minindiv',default='100',help='minimum number of individuals that must be represented in a putatively orthologous sequence set to include in SAM. See sam_from_clust_uniqued.py'+ds)
    parser.add_argument('-ks','--keepseqs',default='200',help='maximum number of unique sequences from a cluster to align for SAM. See sam_from_clust_uniqued.py'+ds)
    
    parser.add_argument('outroot',help='folder to which intermediate files and output will be written')
    parser.add_argument('infiles',nargs='+',help='any number of .uni "uniqued" files created by preprocess_radtag_lane.py')

    opts = parser.parse_args()

    match_engine = opts.match_engine

    ############
    # match_engine specific configuration. Edits for bfast, mummer, etc would go here
    
    # minlen serves as:
    #   blat:   -tileSize
    minlen = opts.minlen

    if match_engine == 'blat':
        blatargstr = opts.blatargstr
        append_argstr = blatargstr.replace('=','').replace(' ','')
    # other match_engine options should summarize run options in a single string (with no special characters)
    # as append_argstr in elif blocks here
    else:
        append_argstr = ''

    #
    ############


    set_mincycles = opts.set_mincycles #set to arbitrarily truncate reads at set length

    cutsite = opts.cutsite
    
    nticks = 20

    nparts = opts.nparts

    outroot = opts.outroot
    infiles = opts.infiles

    ############
    # set up output file names
    
    outprefix = os.path.join(outroot,'rtd')
    if set_mincycles:
        outprefix += '_%sbp' % set_mincycles

    subject = outprefix+'_subj.fa'
    mIDfile = outprefix+'.mIDs'
    queries = ['%s_%04dof%04d_query.fa' % (outprefix,i+1,nparts) for i in range(nparts)]

    sources_file = outprefix+'_sources.txt'

    outprefix += '_%s' % (match_engine + append_argstr)
    outprefix += '_l%s' % minlen

    labelfile = outprefix + '_all.label'
    
    matfile = outprefix + '.mat'
    tabfile = outprefix + '.tab'

    outprefix += '_l%s_mcl_I%0.1f' % (minlen,opts.mclradius)
    
    grfile = outprefix + '.gr'

    clunifile = outprefix + '.cluni'
    sambase = outprefix + '_%sdirt_%sindiv_%sseq' % (opts.clustdirt,opts.minindiv,opts.keepseqs)

    
    #
    ############


    ############
    # create top level output directory
    try:
        os.makedirs(outroot)
        print >>sys.stderr, 'created output directory %s' % outroot
    except OSError:
        print >> sys.stderr, 'output directory %s exists, using' % outroot
    #
    ############

    # skip loads if relevant files already exist
    if not os.path.exists(clunifile):
        if not os.path.exists(grfile):
            if not (os.path.exists(matfile) and os.path.exists(tabfile)):
                if not os.path.exists(labelfile):
                    if not (os.path.exists(mIDfile) and os.path.exists(subject) and all([os.path.exists(query) for query in queries])):
                        # make similarity calc inputs
                        make_similarity_calc_inputs(infiles,set_mincycles,nticks,cutsite,mIDfile,subject,queries)
                    else:
                        print >> sys.stderr, '%s, %s and %s inputs exist, using' % (outprefix+'.mIDs',subject,len(queries))
                    # make similarity triples
                    print >> sys.stderr, 'generate labels %s by %s' % (labelfile,match_engine)
                    labf = run_match(opts.match_engine,opts.parallel_engine,opts.ncores,subject,queries,minlen,blatargstr)
                    print >> sys.stderr, 'match complete, concatenate...',
                    os.system('cat %s > %s' % (' '.join(labf), labelfile))
                    print >> sys.stderr, 'done'
                else:
                    print >> sys.stderr, 'using labels %s' % labelfile
                # make matrix
                print >> sys.stderr, 'Write matrix for clustering'
                os.system('mcxload -abc %s -o %s -write-tab %s' % (labelfile,matfile,tabfile))
            else:
                print >> sys.stderr, '%s and %s present, using' % (matfile,tabfile)
            # make graph
            print >> sys.stderr, 'perform graph clustering'
            os.system('mcl %s -I %s -iliad094o %s' % (matfile,opts.mclradius,grfile))
        else:
            print >> sys.stderr, '%s present, using' % (grfile)
        # make cluni
        print >> sys.stderr, 'merge uniqued and clustering data'
        os.system('%sget_uniqued_lines_by_cluster.py %s %s %s | sort -n > %s' % (radtag_denovo,grfile,tabfile,' '.join(infiles),clunifile))
    else:
        print >> sys.stderr, '%s present, using' % (clunifile)
    # make samiliad094
    if opts.truncate_seqs:
        readlen = get_shortest_readlen(infiles)
        sambase += '_%sbp' % readlen
    else:
        readlen = 0
    cmd = '%ssam_from_clust_uniqued.py -d %s -i %s -k %s -liliad094 %s %s %s' % (radtag_denovo,opts.clustdirt,opts.minindiv,opts.keepseqs,readlen,clunifile,sambase)
    print cmd
    os.system(cmd)
iliad094iliad094iliad094
