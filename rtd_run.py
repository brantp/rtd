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
from preprocess_radtag_lane import get_baseQ,smartopen,get_read_count

from collections import defaultdict
from itertools import groupby
from operator import itemgetter

from subprocess import Popen, PIPE
from editdist import distance
from glob import glob
from copy import deepcopy

import os, sys, re, numpy, gzip
import gdata.spreadsheet.service

def cat(filelist,targetfile):
    '''cats an arbitrarily large filelist to targetfile'''
    fh = smartopen(targetfile,'w')
    print >> sys.stderr, '\n'
    for i,f in enumerate(filelist):
        print >> sys.stderr, '\r%s / %s' % (i,len(filelist)),
        for l in open(f):
            fh.write(l)
    fh.close()

def load_uniqued(all_quality,uniqued,readlen=None,nticks=20,baseQ=None):
    '''given a .uniqued file produced by preprocess_radtag_lane.py

    loads data into all_quality, ensuring sequences remain unique

    all_quality per 20101114 - UPDATE below    
    '''

    nreads = get_read_count(uniqued)
    
    qfh = smartopen(uniqued)
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


    for i,line in enumerate(smartopen(uniqued)):
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

def preprocess_sequence_for_match(all_quality, cutsite, mIDfile, subjects, queries, minlen=20):
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
    import random

    mID_fh = smartopen(mIDfile,'w')

    
    if len(subjects) == 1: #write all subjects to single file
        this_subj_outfile = subjects[0]
        this_subj_fh = smartopen(this_subj_outfile,'w')
        print >> sys.stderr, this_subj_outfile
    else: #write to multiple subject files for parallel execution
        this_subj_outfile = None
        subj_break_at = int(len(all_quality)/(len(subjects)))
        scopy = deepcopy(subjects)

    print >> sys.stderr, 'write sequences'

    gen_queries = []
    gen_subjects = []
    
    if len(queries) == 1: #write all queries to single file
        this_outfile = queries[0]
        this_query_fh = smartopen(this_outfile,'w')
        print >> sys.stderr, this_outfile
    else: #write to multiple query files for parallel execution
        this_outfile = None
        break_at = int(len(all_quality)/(len(queries)))
        qcopy = deepcopy(queries)

    aqkeys = all_quality.keys()
    random.shuffle(aqkeys)
    for i,s in enumerate(aqkeys):
        c = all_quality[s]['tot']
        qsum = all_quality[s]['sum_quality']
        q = qsum / c

        if len(queries) > 1 and i%break_at==0 and len(qcopy) > 0: #move to the next query chunk
            if this_outfile:
                gen_queries.append(this_outfile)
                this_query_fh.close()
            this_outfile = qcopy.pop(0)
            print >> sys.stderr, i,this_outfile
            this_query_fh = smartopen(this_outfile,'w')

        if len(subjects) > 1 and i%subj_break_at==0 and len(scopy) > 0: #move to the next query chunk
            if this_subj_outfile:
                gen_subjects.append(this_subj_outfile)
                this_subj_fh.close()
            this_subj_outfile = scopy.pop(0)
            print >> sys.stderr, i,this_subj_outfile
            this_subj_fh = smartopen(this_subj_outfile,'w')

        if 2 in q:
            first2 = numpy.arange(len(q))[q==2][0]
        else:
            first2 = len(q)

        if first2 > minlen:
            header = '%s.%s.%s.%s' % (i,c,s[:first2],''.join([chr(int(n)+64) for n in q[:first2]]))
            this_query_fh.write('>%s\n%s\n' % (header,s[:first2]))
            mID_fh.write(header+'\t'+('\t'.join(all_quality[s]['mIDs']))+'\n')
            if s.startswith(cutsite) and c > 1:
                this_subj_fh.write('>%s\n%s\n' % (header,s[:first2]))

    gen_queries.append(this_outfile)
    this_query_fh.close()
    gen_subjects.append(this_subj_outfile)
    this_subj_fh.close()
    return gen_subjects, gen_queries

def get_shortest_readlen(unifiles):
    readlen = numpy.inf
    for uniqued in unifiles:
        rl = len(smartopen(uniqued).readline().strip().split()[0])
        if rl < readlen:
            readlen = rl
    return readlen

def make_similarity_calc_inputs(unifiles,set_mincycles,nticks,cutsite,mIDfile,subjects,queries):
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

    gen_subjects,gen_queries = preprocess_sequence_for_match(all_quality,cutsite,mIDfile,subjects,queries)
    if set(gen_subjects) == set(subjects) and set(gen_queries) == set(queries):
        print >> sys.stderr, 'similarity calculation inputs successfully created (%s subjects and %s queries)' % (len(subjects), len(queries))
    else:
        raise ValueError, 'failure to match subjects: given %s generated %s, queries: (given - generated) %s, (generated - given) %s' % (set(subjects),set(gen_subjects),set(queries)-set(gen_queries),set(gen_queries)-set(queries))

    #drop sequences to save on ram for graph processing
    del(all_quality)


def run_lsf_blat(subjects,queries,blattile,blatargstr='',num_batches=100,queue='normal_serial'):
    '''submits mcl_id_triples_by_blat.py jobs to LSF

    intended as an example of parallelization over a compute grid;
    uses a module LSF.py for interaction with scheduler

    '''
    import LSF,run_safe
    
    blatargstr += ' -tileSize=%s' % blattile
    blatargstr += ' -stepSize=%s' % (int(blattile)/2)

    #cmds = []
    labf = []
    to_run_dict = {}
    for q in queries:
        for subject in subjects:
            subjname = os.path.basename(subject).rstrip('.fa').rstrip('_subj')
            outbase = q.rstrip('.fa').rstrip('_query')+'_blat'+'-subj'+subjname+blatargstr.replace('=','').replace(' ','')
            labf.append(outbase+'.label.gz')
            # ESCAPES UNNECESSARY WITH safe_script  
            #cmds.append('%smcl_id_triples_by_blat.py %s %s \\"%s\\" %s' % (radtag_denovo,subject,q,blatargstr,outbase))
            cmd = '%smcl_id_triples_by_blat.py %s %s "%s" %s' % (radtag_denovo,subject,q,blatargstr,outbase)
            to_run_dict[outbase] = run_safe.safe_script(cmd,outbase)


    logfile = os.path.join(os.path.dirname(subjects[0]),'blat-log/blat-log')
    LSF.lsf_run_until_done(to_run_dict, logfile, queue, '-R "select[mem > 20000]"', 'blat2mat', num_batches, 3)

    # REPLACED BY lsf_run_until_done ABOVE
    #logfiles = glob(logfile+'*.lsflog')
    #for lf in logfiles:
    #    try:
    #        os.unlink(lf)
    #    except:
    #        pass
    #print >> sys.stderr, 'LSF %s\nlog: %s' % (cmds,logfile)
    #import time
    #while len(cmds) > 0:
    #    jobids,namedict = LSF.lsf_jobs_submit(cmds,logfile,'normal_serial',bsub_flags='-R "select[mem > 20000]"',jobname_base='blat2mat',num_batches=num_batches)
    #    time.sleep(20)
    #    LSF.lsf_wait_for_jobs(jobids,logfile,namedict=namedict)
    #    logfiles = glob(logfile+'*.lsflog')
    #    cmds = reduce(lambda x,y:x+y, [LSF.lsf_no_success_from_log(lf) for lf in logfiles])

    if not all([os.path.exists(f) for f in labf]):
        raise OSError, 'blat failed'

    return labf

def run_local_blat(subjects,queries,blattile,blatargstr='',num_cores=1):
    '''
    runs blat commands using os.system()
    runs all jobs as a single batch, to run on multiple cores/computers, consider run_parallel_blat()
    '''

    blatargstr += ' -tileSize=%s' % blattile
    blatargstr += ' -stepSize=%s' % (int(blattile)/2)

    cmds = []
    labf = []
    for q in queries:
        for subject in subjects:
            subjname = os.path.basename(subject).rstrip('.fa').rstrip('_subj')
            outbase = q.rstrip('.fa').rstrip('_query')+'_blat'+'-subj'+subjname+blatargstr.replace('=','').replace(' ','')
            labf.append(outbase+'.label.gz')
            cmd = '%smcl_id_triples_by_blat.py %s %s "%s" %s' % (radtag_denovo,subject,q,blatargstr,outbase)
            cmds.append(run_safe.safe_script(cmd,outbase))

    shscr = os.path.join(os.path.dirname(subjects[0]) , 'runblat.sh')
    smartopen(shscr, 'w').writelines([cmd+';\n' for cmd in cmds])
    os.system('chmod +x '+shscr)
    ret = os.system(shscr)
    if ret != 0 or not all([os.path.exists(f) for f in labf]):
        raise OSError, 'blat failed with code %s' % ret
    return labf

def run_parallel_blat(subjects,queries,blattile,blatargstr='',num_cores='+0'):
    '''
    runs blat commands using GUN parallel.

    '''

    blatargstr += ' -tileSize=%s' % blattile
    blatargstr += ' -stepSize=%s' % (int(blattile)/2)

    cmds = []
    labf = []
    for q in queries:
        for subject in subjects:
            subjname = os.path.basename(subject).rstrip('.fa').rstrip('_subj')
            outbase = q.rstrip('.fa').rstrip('_query')+'_blat'+'-subj'+subjname+blatargstr.replace('=','').replace(' ','')
            labf.append(outbase+'.label.gz')
            cmd = '%smcl_id_triples_by_blat.py %s %s "%s" %s' % (radtag_denovo,subject,q,blatargstr,outbase)
            cmds.append(run_safe.safe_script(cmd,outbase))

    shscr = os.path.join(os.path.dirname(subjects[0]) , 'runblat.sh')
    smartopen(shscr, 'w').writelines([cmd+';\n' for cmd in cmds])
    os.system('chmod +x '+shscr)
    ret = os.system('parallel --progress -j %s < %s' % (num_cores,shscr))
    if ret != 0 or not all([os.path.exists(f) for f in labf]):
        raise OSError, 'blat failed with code %s' % ret

    return labf


def run_match(match_engine,parallel_engine,ncores,subjects,queries,minlen,other_argstr):
    if match_engine == 'blat':
        if parallel_engine == 'local':
            return run_local_blat(subjects,queries,minlen,other_argstr,ncores)
        elif parallel_engine == 'lsf':
            return run_lsf_blat(subjects,queries,minlen,other_argstr,ncores)
        elif parallel_engine == 'parallel':
            return run_parallel_blat(subjects,queries,minlen,other_argstr,ncores)

def readlen_from_uniqued(uniqued):
    return len(smartopen(uniqued).readline().strip().split()[0])

def get_uniqued_error(infiles,cdest_searchbase):
    from glob import glob
    print >> sys.stderr, '\nset cluster dirt threshold from per-lane error estimates'
    err_by_uni = {}
    for uniqued in infiles:
        rl = readlen_from_uniqued(uniqued)
        cdest_search = uniqued.rstrip('.gz')+'-rtd/'+cdest_searchbase
        cdests = glob(cdest_search)
        if len(cdests) != 1:
            raise ValueError, 'search string %s did not result in a single .cdest file %s' % (cdest_search,cdests)
        else:
            cd = float(smartopen(cdests[0]).read())
        print >> sys.stderr, '%s: found cluster dirt %s for read length %s. Estimated error: %s' % (uniqued,cd,rl,cd/rl)
        err_by_uni[uniqued] = cd/rl
    return err_by_uni

    
if __name__ == '__main__':

    import argparse

    ds =  ' [%(default)s]'
    #create command line parser
    parser = argparse.ArgumentParser(description='performs distance matrix calculation, graph clustering and SAM generation')

    parser.add_argument('-me','--match_engine',default='blat',choices=['blat',],help='specifies match engine to use for distance matrix calculation'+ds)
    parser.add_argument('-pe','--parallel_engine',default='local',choices=['local','parallel','lsf'],help='specifies parallelization engine to use for distance matrix calculation'+ds)
    parser.add_argument('-l','--minlen',default=8,type=int,help='minimum tile/string match to calculate alignment in distance matrix calculation'+ds)
    parser.add_argument('--blatargstr',default='-minScore=20 -minMatch=2 -repMatch=1000000',help='additional arguments passed to blat if match_engine == "blat"'+ds)
    parser.add_argument('-np','--nparts',default=40,type=int,help='number of match_engine queries to split distance matrix calculation into'+ds)
    parser.add_argument('-ns','--nsubj',default=1,type=int,help='number of match_engine subjects to split distance matrix calculation into'+ds)
    parser.add_argument('-nc','--ncores',default=4,type=int,help='number of cores to run distance matrix sub-parts on'+ds)
     
    parser.add_argument('-I','--mclradius',default=2,type=float,help='radius term for mcl clustering (given as -I term to mcl; see mcl documention)'+ds)
    parser.add_argument('-te','--mclthreads',default=1,type=float,help='expansion threads term for mcl clustering (given as -te term to mcl for values > 1; see mcl documention)'+ds)

    parser.add_argument('-c','--set_mincycles',default=0,type=int,help='arbitrarily truncate reads at set length, rather than shortest read in infiles'+ds)
    parser.add_argument('-ts','--truncate_seqs',action='store_true',help='only output final sequences to the length of cluster seed (as set by --set_mincycles or automatically determined from input data)'+ds)
    parser.add_argument('-s','--cutsite',default='AATTC',help='sequence left behind by restriction enzyme at read1 end of library NOT NECESSARILY FULL R.E. SITE'+ds)
    
    parser.add_argument('-cd','--clustdirt',default='0.05',help='threshold fraction of invalid reads in a putatively orthologous sequence set to permit before removing sequence set. See sam_from_clust_uniqued.py. Set to "None" to invoke dirt estimation from individual preprocess runs'+ds)
    parser.add_argument('-mi','--minindiv',default='100',help='minimum number of individuals that must be represented in a putatively orthologous sequence set to include in SAM. See sam_from_clust_uniqued.py'+ds)
    parser.add_argument('-ks','--keepseqs',default='200',help='maximum number of unique sequences from a cluster to align for SAM. See sam_from_clust_uniqued.py'+ds)
    parser.add_argument('-cs','--cluster_stats_only',action='store_true',help='calculate cluster statistics at supplied thresholds; do not generate alignments'+ds)
    parser.add_argument('-se','--skip_errors',action='store_true',help='ignore errors in the alignment step of bam file generation (skip the offending cluster entirely)'+ds)
    
    parser.add_argument('--cleanup',action='store_true',help='remove intermediate files as each step completes'+ds)
    
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
    nsubj = opts.nsubj

    outroot = opts.outroot
    infiles = opts.infiles

    ############
    # set up output file names
    
    outprefix = os.path.join(outroot,'rtd')
    if set_mincycles:
        outprefix += '_%sbp' % set_mincycles


    #subject = outprefix+'_subj.fa'
    subjects = ['%s_%04dof%04d_subj.fa' % (outprefix,i+1,nsubj) for i in range(nsubj)]
    mIDfile = outprefix+'.mIDs'
    queries = ['%s_%04dof%04d_query.fa' % (outprefix,i+1,nparts) for i in range(nparts)]

    sources_file = outprefix+'_sources.txt'

    outprefix += '_%s' % (match_engine + append_argstr)
    outprefix += '_l%s' % minlen

    labelfile = outprefix + '_all.label'
    labdonefile = labelfile + '.done'
    
    matfile = outprefix + '.mat'
    tabfile = outprefix + '.tab'
    matdonefile = matfile + '.done'

    outprefix += '_l%s_mcl_I%0.1f' % (minlen,opts.mclradius)
    
    grfile = outprefix + '.gr'

    cdest_searchbase = os.path.basename(outprefix) + '*.clstats.cdest'

    clunifile = outprefix + '.cluni'

    if opts.clustdirt == 'None': #set dirt cutoff based on estimated error rate and read lengths
        err_by_uni = get_uniqued_error(infiles,cdest_searchbase)
        mean_err = numpy.mean(err_by_uni.values())
        if set_mincycles:
            set_dirt = mean_err * int(set_mincycles)
            print >> sys.stderr, 'dirt set to mean of error (%s) * specified min_cycles (%s) = %s' % (mean_err, int(set_mincycles), set_dirt)
            opts.clustdirt = set_dirt
        else:
            set_dirt = mean_err * get_shortest_readlen(infiles)
            print >> sys.stderr, 'dirt set to mean of error (%s) * shortest read length (%s) = %s' % (mean_err, get_shortest_readlen(infiles), set_dirt)
            opts.clustdirt = set_dirt


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
            if not (os.path.exists(matfile) and os.path.exists(tabfile) and os.path.exists(matdonefile)):
                if not (os.path.exists(labelfile) and os.path.exists(labdonefile)): #back to concat; testing 20130103
                #if 1: # concatenated label file no longer generated; hardcoded for now
                    if not (os.path.exists(mIDfile) and all([os.path.exists(subj) for subj in subjects]) and all([os.path.exists(query) for query in queries])):
                        # make similarity calc inputs
                        make_similarity_calc_inputs(infiles,set_mincycles,nticks,cutsite,mIDfile,subjects,queries)
                    else:
                        print >> sys.stderr, '%s, %s subjects and %s queries inputs exist, using' % (outprefix+'.mIDs',len(subjects),len(queries))
                    # make similarity triples
                    print >> sys.stderr, 'generate labels %s by %s' % (labelfile,match_engine)
                    labf = run_match(opts.match_engine,opts.parallel_engine,opts.ncores,subjects,queries,minlen,blatargstr)

                    # BACK TO CONCAT LABELS # concatenate replaced by 
                    print >> sys.stderr, 'match complete, concatenate...',
                    cat(labf,labelfile)
                    
                    print >> sys.stderr, 'labels done'
                    os.system('touch %s' % labdonefile)
                else:
                    print >> sys.stderr, 'using labels %s' % labelfile
                # make matrix
                print >> sys.stderr, 'Write matrix for clustering'
                
                # BACK TO CONCAT LABELS # no longer using concatenated label file
                ret = os.system('mcxload -abc %s -o %s -write-tab %s' % (labelfile,matfile,tabfile))                

                # BACK TO CONCAT LABELS 20120103
                #mcxload_ps = Popen('mcxload -abc - -o %s -write-tab %s' % (matfile,tabfile),shell=True,bufsize=1,stdin=PIPE)
                #for lf in labf:
                #    lfh = gzip.open(lf)
                #    for l in lfh:
                #        mcxload_ps.stdin.write(l)
                #mcxload_ps.communicate()
                #ret = mcxload_ps.returncode
                # /NO LABELS
                
                if ret != 0:
                    raise OSError, 'mcxload failed with code %s' % ret
                
                ret = os.system('touch %s' % matdonefile)

            else:
                print >> sys.stderr, '%s and %s present, using' % (matfile,tabfile)
            # make graph
            if opts.mclthreads > 1:
                print >> sys.stderr, 'perform multithreaded (%s) graph clustering' % opts.mclthreads
                ret = os.system('mcl %s -I %s -te %s -o %s' % (matfile,opts.mclradius,opts.mclthreads,grfile))
                if ret != 0:
                    raise OSError, 'mcl failed with code %s' % ret

            else:
                print >> sys.stderr, 'perform graph clustering'
                ret = os.system('mcl %s -I %s -o %s' % (matfile,opts.mclradius,grfile))
                if ret != 0:
                    raise OSError, 'mcl failed with code %s' % ret
        else:
            print >> sys.stderr, '%s present, using' % (grfile)
        # make cluni
        print >> sys.stderr, 'merge uniqued and clustering data'
        ret = os.system('%sget_uniqued_lines_by_cluster.py %s %s %s | sort -n > %s' % (radtag_denovo,grfile,tabfile,' '.join(infiles),clunifile))
        if ret != 0:
            raise OSError, 'get_uniqued_lines_by_cluster.py failed with code %s' % ret
    else:
        print >> sys.stderr, '%s present, using' % (clunifile)

    #cleanup if invoked
    if opts.cleanup and os.path.exists(clunifile) and os.path.getsize(clunifile) > 0:
        print >> sys.stderr, 'file cleanup invoked; remove:'
        print >> sys.stderr, '\tsimilarity calculation mID',mIDfile
        if os.path.exists(mIDfile): os.unlink(mIDfile)
        print >> sys.stderr, '\tsimilarity calculation subjects [%s files]' % len(subjects)
        for subject in subjects:
            if os.path.exists(subject): os.unlink(subject)

        print >> sys.stderr, '\tsimilarity calculation queries [%s files]' % len(queries)
        for query in queries:
            if os.path.exists(query): os.unlink(query)
        try:
            print >> sys.stderr, '\tMCL input label file parts [%s files]' % len(labf)
            for lf in labf:
                if os.path.exists(lf): os.unlink(lf)
        except:
            print >> sys.stderr, '\tMCL input label file parts not defined; skip'
        print >> sys.stderr, '\tMCL input label file',labelfile
        if os.path.exists(labelfile): os.unlink(labelfile)
        print >> sys.stderr, '\tMCL input matrix file',matfile
        if os.path.exists(matfile): os.unlink(matfile)
        print >> sys.stderr, '\tMCL input table file',tabfile
        if os.path.exists(tabfile): os.unlink(tabfile)
        print >> sys.stderr, '\tMCL result graph file',grfile
        if os.path.exists(grfile): os.unlink(grfile)
        print >> sys.stderr, 'done'
        
    
    # make sam
    if opts.truncate_seqs:
        readlen = get_shortest_readlen(infiles)
        sambase += '_%sbp' % readlen
    else:
        readlen = 0
    if opts.cluster_stats_only:
        cmd = '%ssam_from_clust_uniqued.py -d %s -i %s -k %s -l %s %s -cs %s %s' % (radtag_denovo,opts.clustdirt,opts.minindiv,opts.keepseqs,readlen,(opts.skip_errors and '-s' or ''),clunifile,sambase)
    else:
        cmd = '%ssam_from_clust_uniqued.py -d %s -i %s -k %s -l %s %s %s %s' % (radtag_denovo,opts.clustdirt,opts.minindiv,opts.keepseqs,readlen,(opts.skip_errors and '-s' or ''),clunifile,sambase)
    print cmd
    os.system(cmd)
