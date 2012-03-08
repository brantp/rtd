#!/usr/bin/env python

'''
given a SORTED file containing uniqued lines (see preprocess_radtag_lane.py) with cluster and node label prepended
computes multiple alignments across all cluster sequences and outputs SAM formatted alignments taking the most prevalent longest sequence as reference.
'''

import os, sys, re

import musclemap

from collections import defaultdict
from config import RTDROOT

def samline_from_alnpair(rname,raln,qname,qaln,qqual):
    if set(qqual) == set(['#']):
        return None
    
    leader,qseq = re.search('^(-*)(.*?)$',qaln).groups()

    pos = len(leader)+1

    cigar = []
    nm = 0
    md = []
    qi = 0
    for r,q in zip(raln[len(leader):].upper(),qseq.rstrip('-').upper()):
        if q != '-':
            qq = qqual[qi]
            qi += 1
        else:
            qq = None
        if qq == '#' or q == 'N' or r =='N':
            cigar.append('S')
        elif r in ['A','C','G','T'] and q in ['A','C','G','T']:
            cigar.append('M')
            if r != q:
                nm += 1
                md.append(r)
            else:
                md.append(1)
        elif r == '-' and q == '-':
            cigar.append('P')
        elif r == '-':
            cigar.append('I')
            nm += 1
        elif q == '-':
            cigar.append('D')
            nm += 1
            md.append('^'+r)

    #print ''.join(cigar)

    if 'S' in ''.join(cigar).strip('S'):
        return None

    #figure out cigar
    ccnt = 1
    cli = []
    cstate = None
    for c in cigar:
        if cstate == c:
            ccnt += 1
        else:
            if cstate is not None:
                cli.append('%d%s' % (ccnt,cstate))
            cstate = c
            ccnt = 1

       

    cli.append('%d%s' % (ccnt,cstate))
    cstr = ''.join(cli)


    #figure out md
    mdli = []
    mddel = []
    mdcnt = 0
    for c in md+['A']:
        if isinstance(c,int):
            mdcnt += c
            if len(mddel) > 0:
                mdli.append('^'+(''.join(mddel)))
                mddel = []
        else:
            if mdcnt:
                mdli.append(str(mdcnt))
                mdcnt = 0
            if c.startswith('^'):
                mddel.append(c[1:])
            else:
                if len(mddel) > 0:
                    mdli.append('^'+(''.join(mddel)))
                    mddel = []
                mdli.append(c)
           

    mdstr = ''.join(mdli[:-1])
    if mdstr == '':
        mdstr = '0'

    return '\t'.join([qname,'0',rname,str(pos),'30',cstr,'*','0','0',qaln.replace('-',''), qqual, 'NM:i:%s\tMD:Z:%s' % (nm,mdstr)])



def ref_seq_from_clust(clname,cl_aln):
    
    ref_seq = cl_aln[0][1].replace('-','')
    fa_str = '>%s\n%s\n' % (clname,ref_seq)

    return fa_str

def indiv_in_clust(cl_lines,rep_cut = 0):

    if isinstance(cl_lines[0],str):
        cl_lines = [l.strip().split() for l in cl_lines]
        
    ind_cts = defaultdict(int)
    for l in cl_lines:
        for ind,ct in zip( l[5].split(','), [int(i) for i in l[6].split(',')] ):
            if ct >= rep_cut:
                ind_cts[ind] += ct

    return ind_cts
    

def aln_from_clust(clname,cl_lines,keep_seqs=None,seq_len=0):

    if isinstance(cl_lines[0],str):
        cl_lines = [l.strip().split() for l in cl_lines]

    
    if keep_seqs is not None and len(cl_lines) > keep_seqs:
        orig_ind_ct = indiv_in_clust(cl_lines)
        orig_ind = len(indiv_in_clust(cl_lines))
        orig_len = len(cl_lines)
        cl_lines.sort(key = lambda l: (len(l[5].split(',')),sum([int(i) for i in l[6].split(',')]), len(l[2])),reverse=True)
        cl_lines = cl_lines[:keep_seqs]
        now_ind = len(indiv_in_clust(cl_lines))
        now_len = len(cl_lines)
        drop_indiv = set(orig_ind_ct.keys()) - set(indiv_in_clust(cl_lines).keys())
        #summarize!
        print >> sys.stderr, '\tcluster %s abbreviated: orig %s lines, %s indiv now %s lines, %s indiv (dropped: %s)' % \
              (clname, orig_len, orig_ind, now_len, now_ind,[(ind,orig_ind_ct[ind]) for ind in drop_indiv])

    cl_seqs = [l[2] for l in cl_lines]
    cl_nodes = [l[1] for l in cl_lines]
    #20110919 qscore translation functionality moved to get_uniqued_lines_by_cluster.py
    cl_quals = [l[4] for l in cl_lines]

    if seq_len != 0: #truncate sequences
        cl_seqs  = [s[:seq_len] for s in cl_seqs]
        cl_quals = [s[:seq_len] for s in cl_quals]

    lastnode = None
    cl_node_ids = []
    for node in cl_nodes:
        if node != lastnode:
            ct = 0
            lastnode = node
        else:
            ct += 1
        cl_node_ids.append('%s.%03d' % (node,ct))

    cl_aln = sorted( zip( cl_node_ids, \
                          musclemap.muscle(cl_seqs,1), \
                          cl_quals, \
                          [zip( l[5].split(','), [int(i) for i in l[6].split(',')] ) for l in cl_lines] ) , \
                     key=lambda x: (len(x[1].replace('-','').replace('N','')),len(x[3]),len(x[2].replace('#',''))),reverse=True)

    return cl_aln


def write_sam_from_aln(clname,cl_aln,rg_dict,samheader_fh,sambody_fh,ref_fh):

    raln = cl_aln[0][1]

    #sbfh = open(samfile+'.body','w')
    #rofh = open(ref_fasta_file,'w')

    rseq = ref_seq_from_clust(clname,cl_aln)
    ref_fh.write(rseq)

    #headers (@SQ lines)
    headline = '@SQ\tSN:%s\tLN:%s\n' % (clname,len(cl_aln[0][2]))
    samheader_fh.write(headline)

    #body
    for qname,qaln,qqual,inds_cts in cl_aln:

        samline = samline_from_alnpair(clname,raln,qname,qaln,qqual)
        if samline is None: continue
        samfields = samline.split()
        rg_lane = qname.split('.')[1]
        #try:
        #    if any([len(el) != 2 for el in inds_cts]):
        #        print inds_cts
        #except:
        #    print cl_aln
        for ind,ct in inds_cts:
            rg = '%s_%s' % (ind,rg_lane)
            rg_dict[rg] = ind
            for i in range(ct):
                this_samline = '\t'.join([samfields[0]+'.%s.%04d' % (ind,i)] + samfields[1:])
                sambody_fh.write('%s\tRG:Z:%s\n' % (this_samline,rg))

def calc_cluster_dirt(cl_lines):

    cl_ind_ct = defaultdict(list)
    for l in cl_lines:
        f = l.split()
        for ind,ct in zip(f[5].split(','),f[6].split(',')):
            cl_ind_ct[(ind,f[1].split('.')[1])].append(int(ct))
    
    totct = sum([sum(v) for v in cl_ind_ct.values()])
    dirtct = sum([sum(sorted(v,reverse=True)[2:]) for v in cl_ind_ct.values()])
    ctdirt = dirtct/float(totct)

    return ctdirt

if __name__ == '__main__':

    import argparse

    ds =  ' [%(default)s]'
    #create command line parser
    parser = argparse.ArgumentParser(description='generates SAM/BAM by multiple alignment within graph clusters')

    parser.add_argument('-d','--clust_dirt_max',default=0.10,type=float,help='cluster "dirt" threshold for processing (see documentation)'+ds)
    parser.add_argument('-i','--min_indiv',default=20,type=int,help='minimum number of individuals with at least one sequence in a cluster to include cluster'+ds)
    parser.add_argument('-k','--keep_seqs',default=100,type=int,help='only retain this many sequences for processing'+ds)
    parser.add_argument('-l','--seq_len',default=0,type=int,help='arbitrarily truncate sequences in SAM/BAM output at this length if not 0'+ds)

    parser.add_argument('-cs','--calc_only',action='store_true',help='calculate cluster statistics at supplied thresholds; do not generate alignments'+ds)
    
    parser.add_argument('cluniq',help='sorted .cluniq file containing cluster-associated unique sequences')
    parser.add_argument('fbase',help='basename for output files')

    opts = parser.parse_args()

    cluniq = opts.cluniq
    fbase = opts.fbase
    clust_dirt_max = opts.clust_dirt_max
    min_indiv = opts.min_indiv
    keep_seqs = opts.keep_seqs
    seq_len = opts.seq_len
    
    fdir = os.path.dirname(fbase)

    try:
        os.makedirs(fdir)
    except:
        pass

    fh = open(cluniq)
    if not opts.calc_only:
        samheader_fh = open(fbase+'.sam.header','w')
        sambody_fh = open(fbase+'.sam.body','w')
        ref_fh = open(fbase+'.fa','w')
    clstats_fh = open(fbase+'.clstats','w')

    rg_dict = {}

    this_cl = None
    cl_lines = []

    cl_on = 0
    for l in fh:
        if l.split()[0] != this_cl:
            if this_cl is not None:
                cl_dirt = calc_cluster_dirt(cl_lines)
                cl_indiv = len(indiv_in_clust(cl_lines))
                clstats_fh.write('%s\t%s\t%s\t%s\n' % (this_cl,len(cl_lines),cl_indiv,cl_dirt))
                if cl_on % 100 == 0: print >> sys.stderr, '%s\tcluster: %s\tunique seqs: %s\tindiv: %s\tdirt: %s' % (cl_on,this_cl,len(cl_lines),cl_indiv,cl_dirt)
                if not opts.calc_only and cl_dirt < clust_dirt_max and cl_indiv >= min_indiv: 
                    cl_aln = aln_from_clust(this_cl,cl_lines,keep_seqs,seq_len)
                    write_sam_from_aln(this_cl,cl_aln,rg_dict,samheader_fh,sambody_fh,ref_fh)

            cl_on += 1
            this_cl = l.split()[0]
            cl_lines = []
        cl_lines.append(l)

    clstats_fh.write('%s\t%s\t%s\t%s\n' % (this_cl,len(cl_lines),cl_indiv,cl_dirt))
    if not opts.calc_only:
        cl_aln = aln_from_clust(this_cl,cl_lines,keep_seqs)
        if calc_cluster_dirt(cl_lines) < clust_dirt_max and len(indiv_in_clust(cl_lines)) >= min_indiv:
            write_sam_from_aln(this_cl,cl_aln,rg_dict,samheader_fh,sambody_fh,ref_fh)

    clstats_fh.close()
    os.system(os.path.join(RTDROOT,'plot_error.py %s > %s' % (fbase+'.clstats',fbase+'.clstats.cdest' )))
    
    #finish headers (@RG lines)
    if not opts.calc_only:
        if len(rg_dict) == 0:
            print >> sys.stderr, 'readgroup dict is empty; no individuals included in final dataset. Check number of individuals and cluster dirt cutoffs and re-run'
            print >> sys.stderr, 'close output files ...',
            samheader_fh.close()
            sambody_fh.close()
            ref_fh.close()
            print >> sys.stderr, 'done.\nremove output files ...',
            os.unlink(samheader_fh.name)
            os.unlink(sambody_fh.name)
            os.unlink(ref_fh.name)
            print >> sys.stderr, 'done'    
            sys.exit(1)


        for rg in rg_dict:
            headline = '@RG\tID:%s\tPL:Illumina\tLB:%s\tSM:%s\n' % (rg,rg_dict[rg],rg_dict[rg])
            samheader_fh.write(headline)

        samheader_fh.close()
        sambody_fh.close()
        ref_fh.close()

        print >> sys.stderr, 'index reference'
        os.system('samtools faidx %s.fa' % (fbase))
        print >> sys.stderr, 'add headers and sort'
        os.system('cat %s.sam.header %s.sam.body | samtools view -bS - | samtools sort - %s' % (fbase,fbase,fbase))
        print >> sys.stderr, 'index bam'
        os.system('samtools index %s.bam' % (fbase))
        print >> sys.stderr, 'done'
