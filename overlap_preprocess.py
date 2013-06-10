#!/usr/bin/env python

import os,sys
import run_safe
import preprocess_radtag_lane

def idxstr_from_idx(x):
    return x and '_index%s' % x or ''

def get_idxseq(table='DB_multiplex_indices'):
    return dict([(d['idx'],d['seq']) for d in preprocess_radtag_lane.get_table_as_dict(table,suppress_fc_check=True)])

def get_adaptseq(table='DB_adapt_trim_seqs'):
    return dict([(d['adapterstype'],{'r1':d['r1'],'r2':d['r2']}) for d in preprocess_radtag_lane.get_table_as_dict(table,suppress_fc_check=True)])

def get_adapterstype(flowcell,lane,index,table='DB_library_data'):
    if index:
        td = preprocess_radtag_lane.get_table_as_dict(table,sq='flowcell="%s" and lane="%s" and index="%s"' % (flowcell,lane,index))
    else:
        td = preprocess_radtag_lane.get_table_as_dict(table,sq='flowcell="%s" and lane="%s"' % (flowcell,lane))
    adapterstypes = list(set([d.get('adapterstype','') for d in td]))
    if len(adapterstypes) == 1:
        return adapterstypes[0]
    else:
        errstr = 'invalid number of matches: %s' % adapterstypes
        raise ValueError, errstr

def qs_to_q(qs,baseQ):
    return [ord(c)-baseQ for c in qs]

def fq_splitext(fq):
    if fq.endswith('.gz'):
        parts = fq.rsplit('.',2)
        return parts[0], '.'.join(parts[1:])
    else:
        return fq.rsplit('.',1)

def convert_fastq(fq,ofq,out_lnum=4,out_baseQ=33,tickon = 10000):
    nreads = preprocess_radtag_lane.get_read_count(fq)
    lnum,baseQ = preprocess_radtag_lane.get_fastq_properties(fq)
    fh = preprocess_radtag_lane.smartopen(fq)
    ofh = preprocess_radtag_lane.smartopen(ofq,'w')
    for i in xrange(nreads):
        if i%tickon == 0:
            print >> sys.stderr, '\r%s / %s (%0.1f%%)' % (i,nreads,(float(i)/nreads)*100),
        n,s,qs = preprocess_radtag_lane.next_read_from_fh(fh, lnum)
        ofh.write(preprocess_radtag_lane.as_fq_line(n,s,qs_to_q(qs,baseQ),out_baseQ,out_lnum))
    print >> sys.stderr,'\n'
        
def save_previous_and_covert(prev_fq,fq,rm_ext=['.md5','.rc.cache'],**kwargs):
    print >> sys.stderr, 'move %s --> %s' % (fq,prev_fq)
    ret = os.system('mv %s %s' % (fq,prev_fq))
    if ret != 0:
        raise OSError, 'move failed'
    for ext in rm_ext:
        rmtarget = 'rm -f %s%s' % (fq,ext)
        ret = os.system(rmtarget)
        if ret == 0:
            print >> sys.stderr, 'removed %s' % rmtarget
        else:
            print >> sys.stderr, 'failed to remove %s' % rmtarget
    convert_fastq(prev_fq,fq,**kwargs)

def overlap_by_seqprep(r1,r2,outbase,pct_id=0.8,min_ol=10,adaptA='GATCGGAAGAGCACACG',adaptB='AGATCGGAAGAGCGTCGT'):
    '''adaptA is the adapter 1 sequence AS IT APPEARS IN READ 1,
    likewise adaptB is adapter 2 sequence AS IT APPEARS IN READ 2.
    in other words, DB_adapt_trim_seqs r1 (read 1 adapter read-through) is A, r2 is B
    '''
    trim1 = outbase+'.R1.trim.fastq.gz'
    trim2 = outbase+'.R2.trim.fastq.gz'
    drop1 = outbase+'.R1.drop.fastq.gz'
    drop2 = outbase+'.R2.drop.fastq.gz'
    merge = outbase+'.merge.fastq.gz'
    aln = outbase+'.merge.aln.gz'
    cmd = 'SeqPrep -f %s -r %s -1 %s -2 %s -3 %s -4 %s -A %s -B %s -s %s -E %s -o %s -m %s -n %s' % (r1,r2,trim1,trim2,drop1,drop2,adaptA,adaptB,merge,aln,min_ol,1-pct_id,pct_id)
    ss = run_safe.safe_script(cmd,outbase+'-seqprep',force_write=True)
    ret = os.system(ss)
    if ret != 0:
        raise OSError, 'seqprep run failed'
    return merge,trim1,trim2

if __name__ == "__main__":
    import argparse

    ds =  ' [%(default)s]'
    #command parser
    parser = argparse.ArgumentParser(description='runs PE read overlapping (and adapter trimming for either PE or SR) by SeqPrep, followed by preprocess')

    parser.add_argument('-fc','--flowcell',default=None,type=str,help='flowcell name REQUIRED'+ds)
    parser.add_argument('-l','--lane',default=None,type=str,help='lane REQUIRED'+ds)
    parser.add_argument('-idx','--index',default=None,type=str,help='multiplex index'+ds)
    
    parser.add_argument('-pp','--preprocess_argstr',default='',type=str,help='argument string passed to preprocess_radtag_lane.py (e.g. "-w -s CATG" for reads-by-individual and a CATG cutsite) DO NOT INCLUDE --lane --flowcell --index'+ds)
    
    parser.add_argument('-id','--percent_id',default=0.8,type=float,help='% identity required for overlapping and adapter trimming in SeqPrep'+ds)
    parser.add_argument('-ol','--overlap_length',default=10,type=int,help='minimum overlap required for merge/trim in SeqPrep'+ds)

    parser.add_argument('-sb','--seqprep_base',default=None,type=str,help='file basename for seqprep output (file will be in same directory as infiles) if not supplied, will be Sample_lane{lane}_{index}'+ds)

    parser.add_argument('infiles',nargs='+',help='2 fastq files corresponding to reads from a single lane/index, and optionally read 2 sequences for that lane/index')

    opts = parser.parse_args()

    if opts.flowcell is None or opts.lane is None:
        raise ValueError, '--flowcell and --lane (and --index as appropriate) must be specified'

    if len(opts.infiles) != 2: #PE
        errstr = '2 input files must be specified; got %s ' % len(opts.infiles)
        raise ValueError, errstr

    #check fq4-33
    for fq in opts.infiles:
        print >> sys.stderr, '\nfile: %s' % fq
        lnum,baseQ = preprocess_radtag_lane.get_fastq_properties(fq)
        print >> sys.stderr, 'lnum: %s\nbaseQ: %s' % (lnum,baseQ)
        if not (lnum == 4 and baseQ == 33):
            fqbase,fqext = fq_splitext(fq)
            prev_fq = '%s.fq%s-%s%s' % (fqbase,lnum,baseQ,fqext)
            print >> sys.stderr, 'must be 4-line, base 33 fastq to proceed; convert\nnew file will be %s\noriginal kept as %s\n' % (fq,prev_fq)
            save_previous_and_covert(prev_fq,fq)

    adapterstype = get_adapterstype(opts.flowcell,opts.lane,opts.index)
    adaptseq = get_adaptseq()
    adaptA,adaptB = adaptseq[adapterstype]['r1'],adaptseq[adapterstype]['r2']
    print >> sys.stderr, 'use adapterstype: %s\nadaptA: %s\nadaptB: %s' % (adapterstype,adaptA,adaptB)

    #run seqprep
    if opts.seqprep_base:
        sp_base = opts.seqprep_base
    else:
        sp_base = 'Sample_lane%s_%s' % (opts.lane, opts.index and opts.index or 'noidx')
    sp_fullbase = os.path.join(os.path.dirname(opts.infiles[0]),sp_base)

    merge,trim1,trim2 = overlap_by_seqprep(opts.infiles[0],opts.infiles[1],sp_fullbase,pct_id=opts.percent_id,min_ol=opts.overlap_length,adaptA=adaptA,adaptB=adaptB)

    cmd = 'preprocess_radtag_lane.py -iq 33 -suf merge %s -fc %s -l %s %s %s' % (opts.preprocess_argstr,opts.flowcell,opts.lane,(opts.index and '-idx %s' % opts.index or ''),merge)
    print >> sys.stderr, cmd
    ss = run_safe.safe_script(cmd,sp_fullbase+'-preprocess_merge',force_write=True)
    ret = os.system(ss)
    if ret != 0 or not os.path.exists(sp_fullbase+'-preprocess_merge.done'):
        raise OSError, 'merge preprocess failed'
    cmd = 'preprocess_radtag_lane.py -iq 33 -suf trim %s -fc %s -l %s %s %s %s' % (opts.preprocess_argstr,opts.flowcell,opts.lane,(opts.index and '-idx %s' % opts.index or ''),trim1,trim2)
    print >> sys.stderr, cmd
    ss = run_safe.safe_script(cmd,sp_fullbase+'-preprocess_trim',force_write=True)
    ret = os.system(ss)
    if ret != 0 or not os.path.exists(sp_fullbase+'-preprocess_merge.done'):
        raise OSError, 'trim preprocess failed'
