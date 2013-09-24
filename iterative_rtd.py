#!/usr/bin/env python
'''
iterative mapping-assembly by rtd over fragment size
'''

import os,sys,re
import run_safe
from subprocess import Popen,PIPE
from preprocess_radtag_lane import smartopen,get_read_count
from rtd_run import load_uniqued
from collections import defaultdict
from glob import glob

def remove_ext(fname,psplit=False):
    if fname.endswith('.gz'):
        fname = fname[:-3]
    if psplit:
        return os.path.split(os.path.splitext(fname)[0])
    else:
        return os.path.splitext(fname)[0]

def make_bwa_ref_idx(reference, force_index=False):
    idx_suff = ['.sa',
                '.amb',
                '.ann',
                '.pac',
                '.bwt']
    refidxs = [reference+suff for suff in idx_suff]
    if all([os.path.exists(f) for f in refidxs]) and not force_index:
        print >> sys.stderr, 'all bwa reference index files exist (%s)' % refidxs
    else:
        ret = os.system('bwa index %s' % reference)
        if ret == 0 and all([os.path.exists(f) for f in refidxs]):
            print >> sys.stderr, 'all bwa reference index files created (%s)' % refidxs
        else:
            raise OSError, 'reference index creation failed'

def make_stampy_ref_idx(reference, force_index=False):
    idx_suff = ['.sthash','.stidx']
    refidxs = [reference+suff for suff in idx_suff]
    if all([os.path.exists(f) for f in refidxs]) and not force_index:
        print >> sys.stderr, 'all stampy reference index files exist (%s)' % refidxs
    else:
        ret = os.system('stampy.py --overwrite -G %s %s; stampy.py --overwrite -g %s -H %s' % tuple([reference]*4))
        if ret == 0 and all([os.path.exists(f) for f in refidxs]):
            print >> sys.stderr, 'all stampy reference index files created (%s)' % refidxs
        else:
            raise OSError, 'reference index creation failed'

def map_by_bwa(reads,reference,mapped,bwa_args='',make_index=True,force_index=False):
    if make_index: make_bwa_ref_idx(reference,force_index=force_index)
    if reads.endswith('.bam'):
        cmd = 'samtools bam2fq %s | bwa mem %s %s - | samtools view -bS - > %s.bam' % (reads,bwa_args, reference, mapped)
    else:
        cmd = 'bwa mem %s %s %s | samtools view -bS - > %s.bam' % (bwa_args, reference, reads, mapped)

    ss = run_safe.safe_script(cmd,mapped,force_write=True)
    ret = os.system(ss)
    if ret == 0 and os.path.exists(mapped+'.bam'):
        print >> sys.stderr, '%s.bam created' % mapped
    else:
        errstr = 'mapping %s to %s failed' % (reads, reference)
        raise OSError, errstr

    return mapped+'.bam'

def map_by_stampy(reads,reference,mapped,stampy_args='--maxbasequal=60 --bwamark',make_index=True,force_index=False):
    if make_index: make_stampy_ref_idx(reference,force_index=force_index)
    cmd = 'stampy.py --overwrite %s -h %s -g %s -M %s -o %s.sam; samtools view -bS %s.sam > %s.bam' % (stampy_args,reference,reference, reads, mapped,mapped,mapped)
    ss = run_safe.safe_script(cmd,mapped,force_write=True)
    ret = os.system(ss)
    if ret == 0 and os.path.exists(mapped+'.bam'):
        print >> sys.stderr, '%s.bam created' % mapped
    else:
        errstr = 'mapping %s to %s failed' % (reads, reference)
        raise OSError, errstr

    return mapped+'.bam'

def get_unmapped_reads(bam,readnames_only=True,suffix='.unmapped.bam'):
    '''if readnames_only is False:
    generate a bam file of unmapped reads (otherwise, only returns ids of unmapped reads)
    '''
    head_out = os.path.splitext(bam)[0]+'head.sam'
    unmapped_out = os.path.splitext(bam)[0]+'unmapped.sam'
    all_out = os.path.splitext(bam)[0]+suffix
    if readnames_only:
        readnames = []
        ph = Popen('samtools view -f 4 %s' % bam,shell=True,stdout=PIPE)
        for line in ph.stdout:
            readnames.append(line.strip().split()[0])
        return readnames
    #otherwise:
    #header
    print >> sys.stderr, 'get header'
    ret = os.system('samtools view -H %s > %s' % (bam,head_out))
    if ret != 0:
        raise OSError, 'header get failed'
    #body
    print >> sys.stderr, 'get body'
    ret = os.system('samtools view -f 4 %s > %s' % (bam,unmapped_out))
    if ret != 0:
        raise OSError, 'reads get failed'
    #cat; bam
    print >> sys.stderr, 'compress'
    ret = os.system('cat %s %s | samtools view -bS - > %s' % (head_out,unmapped_out,all_out))
    if ret != 0:
        raise OSError, 'cat/bam failed'

    return all_out

from copy import copy
def subtractive_map(reads,reference,bwa=True,stampy=True,make_index=True,force_index=False,readnames_only=True):
    reads_dir,reads_base = remove_ext(reads,psplit=True)
    ref_dir, ref_base = remove_ext(reference,psplit=True)
    mapped_base= reads_base+'-'+ref_base
    input_reads = copy(reads)
    if bwa:
        mapped_base += '-bwa'
        mapped = os.path.join(reads_dir,mapped_base)
        reads = map_by_bwa(reads,reference,mapped,make_index=make_index,force_index=force_index)

    if stampy:
        mapped_base += '-stampy'
        mapped = os.path.join(reads_dir,mapped_base)
        reads = map_by_stampy(reads,reference,mapped,make_index=make_index,force_index=force_index)

    return get_unmapped_reads(reads,readnames_only=readnames_only)

def uniqued_to_fastq(uniqued,id_prefix=''):
    if uniqued.endswith('gz'):
        len_uni = int(Popen('zcat %s | wc -l' % uniqued,shell=True,stdout=PIPE).stdout.read().strip())
    else:
        len_uni = int(Popen('cat %s | wc -l' % uniqued,shell=True,stdout=PIPE).stdout.read().strip())
    fh = smartopen(uniqued)
    outname = remove_ext(uniqued)+'-fromuni.fastq.gz'
    if os.path.exists(outname) and get_read_count(outname) == len_uni:
        print >> sys.stderr, 'output %s exists' % outname
        return outname
    ofh = smartopen(outname,'w')
    print >> sys.stderr, 'convert %s to fastq' % uniqued
    for i,l in enumerate(fh):
        fields = l.strip().split()
        fq_line = '@%s%s\n%s\n+\n%s\n' % (id_prefix,i,fields[0],fields[2])
        ofh.write(fq_line)
        if i % 1000 == 0: print >> sys.stderr, '\r\t%s done' % i,
    ofh.close()
    print >> sys.stderr, '%s done' % outname

    return outname
        
def write_uniqued_by_size(all_quality,outbase,baseQ=33):
    outdir = os.path.dirname(outbase)
    if not os.path.exists(outdir): os.makedirs(outdir)

    outfhs = {}
    ofbysize = {}
    for seq,aqd in all_quality.items():
        #s,c,qstr,indivstr,indcnt,r2,r2cnt
        ind_li,cnt_li = zip(*aqd['count_by_ind'].items())
        outl =  '\t'.join((seq, str(aqd['tot']), ''.join([chr(i+baseQ) for i in map(int,aqd['sum_quality']/float(aqd['tot']))]), ','.join(ind_li),','.join(map(str,cnt_li)),'.','.')) + '\n'
        outf = outbase+'-%s.uniqued.gz' % len(seq)
        if not outf in outfhs:
            outfhs[outf] = smartopen(outf,'w')
            ofbysize[len(seq)] = outf
        outfhs[outf].write(outl)

    for outf,ofh in outfhs.items():
        ofh.close()

    return ofbysize

def get_uniqued_by_size(outbase):
    ofbysize = {}
    unis = glob(outbase+'-*.uniqued.gz')
    for u in unis:
        size = re.search('-(\d+)\.uniqued\.gz',u).groups()[0]
        ofbysize[int(size)] = u
    return ofbysize

def filter_uniqued(uniqued,outfile,lines_to_write):
    ofh = smartopen(outfile,'w')
    for i,l in enumerate(smartopen(uniqued)):
        if i in lines_to_write:
            ofh.write(l)
    ofh.close()

def append_to_ref(target_ref,new_ref,id_prefix):
    nfh = smartopen(new_ref)
    tfh = smartopen(target_ref,'a')
    for l in nfh:
        if l.startswith('>'):
            newl = l.replace('>','>%s_' % id_prefix)
            tfh.write(newl)
        else:
            tfh.write(l)
    nfh.close()
    tfh.close()

def ref_len(ref):
    return int(Popen('grep ">" %s | wc -l' % ref, shell=True,stdout=PIPE).stdout.read())

if __name__ == "__main__":
    #options, someday
    contam_fa = sys.argv[1]
    outdir = sys.argv[2]
    uniqueds = sys.argv[3:]

    bysize_dir = os.path.join(outdir,'by_size/uni_len')
    bysize_done = os.path.join(outdir,'by_size.done')
    denovo_ref = os.path.join(outdir,'denovo.fa')

    if os.path.exists(denovo_ref):
        print >> sys.stderr, 'REMOVE REF: %s' % denovo_ref
        os.unlink(denovo_ref)

    if os.path.exists(bysize_done):
        ofbysize = get_uniqued_by_size(bysize_dir)
    else:
        all_quality = defaultdict(dict)
        
        for uniqued in uniqueds:
            load_uniqued(all_quality,uniqued,count_by_ind=True)
            
        print >> sys.stderr, 'LOAD COMPLETE. WRITE BY-SIZE.'
        ofbysize = write_uniqued_by_size(all_quality,bysize_dir)
        del all_quality
        ret = os.system('touch %s' % bysize_done) 
    
    sizes = sorted(ofbysize.keys(),reverse=True)

    for i in sizes:
        print >> sys.stderr, '\nSTART %s' % i
        uni = ofbysize[i]

        ufq = uniqued_to_fastq(uni)
        nreads = get_read_count(ufq)

        if os.path.exists(denovo_ref):
            dn_len = ref_len(denovo_ref)
            noncontam_ubam = subtractive_map(ufq,contam_fa,stampy=False,readnames_only=False)
            unmapped = subtractive_map(noncontam_ubam,denovo_ref,force_index=True)
        else:
            dn_len = 0
            unmapped = subtractive_map(ufq,contam_fa,stampy=False)

        print >> sys.stderr, '\nGET %s UNMAPPED' % len(unmapped)
        funi = os.path.splitext(uni)[0]+'.filtered.gz'
        filter_uniqued(uni,funi,map(int,unmapped))

        outdir = os.path.splitext(uni)[0]+'-rtd'
        print >> sys.stderr, '\nRTD'
        cmd = "rtd_run.py -pe parallel -np 8 -nc 8 -s CATG -cd 1 -mi 1 --cleanup -te 8 %s %s" % (outdir,funi)
        ret = os.system(cmd)
        if ret != 0:
            raise OSError,cmd
        refadd = glob(outdir+'/*.fa')[0]
        
        append_to_ref(denovo_ref,refadd,i)
        add_len = ref_len(refadd)
        print >> sys.stderr, '\tITERATION %s ADDS %s to DENOVO REF LEN %s\n' % (i,add_len,dn_len)
