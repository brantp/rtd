#!/usr/bin/env python
'''invocation:
vcf_to_rqtl.py path_to_vcf.vcf "P" Q G
where:
P = comma-separated pair of prefixes identifying cross founder species/strains
    for instance, if founding parents of strain "BW" are BW01, BW02, BW03
    and founding parents of strain "PO" are PO01,PO02,PO03,
    this paramter would be "BW,PO"
Q = GATK "QD" score threshold for inclusion
G = individual genotype call quality threshold

'''
# functionality from short_read now copied here
#from short_read_analysis import variant_detection,extract_genotypes_from_mclgr

import os,sys,re,numpy
from collections import defaultdict
import HTSeq,random,sys
from collections import defaultdict
from rtd import preprocess_radtag_lane

def sample_data_from_DB(sampleids, mousedb = 'Hoekstra lab mouse database'):
    td = preprocess_radtag_lane.no_net_get_table_as_dict(mousedb)
    ped = dict([ (d['id'], (d['damid'],d['sireid'])) \
                 for d in td \
                 if d.get('id','') in sampleids \
                 and d.has_key('sireid') \
                 and d.has_key('damid')])
    ped_parents = reduce(lambda x,y:x+y,ped.values())
    ped.update(dict([ (d['id'],(d['damid'],d['sireid'])) \
                      for d in td \
                      if d.get('id','') in ped_parents \
                      and d.has_key('sireid') \
                      and d.has_key('damid')]))
    
    recombinants = [d['id'] for d in td if d.get('id','') in ped.keys() \
                    and ',' in d['damstrain'] and ',' in d['sirestrain']]
    
    parents = []
    for f2 in recombinants:
        for f1 in ped[f2]:
            for g0 in ped[f1]:
                parents.append(g0)
    parents = list(set(parents))
    
    parents_spp = dict([(d['id'],d['damstrain']) for d in td if d.get('id','') in parents])
    parents_spp
    
    return ped, recombinants, parents, parents_spp

def species_tests_by_family(ped, recombinants, parents_spp):
    tests = defaultdict(list)
    for f2 in recombinants:
        g0_by_spp = defaultdict(list)
        for f1 in ped[f2]:
            for g0 in ped[f1]:
                g0_by_spp[parents_spp[g0]].append(g0)
        k = tuple([tuple(set(g0_by_spp[spp])) for spp in sorted(g0_by_spp.keys())])
        tests[k].append(f2)

    return tests

def default_cut_fn(vc):
    #consider adding ReadPosRankSum > -9
    accept = len(vc.alt) == 1 \
             and vc.info.get('QD',0) >= 5 

    return accept

def cross_genotypes_from_htseq_vcf(vcfr, tests, cut_fn=default_cut_fn, gq_cut=20):
    genotypes = defaultdict(dict)
    loci = []

    for i,vc in enumerate(vcfr):
        #THIS LIKELY FAILS IF make_info_dict() hasn't been run...
        vc.unpack_info(vcfr.infodict)
        if i % 10000 == 0: print >> sys.stderr, '\r',i,len(loci),'found',
        if not cut_fn(vc): continue
        loc = '%s.%s' % (vc.pos.chrom,vc.pos.start)
        for (pA,pB),inds in tests.items():
            pAa = set()
            pBa = set()
            if any([float(vc.samples[g0].get('GQ',0)) <  gq_cut for g0 in pA + pB]): continue
            for pAi in pA:
                alleles = vc.samples[pAi]['GT'].split('/')
                pAa = pAa.union(set(alleles))
            for pBi in pB:
                alleles = vc.samples[pBi]['GT'].split('/')
                pBa = pBa.union(set(alleles))
            if not '.' in pAa and not '.' in pBa and len(pAa) == 1 and len(pBa) == 1 and len(pAa.intersection(pBa)) == 0:
                gt_lookup = {list(pAa)[0]: 'A', list(pBa)[0]:'B'}
                for f2 in inds:
                    if float(vc.samples[f2].get('GQ',0)) < gq_cut: continue
                    gt = ''.join([gt_lookup[a] for a in vc.samples[f2]['GT'].split('/')])
                    genotypes[f2][loc] = gt
                if len(loci) and loci[-1] == loc:
                    pass
                else:
                    loci.append(loc)

    return loci, genotypes

def output_cross_radtag_genotypes(loci,genotypes,filename,lg0='X'):
    '''Given list loci and dictionary genotype per genotypes_by_parent, writes file <filename>
    suitable for RQTL

    overloads 20101202:
    - if loci is a dict per maploci from load_cross_radtag_genotypes below, sort by map position in output
    - if filename is not string, use as filehandle (permits passing sys.stdout, for instance)
    '''

    def sortkey(x):
        if x == '':
            return 0
        else:
            return x

    if isinstance(loci,list):
        locnames = loci
        lgs = ['1' for i in range(len(loci))]
        mps = [str(i+1) for i in range(len(loci))]
    elif isinstance(loci,dict):
        for k,v in loci.items():
            if v[0] == 0:
                loci[k] = (lg0,v[1])
        locnames,lgs,mps = zip(*[(loc,str(lg),str(mp)) for loc,(lg,mp) in sorted(loci.items(),key=lambda x:[sortkey(v) for v in x[1]])])

    mID_lookup = dict([(m,str(i)) for i,m in enumerate(sorted(genotypes.keys()))])

    if isinstance(filename,str):
        fh = open(filename ,'w')
        #open(filename+'.mIDlookup','w').write('\n'.join(['%s\t%s' % (i,m) for m,i in sorted(mID_lookup.items())]))
    else:
        fh = filename
        #open(filename.name+'.mIDlookup','w').write('\n'.join(['%s\t%s' % (i,m) for m,i in sorted(mID_lookup.items())]))
        
    fh.write('ID,')
    fh.write(','.join(['%sr' % l for l in locnames]))
    fh.write('\n')
    fh.write(',')
    fh.write(','.join(lgs))
    fh.write('\n')
    fh.write(',')
    fh.write(','.join(mps))
    fh.write('\n')

    out_geno = {}
    
    for mID in genotypes.keys():
        #fh.write(mID_lookup[mID]+',')
        fh.write(mID+',')
        out_geno[mID] = dict([(mkr,genotypes[mID][mkr]) for mkr in locnames if genotypes[mID].has_key(mkr)])
        fh.write(','.join([genotypes[mID].get(mkr,'-') for mkr in locnames]))
        fh.write('\n')

    fh.close()

    return out_geno,mID_lookup

if __name__ == "__main__":

    #vcfn,qd,gq,chi2crit = sys.argv[1:]
    vcfn,outbase,gq,fract_max = sys.argv[1:5]
    gq = float(gq)
    fract_max = float(fract_max)
    #outbase = os.path.splitext(vcfn)[0]

    vcfr = HTSeq.VCF_Reader(vcfn)
    vcfr.parse_meta()
    vcfr.make_info_dict()

    ped, recombinants, parents, parents_spp = sample_data_from_DB(vcfr.sampleids)
    tests = species_tests_by_family(ped, recombinants, parents_spp)
    
    polarized_loci,polarized_geno = cross_genotypes_from_htseq_vcf(vcfr, tests, gq_cut=gq)
    loc_counts = dict([(loc,sum([polarized_geno[ind].has_key(loc) for ind in recombinants])) for loc in polarized_loci])
    mct = max(loc_counts.values())
    keep_sites = [k for k,v in loc_counts.items() if v > mct*fract_max]

    print >> sys.stderr, '\n\nfound %s sites with max count %s individuals. Require %s, keeping %s' % (len(polarized_loci), mct, mct*fract_max, len(keep_sites))
    #chi2-free output:
    ret = output_cross_radtag_genotypes(keep_sites,polarized_geno,'%s_GQ%s_%sind.csv' % (outbase,gq,fract_max))
    

    """ #ditch chi2
    print >> sys.stderr, 'filter X linked, chi2 critical %s' % chi2crit
    xsites,autsites = extract_genotypes_from_mclgr.filter_Xlinked_loci(polarized_loci, polarized_geno,float(chi2crit))
    print >> sys.stderr, '%s X linked, %s autosomal' % (len(xsites),len(autsites))

    print >> sys.stderr, 'write output'
    ret = extract_genotypes_from_mclgr.output_cross_radtag_genotypes(xsites,polarized_geno,'%s_QD%s-GQ%s_%sbp_Xchi%s.csv' % (outbase,qd,gq,site_before,chi2crit))
    ret = extract_genotypes_from_mclgr.output_cross_radtag_genotypes(autsites,polarized_geno,'%s_QD%s-GQ%s_%sbp_autchi%s.csv' % (outbase,qd,gq,site_before,chi2crit))
    print >> sys.stderr, 'wrote:'
    print >> sys.stderr, '%s_QD%s-GQ%s_%sbp_Xchi%s.csv' % (outbase,qd,gq,site_before,chi2crit)
    print >> sys.stderr, '%s_QD%s-GQ%s_%sbp_autchi%s.csv' % (outbase,qd,gq,site_before,chi2crit)
    print >> sys.stderr, 'done'
    """
