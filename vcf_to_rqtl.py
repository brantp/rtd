#!/usr/bin/env python

from short_read_analysis import variant_detection,extract_genotypes_from_mclgr
import os,sys

if __name__ == "__main__":

    parent_str = 'BW,PO'
    #qd = 6
    #gq = 20
    min_indiv = 50
    fh = 0.7
    site_before = 32 #polymorphism must occur before this base in a fragment
    #chi2crit = 30
    
    vcfn,qd,gq,chi2crit = sys.argv[1:]
    
    
    outbase = os.path.splitext(vcfn)[0]

    cut_fn = lambda sd: sd.has_key('QD') and float(sd['QD']) >= float(qd) and len(sd['indiv_gt']) >= min_indiv and sd['fh'] < fh


    print >> sys.stderr, 'loading vcf',vcfn
    vcf = variant_detection.load_vcf(vcfn,cutoff_fn=cut_fn,indiv_gt_phred_cut=float(gq))

    print >> sys.stderr, 'convert to pm/gt matrices'
    pm,gt = extract_genotypes_from_mclgr.genotypes_from_vcf_obj(vcf)

    parents_prefixes = dict(zip(['A', 'B'],parent_str.split(',')))
    parents = dict([(l,[k for k in gt.keys() if k.startswith(p)]) for l,p in parents_prefixes.items()])

    polarized_loci,polarized_geno = extract_genotypes_from_mclgr.genotypes_by_parent(dict([(k,v) for k,v in pm.items() if int(k.split('.')[1]) < site_before]),gt,parents,remove_targets=reduce(lambda x,y: x+y,parents.values()))

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
