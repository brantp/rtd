#!/usr/bin/env python

'''generates R/qtl .csv file from a .vcf based on a pre-existing map, assigning parental allele classes based on .vcf from which the source map was constructed
'''

import os,sys, Util
from collections import defaultdict

gq = 20

from short_read_analysis import variant_detection
from short_read_analysis import extract_genotypes_from_mclgr

def load_vcf(vcf,allele_map,indiv_gt_phred_cut=None,ding_on=100000,return_map=False):
	'''processes a vcf file, adding genotypes satisfying GQ cutoff indiv_gt_phred_cut to a returned cross genotype object
	sites corresponding to keys in allele_map are retained
	'''

	if return_map:
		new_map = defaultdict(dict)
	else:
		vcf_data = {}
	
	i = 0
	for line in open(vcf):
		if i % ding_on == 0: print >> sys.stderr, 'reading',i
		i += 1

		if line.startswith('#CHROM'):
			headers = line[1:].split()
			exp_elements = len(line.split())
			FORMAT = headers.index('FORMAT')
		elif line.startswith('#'):
			continue
		else:
			#extract site stats
			fields = line.split()
			if len(fields) != exp_elements:
				print >>sys.stderr, 'unexpected length, line %s (exp %s obs %s)' % (i,exp_elements,len(fields))
				continue

			#populate site metrics
			sd = dict(zip(headers[:FORMAT],fields[:FORMAT]))
			loc = '%s.%s' % (sd['CHROM'],sd['POS'])
			key = (sd['CHROM'],sd['POS'])

			if not loc in allele_map.keys(): #not interested; skip!
				continue
				

			#temp hack for multiallelic sites
			if ',' in sd['ALT']:
				print >> sys.stderr, '!MULTIALLELIC SITE AT %s' % (key,)
				continue
			#temp hack for GQ-absent sites
			if not 'GQ' in fields[FORMAT]:
				print >> sys.stderr, '!GQ NOT CALCULATED AT %s' % (key,)
				continue

			try:
				infostr = sd.pop('INFO')
				sd.update(dict([el.split('=') for el in infostr.split(';') if '=' in el]))
			except KeyError:
				pass

			print >> sys.stderr, '%s found ...' % loc,
			#populate individual genotype metrics provided each GQ >= indiv_gt_phred_cut if defined
			sd['indiv_gt'] = {}
			for ind,gt in zip(headers[FORMAT+1:],fields[FORMAT+1:]):
				if ':' in gt:
					this_gt = dict(zip(fields[FORMAT].split(':'),gt.split(':')))
					if indiv_gt_phred_cut is None or float(this_gt['GQ']) >= indiv_gt_phred_cut:
						sd['indiv_gt'][ind] = this_gt
						if return_map:
							new_map[ind].update({loc:''.join([allele_map[loc][n] for n in sd['indiv_gt'][ind]['GT'].split('/')])})
			if not return_map:
				vcf_data[key] = sd
			print >> sys.stderr, '%s individuals processed' % len(sd['indiv_gt'])

	if return_map:
		return new_map
	else:
		return vcf_data

source_map_f, source_vcf_f, new_vcf_f, id_header = sys.argv[1:]

print >> sys.stderr, 'load source map:', source_map_f
loci,geno = extract_genotypes_from_mclgr.load_cross_radtag_genotypes(source_map_f,mIDlookup=False,id_header=id_header)
print >> sys.stderr, '%s loci in %s individuals loaded from source map' % (len(loci),len(geno))


print >> sys.stderr, 'load source vcf:', source_vcf_f
source_vcf = load_vcf(source_vcf_f,loci,indiv_gt_phred_cut=gq)
print >> sys.stderr, '%s loci loaded from source vcf' % (len(source_vcf))

allele_map = {}

for loc in loci.keys():
	if not source_vcf.has_key(tuple(loc.split('.'))):
		print >> sys.stderr, 'no key %s for site %s found in source vcf!' % (tuple(loc.split('.')),loc)
		continue
	vcf_loc = source_vcf[tuple(loc.split('.'))]
	AA_ind  = [k for k,v in geno.items() if v.get(loc,'') == 'AA']
	AA_gt = set([vcf_loc['indiv_gt'][ind]['GT'] for ind in AA_ind if ind in vcf_loc['indiv_gt'].keys()])
	#BB_ind  = [k for k,v in geno.items() if v.get(loc,'') == 'BB']
	#BB_gt = set([vcf_loc['indiv_gt'][ind]['GT'] for ind in BB_ind if ind in vcf_loc['indiv_gt'].keys()])
	if len(AA_gt) != 1: #or len(BB_gt) != 1:
		AA_ctd = Util.countdict([vcf_loc['indiv_gt'][ind]['GT'] for ind in AA_ind if ind in vcf_loc['indiv_gt']])
		if len(AA_ctd) == 2 and min(AA_ctd.values()) == 1:
			print >> sys.stderr, 'ignoring 1 invalid AA genotype from vcf'
		else:
			print >> sys.stderr, '%s invalid homozygotes (AA: %s) ' % (loc,AA_ctd)
			continue
	AA_gt = list(AA_gt)[0]
	#BB_gt = list(BB_gt)[0]
	A = set(AA_gt.split('/'))
	#B = set(BB_gt.split('/'))
	if len(A) != 1: #or len(B) != 1:
		print >> sys.stderr, '%s invalid allele mapping (A: %s B: %s)' % (loc,A,B)
		continue
	A = list(A)[0]
	#B = list(B)[0]
	B = int(A) and '0' or '1'
	allele_map[loc] = {A:'A',B:'B'}
	#print >> sys.stderr, '%s %s' % (loc,allele_map[loc])

if len(allele_map) == 0:
	raise ValueError, 'no loci to load!'
print >> sys.stderr, 'load %s loci from %s' % (len(allele_map), new_vcf_f)
new_geno = load_vcf(new_vcf_f,allele_map,gq,return_map=True)
#print >> sys.stderr, new_geno
extract_genotypes_from_mclgr.output_cross_radtag_genotypes(loci, new_geno, sys.stdout)
