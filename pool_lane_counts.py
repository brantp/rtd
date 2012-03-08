import Util
from collections import defaultdict
import gdata_tools
import matplotlib
import pylab
import re

matplotlib.use('Agg')

'''build two dictionaries from one cluni file, 1 for pickcluststokeep fn and 1 for scoring for boxplot construction'''


def clunidict(clunifile):
    '''creates two dictionaries, indbyclust of form {clust1: {i1:r1, i2:r2}, clust2: {i1:r1, i5:r5}}, and fc_ln_idx_indiv_byclust of form {clust1: {(flowcell, lane, index, ind): reads} etc for each cluster in a .cluni file'''
    fh = open(clunifile)
    indbyclust = {}
    fc_ln_idx_indiv_byclust = {}
    for line in fh:
        clustID = line.strip().split()[0]
        if not clustID in indbyclust:
            indbyclust[clustID] = defaultdict(int)
            fc_ln_idx_indiv_byclust[clustID] = defaultdict(int)
        mtch = re.search(r'^.+?\.(.+?)_lane(\d+).+?(?:(?:_index(\d+))|$)',line.strip().split()[1]) #splits out fc, ln, idx info per indiv
        flowcell, lane, index = mtch.groups()
        inds = line.strip().split()[5].split(',') #add parse for flowcell
        reads = line.strip().split()[6].split(',')
        coupled = zip(inds, reads)
        for ind,reads in coupled:
            indbyclust[clustID][ind] += int(reads) 
            fc_ln_idx_indiv_byclust[clustID][(flowcell, lane, index, ind)] += int(reads) 
    return (indbyclust, fc_ln_idx_indiv_byclust)

def pickcluststokeep(indbyclust, fcbyclust, reqreads=4):
    '''takes an indbyclust dictionary and a fcbyclust dictionary and filters the fcbyclust to remove inds with reads < thresh, and clusts with fewer than half the max # of inds seen in a cluster'''
    numinds = []
    for clustID, indcts in indbyclust.items(): #first remove inds that don't have a requisite no. of reads
        for ind, reads in indcts.items(): 
            if reads < reqreads:
                indbyclust[clustID].pop(ind)
    thresh = max(map(len, indbyclust.values()))*0.5
    for clustID, indcts in indbyclust.items():
        if len(indcts) < thresh:
            fcbyclust.pop(clustID)
        else:
            dropk = []
            for k in fcbyclust[clustID]:
                if not k[3] in indbyclust[clustID]: #k[3] designates ind from the key of (flowcell, lane, index, ind)
                    dropk.append(k)
            for k in dropk:
                fcbyclust[clustID].pop(k)          
    return fcbyclust 


def makeboxplot(filteredclusts, dblibrary, figname, pool=False):
    '''takes a filtered dict of clusts worth keeping and creates a boxplot of either by lane (default) or pool'''
    indiv_cluster_count = defaultdict(int) 
    for clust, inddict in filteredclusts.items():
        for ind, reads in inddict.items():
            if ind in indiv_cluster_count.keys():
                indiv_cluster_count[ind]+=1
            else:
                indiv_cluster_count[ind]+=1 
    
    t = gdata_tools.get_table_as_dict(dblibrary)
    db_ind_countd = Util.countdict([d['sampleid'] for d in t if d['sampleid'] in indiv_cluster_count.keys()[3]]) #creates a table of individual dicts from google spreadsheet
    indiv_by_group = defaultdict(list)
    for d in t:
        if 'pool' in d:
            indkey = (d.get('flowcell',None),d.get('lane',None),d.get('index',None),d.get('sampleid',None))
            if indkey in indiv_cluster_count:
                if pool == True:
                    indiv_by_group[(d['flowcell'],d['lane'],d.get('index',None),d['pool'])].append(indiv_cluster_count[indkey]) 
                else:
                    indiv_by_group[(d['flowcell'],d['lane'],d.get('index',None))].append(indiv_cluster_count[indkey])
    
    boxes = []
    labels = []
    for group,indcounts in indiv_by_group.items():
        boxes.append(indcounts)
        labels.append(group)
    boxplt = pylab.figure(1)
    pylab.boxplot(boxes)
    pylab.xticks(arange(1,(len(labels)+1)),labels,fontsize='small') #legend with best location (0) if pools
    boxplt.savefig(figname)
    


#####FROM VCF FILE#####
#def getindsfromvcf(vcffilename):
#    vcf = variant_detection.load_vcf('vcffilename')
#    for sd in vcf.values():
#        for ind in sd['indiv_gt']:
#            indivs.append(ind)
#    indivs = set(indivs)
#    return indivs

#indivs = getindsfromvcf(vcffilename)



          
