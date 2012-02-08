#!/usr/bin/env python

import os,sys
from radtag_denovo import preprocess_radtag_lane
from subprocess import Popen,PIPE

bam,lane,outroot = sys.argv[1:]

idx_len = 6
idx_table = 'DB_multiplex_indices'
idx_field = 'RT'


print >> sys.stderr, 'get index lookup ...',
idx_d = dict([(d['seq'],d['idx']) for d in preprocess_radtag_lane.get_table_as_dict(idx_table)])
print >> sys.stderr, '%s indices' % len(idx_d)

bhandle = Popen('samtools view %s' % bam,shell=True,stdout=PIPE).stdout

try:
    os.makedirs(outroot)
except:
    pass

outfhs = {}
for idx in idx_d.values():
    outf = os.path.join(outroot,'s_%s_sequence_index%s.txt' % (lane,idx))
    outfhs[idx] = open(outf,'w')

outfhs['pass'] = open(os.path.join(outroot,'s_%s_sequence_index%s.txt' % (lane,'PASS')),'w')

print >> sys.stderr, 'process bam'
for samline in bhandle:
    fqstr,idx = preprocess_radtag_lane.sam_line_to_fastq(samline,idx_field,idx_d,idx_len)
    if idx is None:
        outfhs['pass'].write(fqstr)
    else:
        outfhs[idx].write(fqstr)

for outfh in outfhs.values():
    outfh.close()
