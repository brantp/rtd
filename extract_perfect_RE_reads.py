#!/usr/bin/env python

'''
for single read, argv:
cutsite,fq,outfile

for paired end, argv:
cutsite,fq1,fq2,outfile1,outfile2

'''

import preprocess_radtag_lane
import os,sys

barcode_len = 5
tick = 10000 #update progress every this-many reads

if __name__ == "__main__":
    if len(sys.argv) == 4:
        cutsite,fq,outfile = sys.argv[1:]
        rc = preprocess_radtag_lane.get_read_count(fq)
        lnum,baseQ = preprocess_radtag_lane.get_fastq_properties(fq)

        fh = preprocess_radtag_lane.smartopen(fq)
        ofh = preprocess_radtag_lane.smartopen(outfile,'w')

        found = 0
        for i in range(rc):
            if i>0 and i % tick == 0:
                print >> sys.stderr, '\r%s / %s (%0.1f%%) found %s (%0.1f%%)' % \
                      (i,rc,(float(i)/rc)*100,found,(float(found)/i)*100),
            n,s,q = preprocess_radtag_lane.next_read_from_fh(fh,lnum)
            if s[barcode_len:barcode_len+len(cutsite)] == cutsite:
                line = preprocess_radtag_lane.as_fq_line(n,s,q,None,lnum)
                ofh.write(line)
                found += 1
        ofh.close()
    elif len(sys.argv) == 6:
        cutsite,fq1,fq2,outfile1,outfile2 = sys.argv[1:]
        rc1 = preprocess_radtag_lane.get_read_count(fq1)
        rc2 = preprocess_radtag_lane.get_read_count(fq2)
        if rc1 != rc2:
            errstr = 'read count for %s = %s; %s = %s. counts must match' % (fq1,rc1,fq2.rc2)
            raise ValueError, errstr
        lnum,baseQ = preprocess_radtag_lane.get_fastq_properties(fq1)
        
        fh1 = preprocess_radtag_lane.smartopen(fq1)
        ofh1 = preprocess_radtag_lane.smartopen(outfile1,'w')
        fh2 = preprocess_radtag_lane.smartopen(fq2)
        ofh2 = preprocess_radtag_lane.smartopen(outfile2,'w')

        found = 0
        for i in range(rc1):
            if i>0 and i % tick == 0:
                print >> sys.stderr, '\r%s / %s (%0.1f%%) found %s (%0.1f%%)' % \
                      (i,rc1,(float(i)/rc1)*100,found,(float(found)/i)*100),
            n1,s1,q1 = preprocess_radtag_lane.next_read_from_fh(fh1,lnum)
            n2,s2,q2 = preprocess_radtag_lane.next_read_from_fh(fh2,lnum)
            if s1[barcode_len:barcode_len+len(cutsite)] == cutsite:
                line = preprocess_radtag_lane.as_fq_line(n1,s1,q1,None,lnum)
                ofh1.write(line)
                line = preprocess_radtag_lane.as_fq_line(n2,s2,q2,None,lnum)
                ofh2.write(line)
                found += 1
        
        ofh1.close()
        ofh2.close()
    print

        
