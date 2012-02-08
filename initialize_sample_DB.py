#!/usr/bin/env python

'''initializes an empty sample database table given data in config.py
uses example index data in "DB_index_by_well.csv" in this directory
'''

import preprocess_radtag_lane
import config
import sys,os

adapt_data = os.path.join(config.RTDROOT,'DB_index_by_well.csv')

try:
    key, gd_client = preprocess_radtag_lane.get_spreadsheet_key(config.LIBRARY_DATA)
    print >> sys.stderr, 'table %s (%s) exists, skip' % (config.LIBRARY_DATA, key)
except:
    headers = ['sampleID','pool','adapter','adaptersversion','flowcell','lane','index']
    preprocess_radtag_lane.create_empty_table(config.LIBRARY_DATA)
    key, gd_client = preprocess_radtag_lane.get_spreadsheet_key(config.LIBRARY_DATA)
    for i,col in enumerate(headers):
        entry = gd_client.UpdateCell(1,2+i,col,key)


try:
    key, gd_client = preprocess_radtag_lane.get_spreadsheet_key(config.ADAPTER_DATA)
    print >> sys.stderr, 'table %s (%s) exists, skip' % (config.ADAPTER_DATA, key)
except:
    preprocess_radtag_lane.create_empty_table(config.ADAPTER_DATA)
    key, gd_client = preprocess_radtag_lane.get_spreadsheet_key(config.ADAPTER_DATA)
    fh = open(adapt_data)
    headers = fh.readline().strip().split(',')[1:]
    for i,col in enumerate(headers):
        entry = gd_client.UpdateCell(1,2+i,col,key)
    for l in fh:
        d = dict(zip(headers,l.strip().split(',')[1:]))
        entry = gd_client.InsertRow(d,key)
