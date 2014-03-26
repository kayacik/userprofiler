#!/usr/bin/env python

import operator
import sys
from collections import Counter

#####################
#
# COPIED FROM THE GLOBAL LOCATION ON MBP TO HERE
# TO THAT PROFILER WORKS WITH ALL UNISON MACHINES
#
from optparse import OptionParser
usage = "usage: %prog [args]"
parser = OptionParser(usage=usage)
parser.add_option("-b", "--bc", action="store", type="int", dest="bc", help="Bin count. 100 (percentile) by default, change to 4 for quartiles.", default=100)
(opts, args) = parser.parse_args()

c = Counter()
sep = '\t'
fh = sys.stdin
line = fh.readline().rstrip('\n')
linecount = 0
while line:
	linecount +=1
	try:
		key = float( line )
		c[key] += 1
	except Exception,e:
		pass # junk line
	line = fh.readline().rstrip('\n')
total = sum( c.values() )

for i in xrange(1, opts.bc): 
		curtotal = 0L
		for key in sorted(c.iterkeys()):
			curtotal += c[key]
	 		if  curtotal > (total*i)/float(opts.bc):
	 			print str(i) + sep + str( key )
	 			break 