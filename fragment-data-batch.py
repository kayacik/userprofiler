#!/usr/bin/env python
#

import os
import sys
import ast
import time
import datetime
import operator

TSI = 0

if len(sys.argv) != 3:
	print 'Usage: ' + sys.argv[0] + ' [hour|day|week] <dest>'
	exit(1)


TIMEFRAME = None
output = sys.argv[2]
if sys.argv[1] == 'hour':
	 TIMEFRAME = 60*60
elif sys.argv[1] == 'day':
	TIMEFRAME = 24*60*60
elif sys.argv[1] == 'week':
	TIMEFRAME = 7*24*60*60
else:
	print 'Usage: ' + sys.argv[0] + ' [hour|day|week] <dest>'
	exit(1)

if not os.path.exists(output):
    os.makedirs(output) 

start = None
sep = '\t'
fh = sys.stdin
line = fh.readline().rstrip('\n')
data = []
p_frameidx = None
writtensofar = 0
while line:
	token = line.split(sep)
	tm = long( token[TSI] )
	year = datetime.datetime.fromtimestamp(tm).year 
	if year > 2003 and year < 2015: # Filter junk lines
		if start is None:
			start = tm
	
		if tm - start < TIMEFRAME:
			data.append(line)
		else:
			start = tm	
			if len(data) > 0:
				fname = output + '/' + sys.argv[1] + '%0*d' % (3, writtensofar+1)
				outfh = open(fname, 'w')		
				for line in data:
					outfh.write(line + '\n')
				outfh.close()
				writtensofar +=1
				data = []
	line = fh.readline().rstrip('\n')
	
