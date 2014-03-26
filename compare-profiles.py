#!/usr/bin/env python2.7
#

import os
import sys
import ast
import math
import time
import datetime
import operator
from scipy.stats import kde

#
#
# brute forced
def get_alwayson_apps(T, threshold=23):
	d = {}
	aon = []
	for hour in T.keys():
		if 'app' in T[hour]:
			for app in T[hour]['app'].keys():
				if app not in d:
					d[app] = set()
				d[app].add(hour)
	for app in d.keys():
		if len( d[app] ) >= threshold:
			aon.append(app)
	return aon

#
#
#
def remove_alwayson_apps(D, bl):
	for k in D.keys():
		if 'app' in D[k]: # not all locations in S may have an app
			for item in bl:
				D[k]['app'].pop( item, None )
	return D

#
#
#
def read_model_dict(filename):
	strdata = open( filename ).read()
	return ast.literal_eval( strdata )

#
#
#
def prepare_model( M, FRQ_SENSORS, KDE_SENSORS, cutoff=100 ):
	for key in M.keys():
		for sensor in M[key].keys():
			if sensor in FRQ_SENSORS:
				sorted_d = sorted( M[key][sensor].iteritems(), key=operator.itemgetter(1), reverse=True  )[:cutoff]
				items = [x for x,y in sorted_d]
				M[key][sensor] = items
			elif sensor in KDE_SENSORS:
				perc = 1
				total = sum( M[key][sensor].values() )
				li = []
				for i in xrange(1, 100):
					curct = 0.0
					for (val, cnt) in sorted( M[key][sensor].iteritems() ):
						curct += cnt
						if curct > (i*total)/100.00:
							li.append(val)
							break
				M[key][sensor] = li					
	return M

#
#
#
def levenshtein(list1, list2):
    oneago = None
    thisrow = range(1, len(list2) + 1) + [0]
    for x in xrange(len(list1)):
        twoago, oneago, thisrow = oneago, thisrow, [0] * len(list2) + [x + 1]
        for y in xrange(len(list2)):
            delcost = oneago[y] + 1
            addcost = thisrow[y - 1] + 1
            subcost = oneago[y - 1] + (list1[x] != list2[y])
            thisrow[y] = min(delcost, addcost, subcost)
    return thisrow[len(list2) - 1]

#
# http://en.wikipedia.org/wiki/Jaccard_index
#
def compute_jaccard_distance(list1, list2):
	mset1 =  set(list1)
	mset2 =  set(list2)
	AiB = float( len( mset1.intersection(mset2) ) )
	AuB = float( len( mset1.union(mset2) ) )
	try:
		return (AuB - AiB) / AuB
 	except ZeroDivisionError:
 		return 0.0

def normalized_euclidean(list1, list2):
	if len(list1) != len(list2):
		print "[ERROR] Euclidean distance, list lengths do not match"
		print "[ERROR] l1 = " + str(list1)
		print "[ERROR] l2 = " + str(list2)
		return 1.0  
	minval = min( list1 + list2 )
	maxval = max( list1 + list2 )
	if minval == maxval:
		return 0.0 # all the values are identical, dist->0
	for i in xrange(0, len(list1)):
		# normalize lists
		list1[i] = ( list1[i] - minval )  / float( maxval - minval )
		list2[i] = ( list2[i] - minval )  / float( maxval - minval )
	dist = 0.0
	for i in xrange(0, len(list1)):
		dist += (list1[i] - list2[i])**2
	return math.sqrt(dist) / math.sqrt(len(list1)) # normalize the distance to one


#
#
#
def compare_mods(M1, M2, FRQ_SENSORS, KDE_SENSORS):
	D = {}
	# Comparing M2 to M1, hence loops over the hour/loc keys of 2
	for key2 in M2.keys():
		for sensor in M2[key2]:
			dist  = 1.0
			try:
				if sensor in FRQ_SENSORS:
					dist = levenshtein( M1[key2][sensor], M2[key2][sensor] ) / float( len(M2[key2][sensor]) ) # will throw exception if no key
				elif sensor in KDE_SENSORS:
					dist = normalized_euclidean( M1[key2][sensor], M2[key2][sensor] )
			except KeyError:
				dist = 1.0 # M1 either does not have key2 or sensor as keys
			except ZeroDivisionError:
				dist = 1.0 # if M2 has no 'sensor' data but M1 has some. 
			if sensor not in D:
				D[sensor] = []
			D[sensor].append( dist )
	return D			
			
			


if len(sys.argv) != 3:
	print "[ERROR] " + sys.argv[0] + " <profile_folder1> <profile_folder2> " 
	exit(1)


FRQ_SENSORS = ['wifi', 'app', 'call', 'active', 'charge']
KDE_SENSORS = ['usr_load', 'sys_load', 'noise', 'light', 'magx', 'magy', 'magz', 'rotx', 'roty', 'rotz', 'battery'] 


folder1 = sys.argv[1]
folder2 = sys.argv[2]


#print "[INFO] " + str(datetime.datetime.now()) + " Reading T1 model from file."
T1 = read_model_dict(folder1 + '/T.model')
#print "[INFO] " + str(datetime.datetime.now()) + " Reading S1 model from file."
S1 = read_model_dict(folder1 + '/S.model')	
bl_app1 = set(get_alwayson_apps(T1))
T1 = remove_alwayson_apps(T1, bl_app1)
#S1 = remove_alwayson_apps(S1, bl_app1)


#print "[INFO] " + str(datetime.datetime.now()) + " Reading T2 model from file."
T2 = read_model_dict(folder2 + '/T.model')
#print "[INFO] " + str(datetime.datetime.now()) + " Reading S2 model from file."
S2 = read_model_dict(folder2 + '/S.model')	
bl_app2 = set(get_alwayson_apps(T2))
T2 = remove_alwayson_apps(T2, bl_app2)
#S2 = remove_alwayson_apps(S2, bl_app2)


T1 = prepare_model( T1, FRQ_SENSORS, KDE_SENSORS )
S1 = prepare_model( S1, FRQ_SENSORS, KDE_SENSORS )
T2 = prepare_model( T2, FRQ_SENSORS, KDE_SENSORS )
S2 = prepare_model( S2, FRQ_SENSORS, KDE_SENSORS )

dT = compare_mods(T1, T2, FRQ_SENSORS, KDE_SENSORS)
dS = compare_mods(S1, S2, FRQ_SENSORS, KDE_SENSORS)

for item in sorted((FRQ_SENSORS + KDE_SENSORS)):
	print 't_' + item + '\t',
print "t_glob\t",
for item in sorted((FRQ_SENSORS + KDE_SENSORS)):
	print 's_' + item + '\t',
print "s_glob\tx_glob"
tg = 0.0
tc = 0.0
sg = 0.0
sc = 0.0
for item in sorted((FRQ_SENSORS + KDE_SENSORS)):
	if item in dT:
		val = sum(dT[item]) / float( len(dT[item]) )
		tg += val
		tc += 1
	else:
		val = None
	print str(val) + '\t',

try:
	tratio = tg/tc
except:
	tratio = 0.0 	
	
print str(tratio) + '\t',
for item in sorted((FRQ_SENSORS + KDE_SENSORS)):
	if item in dS:
		val = sum(dS[item]) / float( len(dS[item]) )
		sg += val
		sc += 1
	else:
		val = None
	print str(val) + '\t',

try:
	sratio = sg/sc
except:
	sratio = 0.0 	

try:
	gratio = (sg+tg)/(sc+tc)
except:
	gratio = 0.0

print str(sratio) + '\t' + str ( gratio )
#print dT
#print dS