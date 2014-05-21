#!/usr/bin/env python2.7
#

import os
import sys
import ast
import time
import datetime
import operator
from scipy.stats import kde


#
#
#
def reset_score_d(D):
	D['wifi'] = 0.0
	D['app']  = 0.0
	D['usr_load'] = 0.0
	D['sys_load'] = 0.0
	D['noise'] = 0.0
	D['light'] = 0.0
	D['magx'] = 0.0
	D['magy'] = 0.0
	D['magz'] = 0.0
	D['rotx'] = 0.0
	D['roty'] = 0.0
	D['rotz'] = 0.0
	D['accx'] = 0.0
	D['accy'] = 0.0
	D['accz'] = 0.0
	D['call'] = 0.0
	D['active'] = 0.0
	D['charge'] = 0.0
	D['battery'] = 0.0
	return D

#
#
#
def read_model_dict(filename):
	strdata = open( filename ).read()
	return ast.literal_eval( strdata )

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
def compute_comfort(s_score, t_score):
	sensors = {'wifi':1.0, 'app':1.0, 'usr_load':1.0, 'call':1.0,\
					'sys_load':1.0, 'noise':1.0, 'light':1.0, \
					'magx':1.0, 'magy':1.0, 'magz':1.0, \
					'rotx':1.0, 'roty':1.0, 'rotz':1.0, \
					'accx':1.0, 'accy':1.0, 'accz':1.0, \
					'active':1.0, 'battery':1.0, 'charge':1.0}
	wtot_s = 0.0
	cnt_s = 0.0
	wtot_t = 0.0
	cnt_t = 0.0  
	for sensor in sensors.keys():
		if sensor in s_score and s_score[sensor] is not None:
			try:
				cnt_s += sensors[sensor]
				wtot_s += (sensors[sensor] * s_score[sensor])
			except:
				sys.stderr.write( '[ERROR] s_score = ' + str(s_score) + '\n' )
		if sensor in t_score and t_score[sensor] is not None:
			try:
				cnt_t += sensors[sensor]
				wtot_t += (sensors[sensor] * t_score[sensor])
			except:
				sys.stderr.write( '[ERROR] t_score = ' + str(t_score) + '\n' )
	s_comf = 0.0
	if cnt_s != 0:
		s_comf = wtot_s/cnt_s
	t_comf = 0.0
	if cnt_t != 0:
		t_comf = wtot_t/cnt_t
		
	return ( s_comf, t_comf )
	
	
#
#
#
def count_to_freq(D, FRQ_SENSORS, max_count=100):
	for id in D.keys():
		for sensor in D[id]:
			if sensor in FRQ_SENSORS:
				sorted_d = sorted( D[id][sensor].iteritems(),  key=operator.itemgetter(1), reverse=True )
				total =  0.0
				maxval = None
				for (key, count) in sorted_d[0:max_count]:
					total += count
					if maxval is None:
						maxval = count
				D[id][sensor] = {} # remove all 
				for (key, count) in sorted_d[0:max_count]:
					D[id][sensor][key] = (count/total) / (maxval/total) # normalize by maxval so that it provides comfort of 1
	return D		

#
#
#
def count_to_kde(D, KDE_SENSORS):
	problemct = 0
	totalct = 0
	for id in D.keys():
		for sensor in D[id]:
			if sensor in KDE_SENSORS:
				training = []
				maxkey = None
				maxcnt = None
				totvals = float( sum( D[id][sensor].values())  )
				for key in D[id][sensor].keys():
					coef = int( 1000* D[id][sensor][key]/totvals )
					training += ( [key] * coef )
				if len(training) > 0:
					# we have sufficient training data
					try:
						density = kde.gaussian_kde( training )
					except: # Most likely a singular matrix (e.g. [n, n, ... ,n]). Add a little 'noise' to break singularity.
						training += [0.01] * 1
						density = kde.gaussian_kde( training )
					maxval = None
					maxkey = None
					for key in D[id][sensor].keys():
						dk = density(key)[0]
						if maxval is None or maxval < dk:
							maxval = dk
							maxkey = key
					D[id][sensor] = (density, maxval)
				else:
					# we don't have sufficient training data
					D[id][sensor] = None
	return D
	

#
#
# 
def my_float_formatter(fval, points=2):
	formatter = "{0:." + str(points) + "f}"
	return float ( formatter.format(fval) )

#
#
#
def my_truncate(val, scale):
	return int(float(val) / scale) * scale

#
#
#
def write_dict(di, filename):
	fh = open(filename, 'w')
	fh.write( str(di) )
	fh.close()

#
#
#
def update_spatial(S, loc_id_set, tlist, key):
	if len(tlist) > 0:
		for loc_id in loc_id_set: # for each cell id, create a location
			if loc_id not in S:
				S[loc_id] = {}
			if key not in S[loc_id]:
				S[loc_id][key] = {}
			for (ts, li) in tlist:
				for item in li:
					if item not in S[loc_id][key]:
						S[loc_id][key][item] = 0L
					S[loc_id][key][item] += 1
				
#
#
#
def use_spatial(S, loc_id_set, tlist, key, FRQ_SENSORS, KDE_SENSORS):
	score_list = []
	if len(tlist) > 0:
		for loc_id in loc_id_set: # each cell id is a location
			for (ts, li) in tlist:
				if loc_id not in S:
					score_list.append(0.0) # we have not seen this location
				else:
					if key not in S[loc_id]:
						score_list.append(0.0) # we have not seen sensor data for this location
					else:
						for item in li:
							if key in FRQ_SENSORS:
								if item not in S[loc_id][key]:
									score_list.append(0.0) # This item has not been encountered for the location and key
								else:
									score_list.append( S[loc_id][key][item] ) # Return its score
							elif key in KDE_SENSORS:
								if S[loc_id][key] is not None:
									(density, maxval) = S[loc_id][key]
									try:
										score_list.append( float(density(item)) / float(maxval) )
									except:
										pass # possibly /0 don't use kde in case of error
	return my_average( score_list )

#
#
#
def update_temporal(T, tlist, key):
	if len(tlist) > 0:
		for (ts, li) in tlist:
			dt = datetime.datetime.fromtimestamp(ts)
			hourkey = dt.hour
			if hourkey not in T:
				T[hourkey] = {}
			if key not in T[hourkey]:
				T[hourkey][key] = {}
			for item in li:
				if item not in T[hourkey][key]:
					T[hourkey][key][item] = 0L
				T[hourkey][key][item] += 1
		
#
#
#
def use_temporal(T, tlist, key, FRQ_SENSORS, KDE_SENSORS):
	score_list = []
	if len(tlist) > 0:
		for (ts, li) in tlist:
			dt = datetime.datetime.fromtimestamp(ts)
			hourkey = dt.hour
			if hourkey not in T:
				score_list.append(0.0) # we have not seen this hourkey
			else:
				if key not in T[hourkey]:
					score_list.append(0.0) # we have not seen sensor data for this hourkey
				else:
					for item in li:
						if key in FRQ_SENSORS:
							if item not in T[hourkey][key]:
								score_list.append(0.0) # This item has not been encountered for the hourkey and key
							else:
								score_list.append( T[hourkey][key][item] ) # Return its score
						elif key in KDE_SENSORS:
							if T[hourkey][key] is not None:
								(density, maxval) = T[hourkey][key]
								try:
									score_list.append( float(density(item)) / float(maxval) )
								except:
									pass # possibly /0 don't use kde in case of error
	return my_average( score_list )
#
# Average the mag field and rotation values, observed at that second
#
def aggregate_by_ts( li ):
	d = {}
	li2 = []
	for ts, vals in li:
		if ts not in d:
			d[ts] = []
		for val in vals:
			d[ts].append(val)
	for key in d.keys():
		li2.append( (key, [ sum(d[key])/float(len(d[key])) ]) )
	return li2
			

#
#
#
def my_average(li):
	if len(li) == 0:
		return None # Special case which means the sensor did not have any data
	else:
		for i in xrange(0, len(li)):
			if li[i] > 1:
				li[i] = 1.0
		return float( sum(li) ) / float( len(li) )
#
#
#
def time_as_negative(tm, p_tm, timeout=30*60):
	val = 1.0
	print tm-p_tm
	if tm - p_tm < timeout:
		val = float(tm - p_tm) / float( timeout )
	return val

from optparse import OptionParser
parser = OptionParser()
usage = "usage: %prog [options]"
parser = OptionParser(usage=usage)
parser.add_option("-i", "--model-input", dest="input_dir",
                  help="Input model directory.", default=None)
parser.add_option("-o", "--model-output", dest="output_dir",
                  help="Output model directory.", default=None)
parser.add_option("-d", "--detect",
                  action="store_true", dest="detect", default=False,
                  help="Work in detection mode.")
parser.add_option("-t", "--timeout", dest="timeout",
                  help="Number of seconds before score is calculated. Default=600 (10minutes), change to 60 for 1 minute.", default=10*60)
                  
(options, args) = parser.parse_args()


#
# Make sure sane command line arguments are supplied. Else, exit.
#
if options.detect is True and options.input_dir is None:
	print "[ERROR] If in detect mode, must provide input models."
	parser.print_help()
	exit(1)

if options.detect is False and options.output_dir is None:
	print "[ERROR] If in training mode, must provide output model directory."
	parser.print_help()
	exit(1)		

#
# Some constants here
#
SNAME = 1
UNAME = SNAME + 1
sep = '\t'
osep = sep
TIMEOUT = long(options.timeout) # in seconds  
CELLCT = 1 # min num of cells before computing -- we want even buckets for location that describes a place
FRQ_SENSORS = ['wifi', 'app', 'call', 'active', 'charge']
KDE_SENSORS = ['usr_load', 'sys_load', 'noise', 'light', 'magx', 'magy', 'magz', 'rotx', 'roty', 'rotz', 'accx', 'accy', 'accz', 'battery'] 
DEBUG = False
if options.input_dir is not None:
	#
	# Models from file
	#
	print "[INFO] " + str(datetime.datetime.now()) + " Reading T model from file."
	T = read_model_dict(options.input_dir + '/T.model')
	print "[INFO] " + str(datetime.datetime.now()) + " Reading S model from file."
	S = read_model_dict(options.input_dir + '/S.model')	
else:
	#
	# Models from scratch
	#
	T = {}
	for i in xrange(0,24):
		T[i] = {}
	S = {}

if options.detect is True:
	#
	# Prepare models for it
	#
	#
	# Blacklists
	#
	bl_app = set(get_alwayson_apps(T))
	bl_wifi = set()
	
	# remove the 'background' apps
	print "[INFO] " + str(datetime.datetime.now()) + " Removing background applications."
	T = remove_alwayson_apps(T, bl_app)
	S = remove_alwayson_apps(S, bl_app)
	
	#convert counts to frequencies
	print "[INFO] " + str(datetime.datetime.now()) + " Converting counts to frequencies."
	T = count_to_freq(T, FRQ_SENSORS)
	S = count_to_freq(S, FRQ_SENSORS)
	
	#convert counts to kernel density estimation models
	print "[INFO] " + str(datetime.datetime.now()) + " Building kernel density estimators for T."
	T = count_to_kde(T, KDE_SENSORS)
	print "[INFO] " + str(datetime.datetime.now()) + " Building kernel density estimators for S."
	S = count_to_kde(S, KDE_SENSORS)
else:	
	#
	# Blacklists
	#
	bl_app = set()
	bl_wifi = set()

t_score = reset_score_d( {} )
s_score = reset_score_d( {} )
sensor = {}
sensor['cell'] = []
sensor['wifi'] = []
sensor['app']  = []
sensor['usr_load']  = []
sensor['sys_load']  = []
sensor['noise']  = []
sensor['light']  = []
sensor['magx'] = []
sensor['magy'] = []
sensor['magz'] = []
sensor['rotx'] = []
sensor['roty'] = []
sensor['rotz'] = []
sensor['accx'] = []
sensor['accy'] = []
sensor['accz'] = []
sensor['call'] = []
sensor['active'] = [] #display for rice, device for mit
sensor['charge'] = []
sensor['battery'] = []


fh = sys.stdin
line = fh.readline().rstrip('\n')
p_tm = None
linecount = 0
s_patterncount = 0L
t_patterncount = 0L
writtencount = 0L
# NEW: Global score
# accumulates
combined = 0.0

while line:
	linecount +=1
	token = line.split(sep)
	tm = long( token[SNAME - 1] )
	dt = datetime.datetime.fromtimestamp(tm)
	hour = dt.hour
	try:
		# 1. Update sensor data
		#
		if 'CellProbe' in token[SNAME]: # all data
			cell =  token[SNAME + 4]
			cell_ts = float(tm)
			sensor['cell'].append( (tm, cell) )
		elif 'WifiProbe' in token[SNAME] or 'BluetoothProbe' in token[SNAME]: # all data (mit handled as bluetooth)
			wifi_list = token[SNAME + 4:]
			wifi_list = [i for i in wifi_list if i != '' and i not in bl_wifi]
			wifi_ts = float(tm)
			sensor['wifi'].append( (tm, wifi_list) )
		elif 'RunningApplicationsProbe' in token[SNAME]: # all data
			app_list = token[SNAME + 4:]
			app_list = [i for i in app_list if i != '' and i not in bl_app and i.startswith('aexp.sensors') is False]
			app_ts = float(tm)
			sensor['app'].append( (tm, app_list) )
		elif 'SystemProbe' in token[SNAME]: # gcu and rice
			usr_cpu = int( float( token[SNAME + 4] ) * 100 )
			sys_cpu = int( float( token[SNAME + 5] ) * 100 )
			sensor['usr_load'].append( (tm, [usr_cpu]) )
			sensor['sys_load'].append( (tm, [sys_cpu]) )
		elif 'NoiseProbe' in token[SNAME]: # gcu
			min_noise = int( float( token[SNAME + 4] ) )
			mean_noise = int( float( token[SNAME + 5] ) )
			max_noise = int( float( token[SNAME + 6] ) )
			sensor['noise'].append( (tm, [mean_noise]) )
		elif 'LightProbe'  in token[SNAME]: # gcu
			light1 = int( float( token[SNAME + 4] ) )
			sensor['light'].append( (tm, [light1]) )
		elif 'MagneticProbe' in token[SNAME]: # gcu
			magx = my_float_formatter( float( token[SNAME + 4] ) )
			magy = my_float_formatter( float( token[SNAME + 5] ) )
			magz = my_float_formatter( float( token[SNAME + 6] ) )
			sensor['magx'].append( (tm, [magx]) )
			sensor['magy'].append( (tm, [magy]) )
			sensor['magz'].append( (tm, [magz]) )
		elif 'RotationProbe' in token[SNAME]: # gcu
			rotx = my_float_formatter( float( token[SNAME + 4] ) )
			roty = my_float_formatter( float( token[SNAME + 5] ) )
			rotz = my_float_formatter( float( token[SNAME + 6] ) )
			sensor['rotx'].append( (tm, [rotx]) )
			sensor['roty'].append( (tm, [roty]) )
			sensor['rotz'].append( (tm, [rotz]) )
		elif 'AccelProbe' in token[SNAME]: # gcu
			accx = my_float_formatter( float( token[SNAME + 4] ) )
			accy = my_float_formatter( float( token[SNAME + 5] ) )
			accz = my_float_formatter( float( token[SNAME + 6] ) )
			sensor['accx'].append( (tm, [accx]) )
			sensor['accy'].append( (tm, [accy]) )
			sensor['accz'].append( (tm, [accz]) )
		elif 'Mit.CommProbe' in token[SNAME] or 'Rice.CallProbe' in token[SNAME]: # MIT and Rice
			callinfo = ' '.join(token[SNAME + 4:]) # Mit has two cols, rice has one
			sensor['call'].append( (tm, [callinfo]) )
		elif 'Mit.ActiveProbe' in token[SNAME] or 'Rice.DisplayProbe' in token[SNAME]:
			val = token[SNAME + 4]
			if val.lower() == 'true':
				val = 1
			elif val.lower() == 'false':
				val = 0
			else:
				val = int(val)
			sensor['active'].append( (tm, [val]) )
		elif 'Mit.ChargeProbe' in token[SNAME]:
			charge = token[SNAME + 4]
			sensor['charge'].append( (tm, [charge]) )
		elif 'Rice.BatteryProbe' in token[SNAME]:
			batt = int( token[SNAME + 4] )
			sensor['battery'].append( (tm, [batt]) )
	except Exception, e:
		sys.stderr.write( '[ERROR] Error in line ' + str(linecount) + ': ' + str(e) + '. Skipping...\n' )
		sys.stderr.write( '[ERROR] >>>' + line + '<<<\n')	
	# 2. Check if it is time to compute
	#
	if p_tm is None:
		p_tm = tm
	tdiff = tm - p_tm
	# We have seen plenty of sensor data
	if tdiff > TIMEOUT:
		t_patterncount +=1
		loc_id = 'Unknown'
		if len(sensor['cell']) >= CELLCT: # If there is minimum number of cell ids to compute loc_id
			s_patterncount +=1
			# sequence is too restrictive as there may be extra cells or the cell locations may change, use set for now.
			# loc_id = str( sorted( [y for (x,y) in sensor['cell']] ) )
			loc_id = set( [y for (x,y) in sensor['cell']] )		
			# We have spatial data at hand
			# Update or use.
			if options.detect is False:
				# Update the spatial model
				update_spatial(S, loc_id, sensor['wifi'], 'wifi')
				update_spatial(S, loc_id, sensor['app'], 'app')
				update_spatial(S, loc_id, sensor['usr_load'], 'usr_load')
				update_spatial(S, loc_id, sensor['sys_load'], 'sys_load')
				update_spatial(S, loc_id, sensor['noise'], 'noise')
				update_spatial(S, loc_id, sensor['light'], 'light')
				update_spatial(S, loc_id, aggregate_by_ts( sensor['magx'] ), 'magx')
				update_spatial(S, loc_id, aggregate_by_ts( sensor['magy'] ), 'magy')
				update_spatial(S, loc_id, aggregate_by_ts( sensor['magz'] ), 'magz')
				update_spatial(S, loc_id, aggregate_by_ts( sensor['rotx'] ), 'rotx')
				update_spatial(S, loc_id, aggregate_by_ts( sensor['roty'] ), 'roty')
				update_spatial(S, loc_id, aggregate_by_ts( sensor['rotz'] ), 'rotz')
				update_spatial(S, loc_id, aggregate_by_ts( sensor['accx'] ), 'accx')
				update_spatial(S, loc_id, aggregate_by_ts( sensor['accy'] ), 'accy')
				update_spatial(S, loc_id, aggregate_by_ts( sensor['accz'] ), 'accz')
				update_spatial(S, loc_id, sensor['call'], 'call')
				update_spatial(S, loc_id, sensor['active'], 'active')
				update_spatial(S, loc_id, sensor['charge'], 'charge')
				update_spatial(S, loc_id, sensor['battery'], 'battery')
			else:
				# Use the spatial model
				s_score['wifi'] = use_spatial(S, loc_id, sensor['wifi'], 'wifi', FRQ_SENSORS, KDE_SENSORS)
				s_score['app'] = use_spatial(S, loc_id, sensor['app'], 'app', FRQ_SENSORS, KDE_SENSORS)
				s_score['usr_load'] = use_spatial(S, loc_id, sensor['usr_load'], 'usr_load', FRQ_SENSORS, KDE_SENSORS)
				s_score['sys_load'] = use_spatial(S, loc_id, sensor['sys_load'], 'sys_load', FRQ_SENSORS, KDE_SENSORS)
				s_score['noise'] = use_spatial(S, loc_id, sensor['noise'], 'noise', FRQ_SENSORS, KDE_SENSORS)
				s_score['light'] = use_spatial(S, loc_id, sensor['light'], 'light', FRQ_SENSORS, KDE_SENSORS)
				s_score['magx'] = use_spatial(S, loc_id, aggregate_by_ts( sensor['magx'] ), 'magx', FRQ_SENSORS, KDE_SENSORS)
				s_score['magy'] = use_spatial(S, loc_id, aggregate_by_ts( sensor['magy'] ), 'magy', FRQ_SENSORS, KDE_SENSORS)
				s_score['magz'] = use_spatial(S, loc_id, aggregate_by_ts( sensor['magz'] ), 'magz', FRQ_SENSORS, KDE_SENSORS)
				s_score['rotx'] = use_spatial(S, loc_id, aggregate_by_ts( sensor['rotx'] ), 'rotx', FRQ_SENSORS, KDE_SENSORS)
				s_score['roty'] = use_spatial(S, loc_id, aggregate_by_ts( sensor['roty'] ), 'roty', FRQ_SENSORS, KDE_SENSORS)
				s_score['rotz'] = use_spatial(S, loc_id, aggregate_by_ts( sensor['rotz'] ), 'rotz', FRQ_SENSORS, KDE_SENSORS)
				s_score['accx'] = use_spatial(S, loc_id, aggregate_by_ts( sensor['accx'] ), 'accx', FRQ_SENSORS, KDE_SENSORS)
				s_score['accy'] = use_spatial(S, loc_id, aggregate_by_ts( sensor['accy'] ), 'accy', FRQ_SENSORS, KDE_SENSORS)
				s_score['accz'] = use_spatial(S, loc_id, aggregate_by_ts( sensor['accz'] ), 'accz', FRQ_SENSORS, KDE_SENSORS)
				s_score['call'] = use_spatial(S, loc_id, sensor['call'], 'call', FRQ_SENSORS, KDE_SENSORS)
				s_score['active'] = use_spatial(S, loc_id, sensor['active'], 'active', FRQ_SENSORS, KDE_SENSORS)
				s_score['charge'] = use_spatial(S, loc_id, sensor['charge'], 'charge', FRQ_SENSORS, KDE_SENSORS)
				s_score['battery'] = use_spatial(S, loc_id, sensor['battery'], 'battery', FRQ_SENSORS, KDE_SENSORS)
		# Update temporal models no matter if the spatial model is computed or not.
		# Even if we don't know the location, we can update the temporal models
		if options.detect is False:
			# Update the temporal model
			update_temporal(T, sensor['wifi'], 'wifi')
			update_temporal(T, sensor['app'], 'app')
			update_temporal(T, sensor['usr_load'], 'usr_load')
			update_temporal(T, sensor['sys_load'], 'sys_load')
			update_temporal(T, sensor['noise'], 'noise')
			update_temporal(T, sensor['light'], 'light')
			update_temporal(T, aggregate_by_ts( sensor['magx'] ), 'magx')
			update_temporal(T, aggregate_by_ts( sensor['magy'] ), 'magy')
			update_temporal(T, aggregate_by_ts( sensor['magz'] ), 'magz')
			update_temporal(T, aggregate_by_ts( sensor['rotx'] ), 'rotx')
			update_temporal(T, aggregate_by_ts( sensor['roty'] ), 'roty')
			update_temporal(T, aggregate_by_ts( sensor['rotz'] ), 'rotz')
			update_temporal(T, aggregate_by_ts( sensor['accx'] ), 'accx')
			update_temporal(T, aggregate_by_ts( sensor['accy'] ), 'accy')
			update_temporal(T, aggregate_by_ts( sensor['accz'] ), 'accz')
			update_temporal(T, sensor['call'], 'call')
			update_temporal(T, sensor['active'], 'active')
			update_temporal(T, sensor['charge'], 'charge')
			update_temporal(T, sensor['battery'], 'battery')
		else:
			# Use the temporal model
			t_score['wifi'] = use_temporal(T, sensor['wifi'], 'wifi', FRQ_SENSORS, KDE_SENSORS)
			t_score['app']  = use_temporal(T, sensor['app'], 'app', FRQ_SENSORS, KDE_SENSORS)
			t_score['usr_load'] = use_temporal(T, sensor['usr_load'], 'usr_load', FRQ_SENSORS, KDE_SENSORS)
			t_score['sys_load'] = use_temporal(T, sensor['sys_load'], 'sys_load', FRQ_SENSORS, KDE_SENSORS)
			t_score['noise'] = use_temporal(T, sensor['noise'], 'noise', FRQ_SENSORS, KDE_SENSORS)
			t_score['light'] = use_temporal(T, sensor['light'], 'light', FRQ_SENSORS, KDE_SENSORS)
			t_score['magx'] = use_temporal(T, aggregate_by_ts( sensor['magx'] ), 'magx', FRQ_SENSORS, KDE_SENSORS)
			t_score['magy'] = use_temporal(T, aggregate_by_ts( sensor['magy'] ), 'magy', FRQ_SENSORS, KDE_SENSORS)
			t_score['magz'] = use_temporal(T, aggregate_by_ts( sensor['magz'] ), 'magz', FRQ_SENSORS, KDE_SENSORS)
			t_score['rotx'] = use_temporal(T, aggregate_by_ts( sensor['rotx'] ), 'rotx', FRQ_SENSORS, KDE_SENSORS)
			t_score['roty'] = use_temporal(T, aggregate_by_ts( sensor['roty'] ), 'roty', FRQ_SENSORS, KDE_SENSORS)
			t_score['rotz'] = use_temporal(T, aggregate_by_ts( sensor['rotz'] ), 'rotz', FRQ_SENSORS, KDE_SENSORS)
			t_score['accx'] = use_temporal(T, aggregate_by_ts( sensor['accx'] ), 'accx', FRQ_SENSORS, KDE_SENSORS)
			t_score['accy'] = use_temporal(T, aggregate_by_ts( sensor['accy'] ), 'accy', FRQ_SENSORS, KDE_SENSORS)
			t_score['accz'] = use_temporal(T, aggregate_by_ts( sensor['accz'] ), 'accz', FRQ_SENSORS, KDE_SENSORS)
			t_score['call'] = use_temporal(T, sensor['call'], 'call', FRQ_SENSORS, KDE_SENSORS)
			t_score['active'] = use_temporal(T, sensor['active'], 'active', FRQ_SENSORS, KDE_SENSORS)
			t_score['charge'] = use_temporal(T, sensor['charge'], 'charge', FRQ_SENSORS, KDE_SENSORS)
			t_score['battery'] = use_temporal(T, sensor['battery'], 'battery', FRQ_SENSORS, KDE_SENSORS)
		if options.detect is True:
			if writtencount == 0:
				print osep.join( [ "comb_score", "s_score", "t_score", "time_neg_score", \
				"t_wifi", "t_app", "t_usrload", "t_sysload", "t_noise", \
				"t_light", "t_magx", "t_magy", "t_magz", "t_rotx", "t_roty", "t_rotz", "t_accx", "t_accy", "t_accz", "t_call", "t_active", "t_charge", "t_battery",\
				"s_wifi", "s_app", "s_usrload", "s_sysload", "s_noise", \
				"s_light", "s_magx", "s_magy", "s_magz", "s_rotx", "s_roty", "s_rotz", "s_accx", "s_accy", "s_accz", "s_call", "s_active", "s_charge", "s_battery", \
				"tm", "p_tm", "tm_diff", "loc_id"] )
			(sval, tval) = compute_comfort(s_score, t_score)
			# t_as_neg = time_as_negative( tm, p_tm )
			t_as_neg = 1/10.0 # bypass for experiments because we have variable sampling rate
			combined = ( (sval + tval) / 2.0 ) - t_as_neg
			# make sure combined is within [0, 1]
			if combined < -1.0:
				combined = -1.0
			elif combined > 1:
				combined = 1.0
			outstr =  [str(combined), str(sval), str(tval), str(t_as_neg), \
						str(t_score['wifi']), str(t_score['app']), str(t_score['usr_load']), str(t_score['sys_load']), str(t_score['noise']), \
						str(t_score['light']), str(t_score['magx']), str(t_score['magy']), str(t_score['magz']), str(t_score['rotx']), str(t_score['roty']), str(t_score['rotz']), \
						str(t_score['accx']), str(t_score['accy']), str(t_score['accz']), \
						str(t_score['call']), str(t_score['active']), str(t_score['charge']), str(t_score['battery']), \
						str(s_score['wifi']), str(s_score['app']), str(s_score['usr_load']), str(s_score['sys_load']), str(s_score['noise']), \
						str(s_score['light']), str(s_score['magx']), str(s_score['magy']), str(s_score['magz']), str(s_score['rotx']), str(s_score['roty']), str(s_score['rotz']), \
						str(s_score['accx']), str(s_score['accy']), str(s_score['accz']), \
						str(s_score['call']), str(s_score['active']), str(s_score['charge']), str(s_score['battery']), \
						str(tm), str(p_tm), str(tm-p_tm), str(loc_id).replace(' ','')]
			
			print osep.join(outstr)
			writtencount +=1
		# Flush history, whether the bin was used or not
		t_score = reset_score_d(t_score)
		s_score = reset_score_d(s_score)
		sensor['cell'] = []
		sensor['wifi'] = []
		sensor['app']  = []
		sensor['usr_load']  = []
		sensor['sys_load']  = []
		sensor['noise']  = []
		sensor['light']  = []
		sensor['magx'] = []
		sensor['magy'] = []
		sensor['magz'] = []
		sensor['rotx'] = []
		sensor['roty'] = []
		sensor['rotz'] = []
		sensor['accx'] = []
		sensor['accy'] = []
		sensor['accz'] = []
		sensor['call'] = []
		sensor['active'] = []
		sensor['charge'] = []
		sensor['battery'] = []
		p_tm = tm
			
	line = fh.readline().rstrip('\n')
if options.detect is False:
	if not os.path.exists(options.output_dir):
		os.makedirs(options.output_dir) 
	write_dict( T, options.output_dir + '/T.model' )
	write_dict( S, options.output_dir + '/S.model' )

print "[INFO] Number of temporal patterns : " + str(t_patterncount)
print "[INFO] Number of spatial  patterns : " + str(s_patterncount)