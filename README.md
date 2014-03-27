userprofiler
============

Building user behavior models based on mobile device sensor data. The 'profiler' builds a probability distribution function for each 'anchor' and sensor. An anchor is either a place or hour-of-day as computed from the sensor data.
Will provide citation to the published work soon.

Running
---------
To build user profiles incrementally run:
	
	./batch-profile-gcu.sh

This script takes the sample data files (currently one user file is provided) and builds behavior models for each day. This allows us to observe the learning rate and decide when the profile 'sufficiently' covers the user behavior. Result files comparisons.tsv and thresholds.tsv provides how learning rate and detection threshold change over time. Under profiles/ a directory is created for each user (profiles/user2.data.gz/) and each user directory contains subdirectories for each day (profiles/user2.data.gz/day001). Each day, two profiles are built: spatial and temporal. T.model and S.model files are dictionary dumps where a spatial key is a set of cell tower ids observed and temporal key is time of day. 


Platforms
---------

Works on Mac OS 10.8 but should compile well on other platforms. profiler needs scipy installed. It also requires a Python 2.7 mostly because scipy worked for 2.7 on my platform. You will also need bash, gzip, egrep, in other words a functioning Linux-like shell. 