#!/bin/bash -x 
#
# For each user, determine the number of days active
# and incrementally build profiles for edit distance
# comparison (performed separately).

data="./data"
model="./profiles"
fragmenter="./fragment-data-batch.py"
profiler="./profiler.py"
percentile="./percentile.py" # use the local code
compare="./compare-profiles.py"
threshold=90
maxdaycount=2000


rm -rf $model
mkdir -p $model

for file in `ls -1 $data | egrep '.gz'`; do
	mkdir -p $model"/"$file
	dayfiles="$data"/"$file".days
	rm -rf $dayfiles
	mkdir $dayfiles
	gzip -cd $data/$file | $fragmenter day $dayfiles
	cur=""
	pre=""
	daycount=0
	for dayfile in `ls -1 $dayfiles | sort`; do
		(( daycount += 1 ))
		if [ "$daycount" -gt "$maxdaycount" ]; then
			# Done for this user. Maximum day limit reached.
			continue
		fi 
		cur="$dayfile"
		echo "[INFO] current=$cur"
		rm -rf $model"/"$file"/"$cur
		if [ -z "$pre" ]; then
			cat $dayfiles/$dayfile | $profiler -o $model"/"$file"/"$cur/
		else
			cat $dayfiles/$dayfile | $profiler -i $model"/"$file"/"$pre/ -o $model"/"$file"/"$cur/
			cat $dayfiles/$dayfile | $profiler -d -i $model"/"$file"/"$pre/ | egrep -v '(INFO|ERROR)' > $model"/"$file"/"$cur/comfort.data
			comf=`cat $model"/"$file"/"$cur/comfort.data | awk '{print $1}' | $percentile | egrep "$threshold\s" | awk '{print $2}'`
			echo -e $cur'\t'$comf >>   $model"/"$file"/thresholds.tsv"
			$compare $model"/"$file"/"$pre/ $model"/"$file"/"$cur/ >> $model"/"$file"/comparisons.tsv"
		fi
		echo "[INFO] comfort=$comf"
		pre="$cur"
	 done
	 # Clean-up
	 rm -rf $dayfiles
done