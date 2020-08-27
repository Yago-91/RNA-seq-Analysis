#!/bin/bash

#sra file download trials: aria2c, wget and curl

#------------------------Aria2 section--------------------------#
echo "Aria2 section started" $(date +"%r") > down_trials_log.txt
while read line

	do aria2c --file-allocation=none -c -x 16 -s 16 -d . $line
	
done < acc_url.txt
##Where acc_url.txt is the file with downloading urls

echo "Aria2 finished" $(date +"%r") >> down_trials_log.txt

#---------------------------------------------------------------#

#------------------------wget section---------------------------#
echo "wget section started" $(date +"%r") >> down_trials_log.txt
while read line

	do wget --quiet $line &
	
done < acc_url.txt

wait

echo "Wget finished" $(date +"%r") >> down_trials_log.txt

#---------------------------------------------------------------#

#-------------------------curl section--------------------------#
echo "curl section started" $(date +"%r") >> down_trials_log.txt
while read line

    name=$(echo $line | cut -d ' ' -f1)
    url=$(echo $line | cut -d ' ' -f2)
	do curl $url --output $name &
	
done < acc_url_names.txt
#Where acc_url_names.txt is a file with two space-separated columns being 
#the first the file name and the second the url link for the download.

wait

echo "curl finished" $(date +"%r") >> down_trials_log.txt

#---------------------------------------------------------------#
