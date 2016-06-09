#this script uploads a file to an ftp server without having to supply passwords etc. for each submission.
#use for loop in bash to supply this script with one input file at a time. (for i in ./*.fastq.gz; do; bash automate_ftp.sh; done

#!/bin/sh
HOST='webin.ebi.ac.uk'
USER='Webin-1217'
PASS='j6cJPKJN'
FILE=$1

ftp -n $HOST <<END_SCRIPT
quote USER $USER
quote PASS $PASS
bin
prompt
mput $FILE
bye
END_SCRIPT
exit 0