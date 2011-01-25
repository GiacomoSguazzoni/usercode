#!/bin/bash
#
LEN=$(echo ${#1})
echo $LEN
if [ "$LEN" -eq "0" ]; then
    echo "Usage: ./createAndRfcpToDir.sh <new_castor_dir>"
    exit 1
fi
read -p "Do you want to rfmkdir $1 and rfcp there all files in current directory (y/n)? "
if [ $REPLY != "y" ]; then
    echo "Exiting..."
    exit 1
else
    echo Creating $1 on castor...
    rfmkdir $1
    echo Done.
    for i in `ls -1`; do 
	echo Copying $i to $1 with rfcp... 
	rfcp $i $1/.
    done
    echo "Done. Check with the rfdir output below:"
    rfdir $1
    echo "Bye."
    exit 0
fi


