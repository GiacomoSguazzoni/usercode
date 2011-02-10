#!/bin/bash
#
LEN=$(echo ${#1})
echo $LEN
if [ "$LEN" -eq "0" ]; then
    echo "Usage: ./createAndRfcpToDir.sh <existing_castor_base_dir> <dir_to_copy>"
    exit 1
fi
read -p "Do you want to rfmkdir $1/$2 and rfcp there the contents of $2? (y/n) "
if [ $REPLY != "y" ]; then
    echo "Exiting..."
    exit 1
else
    echo Creating $1 on castor...
    rfmkdir -p $1/$2
    echo Done.
    for i in `ls -1 $2`; do 
	echo Copying $i to $1/$2 with rfcp... 
	rfcp $2/$i $1/$2/.
    done
    echo "Done. Check with the rfdir output below:"
    echo "Done."
    echo "Bye."
    exit 0
fi


