#!/bin/bash
#
LEN=$(echo ${#1})
echo $LEN
if [ "$LEN" -eq "0" ]; then
    echo "Usage: ./renameReplicas.sh <castor_dir> <ReplicaNumber2> <ReplicaNumber1>"
    exit 1
fi
read -p "Do you want to rename all files of type *_$2.* to *_$3.* in dir $1 (y/n)? "
if [ $REPLY != "y" ]; then
    echo "Exiting..."
    exit 1
else
    for i in `rfdir $1 | grep "_$2\." | awk '{ print $9 }'`; do 
	echo $1/$i; 
	rfrename $1/$i `echo $1/$i | sed "s/_$2\./_$3\./"`; 
    done
    echo "Done."
#    echo "Done. Check with the rfdir output below:"
#    rfdir $1
    echo "Bye."
    exit 0
fi


