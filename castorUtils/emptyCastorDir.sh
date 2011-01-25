#!/bin/bash
#
LEN=$(echo ${#1})
echo $LEN
if [ "$LEN" -eq "0" ]; then
    echo "Usage: ./emptyCastorDir.sh <castor_directory_path>"
    exit 1
fi

#
#
read -p "Do you want to remove all files in $1 (y/n)? "
if [ $REPLY != "y" ]; then
    echo "Exiting..."
    exit 1
else
    echo "Go..."
    for i in `rfdir $1 | awk '{ print $9 }'`; do
	echo Removing $1/$i;
        rfrm $1/$i;
    done
    echo "Done. Directory should be empty. Below rfdir output:"
    rfdir $1
    echo "Bye."
    exit 0
fi



