#!/bin/bash
#
LEN=$(echo ${#1})
echo $LEN
if [ "$LEN" -eq "0" ]; then
    echo "Usage: ./removeTrailingGridCode.sh <castor_directory_path>"
    exit 1
fi

#
#
read -p "Do you want to rename all files in $1? (y/n)"
if [ $REPLY != "y" ]; then
    echo "Exiting..."
    exit 1
else
    read -p "Are you sure. It is dangerous if already done once? (y/n)"
    if [ $REPLY != "y" ]; then
	echo "Exiting..."
	exit 1
    else
	echo "Go..."
	for i in `rfdir $1 | awk '{ print $9 }'`; do
	    echo Renaming $i to `echo $i | sed 's/_...\.root$/\.root/'`;
	    rfrename $1/$i `echo $1/$i | sed 's/_...\.root$/\.root/'`;
	done
	echo "Done. Below rfdir output for checking..."
	rfdir $1
	echo "Bye."
	exit 0
    fi
fi


