#!/bin/bash

for completefile in $(cmsLs  $1)
do 

#echo $completefile

if [ `expr "$completefile" : $1` != 0 ]; then
    echo "remove $completefile"
    cmsRm $completefile
#    export file=$(basename $completefile)

#    if [ -f /tmp/venturia/$file ]; then
#	echo "File exist already";
#    else
#	cmsStageIn ${EOSBASEDIR}/${EOSDIR}/$file /tmp/venturia/$file
#    fi
fi

done

cmsRmdir $1

