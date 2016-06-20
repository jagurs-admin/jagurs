#!/bin/bash

outfilename="tgs_station"

if [ ${#} -ne "1" ]; then
   echo 'ERROR! Specify target file name!'
   echo "usage: ${0} [filename]"
   exit 0
fi

file="${1}"

if [ ! -r ${file} ]; then
   echo "ERROR! ${file} cannot be read!"
   exit 0
fi

num=`cat ${file} | grep '\[[0-9]*\].*TGS No.=' | wc -l`

echo '=================================================='
echo "Target file: $file"
echo "Number of station: $num"
echo '=================================================='

list=(`cat ${file} | grep '\[[0-9]*\].*TGS No.=' | cut -c 1-8`)

i="0"
while [ ${i} -lt ${num} ]; do
   key=`echo ${list[$i]} | sed 's/\[/\\\\[/' | sed 's/\]/\\\\]/'`
   suffix=`echo ${list[$i]} | cut -c 2-7`
   echo "Output: ${outfilename}.${suffix} (`expr ${i} + 1`/${num})"
   cat ${file} | grep "${key}" | cut -c 9- > ${outfilename}.${suffix}
   i=`expr $i + 1`
done
