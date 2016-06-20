#!/bin/bash

output=output
inputs=()
i=0

x_size_cm=22

usage()
{
   echo "$0 [-o output-file] file1 file2 ..."
   echo ""
   echo "   file1.xyz, file2.xyz : results of grd2xyz"
   echo "   <output-file>.diff   : result of diff -y --suppress-common-lines"
   echo ""
}

while [ $# -gt 0 ]; do
  case $1 in
   -o ) output=$2; shift 2;;
   *  ) inputs=("${inputs[@]}" $1); shift 1;;
  esac 
done
echo "${inputs[@]}"

if [ ${#inputs[@]} -lt 2 ]; then
  usage
  exit 0
fi

x0_min=`grdinfo ${inputs[0]} | grep "x_min" | sed -e "s/.*x_min: \([0-9][0-9.]*\).*/\1/g"`
x0_max=`grdinfo ${inputs[0]} | grep "x_max" | sed -e "s/.*x_max: \([0-9][0-9.]*\).*/\1/g"`
x0_inc=`grdinfo ${inputs[0]} | grep "x_inc" | sed -e "s/.*x_inc: \([0-9][0-9.]*\).*/\1/g"`
y0_min=`grdinfo ${inputs[0]} | grep "y_min" | sed -e "s/.*y_min: \([0-9][0-9.]*\).*/\1/g"`
y0_max=`grdinfo ${inputs[0]} | grep "y_max" | sed -e "s/.*y_max: \([0-9][0-9.]*\).*/\1/g"`
y0_inc=`grdinfo ${inputs[0]} | grep "y_inc" | sed -e "s/.*y_inc: \([0-9][0-9.]*\).*/\1/g"`

x1_min=`grdinfo ${inputs[1]} | grep "x_min" | sed -e "s/.*x_min: \([0-9][0-9.]*\).*/\1/g"`
x1_max=`grdinfo ${inputs[1]} | grep "x_max" | sed -e "s/.*x_max: \([0-9][0-9.]*\).*/\1/g"`
x1_inc=`grdinfo ${inputs[1]} | grep "x_inc" | sed -e "s/.*x_inc: \([0-9][0-9.]*\).*/\1/g"`
y1_min=`grdinfo ${inputs[1]} | grep "y_min" | sed -e "s/.*y_min: \([0-9][0-9.]*\).*/\1/g"`
y1_max=`grdinfo ${inputs[1]} | grep "y_max" | sed -e "s/.*y_max: \([0-9][0-9.]*\).*/\1/g"`
y1_inc=`grdinfo ${inputs[1]} | grep "y_inc" | sed -e "s/.*y_inc: \([0-9][0-9.]*\).*/\1/g"`

echo "x_min   = $x0_min , $x1_min"
echo "x_max   = $x0_max , $x1_max"
echo "x_inc   = $x0_inc , $x1_inc"
echo "y_min   = $y0_min , $y1_min"
echo "y_max   = $y0_max , $y1_max"
echo "y_inc   = $y0_inc , $y1_inc"

x_min=$x0_min
x_max=$x0_max
x_inc=$x0_inc
y_min=$y0_min
y_max=$y0_max
y_inc=$y0_inc

x_cntr=`echo "${x_min} + ( ${x_max} - ${x_min} ) * 0.5" | bc -l`
x_scl=`echo "${x_size_cm} / ( ${x_max} - ${x_min} )" | bc -l`

region=`printf "%.4f/%.4f/%.4f/%.4f" $x_min $x_max $y_min $y_max`
inc=`printf "%.8f/%.8f" $x_inc $y_inc`
x_center=`printf "%.4f" $x_cntr`
x_scale=`printf "%.6f" $x_scl`

echo "x_min  = $x_min"
echo "x_max  = $x_max"
echo "x_inc  = $x_inc"
echo "x_cntr = $x_cntr"
echo "x_scl  = $x_scl"
echo "y_min  = $y_min"
echo "y_max  = $y_max"
echo "y_inc  = $y_inc"
echo "region = $region"
echo "inc    = $inc"

grdedit ${inputs[0]} -R${region}
grdedit ${inputs[1]} -R${region}
grd2xyz ${inputs[0]} > ${inputs[0]}.xyz
grd2xyz ${inputs[1]} > ${inputs[1]}.xyz

diff -y --suppress-common-line ${inputs[0]}.xyz ${inputs[1]}.xyz > ${output}.diff
#diff -y ${inputs[0]}.xyz ${inputs[1]}.xyz > ${output}.diff

#echo "grdmath ${inputs[0]} ${inputs[1]} SUB =  ${output}.diff"
#grdmath ${inputs[0]} ${inputs[1]} SUB =  ${output}.diff
#echo "grdmath ${output}.diff ${inputs[1]} DIV = ${output}.error"
#grdmath ${output}.diff ${inputs[1]} DIV = ${output}.error
#
#echo "pscoast -R${region} -Ba5f1g1 -G128/128/128 -S255/255/255 -Jq${x_center}/${x_scale}c -P -K > ${output}.diff.eps"
#pscoast -R${region} -Ba5f1g1 -G128/128/128 -S255/255/255 -Jq${x_center}/${x_scale}c -K > ${output}.diff.eps
#echo "grd2cpt ${output}.diff > ${output}.diff.cpt"
#grd2cpt ${output}.diff > ${output}.diff.cpt
#if [ $? -eq 1 ]; then
#palette=0.01
#else
#palette=${output}.diff.cpt
#fi
#echo "grdcontour ${output}.diff -R${region} -Ba5f1g1 -C${palette} -Jq${x_center}/${x_scale}c -P -O >> ${output}.diff.eps"
#grdcontour ${output}.diff -R${region} -Ba5f1g1 -C${palette} -Jq${x_center}/${x_scale}c -O >> ${output}.diff.eps
#echo "convert -rotate 90 ${output}.diff.eps ${output}.diff.png"
#convert -rotate 90 ${output}.diff.eps ${output}.diff.png


