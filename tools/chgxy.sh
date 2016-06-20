#!/bin/bash

files=
px=1
py=1
dx=1
dy=1
eps=0.00000001
ffmt="%.10f"
ifmt="%06d"
xbegin="west"
ybegin="north"

usage()
{
   echo "Modified grdfile x/y_min, x/y_max"
   echo ""
   echo "  input   filename : grdfile"
   echo "  outtput filename : grdfile"
   echo ""
}

while [ $# -gt 0 ]; do
    case $1 in
	-ref ) reffile=$2; shift 2;;
	-mod ) modfile=$2; shift 2;;
    esac
done

if [ "$reffile" = "" -o "$modfile" = "" ]; then
    usage
    exit 0
fi

##########################################################################
OPT="-V"

k=0
suffixk=`printf "${ifmt}" $k`

if [ ! -r $reffile.$suffixk ]; then
    echo "not found ref-file : $reffile.$suffixk"
    exit 0
fi

if [ ! -r $modfile.$suffixk ]; then
    echo "not found mod-file : $modfile.$suffixk"
    exit 0
fi

while [ -r $reffile.$suffixk ]
do

    echo $reffile.$suffixk

    x0=`grdinfo $reffile.$suffixk | grep "x_min:" | sed -e "s/.* x_min:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
    y0=`grdinfo $reffile.$suffixk | grep "y_min:" | sed -e "s/.* y_min:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
    x1=`grdinfo $reffile.$suffixk | grep "x_max:" | sed -e "s/.* x_max:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
    y1=`grdinfo $reffile.$suffixk | grep "y_max:" | sed -e "s/.* y_max:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
    x_inc=`grdinfo $reffile.$suffixk | grep "x_inc:" | sed -e "s/.* x_inc:\s* \([0-9][0-9\.]*\)[^0-9].*/\1/g" `
    y_inc=`grdinfo $reffile.$suffixk | grep "y_inc:" | sed -e "s/.* y_inc:\s* \([0-9][0-9\.]*\)[^0-9].*/\1/g" `

# echo " $x0 $x1 $y0 $y1 $x_inc $y_inc "

#    echo "cp -i $modfile.$suffixk $modfile.$suffixk.org"
#    cp -i $modfile.$suffixk $modfile.$suffixk.org

    echo "grdedit $modfile.$suffixk -R$x0/$x1/$y0/$y1"
    grdedit $modfile.$suffixk -R$x0/$x1/$y0/$y1

    if [ ! $? -eq 0 ];then
	echo "grdedit error "
	echo "  in file : $modfile.$suffixk"
	exit 0
    fi

    k=`expr $k + 1`
    suffixk=`printf "${ifmt}" $k`

done
