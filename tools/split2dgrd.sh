#!/bin/bash

files=
px=1
py=1
dx=1
dy=1
eps=0.0000000000001
ffmt="%.15f"
ifmt="%06d"
xbegin="west"
ybegin="north"
reffile=""
# === No copy ==================================================================
nocopy="0"
# ==============================================================================

# new
#xbegin="west"
#ybegin="north"

usage()
{
   echo "2D Split GRD File Merger"
   echo ""
   echo "$0 [-px <proc-x>] [-py <proc-y>] [-dx <edge-x>] [-dy <edge-y>] [-eps <eps>] [-ffmt <fmt>] [-ifmt <fmt>] grdfile"
   echo ""
   echo "  -px <proc-x>    : number of processes on coordinate x (default:$px)"
   echo "  -py <proc-y>    : number of processes on coordinate y (default:$py)"
   echo "  -dx <edge-x>    : overlap edge size on coordinate x (default:$dx)"
   echo "  -dy <edge-y>    : overlap edge size on coordinate y (default:$dy)"
   echo "  -eps <eps>      : epsilon of lon/lat value (default:$eps)"
   echo "  -ffmt <fmt>     : format string for float numbers (default:$ffmt)"
   echo "  -ifmt <fmt>     : format string for int   numbers (default:$ifmt)"
   echo "  -xbegin <we>    : <we>={west|east|w|e}   (default:$xbegin)"
   echo "  -ybegin <ns>    : <ns>={north|south|n|s} (default:$ybegin)"
   echo "  -ref <ref-file> : reference input file (default:$reffile)"
# === No copy ==================================================================
   echo "  -nocopy         : edit input files directly (default:off)"
# ==============================================================================
   echo "  grdfile         : .grd file basename"
   echo ""
   echo "  using commands  : grdpaste, grdinfo, grfpaste"
   echo ""
   echo "  intput  filename : grdfile.000000, grdfile.000001, ..."
   echo "  outtput filename : grdfile"
   echo ""
}

while [ $# -gt 0 ]; do
   case $1 in
     -px  ) px=$2; shift 2;;
     -py  ) py=$2; shift 2;;
     -dx  ) dx=$2; shift 2;;
     -dy  ) dy=$2; shift 2;;
     -eps ) eps=$2; shift 2;;
     -ffmt) ffmt=$2; shift 2;;
     -ifmt) ifmt=$2; shift 2;;
     -xbegin ) 
          case $2 in
            west | w ) xbegin="west";;
            east | e ) xbegin="east";;
            * ) echo "xbegin : invalid value"; usage; exit 1;;
          esac
          shift 2;;
     -ybegin ) 
          case $2 in
            north | n ) ybegin="north";;
            south | s ) ybegin="south";;
            * ) echo "ybegin : invalid value"; usage; exit 1;;
          esac
          shift 2;;
     -ref ) reffile=$2; shift 2;;
# === No copy ==================================================================
     -nocopy ) nocopy=1; shift 1;;
# ==============================================================================
     *    ) files="$files $1"; shift 1;;
   esac
done

if [ "$files" = "" ]; then
   usage
   exit 0
fi
# === No copy ==================================================================
cpcommand="cp"
if [ "$nocopy" == "1" ]; then
   cpcommand="ln"
fi
# ==============================================================================

for file in $files ; do

file=`basename $file .000000`

if [ -r $file ]; then
    echo ""
    echo " $file exist. "
    echo ""
    usage
    exit 0
fi

if [ ! $reffile == "" ]; then
    if [ ! -r $reffile ]; then
	echo "not found $reffile."
	exit 0
    else
	if [ ! -x ./chgxy.sh ]; then
	    echo "not found script [./chgxy.sh]. "
	    exit 0
	fi
	echo "sh ./chgxy.sh -ref $reffile -mod $file"
	sh ./chgxy.sh -ref $reffile -mod $file
    fi
fi

# === Use GMT_TMPDIR with isolation mode. ======================================
export GMT_TMPDIR=`mktemp -d ./gmt.XXXXXX`
# ==============================================================================
# gmt default
# === Use GMT_TMPDIR with isolation mode. ======================================
#gmtdefaults -L > .gmtdefaults4.save
gmtdefaults -L > ${GMT_TMPDIR}/.gmtdefaults4.save
# ==============================================================================
gmtset D_FORMAT $ffmt

# === Use GMT_TMPDIR with isolation mode. ======================================
#tmpfile=__tmp.grd
tmpfile=${GMT_TMPDIR}/__tmp.grd
# ==============================================================================

##########################################################################
OPT="-V"
#OPT=""

k=0
j=0
while [ $j -lt $py ]; do

    suffixj=`printf "${ifmt}" $j`

    i=0
    while [ $i -lt $px ]; do

	suffixk=`printf "${ifmt}" $k`

	infile=$file.$suffixk
   	jfile=__$file.j.$suffixj

	if [ $px -eq 1 ];then
# === No copy ==================================================================
#	    echo "cp $infile $jfile"
#	    cp $infile $jfile
	    echo "$cpcommand $infile $jfile"
	    $cpcommand $infile $jfile
# ==============================================================================
# === jfile must have write permission... ======================================
            chmod +w $jfile
# ==============================================================================
	    k=`expr $k + 1`
	    break;
	fi

	x0=`grdinfo $infile | grep "x_min:" | sed -e "s/.* x_min:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
	y0=`grdinfo $infile | grep "y_min:" | sed -e "s/.* y_min:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
	x1=`grdinfo $infile | grep "x_max:" | sed -e "s/.* x_max:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
	y1=`grdinfo $infile | grep "y_max:" | sed -e "s/.* y_max:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
	x_inc=`grdinfo $infile | grep "x_inc:" | sed -e "s/.* x_inc:\s* \([0-9][-+e0-9\.]*\)[^-+e0-9].*/\1/" `
	y_inc=`grdinfo $infile | grep "y_inc:" | sed -e "s/.* y_inc:\s* \([0-9][-+e0-9\.]*\)[^-+e0-9].*/\1/" `
	x_inc=`echo $x_inc | sed s/e/*10^/`
	y_inc=`echo $y_inc | sed s/e/*10^/`

# echo "x0, x1 = $x0, $x1"
# echo "y0, y1 = $y0, $y1"
# echo "x_inc, y_inc = $x_inc, $y_inc"

	if [ $i -eq 0 ];then

# === No copy ==================================================================
#	    echo "cp $infile $jfile"
#	    cp $infile $jfile
	    echo "$cpcommand $infile $jfile"
	    $cpcommand $infile $jfile
# ==============================================================================
# === jfile must have write permission... ======================================
            chmod +w $jfile
# ==============================================================================

	else

	    if [ "$xbegin" = "east" ]; then
		if [ $i -eq 0 ]; then
		    maxlon=`echo "$x1 - $eps" | bc -l`
		else
		    maxlon=`echo "$x1 - $x_inc*$dx - $eps" | bc -l`
		fi
		minlon=`echo "$x0 + $eps" | bc -l`
	    else
		if [ $i -eq 0 ]; then
		    minlon=`echo "$x0 + $eps" | bc -l`
		else
  		    minlon=`echo "$x0 + $x_inc*$dx" | bc -l`
		fi
		maxlon=`echo "$x1 - $eps" | bc -l`
	    fi

	    minlat=`echo "$y0 + $eps" | bc -l`
	    maxlat=`echo "$y1 - $eps" | bc -l`

	    minlon=`printf "${ffmt}" $minlon`
	    maxlon=`printf "${ffmt}" $maxlon`
	    minlat=`printf "${ffmt}" $minlat`
	    maxlat=`printf "${ffmt}" $maxlat`

	    echo "grdcut $OPT $infile -R$minlon/$maxlon/$minlat/$maxlat -G$tmpfile=cf"
	    grdcut $OPT $infile -R$minlon/$maxlon/$minlat/$maxlat -G$tmpfile=cf

	    if [ ! $? -eq 0 ];then
		echo "grdcut error "
		echo "  in file : $infile"
		exit 0
	    fi

	    if [ "$xbegin" = "east" ]; then
 		if [ $i -ne 0 ]; then
  		    maxlon=`grdinfo $jfile | grep "x_min:" | sed -e "s/.* x_min:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
 		fi
	    else
 		if [ $i -ne 0 ]; then
  		    minlon=`grdinfo $jfile | grep "x_max:" | sed -e "s/.* x_max:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
 		fi
	    fi

 	    echo "grdedit -A $tmpfile -R$minlon/$maxlon/$minlat/$maxlat"
 	    grdedit -A $tmpfile -R$minlon/$maxlon/$minlat/$maxlat

	    if [ ! $? -eq 0 ];then
		echo "grdedit error "
		echo "  in file : $tmpfile"
		exit 0
	    fi

	    echo "grdpaste $OPT $jfile $tmpfile -G${jfile}=cf"
	    grdpaste $OPT $jfile $tmpfile -G${jfile}=cf

	    if [ ! $? -eq 0 ];then
		echo "grdpaste error."
		echo "  in file1 : $jfile"
		echo "  in file2 : $tmpfile"
		echo "  out file : $jfile"
		exit 0
	    else
		echo "rm $tmpfile"
		rm $tmpfile
	    fi

 	    echo "grdedit -A $jfile"
 	    grdedit -A $jfile

	fi

        ##########################################################################

	i=`expr $i + 1`
	k=`expr $k + 1`

    done

    if [ $py -eq 1 ];then
	echo "mv $jfile $file"
	mv $jfile $file
	break;
    fi

    if [ $j -ge 1 ];then

	x0=`grdinfo $jfile | grep "x_min:" | sed -e "s/.* x_min:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
	y0=`grdinfo $jfile | grep "y_min:" | sed -e "s/.* y_min:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
	x1=`grdinfo $jfile | grep "x_max:" | sed -e "s/.* x_max:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
	y1=`grdinfo $jfile | grep "y_max:" | sed -e "s/.* y_max:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
	x_inc=`grdinfo $jfile | grep "x_inc:" | sed -e "s/.* x_inc:\s* \([0-9][0-9\.]*\)[^0-9].*/\1/g" `
	y_inc=`grdinfo $jfile | grep "y_inc:" | sed -e "s/.* y_inc:\s* \([0-9][0-9\.]*\)[^0-9].*/\1/g" `

	if [ "$ybegin" = "south" ]; then
	    minlat=`echo "$y0 + $y_inc*$dy + $eps" | bc -l`
	    maxlat=`echo "$y1 - $eps" | bc -l`
	else
	    minlat=`echo "$y0 + $eps" | bc -l`
	    maxlat=`echo "$y1 - $y_inc*$dy - $eps" | bc -l`
	fi

	minlon=`echo "$x0 + $eps" | bc -l`
	maxlon=`echo "$x1 - $eps" | bc -l`

	echo "grdcut $OPT $jfile -R$minlon/$maxlon/$minlat/$maxlat -G${tmpfile}=cf"
	grdcut $OPT $jfile -R$minlon/$maxlon/$minlat/$maxlat -G${tmpfile}=cf

	if [ ! $? -eq 0 ];then
	    echo "grdcut error."
	    echo "  in file : $jfile"
	    exit 0
	else
	    if [ -r $jfile ];then
		echo "rm $jfile"
		rm $jfile
	    fi
	fi

	if [ $j -eq 1 ];then
	    j1=`expr $j - 1`
	    suffixj1=`printf "${ifmt}" $j1`
	    jfile1=__$file.j.$suffixj1

   	    if [ "$ybegin" = "south" ]; then
  		minlat=`grdinfo $jfile1 | grep "y_max:" | sed -e "s/.* y_max:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
	    else
  		maxlat=`grdinfo $jfile1 | grep "y_min:" | sed -e "s/.* y_min:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
   	    fi

 	    echo "grdedit -A $tmpfile -R$minlon/$maxlon/$minlat/$maxlat"
 	    grdedit -A $tmpfile -R$minlon/$maxlon/$minlat/$maxlat

	    echo "1 grdpaste $OPT $tmpfile $jfile1 -G${file}=cf"
	    grdpaste $OPT $tmpfile $jfile1 -G${file}=cf

	else

  	    if [ "$ybegin" = "south" ]; then
 		minlat=`grdinfo $file | grep "y_max:" | sed -e "s/.* y_max:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
	    else
 		maxlat=`grdinfo $file | grep "y_min:" | sed -e "s/.* y_min:\s* \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
  	    fi

 	    echo "grdedit -A $tmpfile -R$minlon/$maxlon/$minlat/$maxlat"
 	    grdedit -A $tmpfile -R$minlon/$maxlon/$minlat/$maxlat

	    echo "2 grdpaste $OPT $tmpfile $file -G${file}=cf"
	    grdpaste $OPT $tmpfile $file -G${file}=cf
	fi

	if [ ! $? -eq 0 ];then
	    echo "grdpaste error."
	    echo "  in file1 : $jfile"
	    echo "  in file2 : $jfile1"
	    echo "  out file : $jfile"
	    exit 0	    
	else
	    echo "rm  $tmpfile"
	    rm  $tmpfile
	    if [ -r $jfile1 ];then
		echo "rm $jfile1"
		rm $jfile1
	    fi
	fi

	echo "grdedit -A ${file}"
	grdedit -A ${file}

    fi

    j=`expr $j + 1`

done

# revert gmt default set
# === Use GMT_TMPDIR with isolation mode. ======================================
#mv .gmtdefaults4.save  .gmtdefaults4
mv ${GMT_TMPDIR}/.gmtdefaults4.save  ${GMT_TMPDIR}/.gmtdefaults4
# ==============================================================================

done

# === Use GMT_TMPDIR with isolation mode. ======================================
rm -rf ${GMT_TMPDIR}/
# ==============================================================================
