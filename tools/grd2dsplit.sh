#!/bin/bash


files=
px=1
py=1
dx=1
dy=1
eps=0.0000000000001
ffmt="%.10f"
ifmt="%06d"
xbegin="west"
ybegin="north"
# === Error handling for illegal division!!! ===================================
min_bsize="3"
# ==============================================================================

usage()
{
   echo "GRD File 2D Splitter"
   echo ""
   echo "$0 [-px <proc-x>] [-py <proc-y>] [-dx <edge-x>] [-dy <edge-y>] [-eps <eps>] [-ffmt <fmt>] [-ifmt <fmt>] grdfile [grdfile...]"
   echo ""
   echo "  -px <proc-x> : number of processes on coordinate x (default:$px)"
   echo "  -py <proc-y> : number of processes on coordinate y (default:$py)"
   echo "  -dx <edge-x> : overlap edge size on coordinate x (default:$dx)"
   echo "  -dy <edge-y> : overlap edge size on coordinate y (default:$dy)"
   echo "  -eps <eps>   : epsilon of lon/lat value (default:$eps)"
   echo "  -ffmt <fmt>  : format string for float numbers (default:$ffmt)"
   echo "  -ifmt <fmt>  : format string for int   numbers (default:$ifmt)"
   echo "  -xbegin <we> : <we>={west|east|w|e}   (default:$xbegin)"
   echo "  -ybegin <ns> : <ns>={north|south|n|s} (default:$ybegin)"
   echo "  grdfile      : .grd file"
   echo ""
   echo "  using commands  : grdinfo, grfcut"
   echo ""
   echo "  output filename : grdfile.000000, grdfile.000001, ..."
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
     *   ) files="$files $1"; shift 1;;
   esac
done

if [ "$files" = "" ]; then
   usage
   exit 0
fi

echo "process $px, $py"
np=`expr $px \* $py`
for file in $files; do

   x0=`grdinfo $file | grep "x_min:" | sed -e "s/.* x_min: \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
   y0=`grdinfo $file | grep "y_min:" | sed -e "s/.* y_min: \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
   x1=`grdinfo $file | grep "x_max:" | sed -e "s/.* x_max: \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
   y1=`grdinfo $file | grep "y_max:" | sed -e "s/.* y_max: \([0-9\-][0-9\.]*\)[^0-9].*/\1/g" `
   nx=`grdinfo $file | grep "nx:" | sed -e "s/.* nx: \([0-9][0-9]*\)/\1/g" `
   ny=`grdinfo $file | grep "ny:" | sed -e "s/.* ny: \([0-9][0-9]*\)/\1/g" `
# ==============================================================================
#  x_inc=`grdinfo $file | grep "x_inc:" | sed -e "s/.* x_inc: \([0-9][0-9\.]*\)[^0-9].*/\1/g" `
#  y_inc=`grdinfo $file | grep "y_inc:" | sed -e "s/.* y_inc: \([0-9][0-9\.]*\)[^0-9].*/\1/g" `
   x_inc=`grdinfo $file | grep "x_inc:" | sed -e "s/.* x_inc: \([0-9][-+e0-9\.]*\)[^-+e0-9].*/\1/" `
   y_inc=`grdinfo $file | grep "y_inc:" | sed -e "s/.* y_inc: \([0-9][-+e0-9\.]*\)[^-+e0-9].*/\1/" `
   x_inc=`echo $x_inc | sed s/e/*10^/`
   y_inc=`echo $y_inc | sed s/e/*10^/`
# ==============================================================================
   echo "start coord. $x0/$x1/$y0/$y1"
   echo "grid size    $nx, $ny"
   echo "grid disp.   $x_inc, $y_inc"

   nbx=`expr $nx / $px + \( $nx % $px \> 0 \)`
   nby=`expr $ny / $py + \( $ny % $py \> 0 \)`
   echo "block size   $nbx, $nby"
# === Error handling for illegal division!!! ===================================
   last_nbx=`expr $nx - $nbx \* \( $px - 1 \)`
   last_nby=`expr $ny - $nby \* \( $py - 1 \)`
   echo "############################################################"
   echo "# block size of last rank `expr $px \* $py - 1`: $last_nbx x $last_nby"
   if [ "$last_nbx" -lt "$min_bsize" ]; then
      echo "# ERROR!!! Num. of procs. on coordinate x is illegal!"
      echo "#          Block size must be at least $min_bsize."
      echo "#          Please change px!"
      echo "############################################################"
      exit 1
   fi
   if [ "$last_nby" -lt "$min_bsize" ]; then
      echo "# ERROR!!! Num. of procs. on coordinate y is illegal!"
      echo "#          Block size must be at least $min_bsize."
      echo "#          Please change py!"
      echo "############################################################"
      exit 1
   fi
   echo "############################################################"
# ==============================================================================

   k=0
   j=0
   while [ $j -lt $py ]; do
     if [ "$ybegin" = "south" ]; then
       nymin=`expr $j \* $nby - $dy`
       nymax=`expr \( $j + 1 \) \* $nby - 1 + $dy`
     else
       nymax=`expr $ny - \( $j     \) \* $nby + $dy - 1`
       nymin=`expr $ny - \( $j + 1 \) \* $nby - $dy`
     fi
     #echo "y $nymin:$nymax"

     if [ $nymin -lt 0 ]; then
#        nymin=0
        ys=`echo "$y0 + $eps" | bc -l`
        minlat=`printf "${ffmt}" $ys`
     else
        ys=`echo "$y0 + $nymin * $y_inc" | bc -l`
        minlat=`printf "${ffmt}" $ys`
     fi
     if [ $nymax -ge $ny ]; then
        ye=`echo "$y0 + ( $ny - 1 ) * $y_inc - $eps" | bc -l`
        maxlat=`printf "${ffmt}" $ye`
     else
        ye=`echo "$y0 + $nymax * $y_inc" | bc -l`
        maxlat=`printf "${ffmt}" $ye`
     fi
#     minlat=`echo "$y0 + $nymin * $y_inc" | bc -l`
#     maxlat=`echo "$y0 + $nymax * $y_inc" | bc -l`

     i=0
     while [ $i -lt $px ]; do
       if [ "$xbegin" = "east" ]; then
         nxmax=`expr $nx - \( $i     \) \* $nbx + $dx - 1`
         nxmin=`expr $nx - \( $i + 1 \) \* $nbx - $dx`
       else
         nxmin=`expr $i \* $nbx - $dx`
         nxmax=`expr \( $i + 1 \) \* $nbx - 1 + $dx`
       fi
#       echo "x $nxmin:$nxmax"

       if [ $nxmin -lt 0 ]; then
#          nxmin=0
          xs=`echo "$x0 + $eps" | bc -l`
          minlon=`printf "${ffmt}" $xs`
       else
          xs=`echo "$x0 + $nxmin * $x_inc" | bc -l`
          minlon=`printf "${ffmt}" $xs`
       fi
       if [ $nxmax -ge $nx ]; then
          xe=`echo "$x0 + ( $nx - 1 ) * $x_inc - $eps" | bc -l`
          maxlon=`printf "${ffmt}" $xe`
       else
          xe=`echo "$x0 + $nxmax * $x_inc" | bc -l`
          maxlon=`printf "${ffmt}" $xe`
       fi
#       minlon=`echo "$x0 + $nxmin * $x_inc" | bc -l`
#       maxlon=`echo "$x0 + $nxmax * $x_inc" | bc -l`

       suffix=`printf "${ifmt}" $k`
       ofile=`basename $file`.$suffix

#       lb=`echo "( $nxmax - $nxmin - 0.0001 ) * $x_inc " | bc -l`
#       ub=`echo "( $nxmax - $nxmin + 0.0001 ) * $x_inc " | bc -l`
#       bx=`echo "( $maxlon - $minlon )" | bc -l`
#       echo "check range : $lb < $bx < $ub"
#       echo "grdcut $file -R$minlat/$maxlat/$minlon/$maxlon -G$file.$suffix"
#       grdcut $file -R$minlat/$maxlat/$minlon/$maxlon -G$file.$suffix
#       echo "grdcut -V $file -R$minlon/$maxlon/$minlat/$maxlat -G`basename $file`.$suffix=cf"
#       grdcut -V $file -R$minlon/$maxlon/$minlat/$maxlat -G`basename $file`.$suffix=cf
       echo "grdcut -V $file -R$xs/$xe/$ys/$ye -G${ofile}=cf"
       grdcut -V $file -R$xs/$xe/$ys/$ye -G${ofile}=cf
       echo "grdedit $ofile -R$minlon/$maxlon/$minlat/$maxlat"
       grdedit $ofile -R$minlon/$maxlon/$minlat/$maxlat

       i=`expr $i + 1`
       k=`expr $k + 1`
     done 
     j=`expr $j + 1`
   done 
done
