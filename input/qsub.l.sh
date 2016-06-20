#!/bin/sh
#PBS -q L
#PBS -b 16
#PBS -l memsz_job=3gb
#PBS -l elapstim_req=00:10:00
#PBS -l filecap_job=100gb
#PBS -v OMP_NUM_THREADS=4
#PBS -m bea
#PBS -M baba.toshi@tokushima-u.ac.jp
#PBS -w SHOME="/S/data06/G6140/e0564/baba/JAGURS_SAMPLE"

#PBS -v MPIPROGINF=ALL_DETAIL

#PBS -I "${SHOME}/JAGURS-D_V0197/jagurs,ALL:./"
#PBS -I "${SHOME}/tsun.par,ALL:./"
#PBS -I "${SHOME}/gridfile.dat,ALL:./"
#PBS -I "${SHOME}/test_tgs.txt,ALL:./"
#PBS -I "${SHOME}/bathy.SD01.grd.%06j,%j:./"
#PBS -I "${SHOME}/bathy.SD02.grd.%06j,%j:./"
#PBS -I "${SHOME}/bathy.SD03.grd.%06j,%j:./"
#PBS -I "${SHOME}/bathy.SD04.grd.%06j,%j:./"
#PBS -I "${SHOME}/bathy.SD05.grd.%06j,%j:./"
#PBS -I "${SHOME}/disp.SD01.grd.%06j,%j:./"
#PBS -I "${SHOME}/disp.SD02.grd.%06j,%j:./"
#PBS -I "${SHOME}/disp.SD03.grd.%06j,%j:./"
#PBS -I "${SHOME}/disp.SD04.grd.%06j,%j:./"
#PBS -I "${SHOME}/disp.SD05.grd.%06j,%j:./"

#PBS -O "${SHOME}/output.%r/,ALL:./"

mpirun -nnp 1 ./jagurs par=tsun.par
