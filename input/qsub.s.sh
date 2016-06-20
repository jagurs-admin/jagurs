#!/bin/sh
#PBS -q S
#PBS -b 4
#PBS -l memsz_job=3gb
#PBS -l elapstim_req=00:10:00
#PBS -l filecap_job=100gb
#PBS -v OMP_NUM_THREADS=4
#PBS -m bea
#PBS -M baba.toshi@tokushima-u.ac.jp

cd /S/data06/G6140/e0564/baba/JAGURS_SAMPLE
mpirun -nnp 1 ./jagurs par=tsun.par
