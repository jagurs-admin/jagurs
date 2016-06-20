#!/bin/bash
#BSUB -n 128
#BSUB -W 360
#BSUB -R rusage[mem=512]
#BSUB -a ICE
#BSUB -J jagurs
#BSUB -o stdout
#BSUB -e stderr

export OMP_NUM_THREADS=16
hybrid -mpi 8 -omp 16 -mpipn 1 "./JAGURS-D_V0174/jagurs par=tsun.par"
