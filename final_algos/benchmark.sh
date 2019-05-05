#!/bin/bash
#SBATCH -p moh1 # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 12000 # memory pool for all cores
#SBATCH -t 1-00:00 # time (D-HH:MM)
#SBATCH -o test.out # STDOUT
#SBATCH -e test.err # STDERR
./treehash $1 $2 $3 $8 $9 ${10} ${11}
./minhash $1 $4 $5 $8 $9 ${10} ${11}
./lsh $1 $6 $7 $8 $9 ${10} ${11}

