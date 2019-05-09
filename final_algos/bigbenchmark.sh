#!/bin/bash
#SBATCH -p moh1 # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 10 # number of cores
#SBATCH --mem 12000 # memory pool for all cores
#SBATCH -t 1-00:00 # time (D-HH:MM)
#SBATCH -o test.out # STDOUT
#SBATCH -e test.err # STDERR
echo "Running N = 500"
./benchmark.sh "500" $1 $2 $1 $3 $1 $4 > "benches1/N500_b$1.txt"
echo "Running N = 1000"
./benchmark.sh "1000" $1 $2 $1 $3 $1 $4 > "benches1/N1000_b$1.txt"
echo "Running N = 2000"
./benchmark.sh "2000" $1 $2 $1 $3 $1 $4 > "benches1/N2000_b$1.txt"
echo "Running N = 5000"
./benchmark.sh "5000" $1 $2 $1 $3 $1 $4 > "benches1/N5000_b$1.txt"
echo "Running N = 10000"
./benchmark.sh "10000" $1 $2 $1 $3 $1 $4 > "benches1/N10000_b$1.txt"
echo "Running N = 20000"
./benchmark.sh "20000" $1 $2 $1 $3 $1 $4 > "benches1/N20000_b$1.txt"
echo "Running N = 50000"
./benchmark.sh "50000" $1 $2 $1 $3 $1 $4 > "benches1/N50000_b$1.txt"
echo "Running N = 100000"
./benchmark.sh "100000" $1 $2 $1 $3 $1 $4 > "benches1/N100000_b$1.txt"
