#!/bin/sh

# Number of nodes
#SBATCH -N 1
# Number of processor core (32*32=1024, psfc, mit and emiliob nodes have 32 cores per node)
#SBATCH -n 4
# specify how long your job needs. Be HONEST, it affects how long the job may wait for its turn.
#SBATCH --time=12:00:00
# which partition or queue the jobs runs in
#SBATCH -p sched_mit_nse


module load python/3.6.3
python3 run_erl.py