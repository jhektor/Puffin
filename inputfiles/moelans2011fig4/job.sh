#!/bin/sh
# Number of cores on exclusive nodes
#SBATCH --exclusive
#SBATCH -N 4
#SBATCH --tasks-per-node=20

# Walltime hh:mm:ss
#SBATCH -t 01:00:00

# Name the job
#SBATCH -J moelansfig4-2d-test

# Specify outfile and errorfile
#SBATCH -o moelansfig4-2d-test_%j.out
#SBATCH -e moelansfig4-2d-test_%j.err

# Mail when job is finished
#SBATCH --mail-user=johan.hektor@solid.lth.se
#SBATCH --mail-type=ALL


#To run short tests (maximum 1h), it is possible to request extra high priority on Aurora with the help of
#SBATCH --qos=test

# write this script to stdout-file - useful for scripting errors
cat $0

# Load modules
ml load GCC/6.3.0-2.27
ml load OpenMPI/2.0.2
ml load Moose_framework

# Start program
mpirun -bind-to-core ../../puffin-opt -i imcsphere2d.i
