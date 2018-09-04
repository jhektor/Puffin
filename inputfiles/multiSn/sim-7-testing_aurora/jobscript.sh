
#	shell spec.
#!/bin/sh
#	walltime, i.e. how long the job will be on the node, NOT CPU TIME
#SBATCH -t 04:00:00
#	Job name:
#SBATCH -J testAuroraOut
#	Naming output files, %j adds the job number assigned
#SBATCH -o testOutputName_%j.out
#SBATCH -e testOutputName_%j.err
#	Exclusive node access, allocate whole nodes, 20 cores, for the job. 
#SBATCH --exclusive
#	Specify the number of nodes allocated for this job, in this case 4 nodes
#SBATCH	-N 4
#	Specify the number of tasks required per node,commonly equals the number of cores available
#SBATCH --tasks-per-node=20

