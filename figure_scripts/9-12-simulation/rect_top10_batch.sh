#!/bin/bash

# m3 has 22 standard memory nodes, 10 high memory nodes, and some GPU nodes.
# --Standard Memory node
#   Number of nodes: 22
#   Number of cores per node: 24
#   Processor model: 2 x (Intel Xeon CPU E5-2680 v3)
#   Processor frequency: 2.50GHz, with max Turbo frequency 3.30GHz
#   Memory per node: 128GB RAM
#   Partition: m3a
# 
# --Medium Memory node
#   Number of nodes: 13
#   Number of cores per node: 24
#   Processor model: 2 x (Intel Xeon CPU E5-2680 v3)
#   Processor frequency: 2.50GHz, with max Turbo frequency 3.30GHz
#   Memory per node: 256GB RAM
#   Partition:m3d
#        
# --High-Density CPUs node
#   Number of nodes: 52
#   Number of cores per node: 36
#   Processor model: 2 x (Intel Xeon CPU Gold 6150)
#   Processor frequency: 2.70GHz, with max Turbo frequency 3.70GHz
#   Memory per node: 192GB RAM
#   Partition: m3i
#
# --High-Density CPUs with extra high memoery node
#   Number of nodes: 1
#   Number of cores per node: 36
#   Processor model: 2 x (Intel Xeon CPU Gold 6150)
#   Processor frequency: 2.70GHz, with max Turbo frequency 3.70GHz
#   Memory per node: 1TB RAM
#   Partition: m3m
#
#
# --Intel Xeon CPU E5-2680 has 12 cores, and each node has 2 CPU. This means each node has 24 cores in total. 
# --Intel Xeon Gold 6150 has 18 cores, and each node has 2 CPU. This means each node has 36 cores in total.

#SBATCH --account=????                              	# Account
#SBATCH --job-name=rect_top10                       	# Job name
# SBATCH --qos=shortq                                	# Set Quality of Service 
#SBATCH --ntasks=1                                  	# Number of MPI ranks
# SBATCH --ntasks-per-node=1                         	# How many tasks on each CPU or socket
#SBATCH --mem-per-cpu=10GB	                    		#  Memory per processor
#SBATCH --time=4-00                             		# Time limit hrs:min:sec
#SBATCH --output=./massiveOutput/rect_top10_%j.log      # Standard output and error log
#SBATCH --partition=comp                            	# Set partition

echo $channel_id
echo $area_id

echo load matlab
module load matlab/r2019b

echo Start computing 
matlab -nodisplay -r "model_rect_fit_massive_v2($channel_id, $area_id); exit;"

echo finish

