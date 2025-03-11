#!/bin/bash -l

# Set SCC project
#$ -P talbot-lab-data

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m ae

# Give job a name
#$ -N plot_networks

# Combine output and error files into a single file
#$ -j y

# Specify the output file name
#$ -o plot_networks.qlog

# Specify number of cores
#$ -pe omp 8

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $SGE_TASK_ID"
echo "=========================================================="

module load R/4.3.1
Rscript 04_plot_networks.R