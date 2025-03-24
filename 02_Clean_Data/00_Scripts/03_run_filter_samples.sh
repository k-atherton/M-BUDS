#!/bin/bash -l
# Set SCC project
#$ -P microbiome

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m beas

# Give job a name
#$ -N filter_samples

# Combine output and error files into a single file
#$ -j y

# Specify the output file name
#$ -o Bash_Outputs/atherton_filter_samples_bash_output.txt

# Request more than one core
#$ -pe omp 8 #REQUESTS 8 CORES, YOU CAN ONLY REQUEST 1, 4, 8, 16, 28, OR 36

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "=========================================================="

echo "Running on 16S Samples"
module load R
Rscript 03_run_filter_samples.R \
	-a "16S" \
	-n "atherton" \
	-e "N"

echo "Running on ITS samples"
Rscript 03_run_filter_samples.R \
	-a "ITS" \
	-n "atherton" \
	-e "N"

echo "Done."
