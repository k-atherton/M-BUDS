#!/bin/bash -l

# Set SCC project
#$ -P microbiome

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m beas

# Give job a name
#$ -N decontaminate_samples

# Combine output and error files into a single file
#$ -j y

# Specify the output file name
#$ -o Bash_Outputs/atherton_decontaminate_samples_bash_output.txt

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "=========================================================="

echo "Running on 16S Samples"
module load R
Rscript 04_decontaminate_samples.R \
	-a "16S" \
	-n "atherton" \
	-e "N"

echo "Running on ITS Samples"
Rscript 04_decontaminate_samples.R \
	-a "ITS" \
	-n "atherton" \
	-e "N"

echo "Done."
