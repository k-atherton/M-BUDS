#!/bin/bash -l

# Set SCC project
#$ -P microbiome

# Send an email when the job finishes or if it is aborted (by default no email is sent).
#$ -m beas

# Give job a name
#$ -N format_raw_dada2_asv_tables

# Combine output and error files into a single file
#$ -j y

# Specify the output file name
#$ -o Bash_Outputs/atherton_format_raw_dada2_asv_tables_bash_output.txt

# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "=========================================================="

echo "Running on 16S data"
module load R
Rscript 02_format_raw_dada2_asv_tables.R \
	-a "16S" \
	-n "atherton" \
	-e "N"

echo "Running on ITS data"
Rscript 02_format_raw_dada2_asv_tables.R \
	-a "ITS" \
	-n "atherton" \
	-e "N"

echo "Done."
