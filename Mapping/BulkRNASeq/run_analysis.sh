#!/bin/bash
#SBATCH -n 20                               # Request one core
#SBATCH -N 1                               # Request one node (if you request more than one core with -n, also using
                                           # -N 1 means all cores will be on the same node)
#SBATCH -t 0-12:00                         # Runtime in D-HH:MM format
#SBATCH -p short                           # Partition to run in
#SBATCH --mem=70G                          # Memory total in MB (for all cores)
#SBATCH -o hostname_%j.out                 # File to which STDOUT will be written, including job ID
#SBATCH -e hostname_%j.err                 # File to which STDERR will be written, including job ID
#SBATCH --mail-type=ALL                    # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=MohamedNadhir.Djekidel@childrens.harvard.edu   # Email to which notifications will be sent

source ~/.bash_profile
module load gcc/6.2.0
module load star/2.5.2b
module load python/2.7.12
module load samtools/1.3.1
module load trimmomatic/0.36
module load java/jdk-1.8u112
module load R/3.4.1
module load picard/2.8.0
module load deeptools/3.0.0
module load fastqc/0.11.5

source /home/nadhir/jupytervenv/bin/activate

bpipe run -r RNASeq_workflow_SE.groovy

