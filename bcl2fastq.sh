#!/bin/sh
# Move the sequencing output on a folder "foldername", that will be the $1
# sh Bcl2cellranger.sh path/to/foldername runid ScRNAseq1 /path/to/humanreference ScRNAseq2 /path/to/mousereference
#
# In path/to/foldername, create the needed files: index.csv + taghto.csv + tagadt.csv  or index_illu.csv (if not scRNAseq) or index.csv + libraries.csv + featureindex.csv (if using cellranger 4)

runid=$2

module load bcl2fastq2
echo "bcl2fastq --runfolder-dir $1 --output-dir $1/$runid  --sample-sheet=$1/index_illu.csv --mask-short-adapter-reads 0 " | qsub -V -cwd -q max-7d.q
# For some experiments, we had to use --mask-short-adapter-reads 0 or we would get "NNNNN.." in the R1 fastqs ( when multiple library types were present)
