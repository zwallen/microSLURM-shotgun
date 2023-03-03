#!/bin/bash

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# Last updated: 18 Nov 2022                                  #
#                                                            #
# Description: Create report of sequence read numbers and    #
# other data through pipeline.                               #
#                                                            #
# Usage:                                                     #
# ./Pipeline_Report.sh -o output_dir \                       #
#               -p 'commands; to; load; programs' \          #
#               -n node_partition \                          #
#               -t time:in:hours \                           #
#               -k number_of_cores \                         #
#               -m memory_per_cpu \                          #
#               -f notificationEmail@forFailures.edu         #
#                                                            #
# Parameters:                                                #
#     -h    Print the parameter list below then exit.        #
#     -o    (Required) Pipeline output directory.            #
#     -p    (Required) Single quoted string that contains    #
#           commands to load all the necessary programs      #
#           needed to run pipeline steps (e.g. activating    #
#           conda environments, loading modules, adding to   #
#           PATH, etc.).                                     #
#     -n    (Required) Name of the partition to request for  #
#           submitting jobs.                                 #
#     -t    (Required) Time request for running the jobs.    #
#           Specify in hours (e.g. 12:00:00 for 12 hours).   #
#     -k    (Required) Number of cores wanted for jobs.      #
#     -m    (Required) Amount of memory requested for each   #
#           requested core. Specify in Megabytes.            #
#     -f    (Required) E-mail to send notifications to upon  #
#           failure of any jobs.                             #
##############################################################

echo " "
echo "##############################################################"
echo "# Whole Genome Shotgun Metagenomic Processing Pipeline       #"
echo "# Last updated: 18 Nov 2022                                  #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":ho:p:n:t:k:m:f:" opt; do
  case $opt in
    h)
    echo " Description: Create report of sequence read numbers and    "
    echo " other data through pipeline.                               "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " ./Pipeline_Report.sh -o output_dir \                       "
    echo "               -p 'commands; to; load; programs' \          "
    echo "               -n node_partition \                          "
    echo "               -t time:in:hours \                           "
    echo "               -k number_of_cores \                         "
    echo "               -m memory_per_cpu \                          "
    echo "               -f notificationEmail@forFailures.edu         "
    echo "                                                            "
    echo " Parameters:                                                "
    echo "     -h    Print the parameter list below then exit.        "
    echo "     -o    (Required) Pipeline output directory.            "
    echo "     -p    (Required) Single quoted string that contains    "
    echo "           commands to load all the necessary programs      "
    echo "           needed to run pipeline steps (e.g. activating    "
    echo "           conda environments, loading modules, adding to   "
    echo "           PATH, etc.).                                     "
    echo "     -n    (Required) Name of the partition to request for  "
    echo "           submitting jobs.                                 "
    echo "     -t    (Required) Time request for running the jobs.    "
    echo "           Specify in hours (e.g. 12:00:00 for 12 hours).   "
    echo "     -k    (Required) Number of cores wanted for jobs.      "
    echo "     -m    (Required) Amount of memory requested for each   "
    echo "           requested core. Specify in Megabytes.            "
    echo "     -f    (Required) E-mail to send notifications to upon  "
    echo "           failure of any jobs.                             "
    echo " "
    exit 0
    ;;
    o) RESULTS_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    p) PROG_LOAD="$OPTARG"
    ;;
    n) PARTITION="$OPTARG"
    ;;
    t) TIME_REQUEST="$OPTARG"
    ;;
    k) CPU_REQUEST="$OPTARG"
    ;;
    m) MEM_PER_CPU="$OPTARG"
    ;;
    f) FAIL_EMAIL="$OPTARG"
    ;;
    \?) echo "Invalid option: $OPTARG" 1>&2
        exit 1
    ;;
    :) echo "Invalid option: $OPTARG requires an argument" 1>&2
       exit 1
    ;;
  esac
done

# Check that valid arguments were entered

# -o
if [[ -z "$RESULTS_DIR" ]]; then
  echo "ERROR: Argument -o is required, please supply path to pipeline output directory"
  exit 1
fi
if [[ ! -d "$RESULTS_DIR" ]]; then
  echo "ERROR: Argument -o should be a directory, please supply path to pipeline output directory"
  exit 1
fi

# -p
if [[ -z "$PROG_LOAD" ]]; then
  echo "ERROR: Argument -p is required, please supply a single quoted string of commands needed to load required programs (can be an empty string ' ' if none required)"
  exit 1
fi

# -n
if [[ -z "$PARTITION" ]]; then
  echo "ERROR: Argument -n is required, please supply a node partition name to send jobs to."
  exit 1
fi

# -t
if [[ -z "$TIME_REQUEST" ]]; then
  echo "ERROR: Argument -t is required, please supply a max length of time to run each job for."
  exit 1
fi

# -k
if [[ -z "$CPU_REQUEST" ]]; then
  echo "ERROR: Argument -k is required, please supply the number of cores being requested to run each job."
  exit 1
fi

# -m
if [[ -z "$MEM_PER_CPU" ]]; then
  echo "ERROR: Argument -m is required, please supply a memory request for each core of each job."
  exit 1
fi

# -f
if [[ -z "$FAIL_EMAIL" ]]; then
  echo "ERROR: Argument -f is required, please supply an email that can be notified upon failure of any jobs ran during the pipeline"
  exit 1
elif echo $FAIL_EMAIL | grep -q -v '@'; then
  echo "ERROR: Argument -f requires a valid email, please give an email in the form of xxxx@xxxx.xxx"
fi

###### GATHER REPORT DATA #####

SECONDS=0
echo " "
echo "*** Creating pipeline report ***"
echo " "

# Make directory for storing report components
if [ -d "${RESULTS_DIR}/Pipeline_Report" ]
then
	 :
else
	 mkdir ${RESULTS_DIR}/Pipeline_Report
   mkdir ${RESULTS_DIR}/Pipeline_Report/0.ErrorOut
   mkdir ${RESULTS_DIR}/Pipeline_Report/0.Output
fi
DIR=${RESULTS_DIR}/Pipeline_Report

#### For each relevant pipeline step, create files with: ####
#### mean length of reads per sample                     ####
#### number of sequences per sample that pass step       ####
#### adapters/primers detected and removed from reads    ####

#Prep output files
echo "Sample Input Quality_controlled Decontaminated Entropy_filtered" | sed 's/ /\t/g' > ${DIR}/mean_read_length_per_QC_step.txt
echo "Sample Input Quality_controlled Decontaminated Entropy_filtered" | sed 's/ /\t/g' > ${DIR}/remaining_reads_per_step.txt
STATS_FILE=$(ls ${RESULTS_DIR}/3.Quality_Controlled_Sequences/*_stats.txt | sed -n 1p | awk -F'/' '{print $NF}')
echo "Sample Input $(echo $(grep -v '#' ${RESULTS_DIR}/3.Quality_Controlled_Sequences/${STATS_FILE} | awk -F'\t' '{print $1}' | sed 's/ /_/g' | sort))" | sed 's/ /\t/g' > ${DIR}/detected_adapters_primers.txt

#Create script for gathering data and submit
echo '#!/bin/bash' > run.sh
echo "#SBATCH --partition=$PARTITION" >> run.sh
echo "#SBATCH --job-name=report_data" >> run.sh
echo "#SBATCH --error=${DIR}/0.ErrorOut/report_data_%A_%a.err" >> run.sh
echo "#SBATCH --output=${DIR}/0.Output/report_data_%A_%a.out" >> run.sh
echo "#SBATCH --time=$TIME_REQUEST" >> run.sh
echo "#SBATCH --ntasks=1" >> run.sh
echo "#SBATCH --cpus-per-task=$CPU_REQUEST" >> run.sh
echo "#SBATCH --mem-per-cpu=$MEM_PER_CPU" >> run.sh
echo "#SBATCH --mail-type=FAIL" >> run.sh
echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*fastq.gz | wc -l)" >> run.sh
echo "#SBATCH --wait" >> run.sh
echo "$PROG_LOAD" >> run.sh
echo " " >> run.sh
echo "SAMPLE=\$(ls -l ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*fastq.gz | awk '{print \$NF}' | awk -F'/' '{print \$NF}' | awk -F'.fastq.gz' '{print \$1}' | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
echo " " >> run.sh
echo "unzip -qq ${RESULTS_DIR}/1.FastQC_Initial_Reports/\${SAMPLE}_R1_fastqc.zip" >> run.sh
echo "sed -n '/>>Sequence Length Distribution/,/>>END_MODULE/p' \${SAMPLE}_R1_fastqc/fastqc_data.txt | sed '1,2d' | sed '\$d' > ${DIR}/\${SAMPLE}_lengths.txt" >> run.sh
echo "echo \"lens <- read.table('${DIR}/\${SAMPLE}_lengths.txt')\" > ${DIR}/\${SAMPLE}_Rfunc.R" >> run.sh
echo "echo \"cat(paste(round(mean(rep(lens[,1],lens[,2])),1), '±', ifelse(is.na(round(sd(lens[,1],lens[,2]),1)), 0, round(sd(rep(lens[,1],lens[,2])),1)), sep=''), '\n')\" >> ${DIR}/\${SAMPLE}_Rfunc.R" >> run.sh
echo "INPUT_LEN=\$(Rscript --vanilla ${DIR}/\${SAMPLE}_Rfunc.R)" >> run.sh
echo "rm -r \${SAMPLE}_R1_fastqc" >> run.sh
echo " " >> run.sh
echo "echo \"lens <- read.table('${DIR}/\${SAMPLE}_lengths.txt')\" > ${DIR}/\${SAMPLE}_Rfunc.R" >> run.sh
echo "echo \"cat(paste(round(mean(as.numeric(lens[,1])),1), '±', ifelse(is.na(round(sd(lens[,1]),1)), 0, round(sd(lens[,1]),1)), sep=''), '\n')\" >> ${DIR}/\${SAMPLE}_Rfunc.R" >> run.sh
echo "zcat ${RESULTS_DIR}/3.Quality_Controlled_Sequences/\${SAMPLE}_R1.fastq.gz | awk '{if(NR%4==2) print length(\$1)}' > ${DIR}/\${SAMPLE}_lengths.txt" >> run.sh
echo "QC_LEN=\$(Rscript --vanilla ${DIR}/\${SAMPLE}_Rfunc.R)" >> run.sh
echo " " >> run.sh
echo "zcat ${RESULTS_DIR}/4.Decontaminated_Sequences/\${SAMPLE}_R1.fastq.gz | awk '{if(NR%4==2) print length(\$1)}' > ${DIR}/\${SAMPLE}_lengths.txt" >> run.sh
echo "DECONTAM_LEN=\$(Rscript --vanilla ${DIR}/\${SAMPLE}_Rfunc.R)" >> run.sh
echo " " >> run.sh
echo "zcat ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/\${SAMPLE}.fastq.gz | awk '{if(NR%4==2) print length(\$1)}' > ${DIR}/\${SAMPLE}_lengths.txt" >> run.sh
echo "ENTROPY_LEN=\$(Rscript --vanilla ${DIR}/\${SAMPLE}_Rfunc.R)" >> run.sh
echo " " >> run.sh
echo "echo \"\$SAMPLE \$INPUT_LEN \$QC_LEN \$DECONTAM_LEN \$ENTROPY_LEN\" | sed 's/ /\t/g' >> ${DIR}/mean_read_length_per_QC_step.txt" >> run.sh
echo "rm ${DIR}/\${SAMPLE}_Rfunc.R" >> run.sh
echo "rm ${DIR}/\${SAMPLE}_lengths.txt" >> run.sh
echo " " >> run.sh
echo "BEGIN_SEQ=\$(grep 'Input:' ${RESULTS_DIR}/3.Quality_Controlled_Sequences/\${SAMPLE}.log | awk '{print \$2}')" >> run.sh
echo "QC_SEQ=\$(grep 'Result:' ${RESULTS_DIR}/3.Quality_Controlled_Sequences/\${SAMPLE}.log | awk '{print \$2}')" >> run.sh
echo "DECONTAM_SEQ=\$(( \$QC_SEQ - \$(( \$(zcat ${RESULTS_DIR}/4.Decontaminated_Sequences/Extracted_Host_Sequences/\${SAMPLE}.GCA_000001405.28_GRCh38.p13_genomic_contam_1.fastq.gz | wc -l) / 4 * 2 )) ))" >> run.sh
echo "ENTROPY_SEQ=\$(grep 'Result:' ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/\${SAMPLE}.log | awk '{print \$2}')" >> run.sh
echo "echo \"\$SAMPLE \$BEGIN_SEQ \$QC_SEQ \$DECONTAM_SEQ \$ENTROPY_SEQ\" | sed 's/ /\t/g' >> ${DIR}/remaining_reads_per_step.txt" >> run.sh
echo " " >> run.sh
echo "echo \"\$SAMPLE \$BEGIN_SEQ \$(echo \$(grep -v '#' ${RESULTS_DIR}/3.Quality_Controlled_Sequences/\${SAMPLE}_stats.txt | sort -k1,1 | awk -F'\t' '{print \$2}' | sed 's/ /_/g'))\" | sed 's/ /\t/g' >> ${DIR}/detected_adapters_primers.txt" >> run.sh
chmod +x run.sh

sbatch run.sh > /dev/null

rm run.sh

#Sort files by sample ID
awk 'NR > 1' ${DIR}/mean_read_length_per_QC_step.txt | sort -k1,1 | cat <(awk 'NR == 1' ${DIR}/mean_read_length_per_QC_step.txt) - > ${DIR}/temp.txt
mv ${DIR}/temp.txt ${DIR}/mean_read_length_per_QC_step.txt
awk 'NR > 1' ${DIR}/remaining_reads_per_step.txt | sort -k1,1 | cat <(awk 'NR == 1' ${DIR}/remaining_reads_per_step.txt) - > ${DIR}/temp.txt
mv ${DIR}/temp.txt ${DIR}/remaining_reads_per_step.txt
awk 'NR > 1' ${DIR}/detected_adapters_primers.txt | sort -k1,1 | cat <(awk 'NR == 1' ${DIR}/detected_adapters_primers.txt) - > ${DIR}/temp.txt
mv ${DIR}/temp.txt ${DIR}/detected_adapters_primers.txt

#Signal end of report generation
echo " "
echo "Pipeline report complete"
echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
echo " "
#################################################
