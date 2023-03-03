#!/bin/bash

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# Last updated: 18 Nov 2022                                  #
#                                                            #
# Description: Generate FastQC reports for each sequence     #
# read file.                                                 #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
#    FastQC:     For performing quality reports.             #
#    MultiQC:    For summarizing multiple FastQC reports.    #
#                                                            #
# Usage:                                                     #
# ./6.Final_FastQC.sh [-x] \                                 #
#               -o output_dir \                              #
#               -p 'commands; to; load; programs' \          #
#               -n node_partition \                          #
#               -t time:in:hours \                           #
#               -k number_of_cores \                         #
#               -m memory_per_cpu \                          #
#               -f notificationEmail@forFailures.edu         #
#                                                            #
# Parameters:                                                #
#     -h    Print the parameter list below then exit.        #
#     -x    (Required) Are input fastq files merged? Add this#
#           parameter if so. Will cause error otherwise.     #
#           Also add if fastq files were joined before low-  #
#           complexity filtering was performed.              #
#     -o    (Required) Path to pipeline result directory.    #
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
while getopts ":hxo:p:n:t:k:m:f:" opt; do
  case $opt in
    h)
    echo " Description: Generate FastQC reports for each sequence     "
    echo " read file.                                                 "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
    echo "    FastQC:     For performing quality reports.             "
    echo "    MultiQC:    For summarizing multiple FastQC reports.    "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " ./6.Final_FastQC.sh [-x] \                                 "
    echo "               -o output_dir \                              "
    echo "               -p 'commands; to; load; programs' \          "
    echo "               -n node_partition \                          "
    echo "               -t time:in:hours \                           "
    echo "               -k number_of_cores \                         "
    echo "               -m memory_per_cpu \                          "
    echo "               -f notificationEmail@forFailures.edu         "
    echo "                                                            "
    echo " Parameters:                                                "
    echo "     -h    Print the parameter list below then exit.        "
    echo "     -x    (Required) Are input fastq files merged? Add this"
    echo "           parameter if so. Will cause error otherwise.     "
    echo "           Also add if fastq files were joined before low-  "
    echo "           complexity filtering was performed.              "
    echo "     -o    (Required) Path to pipeline result directory.    "
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
    x) MERGE=1
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
  echo "ERROR: Argument -o is required, please supply an output directory"
  exit 1
fi
if [[ ! -d "$RESULTS_DIR" ]]; then
  echo "ERROR: Argument -o should be a directory, please supply an output directory"
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

#Determine if pipeline result directory has been reorganized and edit input/output directories accordingly
if [[ -d "${RESULTS_DIR}/2.Processed_Sequences" ]]; then
  IN_DIR=2.Processed_Sequences
  OUT_DIR=3.FastQC_Final_Reports
else
  IN_DIR=5.Low_Complexity_Filtered_Sequences
  OUT_DIR=6.FastQC_Final_Reports
fi

############# FINAL FASTQC REPORT #############
  SECONDS=0
  echo " "
  echo "*** Running final FastQC report on QCed fastq files ***"
  echo " "

  #Create directory for output
  if [ -d "${RESULTS_DIR}/${OUT_DIR}" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/${OUT_DIR}
	  mkdir ${RESULTS_DIR}/${OUT_DIR}/0.ErrorOut
	  mkdir ${RESULTS_DIR}/${OUT_DIR}/0.Output
  fi

  ##### Run FastQC #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION" >> run.sh
  echo "#SBATCH --job-name=FastQC" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/${OUT_DIR}/0.ErrorOut/FastQC_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/${OUT_DIR}/0.Output/FastQC_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_REQUEST" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_REQUEST" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_PER_CPU" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
    if [[ ! -z "$MERGE" ]]; then
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/${IN_DIR}/*.fastq.gz | wc -l)" >> run.sh
  else
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/${IN_DIR}/*R1.fastq.gz | wc -l)" >> run.sh
  fi
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "FILE=\$(ls ${RESULTS_DIR}/${IN_DIR}/*.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$FILE | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> run.sh
    echo "fastqc \$FILE -d ${RESULTS_DIR}/${OUT_DIR} -o ${RESULTS_DIR}/${OUT_DIR} \\" >> run.sh
  else
    echo "FILE1=\$(ls ${RESULTS_DIR}/${IN_DIR}/*R1.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE2=\$(ls ${RESULTS_DIR}/${IN_DIR}/*R2.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1' '{print \$1}')" >> run.sh
    echo "fastqc \$FILE1 \$FILE2 -d ${RESULTS_DIR}/${OUT_DIR} -o ${RESULTS_DIR}/${OUT_DIR} \\" >> run.sh
  fi
  echo "> ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}.log 2>&1" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION" >> run.sh
  echo "#SBATCH --job-name=MultiQC" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/${OUT_DIR}/0.ErrorOut/MultiQC.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/${OUT_DIR}/0.Output/MultiQC.out" >> run.sh
  echo "#SBATCH --time=$TIME_REQUEST" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_REQUEST" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_PER_CPU" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "multiqc ${RESULTS_DIR}/${OUT_DIR} -o ${RESULTS_DIR}/${OUT_DIR} \\" >> run.sh
  echo "> ${RESULTS_DIR}/${OUT_DIR}/multiqc.log 2>&1" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  #Signal jobs have ended
  echo "Final FastQC reports complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
#################################################
