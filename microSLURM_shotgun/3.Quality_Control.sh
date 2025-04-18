#!/bin/bash

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# Last updated: 18 Nov 2022                                  #
#                                                            #
# Description: Remove adapters, phix sequences using BBDuk.  #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
#    BBDuk:      For adapter and quality trimming of raw wgs #
#                reads, removal of PhiX sequences, and       #
#                removing low complexity sequences.          #
#                                                            #
# Usage:                                                     #
# ./3.Quality_Control.sh [-x] \                              #
#                        -i input_seqs_dir \                 #
#                        -o output_dir \                     #
#                        -p 'commands; to; load; programs' \ #
#                        -n node_partition \                 #
#                        -t time:in:hours \                  #
#                        -k number_of_cores \                #
#                        -m memory_per_cpu \                 #
#                        -f notificationEmail@forFailures.edu#
#                                                            #
# Parameters:                                                #
#     -h    Print the parameter list below then exit.        #
#     -x    (Required) Are input fastq files merged? Add this#
#           parameter if so. Will cause error otherwise.     #
#     -i    (Required) Directory that contains the raw       #
#           fastq files to be processed.                     #
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
while getopts ":hxi:o:p:n:t:k:m:f:" opt; do
  case $opt in
    h)
    echo " Description: Remove adapters, phix sequences using BBDuk.  "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
    echo "    BBDuk:      For adapter and quality trimming of raw wgs "
    echo "                reads, removal of PhiX sequences, and       "
    echo "                removing low complexity sequences.          "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " ./3.Quality_Control.sh [-x] \                              "
    echo "                        -i input_seqs_dir \                 "
    echo "                        -o output_dir \                     "
    echo "                        -p 'commands; to; load; programs' \ "
    echo "                        -n node_partition \                 "
    echo "                        -t time:in:hours \                  "
    echo "                        -k number_of_cores \                "
    echo "                        -m memory_per_cpu \                 "
    echo "                        -f notificationEmail@forFailures.edu"
    echo "                                                            "
    echo " Parameters:                                                "
    echo "     -h    Print the parameter list below then exit.        "
    echo "     -x    (Required) Are input fastq files merged? Add this"
    echo "           parameter if so. Will cause error otherwise.     "
    echo "     -i    (Required) Directory that contains the raw       "
    echo "           fastq files to be processed.                     "
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
    i) SEQ_DIR=$(echo $OPTARG | sed 's#/$##')
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

# -i
if [[ -z "$MERGE" ]]; then
  if [[ -z "$SEQ_DIR" ]]; then
    echo "ERROR: Argument -i is required, please supply a directory with input fastq files"
    exit 1
  fi
  if [[ ! -d "$SEQ_DIR" ]]; then
    echo "ERROR: Argument -i should be a directory, please supply a directory with input fastq files"
    exit 1
  fi
  if ls -l $SEQ_DIR | grep -q ".fastq.gz"; then
    SEQ_EXT=fastq.gz
  elif ls -l $SEQ_DIR | grep -q ".fastq"; then
    SEQ_EXT=fastq
  elif ls -l $SEQ_DIR | grep -q ".fq.gz"; then
    SEQ_EXT=fq.gz
  elif ls -l $SEQ_DIR | grep -q ".fq"; then
    SEQ_EXT=fq
  else
    echo "ERROR: Sequences in input directory should have file extension of either .fastq[.gz] OR .fq[.gz]"
    exit 1
  fi
  found_1=$(ls -l $SEQ_DIR | grep -q "_R1")
  found_2=$(ls -l $SEQ_DIR | grep -q "_R2")
  if [[ -n "$found_1" ]] && [[ -n "$found_2" ]]; then
    echo "ERROR: Sequences in input directory expected to be paired-end Illumina sequence fastq files whose file names contain the strings '_R1' and '_R2'"
    exit 1
  else
    :
  fi
fi

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

# get max memory to use for jobs
MAX_MEM=$(( MEM_PER_CPU * CPU_REQUEST / 1000))

################# QC WITH BBDUK #################

  SECONDS=0
  echo " "
  echo "*** Running BBDuk to perform adapter/quality trimming and filtering on input fastq files ***"
  echo " "

  #Create directory for output
  if [ -d "${RESULTS_DIR}/3.Quality_Controlled_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/3.Quality_Controlled_Sequences
	  mkdir ${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.Output
  fi

  ##### Run BBDuk #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION" >> run.sh
  echo "#SBATCH --job-name=QC" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.ErrorOut/QC_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.Output/QC_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_REQUEST" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_REQUEST" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_PER_CPU" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/2.Merged_Paired_End_Sequences/*.fastq.gz | wc -l)" >> run.sh
  else
    echo "#SBATCH --array=1-$(ls -l ${SEQ_DIR}/*_R1*${SEQ_EXT} | wc -l)" >> run.sh
  fi
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "FILE=\$(ls ${RESULTS_DIR}/2.Merged_Paired_End_Sequences/*.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$FILE | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> run.sh
    echo "bbduk.sh in=\$FILE \\" >> run.sh
    echo "out=${RESULTS_DIR}/3.Quality_Controlled_Sequences/\${FILE_NAME}.fastq.gz \\" >> run.sh
    echo "stats=${RESULTS_DIR}/3.Quality_Controlled_Sequences/\${FILE_NAME}_stats.txt \\" >> run.sh
    echo "ftm=5 qtrim=rl trimq=25 minlen=50 ref=adapters,phix -Xmx${MAX_MEM}g \\" >> run.sh
  else
    echo "FILE1=\$(ls ${SEQ_DIR}/*_R1*${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE2=\$(ls ${SEQ_DIR}/*_R2*${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1' '{print \$1}')" >> run.sh
    echo "bbduk.sh in=\$FILE1 \\" >> run.sh
    echo "in2=\$FILE2 \\" >> run.sh
    echo "out=${RESULTS_DIR}/3.Quality_Controlled_Sequences/\${FILE_NAME}_R1.fastq.gz \\" >> run.sh
    echo "out2=${RESULTS_DIR}/3.Quality_Controlled_Sequences/\${FILE_NAME}_R2.fastq.gz \\" >> run.sh
    echo "stats=${RESULTS_DIR}/3.Quality_Controlled_Sequences/\${FILE_NAME}_stats.txt \\" >> run.sh
    echo "ftm=5 tpe tbo qtrim=rl trimq=25 minlen=50 ref=adapters,phix -Xmx${MAX_MEM}g \\" >> run.sh
  fi
  echo "> ${RESULTS_DIR}/3.Quality_Controlled_Sequences/\${FILE_NAME}.log 2>&1"  >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  #Signal jobs have ended
  echo "Quality control with BBDuk complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "

#################################################
