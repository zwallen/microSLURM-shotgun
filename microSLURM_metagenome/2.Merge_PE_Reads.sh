#!/bin/bash

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# Last updated: 18 Nov 2022                                  #
#                                                            #
# Description: Merge paired-end reads using BBMerge.         #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
#    BBMerge:    For merging paired-end reads.               #
#                                                            #
# Usage:                                                     #
# ./2.Merge_PE_Reads.sh -o output_dir \                      #
#                       -s \                                 #
#                     OR                                     #
#                       -i input_seqs_dir \                  #
#                       -o output_dir \                      #
#                       -p 'commands; to; load; programs' \  #
#                       -a path/to/adapters.fa \             #
#                       -n node_partition \                  #
#                       -t time:in:hours \                   #
#                       -k number_of_cores \                 #
#                       -m memory_per_cpu \                  #
#                       -f notificationEmail@forFailures.edu #
#                                                            #
# Parameters:                                                #
#     -h    Print the parameter list below then exit.        #
#     -o    (Required) Path to pipeline result directory.    #
#     -s    (Required) Skip merging of paired-end reads, just#
#           make empty directory to keep numbering consistent#
#   OR                                                       #
#     -i    (Required) Directory that contains the raw       #
#           fastq files to be processed. Sequences must have #
#           file extensions .fastq OR .fq,                   #
#           and can be gzipped or not.                       #
#     -o    (Required) Path to pipeline result directory.    #
#     -p    (Required) Single quoted string that contains    #
#           commands to load all the necessary programs      #
#           needed to run pipeline steps (e.g. activating    #
#           conda environments, loading modules, adding to   #
#           PATH, etc.).                                     #
#     -a    (Required) Path to adapters.fa file that comes   #
#           packaged with BBMerge and BBDuk.                 #
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
while getopts ":hsi:o:p:a:n:t:k:m:f:" opt; do
  case $opt in
    h)
    echo " Description: Merge paired-end reads using BBMerge.         "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
    echo "    BBMerge:    For merging paired-end reads.               "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " ./2.Merge_PE_Reads.sh -o output_dir \                      "
    echo "                       -s \                                 "
    echo "                     OR                                     "
    echo "                       -i input_seqs_dir \                  "
    echo "                       -o output_dir \                      "
    echo "                       -p 'commands; to; load; programs' \  "
    echo "                       -a path/to/adapters.fa \             "
    echo "                       -n node_partition \                  "
    echo "                       -t time:in:hours \                   "
    echo "                       -k number_of_cores \                 "
    echo "                       -m memory_per_cpu \                  "
    echo "                       -f notificationEmail@forFailures.edu "
    echo "                                                            "
    echo " Parameters:                                                "
    echo "     -h    Print the parameter list below then exit.        "
    echo "     -o    (Required) Path to pipeline result directory.    "
    echo "     -s    (Required) Skip merging of paired-end reads, just"
    echo "           make empty directory to keep numbering consistent"
    echo "   OR                                                       "
    echo "     -i    (Required) Directory that contains the raw       "
    echo "           fastq files to be processed. Sequences must have "
    echo "           file extensions .fastq OR .fq,                   "
    echo "           and can be gzipped or not.                       "
    echo "     -o    (Required) Path to pipeline result directory.    "
    echo "     -p    (Required) Single quoted string that contains    "
    echo "           commands to load all the necessary programs      "
    echo "           needed to run pipeline steps (e.g. activating    "
    echo "           conda environments, loading modules, adding to   "
    echo "           PATH, etc.).                                     "
    echo "     -a    (Required) Path to adapters.fa file that comes   "
    echo "           packaged with BBMerge and BBDuk.                 "
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
    s) SKIP=1
    ;;
    i) SEQ_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    p) PROG_LOAD="$OPTARG"
    ;;
    a) ADAPTERS="$OPTARG"
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

if [[ ! -z "$SKIP" ]]; then
  echo " "
  echo "*** Skipping running of BBMerge for paired read merging ***"
  echo "*** Making empty directory to keep numbering consistent ***"
  echo " "

  #Create directory
  if [ -d "${RESULTS_DIR}/2.Merged_Paired_End_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/2.Merged_Paired_End_Sequences
	  mkdir ${RESULTS_DIR}/2.Merged_Paired_End_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/2.Merged_Paired_End_Sequences/0.Output
  fi
  exit 0
fi

# Check that valid arguments were entered

# -i
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

# -a
if [[ -z "$ADAPTERS" ]]; then
  echo "ERROR: Argument -a is required, please supply path to the adapters.fa file that comes with BBMerge and BBDuk"
  exit 1
fi
if [[ -d "$ADAPTERS" ]]; then
  echo "ERROR: Argument -a should be the path to a single file, not a directory, please supply path to the adapters.fa file that comes with BBMerge and BBDuk"
  exit 1
fi
if echo $ADAPTERS | grep -q -v "adapters\.fa"; then
  echo "ERROR: path given to -o does not contain the file name adapters.fa, please supply the adapters.fa file that comes with BBMerge and BBDuk to this argument"
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

######### MERGE PAIRED READS WITH BBMERGE #######

  SECONDS=0
  echo " "
  echo "*** Running BBMerge to perform paired read merging ***"
  echo " "

  #Create directory for output
  if [ -d "${RESULTS_DIR}/2.Merged_Paired_End_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/2.Merged_Paired_End_Sequences
	  mkdir ${RESULTS_DIR}/2.Merged_Paired_End_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/2.Merged_Paired_End_Sequences/0.Output
  fi

  ##### Run BBMerge #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION" >> run.sh
  echo "#SBATCH --job-name=Merge" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/2.Merged_Paired_End_Sequences/0.ErrorOut/Merge_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/2.Merged_Paired_End_Sequences/0.Output/Merge_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_REQUEST" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_REQUEST" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_PER_CPU" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --array=1-$(ls -l ${SEQ_DIR}/*_R1*${SEQ_EXT} | wc -l)" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "FILE1=\$(ls ${SEQ_DIR}/*_R1*${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
  echo "FILE2=\$(ls ${SEQ_DIR}/*_R2*${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
  echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1' '{print \$1}')" >> run.sh
  echo "bbmerge.sh in1=\$FILE1 \\" >> run.sh
  echo "in2=\$FILE2 \\" >> run.sh
  echo "out=${RESULTS_DIR}/2.Merged_Paired_End_Sequences/\${FILE_NAME}.fastq.gz \\" >> run.sh
  echo "adapters=${ADAPTERS} \\" >> run.sh
  echo "rem iterations=5 extend2=20 ecct t=2 -Xmx${MAX_MEM}g \\" >> run.sh
  echo "> ${RESULTS_DIR}/2.Merged_Paired_End_Sequences/\${FILE_NAME}.log 2>&1"  >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  #Signal jobs have ended
  echo "Merging of paired end reads with BBMerge complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "

#################################################
