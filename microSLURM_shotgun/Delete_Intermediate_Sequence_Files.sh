#!/bin/bash

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# Last updated: 18 Nov 2022                                  #
#                                                            #
# Description: Delete sequence files generated during        #
# intermediate pipeline steps (merging, quality              #
# trimming/filtering, human decontamination) while keeping   #
# potentially important components of each step to keep      #
# record of (log files, extracted host sequences, extracted  #
# low-complexity sequences). Will reorganize output folder if#
# flag is given since the default numbering of folders will  #
# no longer make sense.                                      #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
#                                                            #
# Usage:                                                     #
# ./Delete_Intermediate_Sequence_Files.sh -o output_dir      #
#                                                            #
# Parameters:                                                #
#     -h    Print the parameter list below then exit.        #
#     -i    (Required) Directory that contains the raw       #
#           fastq files that were processed.                 #
#     -o    (Required) Path to pipeline result directory.    #
#     -x    (Optional) Add this flag if paired sequences were#
#           merged during pipeline.                          #
##############################################################

echo " "
echo "##############################################################"
echo "# Whole Genome Shotgun Metagenomic Processing Pipeline       #"
echo "# Last updated: 18 Nov 2022                                  #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":hi:o:x" opt; do
  case $opt in
    h)
    echo " Description: Delete sequence files generated during        "
    echo " intermediate pipeline steps (merging, quality              "
    echo " trimming/filtering, human decontamination) while keeping   "
    echo " potentially important components of each step to keep      "
    echo " record of (log files, extracted host sequences, extracted  "
    echo " low-complexity sequences). Will reorganize output folder if"
    echo " flag is given since the default numbering of folders will  "
    echo " no longer make sense.                                      "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " ./Delete_Intermediate_Sequence_Files.sh -o output_dir      "
    echo "                                                            "
    echo " Parameters:                                                "
    echo "     -h    Print the parameter list below then exit.        "
    echo "     -i    (Required) Directory that contains the raw       "
    echo "           fastq files that were processed.                 "
    echo "     -o    (Required) Path to pipeline result directory.    "
    echo "     -x    (Optional) Add this flag if paired sequences were"
    echo "           merged during pipeline.                          "
    echo " "
    exit 0
    ;;
    i) SEQ_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    o) RESULTS_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    x) MERGE=1
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
  echo "ERROR: Argument -o is required, please supply the pipeline output directory"
  exit 1
fi
if [[ ! -d "$RESULTS_DIR" ]]; then
  echo "ERROR: Argument -o should be a directory, please supply the pipeline output directory"
  exit 1
fi

###### CREATE DIRECTORY FOR PIPELINE OUTPUT #####
echo "*** Removing intermediate sequence files and reorganizing ***"
echo " "
  
#Create new directories needed
mkdir ${RESULTS_DIR}/2.Processed_Sequences
mkdir ${RESULTS_DIR}/2.Processed_Sequences/Log_Files

#Create combined log file for each sample
for file in ${SEQ_DIR}/*_R1*${SEQ_EXT}; do
  FILE_NAME=$(echo $file | awk -F '/' '{print $NF}' | awk -F '_R1' '{print $1}')

  #Initialize log file
  echo "*** Log file for sample ${FILE_NAME} ***" > ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log

  #Grab log file contents
  if [[ ! -z $MERGE ]]; then
    echo "### Log for merging of paired-end reads" >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
    echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
    cat ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log \
        ${RESULTS_DIR}/2.Merged_Paired_End_Sequences/${FILE_NAME}.log \
        > temp
    mv temp ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
    echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  fi

  echo "### Log for quality trimming/filtering of reads" >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  cat ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log \
      ${RESULTS_DIR}/3.Quality_Controlled_Sequences/${FILE_NAME}.log \
      > temp
  mv temp ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  cat ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log \
      ${RESULTS_DIR}/3.Quality_Controlled_Sequences/${FILE_NAME}_stats.txt \
      > temp
  mv temp ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log

  echo "### Log for removing host sequences" >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  cat ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log \
      ${RESULTS_DIR}/4.Decontaminated_Sequences/${FILE_NAME}.log \
      > temp
  mv temp ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  cat ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log \
      ${RESULTS_DIR}/4.Decontaminated_Sequences/${FILE_NAME}_refstats.log \
      > temp
  mv temp ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log

  echo "### Log for removing low-complexity sequences" >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  cat ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log \
      ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/${FILE_NAME}.log \
      > temp
  mv temp ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
  echo " " >> ${RESULTS_DIR}/2.Processed_Sequences/Log_Files/${FILE_NAME}.log
done

#Grab extracted sequences from intermediate directories
cp -r ${RESULTS_DIR}/4.Decontaminated_Sequences/Extracted_Host_Sequences \
      ${RESULTS_DIR}/2.Processed_Sequences/Extracted_Host_Sequences
cp -r ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/Removed_Sequences \
      ${RESULTS_DIR}/2.Processed_Sequences/Extracted_Low_Complexity_Sequences
  
#Grab processed sequences
cp ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*fastq.gz \
   ${RESULTS_DIR}/2.Processed_Sequences/
ERROR=$(diff <(ls ${RESULTS_DIR}/2.Processed_Sequences/*fastq.gz | awk -F'/' '{print $NF}') \
             <(ls ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*fastq.gz | awk -F'/' '{print $NF}'))
if [[ ! -z "$ERROR" ]]; then
  echo "ERROR: Something went wrong when copying processed sequences to new directory, copied sequences do not match sequences original sequences"
  exit 1
fi

#Modify numbering of post sequence processing directories
mv ${RESULTS_DIR}/6.FastQC_Final_Reports ${RESULTS_DIR}/3.FastQC_Final_Reports
mv ${RESULTS_DIR}/7.Taxonomic_Profiling ${RESULTS_DIR}/4.Taxonomic_Profiling
mv ${RESULTS_DIR}/8.Functional_Profiling ${RESULTS_DIR}/5.Functional_Profiling
  
#Remove directories with intermediate sequences
rm -rf ${RESULTS_DIR}/2.Merged_Paired_End_Sequences
rm -rf ${RESULTS_DIR}/3.Quality_Controlled_Sequences
rm -rf ${RESULTS_DIR}/4.Decontaminated_Sequences
rm -rf ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences

echo " "
echo "Intermediate sequences have been removed and pipeline output directory reorganized"
echo " "
