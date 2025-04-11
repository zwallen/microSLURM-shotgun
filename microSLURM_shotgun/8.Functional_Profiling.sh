#!/bin/bash

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# Last updated: 18 Nov 2022                                  #
#                                                            #
# Description: Perform functional profiling using HUMAnN     #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
#    HUMAnN:     For generating taxonomic, gene family, and  #
#                pathway abundances.                         #
#    ChocoPhlAn: Database used for taxonomic profiling. Can  #
#                Downloaded for HUMAnN using                 #
#                    $ humann_databases \\                   #
#                    $ --download chocophlan full \\         #
#                    $ $INSTALL_LOCATION                     #
#    UniRef:     Database used for functional profiling. Can #
#                be any of the UniRef databases downloaded   #
#                using humann_databases utility program.     #
#                                                            #
# Usage:                                                     #
# ./8.Functional_Profiling.sh [-x] \                         #
#                        -o output_dir \                     #
#                        -p 'commands; to; load; programs' \ #
#                        -c path/to/chocophlan/dir \         #
#                        -u path/to/uniref/dir \             #
#                        -n parition \                       #
#                        -t time:in:hours \                  #
#                        -k cores \                          #
#                        -m mem_per_cpu \                    #
#                        -f notificationEmail@forFailures.edu#
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
#     -c    (Required) Path to humann ChocoPhlAn database.   #
#     -u    (Required) Path to humann UniRef90 database.     #
#     -n    (Required) Name of the partition to request for  #
#           submitting jobs. Provide comma separated list of #
#           two names one for submitting HUMAnN job and one  #
#           for submitting additional MetaPhlAn run.         #
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
while getopts ":hxo:p:c:u:n:t:k:m:f:" opt; do
  case $opt in
    h)
    echo " Description: Perform functional profiling using HUMAnN.    "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
    echo "    HUMAnN:     For generating taxonomic, gene family, and  "
    echo "                pathway abundances.                         "
    echo "    ChocoPhlAn: Database used for taxonomic profiling. Can  "
    echo "                Downloaded for HUMAnN using                 "
    echo "                    $ humann_databases \\                   "
    echo "                    $ --download chocophlan full \\         "
    echo "                    $ \$INSTALL_LOCATION                    "
    echo "    UniRef:     Database used for functional profiling. Can "
    echo "                be any of the UniRef databases downloaded   "
    echo "                using humann_databases utility program.     "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " ./8.Functional_Profiling.sh [-x] \                         "
    echo "                        -o output_dir \                     "
    echo "                        -p 'commands; to; load; programs' \ "
    echo "                        -c path/to/chocophlan/dir \         "
    echo "                        -u path/to/uniref/dir \             "
    echo "                        -n parition \                       "
    echo "                        -t time:in:hours \                  "
    echo "                        -k cores \                          "
    echo "                        -m mem_per_cpu \                    "
    echo "                        -f notificationEmail@forFailures.edu"
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
    echo "     -c    (Required) Path to humann ChocoPhlAn database.   "
    echo "     -u    (Required) Path to humann UniRef90 database.     "
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
    c) CHOCO=$(echo $OPTARG | sed 's#/$##')
    ;;
    u) UNIREF=$(echo $OPTARG | sed 's#/$##')
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

# -c
if [[ -z "$CHOCO" ]]; then
  echo "ERROR: Argument -c is required, please supply path to ChocoPhlAn database directory"
  exit 1
fi
if [[ ! -d "$CHOCO" ]]; then
  echo "ERROR: Argument -c should be a directory, please supply path to ChocoPhlAn database directory"
  exit 1
fi

# -u
if [[ -z "$UNIREF" ]]; then
  echo "ERROR: Argument -u is required, please supply path to UniRef database directory"
  exit 1
fi
if [[ ! -d "$UNIREF" ]]; then
  echo "ERROR: Argument -u should be a directory, please supply path to UniRef database directory"
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
  OUT_DIR=5.Functional_Profiling
else
  IN_DIR=5.Low_Complexity_Filtered_Sequences
  OUT_DIR=8.Functional_Profiling
fi

###### FUNCTIONAL PROFILING #######

  SECONDS=0
  echo " "
  echo "*** Running HUMAnN workflow for performing functional profiling ***"
  echo " "

  #Create directory for profiling output
  if [ -d "${RESULTS_DIR}/${OUT_DIR}" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/${OUT_DIR}
	  mkdir ${RESULTS_DIR}/${OUT_DIR}/0.ErrorOut
	  mkdir ${RESULTS_DIR}/${OUT_DIR}/0.Output
  fi

  ##### Run HUMAnN workflow #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION" >> run.sh
  echo "#SBATCH --job-name=Profiling" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/${OUT_DIR}/0.ErrorOut/Functional_profiling_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/${OUT_DIR}/0.Output/Functional_profiling_%A_%a.out" >> run.sh
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
  else
    echo "FILE1=\$(ls ${RESULTS_DIR}/${IN_DIR}/*R1.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE2=\$(ls ${RESULTS_DIR}/${IN_DIR}/*R2.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1' '{print \$1}')" >> run.sh
    echo "zcat \$FILE1 \$FILE2 > ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}.temp.fastq" >> run.sh
    echo "FILE=${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}.temp.fastq" >> run.sh
  fi
  echo "humann --input \$FILE \\" >> run.sh
  echo "--output ${RESULTS_DIR}/${OUT_DIR} \\" >> run.sh
  echo "--output-basename \$FILE_NAME \\" >> run.sh
  echo "--nucleotide-database $CHOCO \\" >> run.sh
  echo "--protein-database $UNIREF \\" >> run.sh
  echo "--metaphlan-options '-t rel_ab' \\" >> run.sh
  echo "--prescreen-threshold 0.01 \\" >> run.sh
  echo "--threads $CPU_REQUEST \\" >> run.sh
  echo "--verbose \\" >> run.sh
  echo "> ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}.log 2>&1" >> run.sh
  echo "rm ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}.temp.fastq" >> run.sh
  echo "rm ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_humann_temp/*_bowtie2_*" >> run.sh
  echo "rm ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_humann_temp/*_diamond_*" >> run.sh
  echo "rm ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_humann_temp/*_custom_chocophlan_database.ffn" >> run.sh
  if [[ -z "$MERGE" ]] && [[ -z "$JOIN" ]]; then
    echo "rm ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}.temp.fastq" >> run.sh
  fi
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  ##### Clean up #####
  echo "Cleaning things up... merging per sample tables and removing unneeded files"
  echo " "
  for dir in ${RESULTS_DIR}/${OUT_DIR}/*humann_temp; do
    FILE_NAME=$(echo $dir | awk -F '/' '{print $NF}' | awk -F '_humann_temp' '{print $1}')

    mkdir ${RESULTS_DIR}/${OUT_DIR}/${FILE_NAME}_Profiles

    mv ${dir}/${FILE_NAME}_metaphlan_bugs_list.tsv \
    ${RESULTS_DIR}/${OUT_DIR}/${FILE_NAME}_Profiles/

    mv ${dir}/${FILE_NAME}_metaphlan_bowtie2.txt \
    ${RESULTS_DIR}/${OUT_DIR}/${FILE_NAME}_Profiles/

    mv ${dir}/${FILE_NAME}.log \
    ${RESULTS_DIR}/${OUT_DIR}/${FILE_NAME}_Profiles/${FILE_NAME}_humann.log

    rm -r ${dir}

    mv ${RESULTS_DIR}/${OUT_DIR}/${FILE_NAME}*.* \
    ${RESULTS_DIR}/${OUT_DIR}/${FILE_NAME}_Profiles/
  done

  ##### Merge tables #####
  mkdir ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables

  for dir in ${RESULTS_DIR}/${OUT_DIR}/*_Profiles; do
    FILE_NAME=$(echo $dir | awk -F '/' '{print $NF}' | awk -F '_Profiles' '{print $1}')

    cp ${dir}/${FILE_NAME}_metaphlan_bugs_list.tsv ${RESULTS_DIR}/${OUT_DIR}/

    cp ${dir}/${FILE_NAME}_genefamilies.tsv ${RESULTS_DIR}/${OUT_DIR}/

    cp ${dir}/${FILE_NAME}_pathabundance.tsv ${RESULTS_DIR}/${OUT_DIR}/

    cp ${dir}/${FILE_NAME}_pathcoverage.tsv ${RESULTS_DIR}/${OUT_DIR}/
  done

  #Create individual scripts for running table merging
  echo '#!/bin/bash' > run1.sh
  echo "$PROG_LOAD" >> run1.sh
  echo "merge_metaphlan_tables.py ${RESULTS_DIR}/${OUT_DIR}/*_metaphlan_bugs_list.tsv \\" >> run1.sh
  echo "> ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables/metaphlan_bugs_list.tsv" >> run1.sh
  echo "sed -i '2s/_metaphlan_bugs_list//g' ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables/metaphlan_bugs_list.tsv" >> run1.sh
  echo "rm ${RESULTS_DIR}/${OUT_DIR}/*_metaphlan_bugs_list.tsv" >> run1.sh
  chmod +x run1.sh

  echo '#!/bin/bash' > run2.sh
  echo "$PROG_LOAD" >> run2.sh
  echo "humann_join_tables --input ${RESULTS_DIR}/${OUT_DIR}/ \\" >> run2.sh
  echo "--output ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables/humann_genefamilies.tsv \\" >> run2.sh
  echo "--file_name genefamilies.tsv" >> run2.sh
  echo "rm ${RESULTS_DIR}/${OUT_DIR}/*genefamilies.tsv" >> run2.sh
  chmod +x run2.sh

  echo '#!/bin/bash' > run3.sh
  echo "$PROG_LOAD" >> run3.sh
  echo "humann_join_tables --input ${RESULTS_DIR}/${OUT_DIR}/ \\" >> run3.sh
  echo "--output ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables/humann_pathabundance.tsv \\" >> run3.sh
  echo "--file_name pathabundance.tsv" >> run3.sh
  echo "rm ${RESULTS_DIR}/${OUT_DIR}/*pathabundance.tsv" >> run3.sh
  chmod +x run3.sh

  echo '#!/bin/bash' > run4.sh
  echo "$PROG_LOAD" >> run4.sh
  echo "humann_join_tables --input ${RESULTS_DIR}/${OUT_DIR}/ \\" >> run4.sh
  echo "--output ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables/humann_pathcoverage.tsv \\" >> run4.sh
  echo "--file_name pathcoverage.tsv" >> run4.sh
  echo "rm ${RESULTS_DIR}/${OUT_DIR}/*pathcoverage.tsv" >> run4.sh
  chmod +x run4.sh

  #Create script for running and submitting individual scripts
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION" >> run.sh
  echo "#SBATCH --job-name=Merge_tables" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/${OUT_DIR}/0.ErrorOut/Merge_tables_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/${OUT_DIR}/0.Output/Merge_tables_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_REQUEST" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_REQUEST" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_PER_CPU" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --array=1-4" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "./run\${SLURM_ARRAY_TASK_ID}.sh" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run*sh

  #Signal workflow has completed
  echo "Running of HUMAnN workflow complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
#################################################
