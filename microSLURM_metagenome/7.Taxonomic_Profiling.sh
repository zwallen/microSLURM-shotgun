#!/bin/bash

##############################################################
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# Last updated: 18 Nov 2022                                  #
#                                                            #
# Description: Perform taxonomic profiling using MetaPhlAn   #
#              or Kraken2/Bracken                            #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
#    MetaPhlAn:  For generating taxonomic abundances.        #
#    ChocoPhlAn: Database used for taxonomic profiling.      #
#                Downloaded for MetaPhlAn using              #
#                    $ metaphlan --install                   #
# Alternatively                                              #
#    Kraken2:    For taxonomic profiling.                    #
#    Bracken:    For calculating taxonomic abundances.       #
#    Kraken2 database:                                       #
#                Database used for taxonomic profiling with  #
#                Kraken2. Should be built using the Kraken2  #
#                command 'kraken2-build', then the kmer      #
#                distribution calculated using 'braken-build'#
#                Alternatively, can download the pre-built   #
#                standard database from the repository:      #
#                https://benlangmead.github.io/aws-indexes/k2#
#    kreport2mpa.py and combine_mpa.py:                      #
#                Scripts from KrakenTools to convert Bracken #
#                output to MetaPhlAn style output and combine#
#                individual sample MetaPhlAn style outputs.  #
#                Make sure they are in your $PATH.           #
#                                                            #
# Usage:                                                     #
# ./7.Taxonomic_Profiling.sh [-x] \                          #
#                        -o output_dir \                     #
#                        -p 'commands; to; load; programs' \ #
#                        -n node_partition \                 #
#                        -t time:in:hours \                  #
#                        -k number_of_cores \                #
#                        -m memory_per_cpu \                 #
#                        -f notificationEmail@forFailures.edu#
# if wanting to use Kraken2/Bracken instead of MetaPhlAn:    #
#                        -b /path/to/kraken2/database        #
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
#     -b    (Optional) Use Kraken2/Bracken to perform        #
#           taxonomic profiling and supply path to kraken2   #
#           database directory here.                         #
##############################################################

echo " "
echo "##############################################################"
echo "# Whole Genome Shotgun Metagenomic Processing Pipeline       #"
echo "# Last updated: 18 Nov 2022                                  #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":hxo:p:n:t:k:m:f:b:" opt; do
  case $opt in
    h)
    echo " Description: Perform taxonomic profiling using MetaPhlAn.  "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
    echo "    MetaPhlAn:  For generating taxonomic abundances.        "
    echo "    ChocoPhlAn: Database used for taxonomic profiling.      "
    echo "                Downloaded for MetaPhlAn using              "
    echo "                    $ metaphlan --install                   "
    echo " Alternatively                                              "
    echo "    Kraken2:    For taxonomic profiling.                    "
    echo "    Bracken:    For calculating taxonomic abundances.       "
    echo "    Kraken2 database:                                       "
    echo "                Database used for taxonomic profiling with  "
    echo "                Kraken2. Should be built using the Kraken2  "
    echo "                command 'kraken2-build', then the kmer      "
    echo "                distribution calculated using 'braken-build'"
    echo "                Alternatively, can download the pre-built   "
    echo "                standard database from the repository:      "
    echo "                https://benlangmead.github.io/aws-indexes/k2"
    echo "    kreport2mpa.py and combine_mpa.py:                      "
    echo "                Scripts from KrakenTools to convert Bracken "
    echo "                output to MetaPhlAn style output and combine"
    echo "                individual sample MetaPhlAn style outputs.  "
    echo "                Make sure they are in your \$PATH.          "
    echo "                                                            "
    echo " Usage:                                                     "
    echo " ./7.Taxonomic_Profiling.sh [-x] \                          "
    echo "                        -o output_dir \                     "
    echo "                        -p 'commands; to; load; programs' \ "
    echo "                        -n node_partition \                 "
    echo "                        -t time:in:hours \                  "
    echo "                        -k number_of_cores \                "
    echo "                        -m memory_per_cpu \                 "
    echo "                        -f notificationEmail@forFailures.edu"
    echo " if wanting to use Kraken2/Bracken instead of MetaPhlAn:    "
    echo "                        -b /path/to/kraken2/database        "
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
    echo "     -b    (Optional) Use Kraken2/Bracken to perform        "
    echo "           taxonomic profiling and supply path to kraken2   "
    echo "           database directory here.                         "
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
    b) BRACKEN="$OPTARG"
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

# -b
if [[ ! -z "$BRACKEN" ]]; then
  if [[ ! -d "$BRACKEN" ]]; then
    echo "ERROR: Argument -b should be a directory, please supply a the Kraken2 database directory"
    exit 1
  fi
  if [[ ! -f "$BRACKEN/hash.k2d" ]]; then
    echo "ERROR: Did not find file 'hash.k2d' in supplied Kraken2 database directory, needed for taxonomic profiling"
    exit 1
  fi
  if [[ ! -f "$BRACKEN/opts.k2d" ]]; then
    echo "ERROR: Did not find file 'opts.k2d' in supplied Kraken2 database directory, needed for taxonomic profiling"
    exit 1
  fi
  if [[ ! -f "$BRACKEN/taxo.k2d" ]]; then
    echo "ERROR: Did not find file 'taxo.k2d' in supplied Kraken2 database directory, needed for taxonomic profiling"
    exit 1
  fi
  if [[ ! -f "$BRACKEN/database50mers.kmer_distrib" ]]; then
    echo "ERROR: Did not find file 'database50mers.kmer_distrib' in supplied Kraken2 database directory, was bracken-build ran with '-l 50'?"
    exit 1
  fi
fi

#Determine if pipeline result directory has been reorganized and edit input/output directories accordingly
if [[ -d "${RESULTS_DIR}/2.Processed_Sequences" ]]; then
  IN_DIR=2.Processed_Sequences
  OUT_DIR=4.Taxonomic_Profiling
else
  IN_DIR=5.Low_Complexity_Filtered_Sequences
  OUT_DIR=7.Taxonomic_Profiling
fi

###### TAXONOMIC PROFILING #######

#Create directory for profiling output
if [ -d "${RESULTS_DIR}/${OUT_DIR}" ]
then
  :
else
  mkdir ${RESULTS_DIR}/${OUT_DIR}
  mkdir ${RESULTS_DIR}/${OUT_DIR}/0.ErrorOut
  mkdir ${RESULTS_DIR}/${OUT_DIR}/0.Output
fi

if [[ -z $BRACKEN ]]; then
  SECONDS=0
  echo " "
  echo "*** Running MetaPhlAn workflow for performing taxonomic profiling ***"
  echo " "

  ##### Run MetaPhlAn workflow #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION" >> run.sh
  echo "#SBATCH --job-name=Profiling" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/${OUT_DIR}/0.ErrorOut/Taxonomic_Profiling_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/${OUT_DIR}/0.Output/Taxonomic_Profiling_%A_%a.out" >> run.sh
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
  echo "if [ -d \"${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles\" ]" >> run.sh
  echo "then" >> run.sh
  echo "  :" >> run.sh
  echo "else" >> run.sh
  echo "  mkdir ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles" >> run.sh
  echo "fi" >> run.sh
  echo "metaphlan \$FILE \\" >> run.sh
  echo "--input_type fastq \\" >> run.sh
  echo "-t rel_ab \\" >> run.sh
  echo "--nproc $CPU_REQUEST \\" >> run.sh
  echo "--bowtie2out ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles/\${FILE_NAME}_metaphlan_bowtie2.txt \\" >> run.sh
  echo "-o ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles/\${FILE_NAME}_metaphlan_rel_ab.tsv" >> run.sh
  if [[ -z "$MERGE" ]] && [[ -z "$JOIN" ]]; then
      echo "rm ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}.temp.fastq" >> run.sh
  fi
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  ##### Calculate taxonomic count data #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION" >> run.sh
  echo "#SBATCH --job-name=get_taxa_counts" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/${OUT_DIR}/0.ErrorOut/get_taxa_counts_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/${OUT_DIR}/0.Output/get_taxa_counts_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_REQUEST" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_REQUEST" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_PER_CPU" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --array=1-$(ls -d -l ${RESULTS_DIR}/${OUT_DIR}/*Profiles | wc -l)" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "DIR=\$(ls -d ${RESULTS_DIR}/${OUT_DIR}/*Profiles | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
  echo "FILE_NAME=\$(echo \$DIR | awk -F '/' '{print \$NF}' | awk -F '_Profiles' '{print \$1}')" >> run.sh
  echo "metaphlan \${DIR}/\${FILE_NAME}_metaphlan_bowtie2.txt \\" >> run.sh
  echo "--input_type bowtie2out \\" >> run.sh
  echo "-t rel_ab \\" >> run.sh
  echo "--unknown_estimation \\" >> run.sh
  echo "--nproc $CPU_REQUEST \\" >> run.sh
  echo "-o \${DIR}/\${FILE_NAME}_metaphlan_bugs_list_w_unknown.tsv" >> run.sh
  echo "NREADS=\$(grep '#nreads' \${DIR}/\${FILE_NAME}_metaphlan_bowtie2.txt | awk '{print \$2}')" >> run.sh
  echo "grep -v '^#' \${DIR}/\${FILE_NAME}_metaphlan_bugs_list_w_unknown.tsv | \\" >> run.sh
  echo "awk -v nreads=\"\$NREADS\" '{OFMT=\"%f\";print \$1,\$2,\$3,(\$3/100)*nreads,\$4}' OFS='\t' | \\" >> run.sh
  echo "cat <(grep '^#' \${DIR}/\${FILE_NAME}_metaphlan_bugs_list_w_unknown.tsv) - | \\" >> run.sh
  echo "sed 's/#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species/#clade_name\tNCBI_tax_id\trelative_abundance\tcounts\tadditional_species/' \\" >> run.sh
  echo "> \${DIR}/\${FILE_NAME}_metaphlan_rel_ab_w_counts.tsv" >> run.sh
  echo "rm \${DIR}/\${FILE_NAME}_metaphlan_bugs_list_w_unknown.tsv" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  ##### Merge tables #####
  #Create directory for merged tables
  if [ -d "${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables" ]
  then
	  :
  else
    mkdir ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables
  fi

  #For each sample get counts and relative abundances
  for dir in ${RESULTS_DIR}/${OUT_DIR}/*_Profiles; do
    FILE_NAME=$(echo $dir | awk -F '/' '{print $NF}' | awk -F '_Profiles' '{print $1}')

    cp ${dir}/${FILE_NAME}_metaphlan_rel_ab.tsv ${RESULTS_DIR}/${OUT_DIR}/

    awk -F'\t' '{print $1,$2,$3}' OFS='\t' ${dir}/${FILE_NAME}_metaphlan_rel_ab_w_counts.tsv \
    > ${RESULTS_DIR}/${OUT_DIR}/${FILE_NAME}_metaphlan_rel_ab_w_unknown.tsv

    awk -F'\t' '{print $1,$2,$4}' OFS='\t' ${dir}/${FILE_NAME}_metaphlan_rel_ab_w_counts.tsv | \
    sed 's/counts/relative_abundance/' \
    > ${RESULTS_DIR}/${OUT_DIR}/${FILE_NAME}_metaphlan_counts.tsv
  done

  #Create individual scripts for running table merging
  echo '#!/bin/bash' > run1.sh
  echo "$PROG_LOAD" >> run1.sh
  echo "merge_metaphlan_tables.py ${RESULTS_DIR}/${OUT_DIR}/*_metaphlan_rel_ab.tsv \\" >> run1.sh
  echo "> ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables/metaphlan_rel_ab.tsv" >> run1.sh
  echo "sed -i '2s/_metaphlan_rel_ab//g' ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables/metaphlan_rel_ab.tsv" >> run1.sh
  echo "rm ${RESULTS_DIR}/${OUT_DIR}/*_metaphlan_rel_ab.tsv" >> run1.sh
  chmod +x run1.sh

  echo '#!/bin/bash' > run2.sh
  echo "$PROG_LOAD" >> run2.sh
  echo "merge_metaphlan_tables.py ${RESULTS_DIR}/${OUT_DIR}/*_metaphlan_rel_ab_w_unknown.tsv \\" >> run2.sh
  echo "> ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables/metaphlan_rel_ab_w_unknown.tsv" >> run2.sh
  echo "sed -i '2s/_metaphlan_rel_ab_w_unknown//g' ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables/metaphlan_rel_ab_w_unknown.tsv" >> run2.sh
  echo "rm ${RESULTS_DIR}/${OUT_DIR}/*_metaphlan_rel_ab_w_unknown.tsv" >> run2.sh
  chmod +x run2.sh

  echo '#!/bin/bash' > run3.sh
  echo "$PROG_LOAD" >> run3.sh
  echo "merge_metaphlan_tables.py ${RESULTS_DIR}/${OUT_DIR}/*_metaphlan_counts.tsv \\" >> run3.sh
  echo "> ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables/metaphlan_counts.tsv" >> run3.sh
  echo "sed -i '2s/_metaphlan_counts//g' ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables/metaphlan_counts.tsv" >> run3.sh
  echo "rm ${RESULTS_DIR}/${OUT_DIR}/*_metaphlan_counts.tsv" >> run3.sh
  chmod +x run3.sh

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
  echo "#SBATCH --array=1-3" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "./run\${SLURM_ARRAY_TASK_ID}.sh" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run*sh

  #Signal workflow has completed
  echo "Running of MetaPhlAn workflow complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
else
  SECONDS=0
  echo " "
  echo "*** Running Kraken2/Bracken workflow for performing taxonomic profiling ***"
  echo " "

  ##### Run Kraken2/Bracken workflow #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION" >> run.sh
  echo "#SBATCH --job-name=Profiling" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/${OUT_DIR}/0.ErrorOut/Taxonomic_Profiling_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/${OUT_DIR}/0.Output/Taxonomic_Profiling_%A_%a.out" >> run.sh
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
  fi
  echo "if [ -d \"${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles\" ]" >> run.sh
  echo "then" >> run.sh
  echo "  :" >> run.sh
  echo "else" >> run.sh
  echo "  mkdir ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles" >> run.sh
  echo "fi" >> run.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "kraken2 \$FILE \\" >> run.sh
  else
    echo "kraken2 \$FILE1 \$FILE2 \\" >> run.sh
    echo "--paired \\" >> run.sh
  fi
  echo "--db $BRACKEN \\" >> run.sh
  echo "--threads $CPU_REQUEST \\" >> run.sh
  echo "--confidence 0.5 \\" >> run.sh
  echo "--minimum-base-quality 0 \\" >> run.sh
  echo "--minimum-hit-groups 3 \\" >> run.sh
  echo "--report ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles/\${FILE_NAME}_kraken2_report.tsv \\" >> run.sh
  echo "--report-zero-counts \\" >> run.sh
  echo "--use-names \\" >> run.sh
  echo "--gzip-compressed \\" >> run.sh
  echo "--unclassified-out ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles/\${FILE_NAME}_kraken2_unclassified#.fastq \\" >> run.sh
  echo "--output ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles/\${FILE_NAME}_kraken2_classified.tsv \\" >> run.sh
  echo "> ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles/\${FILE_NAME}.log 2>&1" >> run.sh
  echo "bracken -d $BRACKEN \\" >> run.sh
  echo "-i ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles/\${FILE_NAME}_kraken2_report.tsv \\" >> run.sh
  echo "-r 50 -l S -t 30 \\" >> run.sh
  echo "-w ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles/\${FILE_NAME}_bracken_report.tsv \\" >> run.sh
  echo "-o ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles/\${FILE_NAME}_bracken_counts.tsv \\" >> run.sh
  echo ">> ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles/\${FILE_NAME}.log 2>&1" >> run.sh
  echo "kreport2mpa.py -r ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles/\${FILE_NAME}_bracken_report.tsv \\" >> run.sh
  echo "-o ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles/\${FILE_NAME}_bracken_counts_mpa.tsv" >> run.sh
  echo "sed -i \"1s/^/#clade_name\t\${FILE_NAME}\n/\" ${RESULTS_DIR}/${OUT_DIR}/\${FILE_NAME}_Profiles/\${FILE_NAME}_bracken_counts_mpa.tsv" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  ##### Merge tables #####
  #Create directory for merged tables
  if [ -d "${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables" ]
  then
	  :
  else
    mkdir ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables
  fi

  #For each sample get counts
  for dir in ${RESULTS_DIR}/${OUT_DIR}/*_Profiles; do
    FILE_NAME=$(echo $dir | awk -F '/' '{print $NF}' | awk -F '_Profiles' '{print $1}')
    cp ${dir}/${FILE_NAME}_bracken_counts_mpa.tsv ${RESULTS_DIR}/${OUT_DIR}/
  done
  
  #Create script for running program and submit
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
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "combine_mpa.py -i ${RESULTS_DIR}/${OUT_DIR}/*_bracken_counts_mpa.tsv \\" >> run.sh
  echo "-o ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables/bracken_counts_mpa.tsv" >> run.sh
  echo "sed -i '1s/#Classification/clade_name/' ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables/bracken_counts_mpa.tsv" >> run.sh
  echo "sed -i \"1s/^/#$(echo ${BRACKEN} | awk -F'/' '{print $NF}')\n/\" ${RESULTS_DIR}/${OUT_DIR}/Merged_Sample_Tables/bracken_counts_mpa.tsv" >> run.sh
  echo "rm ${RESULTS_DIR}/${OUT_DIR}/*_bracken_counts_mpa.tsv" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh
  
  #Signal workflow has completed
  echo "Running of Kraken2/Bracken workflow complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################
