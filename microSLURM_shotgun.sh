#!/bin/bash

##############################################################
# microSLURM_shotgun                                         #
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# Last updated: 18 Nov 2022                                  #
#                                                            #
# Description: This is a wrapper program that wraps various  #
# programs to process raw paired-end whole genome shotgun    #
# metagenomic sequences. The end product of this pipeline    #
# is taxonomic and functional (gene and pathway) abundances  #
# that are ready for further statistical analyses.           #
#                                                            #
# Required programs and databases:                           #
#    SLURM:      Program is designed to work with a SLURM    #
#                high performance computing cluster          #
#                scheduling system.                          #
#    FastQC:     For performing quality reports.             #
#    MultiQC:    For summarizing multiple FastQC reports.    #
#    BBMerge:    For merging paired-end reads.               #
#    BBDuk:      For adapter and quality trimming of raw wgs #
#                reads, removal of PhiX sequences, and       #
#                removing low complexity sequences.          #
#    BBMap/                                                  #
#    BBSplit:    For removing host contamination from wgs    #
#                reads. Requires FASTA genome reference file #
#                to map reads against.                       #
#    MetaPhlAn:  For generating taxonomic abundances.        #
#    HUMAnN:     For generating gene family, and pathway     #
#                abundances.                                 #
#    ChocoPhlAn: Database used for taxonomic profiling.      #
#                Downloaded for MetaPhlAn using              #
#                    $ metaphlan --install                   #
#                and for HUMAnN using                        #
#                    $ humann_databases \                    #
#                    $ --download chocophlan full \          #
#                    $ $INSTALL_LOCATION                     #
#    UniRef:     Database used for functional profiling. Can #
#                be any of the UniRef databases downloaded   #
#                using humann_databases utility program.     #
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
# microSLURM_shotgun.sh -i input_seqs_dir \               #
#                    -o output_dir \                         #
#                    -p 'commands; to; load; programs' \     #
#                    -r path/to/host/ref/file.fa \           #
#                    -c path/to/chocophlan/dir \             #
#                    -u path/to/uniref/dir \                 #
#                    -f notificationEmail@forFailures.edu \  #
#                    -n list,of,node,partitions \            #
#                    -t list,of,time,requests \              #
#                    -k list,of,cpu,number,requests \        #
#                    -m list,of,memory,per,cpu,requests \    #
# if wanting to use Kraken2/Bracken instead of MetaPhlAn:    #
#                    -b /path/to/kraken2/database \          #
#                    [additional options]                    #
#                                                            #
# Parameters:                                                #
#     -h    Print the parameter list below then exit.        #
#                                                            #
# Required analysis parameters                               #
#     -i    (Required) Directory that contains the raw       #
#           fastq files to be processed. Sequences must have #
#           file extensions .fastq OR .fq,                   #
#           and can be gzipped or not.                       #
#     -o    (Required) Directory to put output of pipeline   #
#           into. NOTE: make sure output directory is in an  #
#           area that has plenty of data storage space       #
#           available if processing large datasets.          #
#     -p    (Required) Single quoted string that contains    #
#           commands to load all the necessary programs      #
#           needed to run pipeline steps (e.g. activating    #
#           conda environments, loading modules, adding to   #
#           PATH, etc.).                                     #
#     -r    (Required) Path to reference genome file of host.#
#           Should be in FASTA format, and uncompressed.     #
#     -c    (Required) Path to humann ChocoPhlAn database.   #
#     -u    (Required) Path to humann UniRef90 database.     #
#                                                            #
# Required SLURM parameters                                  #
# Note: For each parameter, a comma separated list of 8      #
# entries must be given (one for each pipeline step). If not #
# performing a certain step (e.g. merging of paired reads)   #
# then you can put a "NA" in the list.                       #
#     -n    (Required) Names of the partitions to request for#
#           submitting jobs.                                 #
#     -t    (Required) Time requests for running jobs.       #
#           Specify in hours (e.g. 12:00:00 for 12 hours).   #
#     -k    (Required) Number of cores wanted for jobs.      #
#     -m    (Required) Amount of memory requested for each   #
#           requested core. Specify in Megabytes.            #
#     -f    (Required) E-mail to send notifications to upon  #
#           failure of any jobs.                             #
#                                                            #
# Optional pipeline parameters                               #
#     -b    (Optional) Use Kraken2/Bracken to perform        #
#           taxonomic profiling and supply path to kraken2   #
#           database directory here.                         #
#     -x    (Optional) Merge paired-end reads before         #
#           performing the pipeline using BBMerge.           #
#     -a    (Optional) Path to adapters.fa file that comes   #
#           packaged with BBMerge and BBDuk. Required when   #
#           merging reads.                                   #
#     -j    (Optional) Join forward and reverse reads        #
#           together into one file (via concatenation) prior #
#           to performing low-complexity sequence filtering. #
#           May want to do this if using MetaPhlAn for       #
#           taxonomic profiling as this will be done anyway  #
#           prior to running MetaPhlAn and low-complexity    #
#           sequences caused by technical artificat are more #
#           prone to be located in reverse reads while the   #
#           forward reads are still good.                    #
#     -s    (Optional) Skip certain steps in the pipeline if #
#           need be. Provide a comma separated list of steps #
#           that you wish to skip in the pipeline. List may  #
#           have the values: fastqc_initial, QC, decontam,   #
#           entropy_filter, fastqc_final,                    #
#           taxonomic_profiling, functional_profiling        #
#     -d    (Optional) Delete sequence files generated during#
#           intermediate pipeline steps (merging, quality    #
#           trimming/filtering, human decontamination) while #
#           keeping potentially important components of each #
#           step to keep record of (log files, extracted host#
#           sequences, extracted low-complexity sequences).  #
#           Will reorganize output folder if flag is given   #
#           since the default numbering of folder will no    #
#           longer make sense.                               #
##############################################################

echo " "
echo "##############################################################"
echo "# microSLURM_shotgun                                         #"
echo "# Whole Genome Shotgun Metagenomic Processing Pipeline       #"
echo "# Last updated: 18 Nov 2022                                  #"
echo "##############################################################"
echo " "

# Argument parsing
while getopts ":hi:o:p:r:c:u:n:t:k:m:f:b:xa:js:d" opt; do
  case $opt in
    h)
    echo " Description: This is a wrapper program that wraps various  "
    echo " programs to process raw paired-end whole genome shotgun    "
    echo " metagenomic sequences. The end product of this pipeline    "
    echo " is taxonomic and functional (gene and pathway) abundances  "
    echo " that are ready for further statistical analyses.           "
    echo "                                                            "
    echo " Required programs and databases:                           "
    echo "    SLURM:      Program is designed to work with a SLURM    "
    echo "                high performance computing cluster          "
    echo "                scheduling system.                          "
    echo "    FastQC:     For performing quality reports.             "
    echo "    MultiQC:    For summarizing multiple FastQC reports.    "
    echo "    BBMerge:    For merging paired-end reads.               "
    echo "    BBDuk:      For adapter and quality trimming of raw wgs "
    echo "                reads, removal of PhiX sequences, and       "
    echo "                removing low complexity sequences.          "
    echo "    BBMap/                                                  "
    echo "    BBSplit:    For removing host contamination from wgs    "
    echo "                reads. Requires FASTA genome reference file "
    echo "                to map reads against.                       "
    echo "    MetaPhlAn:  For generating taxonomic abundances.        "
    echo "    HUMAnN:     For generating gene family, and pathway     "
    echo "                abundances.                                 "
    echo "    ChocoPhlAn: Database used for taxonomic profiling.      "
    echo "                Downloaded for MetaPhlAn using              "
    echo "                    $ metaphlan --install                   "
    echo "                and for HUMAnN using                        "
    echo "                    $ humann_databases \\                   "
    echo "                    $ --download chocophlan full \\         "
    echo "                    $ \$INSTALL_LOCATION                    "
    echo "    UniRef:     Database used for functional profiling. Can "
    echo "                be any of the UniRef databases downloaded   "
    echo "                using humann_databases utility program.     "
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
    echo " microSLURM_shotgun.sh -i input_seqs_dir \                  "
    echo "                    -o output_dir \                         "
    echo "                    -p 'commands; to; load; programs' \     "
    echo "                    -r path/to/host/ref/file.fa \           "
    echo "                    -c path/to/chocophlan/dir \             "
    echo "                    -u path/to/uniref/dir \                 "
    echo "                    -f notificationEmail@forFailures.edu \  "
    echo "                    -n list,of,node,partitions \            "
    echo "                    -t list,of,time,requests \              "
    echo "                    -k list,of,cpu,number,requests \        "
    echo "                    -m list,of,memory,per,cpu,requests \    "
    echo " if wanting to use Kraken2/Bracken instead of MetaPhlAn:    "
    echo "                    -b /path/to/kraken2/database \          "
    echo "                    [additional options]                    "
    echo "                                                            "
    echo " Parameters:                                                "
    echo "     -h    Print the parameter list below then exit.        "
    echo "                                                            "
    echo " Required analysis parameters                               "
    echo "     -i    (Required) Directory that contains the raw       "
    echo "           fastq files to be processed. Sequences must have "
    echo "           file extensions .fastq OR .fq,                   "
    echo "           and can be gzipped or not.                       "
    echo "     -o    (Required) Directory to put output of pipeline   "
    echo "           into. NOTE: make sure output directory is in an  "
    echo "           area that has plenty of data storage space       "
    echo "           available if processing large datasets.          "
    echo "     -p    (Required) Single quoted string that contains    "
    echo "           commands to load all the necessary programs      "
    echo "           needed to run pipeline steps (e.g. activating    "
    echo "           conda environments, loading modules, adding to   "
    echo "           PATH, etc.).                                     "
    echo "     -r    (Required) Path to reference genome file of host."
    echo "           Should be in FASTA format, and uncompressed.     "
    echo "     -c    (Required) Path to humann ChocoPhlAn database.   "
    echo "     -u    (Required) Path to humann UniRef90 database.     "
    echo "                                                            "
    echo " Required SLURM parameters                                  "
    echo " Note: For each parameter, a comma separated list of 8      "
    echo " entries must be given (one for each pipeline step). If not "
    echo " performing a certain step (e.g. merging of paired reads)   "
    echo " then you can put a 'NA' in the list.                       "
    echo "    -n    (Required) Names of the partitions to request for "
    echo "           submitting jobs.                                 "
    echo "     -t    (Required) Time requests for running jobs.       "
    echo "           Specify in hours (e.g. 12:00:00 for 12 hours).   "
    echo "     -k    (Required) Number of cores wanted for jobs.      "
    echo "     -m    (Required) Amount of memory requested for each   "
    echo "           requested core. Specify in Megabytes.            "
    echo "     -f    (Required) E-mail to send notifications to upon  "
    echo "           failure of any jobs.                             "
    echo "                                                            "
    echo " Optional pipeline parameters                               "
    echo "     -b    (Optional) Use Kraken2/Bracken to perform        "
    echo "           taxonomic profiling and supply path to kraken2   "
    echo "           database directory here.                         "
    echo "     -x    (Optional) Merge paired-end reads before         "
    echo "           performing the pipeline using BBMerge.           "
    echo "     -a    (Optional) Path to adapters.fa file that comes   "
    echo "           packaged with BBMerge and BBDuk. Required when   "
    echo "           merging reads.                                   "
    echo "     -j    (Optional) Join forward and reverse reads        "
    echo "           together into one file (via concatenation) prior "
    echo "           to performing low-complexity sequence filtering. "
    echo "           May want to do this if using MetaPhlAn for       "
    echo "           taxonomic profiling as this will be done anyway  "
    echo "           prior to running MetaPhlAn and low-complexity    "
    echo "           sequences caused by technical artificat are more "
    echo "           prone to be located in reverse reads while the   "
    echo "           forward reads are still good.                    "
    echo "     -s    (Optional) Skip certain steps in the pipeline if "
    echo "           need be. Provide a comma separated list of steps "
    echo "           that you wish to skip in the pipeline. List may  "
    echo "           have the values: fastqc_initial, QC, decontam,   "
    echo "           entropy_filter, fastqc_final,                    "
    echo "           taxonomic_profiling, functional_profiling        "
    echo "     -d    (Optional) Delete sequence files generated during"
    echo "           intermediate pipeline steps (merging, quality    "
    echo "           trimming/filtering, human decontamination) while "
    echo "           keeping potentially important components of each "
    echo "           step to keep record of (log files, extracted host"
    echo "           sequences, extracted low-complexity sequences).  "
    echo "           Will reorganize output folder if flag is given   "
    echo "           since the default numbering of folder will no    "
    echo "           longer make sense.                               "
    echo " "
    exit 0
    ;;
    i) SEQ_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    o) OUT_DIR=$(echo $OPTARG | sed 's#/$##')
    ;;
    p) PROG_LOAD="$OPTARG"
    ;;
    r) HOST_REF="$OPTARG"
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
    b) BRACKEN="$OPTARG"
    ;;
    x) MERGE=1
    ;;
    a) ADAPTERS="$OPTARG"
    ;;
    j) JOIN=1
    ;;
    s) SKIP="$OPTARG"
    ;;
    d) DELETE=1
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
if [[ -z "$OUT_DIR" ]]; then
  echo "ERROR: Argument -o is required, please supply an output directory"
  exit 1
fi
if [[ ! -d "$OUT_DIR" ]]; then
  echo "ERROR: Argument -o should be a directory, please supply an output directory"
  exit 1
fi

# -p
if [[ -z "$PROG_LOAD" ]]; then
  echo "ERROR: Argument -p is required, please supply a single quoted string of commands needed to load required programs (can be an empty string ' ' if none required)"
  exit 1
fi

# -r
if [[ -z "$HOST_REF" ]]; then
  echo "ERROR: Argument -r is required, please supply a host genome reference file in FASTA format"
  exit 1
fi
if [[ -d "$HOST_REF" ]]; then
  echo "ERROR: Argument -r should be a single FASTA file, please supply a host genome reference file in FASTA format"
  exit 1
fi
if ls -l $HOST_REF | grep -q ".fa"; then
  :
elif ls -l $HOST_REF | grep -q ".fna"; then
  :
elif ls -l $HOST_REF | grep -q ".fasta"; then
  :
else
  echo "ERROR: Expecting host genome reference file to be in FASTA format with extension '.fa', '.fna', or '.fasta'"
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
PARTITION_1=$(echo $PARTITION | awk -F',' '{print $1}')
PARTITION_2=$(echo $PARTITION | awk -F',' '{print $2}')
PARTITION_3=$(echo $PARTITION | awk -F',' '{print $3}')
PARTITION_4=$(echo $PARTITION | awk -F',' '{print $4}')
PARTITION_5=$(echo $PARTITION | awk -F',' '{print $5}')
PARTITION_6=$(echo $PARTITION | awk -F',' '{print $6}')
PARTITION_7=$(echo $PARTITION | awk -F',' '{print $7}')
PARTITION_8=$(echo $PARTITION | awk -F',' '{print $8}')

# -t
if [[ -z "$TIME_REQUEST" ]]; then
  echo "ERROR: Argument -t is required, please supply a max length of time to run each job for."
  exit 1
fi
TIME_1=$(echo $TIME_REQUEST | awk -F',' '{print $1}')
TIME_2=$(echo $TIME_REQUEST | awk -F',' '{print $2}')
TIME_3=$(echo $TIME_REQUEST | awk -F',' '{print $3}')
TIME_4=$(echo $TIME_REQUEST | awk -F',' '{print $4}')
TIME_5=$(echo $TIME_REQUEST | awk -F',' '{print $5}')
TIME_6=$(echo $TIME_REQUEST | awk -F',' '{print $6}')
TIME_7=$(echo $TIME_REQUEST | awk -F',' '{print $7}')
TIME_8=$(echo $TIME_REQUEST | awk -F',' '{print $8}')

# -k
if [[ -z "$CPU_REQUEST" ]]; then
  echo "ERROR: Argument -k is required, please supply the number of cores being requested to run each job."
  exit 1
fi
CPU_1=$(echo $CPU_REQUEST | awk -F',' '{print $1}')
CPU_2=$(echo $CPU_REQUEST | awk -F',' '{print $2}')
CPU_3=$(echo $CPU_REQUEST | awk -F',' '{print $3}')
CPU_4=$(echo $CPU_REQUEST | awk -F',' '{print $4}')
CPU_5=$(echo $CPU_REQUEST | awk -F',' '{print $5}')
CPU_6=$(echo $CPU_REQUEST | awk -F',' '{print $6}')
CPU_7=$(echo $CPU_REQUEST | awk -F',' '{print $7}')
CPU_8=$(echo $CPU_REQUEST | awk -F',' '{print $8}')

# -m
if [[ -z "$MEM_PER_CPU" ]]; then
  echo "ERROR: Argument -m is required, please supply a memory request for each core of each job."
  exit 1
fi
MEM_1=$(echo $MEM_PER_CPU | awk -F',' '{print $1}')
MEM_2=$(echo $MEM_PER_CPU | awk -F',' '{print $2}')
MEM_3=$(echo $MEM_PER_CPU | awk -F',' '{print $3}')
MEM_4=$(echo $MEM_PER_CPU | awk -F',' '{print $4}')
MEM_5=$(echo $MEM_PER_CPU | awk -F',' '{print $5}')
MEM_6=$(echo $MEM_PER_CPU | awk -F',' '{print $6}')
MEM_7=$(echo $MEM_PER_CPU | awk -F',' '{print $7}')
MEM_8=$(echo $MEM_PER_CPU | awk -F',' '{print $8}')

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

# -m
if [[ ! -z "$MERGE" ]]; then
  if [[ -z "$ADAPTERS" ]]; then
    echo "ERROR: when specifying the -m parameter, the -a parameter must also be specified"
    exit 1
  fi
fi

# -a
if [[ ! -z "$ADAPTERS" ]]; then
  if [[ -z "$MERGE" ]]; then
    echo "ERROR: when specifying the -a parameter, the -m parameter must also be specified"
    exit 1
  fi
  if [[ -d "$ADAPTERS" ]]; then
    echo "ERROR: Argument -a should be the path to a single file, not a directory, please supply path to the adapters.fa file that comes with BBMerge and BBDuk"
    exit 1
  fi
  if echo $ADAPTERS | grep -q -v "adapters\.fa"; then
    echo "ERROR: path given to -a does not contain the file name adapters.fa, please supply the adapters.fa file that comes with BBMerge and BBDuk to this argument"
    exit 1
  fi
fi

# -s
if [[ ! -z "$SKIP" ]]; then
  if echo $SKIP | grep -q "fastqc_initial"; then
    :
  elif echo $SKIP | grep -q "QC"; then
    :
  elif echo $SKIP | grep -q "decontam"; then
    :
  elif echo $SKIP | grep -q "entropy_filter"; then
    :
  elif echo $SKIP | grep -q "fastqc_final"; then
    :
  elif echo $SKIP | grep -q "taxonomic_profiling"; then
    :
  elif echo $SKIP | grep -q "functional_profiling"; then
    :
  else
    echo "ERROR: Invalid argument given to -s, please specify one or more of: fastqc_initial,QC,decontam,entropy_filter,taxonomic_profiling,functional_profiling"
    exit 1
  fi
fi

###### CREATE DIRECTORY FOR PIPELINE OUTPUT #####
DATE=$(date | awk '{print $3"_"$2"_"$6}')
if [ -d "${OUT_DIR}/Metagenomic_Pipeline_${DATE}" ]
then
	:
else
	mkdir ${OUT_DIR}/Metagenomic_Pipeline_${DATE}
fi

echo " "
echo "Directory for shotgun metagenomic pipeline output: Metagenomic_Pipeline_${DATE}"
echo " "
RESULTS_DIR="${OUT_DIR}/Metagenomic_Pipeline_${DATE}"

############# INITIAL FASTQC REPORT #############
if echo $SKIP | grep -q "fastqc_initial"; then
  echo " "
  echo "*** Skipping running of initial FastQC report generation on input fastq files ***"
  echo " "

  #Create directory
  if [ -d "${RESULTS_DIR}/1.FastQC_Initial_Reports" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports/0.ErrorOut
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports/0.Output
  fi
else
  SECONDS=0
  echo " "
  echo "*** Running FastQC on input fastq files ***"
  echo " "

  #Create directory for output
  if [ -d "${RESULTS_DIR}/1.FastQC_Initial_Reports" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports/0.ErrorOut
	  mkdir ${RESULTS_DIR}/1.FastQC_Initial_Reports/0.Output
  fi

  ##### Run FastQC #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_1" >> run.sh
  echo "#SBATCH --job-name=FastQC" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/1.FastQC_Initial_Reports/0.ErrorOut/FastQC_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/1.FastQC_Initial_Reports/0.Output/FastQC_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_1" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_1" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_1" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --array=1-$(ls -l ${SEQ_DIR}/*_R1*${SEQ_EXT} | wc -l)" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "FILE1=\$(ls ${SEQ_DIR}/*_R1*${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
  echo "FILE2=\$(ls ${SEQ_DIR}/*_R2*${SEQ_EXT} | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
  echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1' '{print \$1}')" >> run.sh
  echo "fastqc \$FILE1 \$FILE2 -d ${RESULTS_DIR}/1.FastQC_Initial_Reports -o ${RESULTS_DIR}/1.FastQC_Initial_Reports \\" >> run.sh
  echo "> ${RESULTS_DIR}/1.FastQC_Initial_Reports/\${FILE_NAME}.log 2>&1" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_1" >> run.sh
  echo "#SBATCH --job-name=MultiQC" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/1.FastQC_Initial_Reports/0.ErrorOut/MultiQC.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/1.FastQC_Initial_Reports/0.Output/MultiQC.out" >> run.sh
  echo "#SBATCH --time=$TIME_1" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_1" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_1" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "multiqc ${RESULTS_DIR}/1.FastQC_Initial_Reports -o ${RESULTS_DIR}/1.FastQC_Initial_Reports \\" >> run.sh
  echo "> ${RESULTS_DIR}/1.FastQC_Initial_Reports/multiqc.log 2>&1" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  #Signal jobs have ended
  echo "Initial FastQC reports complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

######### MERGE PAIRED READS WITH BBMERGE #######
if [[ -z "$MERGE" ]]; then
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
else
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

  # get max memory to use for jobs
  MAX_MEM=$(( MEM_2 * CPU_2 / 1000))

  ##### Run BBMerge #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_2" >> run.sh
  echo "#SBATCH --job-name=Merge" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/2.Merged_Paired_End_Sequences/0.ErrorOut/Merge_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/2.Merged_Paired_End_Sequences/0.Output/Merge_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_2" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_2" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_2" >> run.sh
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
  echo "rem iterations=5 extend2=20 ecct t=$CPU_2 -Xmx${MAX_MEM}g \\" >> run.sh
  echo "> ${RESULTS_DIR}/2.Merged_Paired_End_Sequences/\${FILE_NAME}.log 2>&1"  >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  #Signal jobs have ended
  echo "Merging of paired end reads with BBMerge complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

################# QC WITH BBDUK #################
if echo $SKIP | grep -q "QC"; then
  echo " "
  echo "*** Skipping running of BBDuk for adapter/quality trimming and filtering of input fastq files ***"
  echo " "

  #Create directory
  if [ -d "${RESULTS_DIR}/3.Quality_Controlled_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/3.Quality_Controlled_Sequences
	  mkdir ${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.Output
  fi
else
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

  # get max memory to use for jobs
  MAX_MEM=$(( MEM_3 * CPU_3 / 1000))

  ##### Run BBDuk #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_3" >> run.sh
  echo "#SBATCH --job-name=QC" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.ErrorOut/QC_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/3.Quality_Controlled_Sequences/0.Output/QC_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_3" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_3" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_3" >> run.sh
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
fi
#################################################

##### REMOVAL OF HOST READS WITH BBMAP/BBSPLIT ######
if echo $SKIP | grep -q "decontam"; then
  echo " "
  echo "*** Skipping running of BBMap/BBSplit for removal of contaminant host reads ***"
  echo " "

  #Create directory
  if [ -d "${RESULTS_DIR}/4.Decontaminated_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/4.Decontaminated_Sequences
	  mkdir ${RESULTS_DIR}/4.Decontaminated_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/4.Decontaminated_Sequences/0.Output
  fi
else
  SECONDS=0
  echo " "
  echo "*** Running BBMap/BBSplit to perform removal of contaminant host reads ***"
  echo " "

  #Create directory
  if [ -d "${RESULTS_DIR}/4.Decontaminated_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/4.Decontaminated_Sequences
	  mkdir ${RESULTS_DIR}/4.Decontaminated_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/4.Decontaminated_Sequences/0.Output
  fi

  # get max memory to use for jobs
  MAX_MEM=$(( MEM_4 * CPU_4 / 1000))

  ##### Run BBMap/BBSplit #####
  #Index given reference genome file
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_4" >> run.sh
  echo "#SBATCH --job-name=Host_Ref_Indexing" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/4.Decontaminated_Sequences/0.ErrorOut/Host_Ref_Indexing.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/4.Decontaminated_Sequences/0.Output/Host_Ref_Indexing.out" >> run.sh
  echo "#SBATCH --time=$TIME_4" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_4" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_4" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "bbsplit.sh ref=${HOST_REF} path=${RESULTS_DIR}/4.Decontaminated_Sequences \\" >> run.sh
  echo "t=$CPU_4 -Xmx${MAX_MEM}g \\" >> run.sh
  echo "> ${RESULTS_DIR}/4.Decontaminated_Sequences/Host_Ref_Indexing.log 2>&1" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_4" >> run.sh
  echo "#SBATCH --job-name=Decontam" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/4.Decontaminated_Sequences/0.ErrorOut/Decontam_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/4.Decontaminated_Sequences/0.Output/Decontam_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_4" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_4" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_4" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/3.Quality_Controlled_Sequences/*.fastq.gz | wc -l)" >> run.sh
  else
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/3.Quality_Controlled_Sequences/*R1.fastq.gz | wc -l)" >> run.sh
  fi
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "FILE=\$(ls ${RESULTS_DIR}/3.Quality_Controlled_Sequences/*.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$FILE | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> run.sh
    echo "bbsplit.sh in=\${FILE} \\" >> run.sh
    echo "out=${RESULTS_DIR}/4.Decontaminated_Sequences/\${FILE_NAME}.fastq.gz \\" >> run.sh
  else
    echo "FILE1=\$(ls ${RESULTS_DIR}/3.Quality_Controlled_Sequences/*R1.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE2=\$(ls ${RESULTS_DIR}/3.Quality_Controlled_Sequences/*R2.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1' '{print \$1}')" >> run.sh
    echo "bbsplit.sh in1=\${FILE1} \\" >> run.sh
    echo "in2=\${FILE2} \\" >> run.sh
    echo "outu1=${RESULTS_DIR}/4.Decontaminated_Sequences/\${FILE_NAME}_R1.fastq.gz \\" >> run.sh
    echo "outu2=${RESULTS_DIR}/4.Decontaminated_Sequences/\${FILE_NAME}_R2.fastq.gz \\" >> run.sh
  fi
  echo "path=${RESULTS_DIR}/4.Decontaminated_Sequences \\" >> run.sh
  echo "basename=${RESULTS_DIR}/4.Decontaminated_Sequences/\${FILE_NAME}.%_contam_#.fastq.gz \\" >> run.sh
  echo "refstats=${RESULTS_DIR}/4.Decontaminated_Sequences/\${FILE_NAME}_refstats.log \\" >> run.sh
  echo "t=$CPU_4 -Xmx${MAX_MEM}g \\" >> run.sh
  echo "> ${RESULTS_DIR}/4.Decontaminated_Sequences/\${FILE_NAME}.log 2>&1" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh
  rm -r ${RESULTS_DIR}/4.Decontaminated_Sequences/ref

  #Move extracted host sequences to their own directory
  mkdir ${RESULTS_DIR}/4.Decontaminated_Sequences/Extracted_Host_Sequences
  mv ${RESULTS_DIR}/4.Decontaminated_Sequences/*contam* ${RESULTS_DIR}/4.Decontaminated_Sequences/Extracted_Host_Sequences/

  #Signal jobs have ended
  echo "Running of BBMap/BBSplit complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

################# LOW COMPLEXITY SEQUENCE FILTERING WITH BBDUK #################
if echo $SKIP | grep -q "entropy_filter"; then
  echo " "
  echo "*** Skipping running of BBDuk for low complexity sequence filtering ***"
  echo " "

  #Create directory
  if [ -d "${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences
	  mkdir ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/0.Output
  fi
else
  SECONDS=0
  echo " "
  echo "*** Running BBDuk to perform low complexity sequence filtering ***"
  echo " "

  #Create directory for output
  if [ -d "${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences
	  mkdir ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/0.ErrorOut
	  mkdir ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/0.Output
    mkdir ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/Removed_Sequences
  fi

  # get max memory to use for jobs
  MAX_MEM=$(( MEM_5 * CPU_5 / 1000))

  ##### Run BBDuk #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_5" >> run.sh
  echo "#SBATCH --job-name=entropy_filter" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/0.ErrorOut/entropy_filter_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/0.Output/entropy_filter_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_5" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_5" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_5" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/4.Decontaminated_Sequences/*.fastq.gz | wc -l)" >> run.sh
  else
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/4.Decontaminated_Sequences/*R1.fastq.gz | wc -l)" >> run.sh
  fi
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  if [[ ! -z "$MERGE" ]]; then
    echo "FILE=\$(ls ${RESULTS_DIR}/4.Decontaminated_Sequences/*.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$FILE | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> run.sh
    echo "bbduk.sh in=\${FILE} \\" >> run.sh
    echo "out=${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/\${FILE_NAME}.fastq.gz \\" >> run.sh
    echo "outm=${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/Removed_Sequences/\${FILE_NAME}.fastq.gz \\" >> run.sh
    echo "entropy=0.01 entropywindow=50 entropyk=5 -Xmx${MAX_MEM}g \\" >> run.sh
    echo "> ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/\${FILE_NAME}.log 2>&1" >> run.sh
  elif [[ ! -z "$JOIN" ]]; then
    echo "FILE1=\$(ls ${RESULTS_DIR}/4.Decontaminated_Sequences/*R1.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE2=\$(ls ${RESULTS_DIR}/4.Decontaminated_Sequences/*R2.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1' '{print \$1}')" >> run.sh
    echo "zcat \$FILE1 \$FILE2 > ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/\${FILE_NAME}.temp.fastq" >> run.sh
    echo "FILE=${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/\${FILE_NAME}.temp.fastq" >> run.sh
    echo "bbduk.sh in=\${FILE} \\" >> run.sh
    echo "out=${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/\${FILE_NAME}.fastq.gz \\" >> run.sh
    echo "outm=${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/Removed_Sequences/\${FILE_NAME}.fastq.gz \\" >> run.sh
    echo "entropy=0.01 entropywindow=50 entropyk=5 -Xmx${MAX_MEM}g \\" >> run.sh
    echo "> ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/\${FILE_NAME}.log 2>&1" >> run.sh
    echo "rm \$FILE" >> run.sh
  else
    echo "FILE1=\$(ls ${RESULTS_DIR}/4.Decontaminated_Sequences/*R1.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE2=\$(ls ${RESULTS_DIR}/4.Decontaminated_Sequences/*R2.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1' '{print \$1}')" >> run.sh
    echo "bbduk.sh in=\$FILE1 \\" >> run.sh
    echo "in2=\$FILE2 \\" >> run.sh
    echo "out=${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/\${FILE_NAME}_R1.fastq.gz \\" >> run.sh
    echo "out2=${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/\${FILE_NAME}_R2.fastq.gz \\" >> run.sh
    echo "outm=${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/Removed_Sequences/\${FILE_NAME}_R1.fastq.gz \\" >> run.sh
    echo "outm2=${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/Removed_Sequences/\${FILE_NAME}_R2.fastq.gz \\" >> run.sh
    echo "entropy=0.01 entropywindow=50 entropyk=5 -Xmx${MAX_MEM}g \\" >> run.sh
    echo "> ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/\${FILE_NAME}.log 2>&1" >> run.sh
  fi
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  #Signal jobs have ended
  echo "Low complexity sequence filtering with BBDuk complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

############# FINAL FASTQC REPORT #############
if echo $SKIP | grep -q "fastqc_final"; then
  echo " "
  echo "*** Skipping running of final FastQC report generation on quality controlled fastq files ***"
  echo " "

  #Create directory
  if [ -d "${RESULTS_DIR}/6.FastQC_Final_Reports" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/6.FastQC_Final_Reports
	  mkdir ${RESULTS_DIR}/6.FastQC_Final_Reports/0.ErrorOut
	  mkdir ${RESULTS_DIR}/6.FastQC_Final_Reports/0.Output
  fi
else
  SECONDS=0
  echo " "
  echo "*** Running final FastQC report on QCed fastq files ***"
  echo " "

  #Create directory for output
  if [ -d "${RESULTS_DIR}/6.FastQC_Final_Reports" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/6.FastQC_Final_Reports
	  mkdir ${RESULTS_DIR}/6.FastQC_Final_Reports/0.ErrorOut
	  mkdir ${RESULTS_DIR}/6.FastQC_Final_Reports/0.Output
  fi

  ##### Run FastQC #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_6" >> run.sh
  echo "#SBATCH --job-name=FastQC" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/6.FastQC_Final_Reports/0.ErrorOut/FastQC_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/6.FastQC_Final_Reports/0.Output/FastQC_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_6" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_6" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_6" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  if [[ ! -z "$MERGE" ]] || [[ ! -z "$JOIN" ]]; then
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*.fastq.gz | wc -l)" >> run.sh
  else
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*R1.fastq.gz | wc -l)" >> run.sh
  fi
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  if [[ ! -z "$MERGE" ]] || [[ ! -z "$JOIN" ]]; then
    echo "FILE=\$(ls ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$FILE | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> run.sh
    echo "fastqc \$FILE -d ${RESULTS_DIR}/6.FastQC_Final_Reports -o ${RESULTS_DIR}/6.FastQC_Final_Reports \\" >> run.sh
  else
    echo "FILE1=\$(ls ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*R1.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE2=\$(ls ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*R2.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1' '{print \$1}')" >> run.sh
    echo "fastqc \$FILE1 \$FILE2 -d ${RESULTS_DIR}/6.FastQC_Final_Reports -o ${RESULTS_DIR}/6.FastQC_Final_Reports \\" >> run.sh
  fi
  echo "> ${RESULTS_DIR}/6.FastQC_Final_Reports/\${FILE_NAME}.log 2>&1" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_6" >> run.sh
  echo "#SBATCH --job-name=MultiQC" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/6.FastQC_Final_Reports/0.ErrorOut/MultiQC.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/6.FastQC_Final_Reports/0.Output/MultiQC.out" >> run.sh
  echo "#SBATCH --time=$TIME_6" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_6" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_6" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  echo "multiqc ${RESULTS_DIR}/6.FastQC_Final_Reports -o ${RESULTS_DIR}/6.FastQC_Final_Reports \\" >> run.sh
  echo "> ${RESULTS_DIR}/6.FastQC_Final_Reports/multiqc.log 2>&1" >> run.sh
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  #Signal jobs have ended
  echo "Final FastQC reports complete"
  echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
  echo " "
fi
#################################################

###### TAXONOMIC PROFILING #######
if echo $SKIP | grep -q "taxonomic_profiling"; then
  echo " "
  echo "*** Skipping running of taxonomic profiling ***"
  echo " "

  #Create directory for profiling output
  if [ -d "${RESULTS_DIR}/7.Taxonomic_Profiling" ]
  then
    :
  else
    mkdir ${RESULTS_DIR}/7.Taxonomic_Profiling
    mkdir ${RESULTS_DIR}/7.Taxonomic_Profiling/0.ErrorOut
    mkdir ${RESULTS_DIR}/7.Taxonomic_Profiling/0.Output
  fi
else
  #Create directory for profiling output
  if [ -d "${RESULTS_DIR}/7.Taxonomic_Profiling" ]
  then
    :
  else
    mkdir ${RESULTS_DIR}/7.Taxonomic_Profiling
    mkdir ${RESULTS_DIR}/7.Taxonomic_Profiling/0.ErrorOut
    mkdir ${RESULTS_DIR}/7.Taxonomic_Profiling/0.Output
  fi

  if [[ -z $BRACKEN ]]; then
    SECONDS=0
    echo " "
    echo "*** Running MetaPhlAn workflow for performing taxonomic profiling ***"
    echo " "

    ##### Run MetaPhlAn workflow #####
    #Create script for running program and submit
    echo '#!/bin/bash' > run.sh
    echo "#SBATCH --partition=$PARTITION_7" >> run.sh
    echo "#SBATCH --job-name=Profiling" >> run.sh
    echo "#SBATCH --error=${RESULTS_DIR}/7.Taxonomic_Profiling/0.ErrorOut/Taxonomic_Profiling_%A_%a.err" >> run.sh
    echo "#SBATCH --output=${RESULTS_DIR}/7.Taxonomic_Profiling/0.Output/Taxonomic_Profiling_%A_%a.out" >> run.sh
    echo "#SBATCH --time=$TIME_7" >> run.sh
    echo "#SBATCH --ntasks=1" >> run.sh
    echo "#SBATCH --cpus-per-task=$CPU_7" >> run.sh
    echo "#SBATCH --mem-per-cpu=$MEM_7" >> run.sh
    echo "#SBATCH --mail-type=FAIL" >> run.sh
    echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
    if [[ ! -z "$MERGE" ]] || [[ ! -z "$JOIN" ]]; then
      echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*.fastq.gz | wc -l)" >> run.sh
    else
      echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*R1.fastq.gz | wc -l)" >> run.sh
    fi
    echo "#SBATCH --wait" >> run.sh
    echo "$PROG_LOAD" >> run.sh
    if [[ ! -z "$MERGE" ]] || [[ ! -z "$JOIN" ]]; then
      echo "FILE=\$(ls ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
      echo "FILE_NAME=\$(echo \$FILE | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> run.sh
    else
      echo "FILE1=\$(ls ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*R1.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
      echo "FILE2=\$(ls ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*R2.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
      echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1' '{print \$1}')" >> run.sh
      echo "zcat \$FILE1 \$FILE2 > ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}.temp.fastq" >> run.sh
      echo "FILE=${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}.temp.fastq" >> run.sh
    fi
    echo "if [ -d \"${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles\" ]" >> run.sh
    echo "then" >> run.sh
    echo "  :" >> run.sh
    echo "else" >> run.sh
    echo "  mkdir ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles" >> run.sh
    echo "fi" >> run.sh
    echo "metaphlan \$FILE \\" >> run.sh
    echo "--input_type fastq \\" >> run.sh
    echo "-t rel_ab \\" >> run.sh
    echo "--nproc $CPU_7 \\" >> run.sh
    echo "--bowtie2out ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles/\${FILE_NAME}_metaphlan_bowtie2.txt \\" >> run.sh
    echo "-o ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles/\${FILE_NAME}_metaphlan_rel_ab.tsv" >> run.sh
    if [[ -z "$MERGE" ]] && [[ -z "$JOIN" ]]; then
      echo "rm ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}.temp.fastq" >> run.sh
    fi
    chmod +x run.sh

    sbatch run.sh > /dev/null

    rm run.sh

    ##### Calculate taxonomic count data #####
    #Create script for running program and submit
    echo '#!/bin/bash' > run.sh
    echo "#SBATCH --partition=$PARTITION_7" >> run.sh
    echo "#SBATCH --job-name=get_taxa_counts" >> run.sh
    echo "#SBATCH --error=${RESULTS_DIR}/7.Taxonomic_Profiling/0.ErrorOut/get_taxa_counts_%A_%a.err" >> run.sh
    echo "#SBATCH --output=${RESULTS_DIR}/7.Taxonomic_Profiling/0.Output/get_taxa_counts_%A_%a.out" >> run.sh
    echo "#SBATCH --time=$TIME_7" >> run.sh
    echo "#SBATCH --ntasks=1" >> run.sh
    echo "#SBATCH --cpus-per-task=$CPU_7" >> run.sh
    echo "#SBATCH --mem-per-cpu=$MEM_7" >> run.sh
    echo "#SBATCH --mail-type=FAIL" >> run.sh
    echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
    echo "#SBATCH --array=1-$(ls -d -l ${RESULTS_DIR}/7.Taxonomic_Profiling/*Profiles | wc -l)" >> run.sh
    echo "#SBATCH --wait" >> run.sh
    echo "$PROG_LOAD" >> run.sh
    echo "DIR=\$(ls -d ${RESULTS_DIR}/7.Taxonomic_Profiling/*Profiles | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$DIR | awk -F '/' '{print \$NF}' | awk -F '_Profiles' '{print \$1}')" >> run.sh
    echo "metaphlan \${DIR}/\${FILE_NAME}_metaphlan_bowtie2.txt \\" >> run.sh
    echo "--input_type bowtie2out \\" >> run.sh
    echo "-t rel_ab \\" >> run.sh
    echo "--unknown_estimation \\" >> run.sh
    echo "--nproc $CPU_7 \\" >> run.sh
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
    if [ -d "${RESULTS_DIR}/7.Taxonomic_Profiling/Merged_Sample_Tables" ]
    then
  	  :
    else
      mkdir ${RESULTS_DIR}/7.Taxonomic_Profiling/Merged_Sample_Tables
    fi

    #For each sample get counts and relative abundances
    for dir in ${RESULTS_DIR}/7.Taxonomic_Profiling/*_Profiles; do
      FILE_NAME=$(echo $dir | awk -F '/' '{print $NF}' | awk -F '_Profiles' '{print $1}')

      cp ${dir}/${FILE_NAME}_metaphlan_rel_ab.tsv ${RESULTS_DIR}/7.Taxonomic_Profiling/

      awk -F'\t' '{print $1,$2,$3}' OFS='\t' ${dir}/${FILE_NAME}_metaphlan_rel_ab_w_counts.tsv \
      > ${RESULTS_DIR}/7.Taxonomic_Profiling/${FILE_NAME}_metaphlan_rel_ab_w_unknown.tsv

      awk -F'\t' '{print $1,$2,$4}' OFS='\t' ${dir}/${FILE_NAME}_metaphlan_rel_ab_w_counts.tsv | \
      sed 's/counts/relative_abundance/' \
      > ${RESULTS_DIR}/7.Taxonomic_Profiling/${FILE_NAME}_metaphlan_counts.tsv
    done

    #Create individual scripts for running table merging
    echo '#!/bin/bash' > run1.sh
    echo "$PROG_LOAD" >> run1.sh
    echo "merge_metaphlan_tables.py ${RESULTS_DIR}/7.Taxonomic_Profiling/*_metaphlan_rel_ab.tsv \\" >> run1.sh
    echo "> ${RESULTS_DIR}/7.Taxonomic_Profiling/Merged_Sample_Tables/metaphlan_rel_ab.tsv" >> run1.sh
    echo "sed -i '2s/_metaphlan_rel_ab//g' ${RESULTS_DIR}/7.Taxonomic_Profiling/Merged_Sample_Tables/metaphlan_rel_ab.tsv" >> run1.sh
    echo "rm ${RESULTS_DIR}/7.Taxonomic_Profiling/*_metaphlan_rel_ab.tsv" >> run1.sh
    chmod +x run1.sh

    echo '#!/bin/bash' > run2.sh
    echo "$PROG_LOAD" >> run2.sh
    echo "merge_metaphlan_tables.py ${RESULTS_DIR}/7.Taxonomic_Profiling/*_metaphlan_rel_ab_w_unknown.tsv \\" >> run2.sh
    echo "> ${RESULTS_DIR}/7.Taxonomic_Profiling/Merged_Sample_Tables/metaphlan_rel_ab_w_unknown.tsv" >> run2.sh
    echo "sed -i '2s/_metaphlan_rel_ab_w_unknown//g' ${RESULTS_DIR}/7.Taxonomic_Profiling/Merged_Sample_Tables/metaphlan_rel_ab_w_unknown.tsv" >> run2.sh
    echo "rm ${RESULTS_DIR}/7.Taxonomic_Profiling/*_metaphlan_rel_ab_w_unknown.tsv" >> run2.sh
    chmod +x run2.sh

    echo '#!/bin/bash' > run3.sh
    echo "$PROG_LOAD" >> run3.sh
    echo "merge_metaphlan_tables.py ${RESULTS_DIR}/7.Taxonomic_Profiling/*_metaphlan_counts.tsv \\" >> run3.sh
    echo "> ${RESULTS_DIR}/7.Taxonomic_Profiling/Merged_Sample_Tables/metaphlan_counts.tsv" >> run3.sh
    echo "sed -i '2s/_metaphlan_counts//g' ${RESULTS_DIR}/7.Taxonomic_Profiling/Merged_Sample_Tables/metaphlan_counts.tsv" >> run3.sh
    echo "rm ${RESULTS_DIR}/7.Taxonomic_Profiling/*_metaphlan_counts.tsv" >> run3.sh
    chmod +x run3.sh

    #Create script for running and submitting individual scripts
    echo '#!/bin/bash' > run.sh
    echo "#SBATCH --partition=$PARTITION_7" >> run.sh
    echo "#SBATCH --job-name=Merge_tables" >> run.sh
    echo "#SBATCH --error=${RESULTS_DIR}/7.Taxonomic_Profiling/0.ErrorOut/Merge_tables_%A_%a.err" >> run.sh
    echo "#SBATCH --output=${RESULTS_DIR}/7.Taxonomic_Profiling/0.Output/Merge_tables_%A_%a.out" >> run.sh
    echo "#SBATCH --time=$TIME_7" >> run.sh
    echo "#SBATCH --ntasks=1" >> run.sh
    echo "#SBATCH --cpus-per-task=$CPU_7" >> run.sh
    echo "#SBATCH --mem-per-cpu=$MEM_7" >> run.sh
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
    echo "#SBATCH --partition=$PARTITION_7" >> run.sh
    echo "#SBATCH --job-name=Profiling" >> run.sh
    echo "#SBATCH --error=${RESULTS_DIR}/7.Taxonomic_Profiling/0.ErrorOut/Taxonomic_Profiling_%A_%a.err" >> run.sh
    echo "#SBATCH --output=${RESULTS_DIR}/7.Taxonomic_Profiling/0.Output/Taxonomic_Profiling_%A_%a.out" >> run.sh
    echo "#SBATCH --time=$TIME_7" >> run.sh
    echo "#SBATCH --ntasks=1" >> run.sh
    echo "#SBATCH --cpus-per-task=$CPU_7" >> run.sh
    echo "#SBATCH --mem-per-cpu=$MEM_7" >> run.sh
    echo "#SBATCH --mail-type=FAIL" >> run.sh
    echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
    if [[ ! -z "$MERGE" ]] || [[ ! -z "$JOIN" ]]; then
      echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*.fastq.gz | wc -l)" >> run.sh
    else
      echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*R1.fastq.gz | wc -l)" >> run.sh
    fi
    echo "#SBATCH --wait" >> run.sh
    echo "$PROG_LOAD" >> run.sh
    if [[ ! -z "$MERGE" ]] || [[ ! -z "$JOIN" ]]; then
      echo "FILE=\$(ls ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
      echo "FILE_NAME=\$(echo \$FILE | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> run.sh
    else
      echo "FILE1=\$(ls ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*R1.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
      echo "FILE2=\$(ls ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*R2.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
      echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1' '{print \$1}')" >> run.sh
    fi
    echo "if [ -d \"${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles\" ]" >> run.sh
    echo "then" >> run.sh
    echo "  :" >> run.sh
    echo "else" >> run.sh
    echo "  mkdir ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles" >> run.sh
    echo "fi" >> run.sh
    if [[ ! -z "$MERGE" ]] || [[ ! -z "$JOIN" ]]; then
      echo "kraken2 \$FILE \\" >> run.sh
    else
      echo "kraken2 \$FILE1 \$FILE2 \\" >> run.sh
      echo "--paired \\" >> run.sh
    fi
    echo "--db $BRACKEN \\" >> run.sh
    echo "--threads $CPU_7 \\" >> run.sh
    echo "--confidence 0.5 \\" >> run.sh
    echo "--minimum-base-quality 0 \\" >> run.sh
    echo "--minimum-hit-groups 3 \\" >> run.sh
    echo "--report ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles/\${FILE_NAME}_kraken2_report.tsv \\" >> run.sh
    echo "--report-zero-counts \\" >> run.sh
    echo "--use-names \\" >> run.sh
    echo "--gzip-compressed \\" >> run.sh
    echo "--unclassified-out ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles/\${FILE_NAME}_kraken2_unclassified#.fastq \\" >> run.sh
    echo "--output ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles/\${FILE_NAME}_kraken2_classified.tsv \\" >> run.sh
    echo "> ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles/\${FILE_NAME}.log 2>&1" >> run.sh
    echo "bracken -d $BRACKEN \\" >> run.sh
    echo "-i ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles/\${FILE_NAME}_kraken2_report.tsv \\" >> run.sh
    echo "-r 50 -l S -t 30 \\" >> run.sh
    echo "-w ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles/\${FILE_NAME}_bracken_report.tsv \\" >> run.sh
    echo "-o ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles/\${FILE_NAME}_bracken_counts.tsv \\" >> run.sh
    echo ">> ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles/\${FILE_NAME}.log 2>&1" >> run.sh
    echo "kreport2mpa.py -r ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles/\${FILE_NAME}_bracken_report.tsv \\" >> run.sh
    echo "-o ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles/\${FILE_NAME}_bracken_counts_mpa.tsv" >> run.sh
    echo "sed -i \"1s/^/#clade_name\t\${FILE_NAME}\n/\" ${RESULTS_DIR}/7.Taxonomic_Profiling/\${FILE_NAME}_Profiles/\${FILE_NAME}_bracken_counts_mpa.tsv" >> run.sh
    chmod +x run.sh

    sbatch run.sh > /dev/null

    rm run.sh

    ##### Merge tables #####
    #Create directory for merged tables
    if [ -d "${RESULTS_DIR}/7.Taxonomic_Profiling/Merged_Sample_Tables" ]
    then
  	  :
    else
      mkdir ${RESULTS_DIR}/7.Taxonomic_Profiling/Merged_Sample_Tables
    fi

    #For each sample get counts
    for dir in ${RESULTS_DIR}/7.Taxonomic_Profiling/*_Profiles; do
      FILE_NAME=$(echo $dir | awk -F '/' '{print $NF}' | awk -F '_Profiles' '{print $1}')
      cp ${dir}/${FILE_NAME}_bracken_counts_mpa.tsv ${RESULTS_DIR}/7.Taxonomic_Profiling/
    done
    
    #Create script for running program and submit
    echo '#!/bin/bash' > run.sh
    echo "#SBATCH --partition=$PARTITION_7" >> run.sh
    echo "#SBATCH --job-name=Merge_tables" >> run.sh
    echo "#SBATCH --error=${RESULTS_DIR}/7.Taxonomic_Profiling/0.ErrorOut/Merge_tables_%A_%a.err" >> run.sh
    echo "#SBATCH --output=${RESULTS_DIR}/7.Taxonomic_Profiling/0.Output/Merge_tables_%A_%a.out" >> run.sh
    echo "#SBATCH --time=$TIME_7" >> run.sh
    echo "#SBATCH --ntasks=1" >> run.sh
    echo "#SBATCH --cpus-per-task=$CPU_7" >> run.sh
    echo "#SBATCH --mem-per-cpu=$MEM_7" >> run.sh
    echo "#SBATCH --mail-type=FAIL" >> run.sh
    echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
    echo "#SBATCH --wait" >> run.sh
    echo "$PROG_LOAD" >> run.sh
    echo "combine_mpa.py -i ${RESULTS_DIR}/7.Taxonomic_Profiling/*_bracken_counts_mpa.tsv \\" >> run.sh
    echo "-o ${RESULTS_DIR}/7.Taxonomic_Profiling/Merged_Sample_Tables/bracken_counts_mpa.tsv" >> run.sh
    echo "sed -i '1s/#Classification/clade_name/' ${RESULTS_DIR}/7.Taxonomic_Profiling/Merged_Sample_Tables/bracken_counts_mpa.tsv" >> run.sh
    echo "sed -i \"1s/^/#$(echo ${BRACKEN} | awk -F'/' '{print $NF}')\n/\" ${RESULTS_DIR}/7.Taxonomic_Profiling/Merged_Sample_Tables/bracken_counts_mpa.tsv" >> run.sh
    echo "rm ${RESULTS_DIR}/7.Taxonomic_Profiling/*_bracken_counts_mpa.tsv" >> run.sh
    chmod +x run.sh

    sbatch run.sh > /dev/null

    rm run.sh

    #Signal workflow has completed
    echo "Running of Kraken2/Bracken workflow complete"
    echo "Elapsed time: $(($SECONDS / 3600)) hr : $(($(($SECONDS % 3600)) / 60)) min : $(($SECONDS % 60)) sec"
    echo " "
  fi
fi
#################################################

###### FUNCTIONAL PROFILING #######
if echo $SKIP | grep -q "functional_profiling"; then
  echo " "
  echo "*** Skipping running of HUMAnN workflow for performing functional profiling ***"
  echo " "
else
  SECONDS=0
  echo " "
  echo "*** Running HUMAnN workflow for performing functional profiling ***"
  echo " "

  #Create directory for profiling output
  if [ -d "${RESULTS_DIR}/8.Functional_Profiling" ]
  then
	  :
  else
	  mkdir ${RESULTS_DIR}/8.Functional_Profiling
	  mkdir ${RESULTS_DIR}/8.Functional_Profiling/0.ErrorOut
	  mkdir ${RESULTS_DIR}/8.Functional_Profiling/0.Output
  fi

  ##### Run HUMAnN workflow #####
  #Create script for running program and submit
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_8" >> run.sh
  echo "#SBATCH --job-name=Profiling" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/8.Functional_Profiling/0.ErrorOut/Functional_profiling_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/8.Functional_Profiling/0.Output/Functional_profiling_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_8" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_8" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_8" >> run.sh
  echo "#SBATCH --mail-type=FAIL" >> run.sh
  echo "#SBATCH --mail-user=${FAIL_EMAIL}" >> run.sh
  if [[ ! -z "$MERGE" ]] || [[ ! -z "$JOIN" ]]; then
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*.fastq.gz | wc -l)" >> run.sh
  else
    echo "#SBATCH --array=1-$(ls -l ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*R1.fastq.gz | wc -l)" >> run.sh
  fi
  echo "#SBATCH --wait" >> run.sh
  echo "$PROG_LOAD" >> run.sh
  if [[ ! -z "$MERGE" ]] || [[ ! -z "$JOIN" ]]; then
    echo "FILE=\$(ls ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$FILE | awk -F '/' '{print \$NF}' | awk -F '.fastq.gz' '{print \$1}')" >> run.sh
  else
    echo "FILE1=\$(ls ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*R1.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE2=\$(ls ${RESULTS_DIR}/5.Low_Complexity_Filtered_Sequences/*R2.fastq.gz | sed -n \${SLURM_ARRAY_TASK_ID}p)" >> run.sh
    echo "FILE_NAME=\$(echo \$FILE1 | awk -F '/' '{print \$NF}' | awk -F '_R1' '{print \$1}')" >> run.sh
    echo "zcat \$FILE1 \$FILE2 > ${RESULTS_DIR}/8.Functional_Profiling/\${FILE_NAME}.temp.fastq" >> run.sh
    echo "FILE=${RESULTS_DIR}/8.Functional_Profiling/\${FILE_NAME}.temp.fastq" >> run.sh
  fi
  echo "humann --input \$FILE \\" >> run.sh
  echo "--output ${RESULTS_DIR}/8.Functional_Profiling \\" >> run.sh
  echo "--output-basename \$FILE_NAME \\" >> run.sh
  echo "--nucleotide-database $CHOCO \\" >> run.sh
  echo "--protein-database $UNIREF \\" >> run.sh
  echo "--metaphlan-options '-t rel_ab' \\" >> run.sh
  echo "--prescreen-threshold 0.01 \\" >> run.sh
  echo "--threads $CPU_8 \\" >> run.sh
  echo "--verbose \\" >> run.sh
  echo "> ${RESULTS_DIR}/8.Functional_Profiling/\${FILE_NAME}.log 2>&1" >> run.sh
  echo "rm ${RESULTS_DIR}/8.Functional_Profiling/\${FILE_NAME}_humann_temp/*_bowtie2_*" >> run.sh
  echo "rm ${RESULTS_DIR}/8.Functional_Profiling/\${FILE_NAME}_humann_temp/*_diamond_*" >> run.sh
  echo "rm ${RESULTS_DIR}/8.Functional_Profiling/\${FILE_NAME}_humann_temp/*_custom_chocophlan_database.ffn" >> run.sh
  if [[ -z "$MERGE" ]] && [[ -z "$JOIN" ]]; then
    echo "rm ${RESULTS_DIR}/8.Functional_Profiling/\${FILE_NAME}.temp.fastq" >> run.sh
  fi
  chmod +x run.sh

  sbatch run.sh > /dev/null

  rm run.sh

  ##### Clean up #####
  echo "Cleaning things up... merging per sample tables and removing unneeded files"
  echo " "
  for dir in ${RESULTS_DIR}/8.Functional_Profiling/*humann_temp; do
    FILE_NAME=$(echo $dir | awk -F '/' '{print $NF}' | awk -F '_humann_temp' '{print $1}')

    mkdir ${RESULTS_DIR}/8.Functional_Profiling/${FILE_NAME}_Profiles

    mv ${dir}/${FILE_NAME}_metaphlan_bugs_list.tsv \
    ${RESULTS_DIR}/8.Functional_Profiling/${FILE_NAME}_Profiles/

    mv ${dir}/${FILE_NAME}_metaphlan_bowtie2.txt \
    ${RESULTS_DIR}/8.Functional_Profiling/${FILE_NAME}_Profiles/

    mv ${dir}/${FILE_NAME}.log \
    ${RESULTS_DIR}/8.Functional_Profiling/${FILE_NAME}_Profiles/${FILE_NAME}_humann.log

    rm -r ${dir}

    mv ${RESULTS_DIR}/8.Functional_Profiling/${FILE_NAME}*.* \
    ${RESULTS_DIR}/8.Functional_Profiling/${FILE_NAME}_Profiles/
  done

  ##### Merge tables #####
  mkdir ${RESULTS_DIR}/8.Functional_Profiling/Merged_Sample_Tables

  for dir in ${RESULTS_DIR}/8.Functional_Profiling/*_Profiles; do
    FILE_NAME=$(echo $dir | awk -F '/' '{print $NF}' | awk -F '_Profiles' '{print $1}')

    cp ${dir}/${FILE_NAME}_metaphlan_bugs_list.tsv ${RESULTS_DIR}/8.Functional_Profiling/

    cp ${dir}/${FILE_NAME}_genefamilies.tsv ${RESULTS_DIR}/8.Functional_Profiling/

    cp ${dir}/${FILE_NAME}_pathabundance.tsv ${RESULTS_DIR}/8.Functional_Profiling/

    cp ${dir}/${FILE_NAME}_pathcoverage.tsv ${RESULTS_DIR}/8.Functional_Profiling/
  done

  #Create individual scripts for running table merging
  echo '#!/bin/bash' > run1.sh
  echo "$PROG_LOAD" >> run1.sh
  echo "merge_metaphlan_tables.py ${RESULTS_DIR}/8.Functional_Profiling/*_metaphlan_bugs_list.tsv \\" >> run1.sh
  echo "> ${RESULTS_DIR}/8.Functional_Profiling/Merged_Sample_Tables/metaphlan_bugs_list.tsv" >> run1.sh
  echo "sed -i '2s/_metaphlan_bugs_list//g' ${RESULTS_DIR}/8.Functional_Profiling/Merged_Sample_Tables/metaphlan_bugs_list.tsv" >> run1.sh
  echo "rm ${RESULTS_DIR}/8.Functional_Profiling/*_metaphlan_bugs_list.tsv" >> run1.sh
  chmod +x run1.sh

  echo '#!/bin/bash' > run2.sh
  echo "$PROG_LOAD" >> run2.sh
  echo "humann_join_tables --input ${RESULTS_DIR}/8.Functional_Profiling/ \\" >> run2.sh
  echo "--output ${RESULTS_DIR}/8.Functional_Profiling/Merged_Sample_Tables/humann_genefamilies.tsv \\" >> run2.sh
  echo "--file_name genefamilies.tsv" >> run2.sh
  echo "rm ${RESULTS_DIR}/8.Functional_Profiling/*genefamilies.tsv" >> run2.sh
  chmod +x run2.sh

  echo '#!/bin/bash' > run3.sh
  echo "$PROG_LOAD" >> run3.sh
  echo "humann_join_tables --input ${RESULTS_DIR}/8.Functional_Profiling/ \\" >> run3.sh
  echo "--output ${RESULTS_DIR}/8.Functional_Profiling/Merged_Sample_Tables/humann_pathabundance.tsv \\" >> run3.sh
  echo "--file_name pathabundance.tsv" >> run3.sh
  echo "rm ${RESULTS_DIR}/8.Functional_Profiling/*pathabundance.tsv" >> run3.sh
  chmod +x run3.sh

  echo '#!/bin/bash' > run4.sh
  echo "$PROG_LOAD" >> run4.sh
  echo "humann_join_tables --input ${RESULTS_DIR}/8.Functional_Profiling/ \\" >> run4.sh
  echo "--output ${RESULTS_DIR}/8.Functional_Profiling/Merged_Sample_Tables/humann_pathcoverage.tsv \\" >> run4.sh
  echo "--file_name pathcoverage.tsv" >> run4.sh
  echo "rm ${RESULTS_DIR}/8.Functional_Profiling/*pathcoverage.tsv" >> run4.sh
  chmod +x run4.sh

  #Create script for running and submitting individual scripts
  echo '#!/bin/bash' > run.sh
  echo "#SBATCH --partition=$PARTITION_8" >> run.sh
  echo "#SBATCH --job-name=Merge_tables" >> run.sh
  echo "#SBATCH --error=${RESULTS_DIR}/8.Functional_Profiling/0.ErrorOut/Merge_tables_%A_%a.err" >> run.sh
  echo "#SBATCH --output=${RESULTS_DIR}/8.Functional_Profiling/0.Output/Merge_tables_%A_%a.out" >> run.sh
  echo "#SBATCH --time=$TIME_8" >> run.sh
  echo "#SBATCH --ntasks=1" >> run.sh
  echo "#SBATCH --cpus-per-task=$CPU_8" >> run.sh
  echo "#SBATCH --mem-per-cpu=$MEM_8" >> run.sh
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
fi
#################################################

###### REMOVE INTERMEDIATE SEQUENCE FILES #######

if [[ ! -z $DELETE ]]; then
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
fi
#################################################

echo "*** Metagenomics pipeline complete ***"
