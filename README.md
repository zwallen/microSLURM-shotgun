This GitHub repository houses a wrapper program (`microSLURM_metagenome.sh`) for processing shotgun metagenomic sequences (derived from Illumina paired-end whole genome shotgun sequencing) to taxonomic (relative abundances and abundance counts) and functional (gene and pathway) profiles on a computing cluster using a SLURM scheduling system. The overall pipeline includes performing an initial quality assessment on raw sequences using FastQC [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/] and MultiQC [https://multiqc.info/], sequence processing using BBMap suite of tools [https://sourceforge.net/projects/bbmap/] including merging of paired-end reads using BBMerge (optional), adapter removal and quality trimming/filtering using BBDuk, removal of host contaminant reads using BBSplit/BBMap, removal of low complexity sequences using BBDuk, a final quality assessment of processed sequences, and lastly taxonomic and functional profiling using MetaPhlAn [https://huttenhower.sph.harvard.edu/metaphlan/] and HUMAnN [https://huttenhower.sph.harvard.edu/humann/] workflows. Alternatively, one can use Kraken2 [https://github.com/DerrickWood/kraken2/wiki] and Bracken [https://ccb.jhu.edu/software/bracken/] instead of MetaPhlAn for taxonomic profiling.

The following gives an overview of the overall structure of the repository:

## Directory tree for repository
```
microSLURM_metagenome
|
|-- Environment -- Directory that contains a .yml file to create a conda environment 
|                  with all the necessary packages for running the pipeline scripts.
|                  Optional to build this environment, provided for convenience to 
|                  ensure proper programs and versions are downloaded.
|
|-- Reference_Files -- Directory that contains the shell script 0.get_reference_files.sh.
|                      Running 0.get_reference_files.sh will download all the necessary
|                      and most up to date reference files and databases used in the pipeline.
|
|-- microSLURM_metagenome -- Directory that contains separate shell scripts 
|                            for each step of the pipeline. Instead of running 
|                            the whole pipeline using the main wrapper script
|                            script, one can run it in chunks using these scripts. 
|                            Useful for when the dataset is large, and will not 
|                            finish running with using the main wrapper script, 
|                            or if the main wrapper script failed or had to be 
|                            stopped at a certain step. Each pipeline chunk comes 
|                            with its own example job script for submitting it to 
|                            a SLURM scheduler.
|
|-- microSLURM_metagenome.job -- An example sbatch job script for submitting the 
|                                microSLURM_metagenome.sh to a SLURM 
|                                scheduling system on a high performance computing cluster.
|
|-- microSLURM_metagenome.sh -- The wrapper program that runs the metagenomic pipeline.

```
## Important notes about the pipeline program

### Submission of internal sbatch jobs
The `microSLURM_metagenome` scripts will internally submit jobs for each step of the pipeline. For each job step, a `run.sh` file is created in the current directory that contains the code for the currently running step, and is deleted once the step completes. Partitions, time limits, number of cores (cpus), and memory per cpu can be requested through the wrapper scripts to meet the demands of a particular SLURM scheduler, dataset, etc.

### Required programs/databases and parameter descriptions
For descriptions of required programs/databases and parameters for `microSLURM_metagenome.sh` script, run the respective script with parameter `-h`. As stated in the directory tree above, the directory `Environment` contains a `.yml` file that can be used to build a conda environment with all the necessary programs for running the `microSLURM_metagenome.sh` script.

### Separate shell scripts for individual pipeline steps
The directory `microSLURM_metagenome/` contains separate shell scripts that can perform each step of the pipeline individually, and are numbered in the sequence they are performed in the wrapper script `microSLURM_metagenome.sh`. These have been provided as an alternative to running the whole pipeline in one shot. Using these can be useful if the dataset being processed is very large, and will not complete in time using the one shot `microSLURM_metagenome.sh` wrapper script. These are also useful if the main wrapper script has failed, or stopped, at a certain step, and you want to continue with the pipeline without having to re-run the `microSLURM_metagenome.sh` script.

### Removing intermediate sequence files generated during pipeline
During running of the `microSLURM_metagenome.sh` script, intermediate sequence files are generated after each step of sequence processing (merging paired-reads if specified, quality trimming/filtering, decontamination) with the final quality controlled sequences being generated after the low-complexity sequence filtering step. This may take up a lot of storage space depending on the dataset size, so if you would like intermediate files removed at the end of the pipeline, keeping only the processed, quality controlled sequences, use the flag `-d` to remove them and their corresponding directories. Sequences that have been extracted from metagenomic reads during sequence processing (host sequences, low-complexity sequences) will be stored for record purposes, along with a combined log file for each sample joining all logs from each step. As the numbering of directories no longer makes sense after removing the intermediate sequence directories, the pipeline output directory is reorganized to include a directory for initial FastQC reports, the processed sequences, final FastQC reports, and results for taxonomic and functional profiling all with new numbering. The `microSLURM_metagenome/` directory contains a stand-a-lone script for performing this deletion and reorganization if wanting to do at a later time, or if running the pipeline in chunks.

### Generating a pipeline report for tracking number of sequences (and other metrics) through pipeline
The `microSLURM_metagenome/` directory that houses separate shell scripts for each step of the pipeline also houses a script that can be used to generate a pipeline report detailing how many sequences were outputted from each step of the pipeline which can be useful for assessing sequence loss through the pipeline (`Pipeline_Report.sh`). Currently, this script only works if the pipeline was run without the `-d` flag (no removal of intermediate sequence files).

```
./microSLURM_metagenome.sh -h

##############################################################
# microSLURM_metagenome                                      #
# Whole Genome Shotgun Metagenomic Processing Pipeline       #
# Last updated: 19 Oct 2022                                  #
##############################################################
 
 Description: This is a wrapper program that wraps various  
 programs to process raw paired-end whole genome shotgun    
 metagenomic sequences. The end product of this pipeline    
 is taxonomic and functional (gene and pathway) abundances  
 that are ready for further statistical analyses.           
                                                            
 Required programs and databases:                           
    SLURM:      Program is designed to work with a SLURM    
                high performance computing cluster          
                scheduling system.                          
    R base:     For performing functions in pipeline script.
    FastQC:     For performing initial quality reports.     
    MultiQC:    For summarizing multiple FastQC reports.    
    BBMerge:    For merging paired-end reads.               
    BBDuk:      For adapter and quality trimming of raw wgs 
                reads, removal of PhiX sequences, and       
                removing low complexity sequences.          
    BBMap/                                                  
    BBSplit:    For removing host contamination from wgs    
                reads. Requires FASTA genome reference file 
                to map reads against.                       
    MetaPhlAn:  For generating taxonomic abundances.        
    HUMAnN:     For generating gene family, and pathway     
                abundances.                                 
    ChocoPhlAn: Database used for taxonomic profiling.      
                Downloaded for MetaPhlAn using              
                    $ metaphlan --install                   
                and for HUMAnN using                        
                    $ humann_databases \                   
                    $ --download chocophlan full \         
                    $ $INSTALL_LOCATION                    
    UniRef:     Database used for functional profiling. Can 
                be any of the UniRef databases downloaded   
                using humann_databases utility program.     
 Alternatively                                              
    Kraken2:    For taxonomic profiling.                    
    Bracken:    For calculating taxonomic abundances.       
    Kraken2 database:                                       
                Database used for taxonomic profiling with  
                Kraken2. Should be built using the Kraken2  
                command 'kraken2-build', then the kmer      
                distribution calculated using 'braken-build'
                Alternatively, can download the pre-built   
                standard database from the repository:      
                https://benlangmead.github.io/aws-indexes/k2
    kreport2mpa.py and combine_mpa.py:                      
                Scripts from KrakenTools to convert Bracken 
                output to MetaPhlAn style output and combine
                individual sample MetaPhlAn style outputs.  
                Make sure they are in your $PATH.          
                                                            
 Usage:                                                     
 microSLURM_metagenome.sh -i input_seqs_dir \  
                    -o output_dir \                         
                    -p 'commands; to; load; programs' \     
                    -r path/to/host/ref/file.fa \           
                    -c path/to/chocophlan/dir \             
                    -u path/to/uniref/dir \                 
                    -f notificationEmail@forFailures.edu \  
                    -n list,of,node,partitions \            
                    -t list,of,time,requests \              
                    -k list,of,cpu,number,requests \        
                    -m list,of,memory,per,cpu,requests \    
 if wanting to use Kraken2/Bracken instead of MetaPhlAn:    
                    -b /path/to/kraken2/database \          
                    [additional options]                    
                                                            
 Parameters:                                                
     -h    Print the parameter list below then exit.        
                                                            
 Required analysis parameters                               
     -i    (Required) Directory that contains the raw       
           fastq files to be processed. Sequences must have 
           file extensions .fastq OR .fq,                   
           and can be gzipped or not.                       
     -o    (Required) Directory to put output of pipeline   
           into. NOTE: make sure output directory is in an  
           area that has plenty of data storage space       
           available if processing large datasets.          
     -p    (Required) Single quoted string that contains    
           commands to load all the necessary programs      
           needed to run pipeline steps (e.g. activating    
           conda environments, loading modules, adding to   
           PATH, etc.).                                     
     -r    (Required) Path to reference genome file of host.
           Should be in FASTA format, and uncompressed.     
     -c    (Required) Path to humann ChocoPhlAn database.   
     -u    (Required) Path to humann UniRef90 database.     
                                                            
 Required SLURM parameters                                  
 Note: For each parameter, a comma separated list of 8      
 entries must be given (one for each pipeline step). If not 
 performing a certain step (e.g. merging of paired reads)   
 then you can put a 'NA' in the list.                       
     -n    (Required) Names of the partitions to request for 
           submitting jobs.                                 
     -t    (Required) Time requests for running jobs.       
           Specify in hours (e.g. 12:00:00 for 12 hours).   
     -k    (Required) Number of cores wanted for jobs.      
     -m    (Required) Amount of memory requested for each   
           requested core. Specify in Megabytes.            
     -f    (Required) E-mail to send notifications to upon  
           failure of any jobs.                             
                                                            
 Optional pipeline parameters                               
     -b    (Optional) Use Kraken2/Bracken to perform        
           taxonomic profiling and supply path to kraken2   
           database directory here.                         
     -x    (Optional) Merge paired-end reads before         
           performing the pipeline using BBMerge.           
     -a    (Optional) Path to adapters.fa file that comes   
           packaged with BBMerge and BBDuk. Required when   
           merging reads.                                   
     -j    (Optional) Join forward and reverse reads        
           together into one file (via concatenation) prior 
           to performing low-complexity sequence filtering. 
           May want to do this if using MetaPhlAn for       
           taxonomic profiling as this will be done anyway  
           prior to running MetaPhlAn and low-complexity    
           sequences caused by technical artificat are more 
           prone to be located in reverse reads while the   
           forward reads are still good.                    
     -s    (Optional) Skip certain steps in the pipeline if 
           need be. Provide a comma separated list of steps 
           that you wish to skip in the pipeline. List may  
           have the values: fastqc_initial, QC, decontam,   
           entropy_filter, fastqc_final,                    
           taxonomic_profiling, functional_profiling        
     -d    (Optional) Delete sequence files generated during
           intermediate pipeline steps (merging, quality    
           trimming/filtering, human decontamination) while 
           keeping potentially important components of each 
           step to keep record of (log files, extracted host
           sequences, extracted low-complexity sequences).  
           Will reorganize output folder if flag is given   
           since the default numbering of folder will no    
           longer make sense.
```
