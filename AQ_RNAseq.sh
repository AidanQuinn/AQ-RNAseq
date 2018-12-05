#!/bin/bash


# Functions
welcome() {
    clear
    echo "

#################################################################################################################
#                                                                                                               #
#                                    Welome to AQ_RNAseq pipeline!                                              #
#                                      · v.0.1.0 · 2018.12.05 ·                                                 #
#                                                                                                               #
#################################################################################################################
"
}

usage(){
    echo "
Usage: AQ_RNAseq.sh [options]
                                                                                                               
     Required Arguments:                                                                                       
        -g | --genome_index            Can be one of [ GRCh38 or mm10 ]. The full path to the STAR index       
                                       must be specified in AQ_RNAseq_settings.txt as genome_index_dir.        
                                                                                                               
        -s | --samples_file            Sample file should be formated as tsv txt file.
        
        -p | --parameters_file         File containing AQ_RNAseq run basic paramerts. Tab separated key 
                                       value pairs. For example see default_params.txt included.
               
     Optional Arguments:                                                                                       
        -m | --align_mode              Can be one of [ sinlge or double ]. Default, single. Passed to 
                                       STAR for alignment mode. Double pass alignment is more sensitive
                                       for discovery of novel splice junctions, but requires more time.

        -p | --is_paired               [ yes or no ] Default: no. Indicates that reads are paired-end. If
                                       is yes, fastq files should end in {sample}_1.fastq.gz and 
                                       {sample}_2.fastq.gz for mate-pair 1 and 2, respectively.
    
        -t | --threads_max             Maximum number of threads to use for all steps.

        -h | --help                    Prints this message.

##################################################################################################################
"
}

# Print weclome message
welcome

# Collect options or print usage statement
POSITIONAL=()
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -g | --genome_index )
            genome_index="$2"
            shift # past argument
            shift # past value
            ;;
        -s | --samples_file )
            samples_file="$2"
            shift
            shift
            ;;
        -p | --parameters_file )
            parameters_file="$2"
            shift
            shift
            ;;
        -m | --align_mode )
            align_mode="$2"
            shift
            shift
            ;;
        -p | --is_paired )
            is_paired="$2"
            shift
            shift
            ;;
        -t | --threads_max )
            threads_max="$2"
            shift
            shift
            ;;
        -v | --verbose )
            verbose="TRUE"
            shift
            shift
            ;;
        -h | --help )
            usage
            exit 1
    esac
done

## Debugging
printf "Params:\n 
--genome_index=${genome_index}\n
--samples_file=${samples_file}\n"

##################################################################################################################
# RUN PARAMETERS
##################################################################################################################
# Read parameters from parameters_file
if [ "${verbose}" == "TRUE" ]; then 
    printf "Reading parameters from $parameters_file \n" 
fi

source $parameters_file


##################################################################################################################
# INPUT FILES
##################################################################################################################
# Check that samples_file exists
if [ ! -f "${samples_file}" ]; then
    printf "[ERROR] Samples File not found. Looked here: ${samples_file} \n \n"
fi

while read col1 col2 col3 ; do 
    samples+=($col1)
    is_paired+=($col2) 
    read_len+=($col3)
    #ar4+=($col4) # for future use, if want more info
done < <(sed 's/;/\t/g' $samples_file)

## Debugging
printf 'Found %s samples: \n' "${#samples[@]}"
printf '%s\n' "${samples[@]}"
echo ""

##################################################################################################################
# ALIGNMENT
##################################################################################################################

# Dynamic read_len support
if [ "${#read_len[@]}" -gt 1 ]; then
    echo "Dynamic Read Length not yet supported. Feature comming soon."
    exit 1
fi

# Collect and validate required params
star_index_path="${genome_index_dir}/${genome_index}_${read_len}"
if [ ! -d "${star_index_path}" ]; then
    printf "[ERROR] Could not find genome index. Looked here: ${star_index_path} \n\n"
    exit 1
fi







## Debugging
#echo "--genome_index = ${genome_index}"
#echo "--samples_file = ${samples_file}"
#echo "--genome_index_dir = ${genome_index_dir}"








