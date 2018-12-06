#!/bin/bash


# Functions
welcome() {
    clear
    echo "#################################################################################################################
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

        -f | --force_overwrite         Danger! This option will allow the pipeline to overwrite files but
                                       will ensures a clean output. Default: FALSE [ TRUE or FALSE]

        -h | --help                    Prints this message.

##################################################################################################################
"
}

# Print weclome message
welcome

# Set defaults
threads_max=1
force_overwrite="FALSE"

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
        -f | --force_overwrite )
            force_overwrite="$2"
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

# Warn about overwrite
if [ "$force_overwrite" == "TRUE" ]; then
    printf "[WARNING] Set to overwrite files! Ctl + c now if you didn't mean to overwrite..."
    for i in {5..1}; do
        printf '%s...' "$i"
        sleep 1
    done
    printf "\n"
fi

## Debugging
#printf "Params:\n 
#--genome_index=${genome_index}\n
#--samples_file=${samples_file}\n"

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

# Read the columns of samples_file as separate arrays
while read col1 col2 col3 col4; do 
    samples+=($col1)
    read_lens+=($col2)
    fastq_paths_1+=($col3)
    fastq_paths_2+=($col4)
done < <(sed 's/;/\t/g' $samples_file)

# Check that we have a read len for every sample
if [ ! "${#samples[@]}" -eq "${#read_lens[@]}" ]; then
    echo "[ERROR] Check your samples file. number of sample != number of read length.\n"
    exit 1
fi

# Debugging
printf 'Found %s samples: \n' "${#samples[@]}"
printf '%s\n' "${samples[@]}"
echo ""
printf 'Read lengths are:\n'
printf '%s\n' "${read_lens[@]}"

# Get status of status of mate-pairs
for i in "${!samples[@]}"; do
    if [ -f "${fastq_paths_2[${i}]}" ]; then
        is_paireds+=("TRUE")
    else
        is_paireds+=("FALSE")
    fi 
done

# Debugging
printf 'Paired statuses are:\n'
printf '%s\n' "${is_paireds[@]}"

# Validate fastq files
for i in "${!samples[@]}"; do
    if [ ! -f "$fastq_paths_1[$i]" ]; then
        printf '[ERROR] Could not find fastq files for %s \n' "${samples[$i]}"
        exit 1
    fi
done


##################################################################################################################
# ALIGNMENT
##################################################################################################################

# NOTES:
# · Add dynamic read_len support
# · Check that STAR is in $PATH

# Collect and validate additional required params
# star index:
star_index_path="${genome_index_dir}/${genome_index}_${read_len}"
if [ ! -d "${star_index_path}" ]; then
    printf "[ERROR] Could not find genome index. Looked here: ${star_index_path} \n\n"
    exit 1
fi


# Make aligned directory in current working dir
if [ "${force_overwrite}" == "TRUE" ]; then
    rm -rf ./aligned; mkdir ./aligned
else
    mkdir -p ./aligned
fi


# Begin Alignment
printf 'Beginning STAR alignemnt on %s threads...\n' "${threads_max}"
for i in "${!samples[@]}"; do
    printf 'Processing %s' "${samples[$i]}"
    STAR \
        --runThreadN "${threads_max}" \
        --genomeDir "${star_index_path}" \
        --readFilesIn "${fastq_paths_1[$i]}" "${fastq_paths_2[$i]}" \
        --readFilesCommand zcat \
        --outFileNamePrefix "./aligned/${samples[$i]}_${genome_index}_" \
        --outSAMtype BAM Unsorted
done

## Debugging
#echo "--genome_index = ${genome_index}"
#echo "--samples_file = ${samples_file}"
#echo "--genome_index_dir = ${genome_index_dir}"








