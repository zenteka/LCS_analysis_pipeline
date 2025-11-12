#!/bin/bash
# LCS PIPELINE - Loss of Chromatin Segregation analysis from the WGD paper.
# AUTHOR: Karol Piera and Luca Nanni
# Email: karol.piera@unil.ch

# LCS is to be run on hic pipeline output. Desired comparison should be provided in a tsv, csv or xlsx format.
# It is a module-based analysis. Ouptut of module 1 to 4 is required to run the main LCS analysis.

# TODO:
## Flip names for plotting
## Consolidate included / excluded chromosomes across modules for maximum compatibility. Include handling of non-standard chromosomes.

SCRIPT_DIR=$(dirname $0)

usage () {
    	  echo "|---------------------------------------------------------------------------------------------------|"  
	  echo "| Run Loss of chromatin segregation analysis.                                                       |"
	  echo "| Author: Karol Piera and Luca Nanni                                                                |"
	  echo "| This script runs the full LCS pipeline including. It consists of 5 steps.                         |"
	  echo "| The pipeline consists of 5 steps.                						    |"
	  echo "| The steps are:                                                                                    |"
	  echo "| 1) Calder - Compute compartments and subcompartments. (Can be skipped if output exists)           |"
	  echo "| 2) Dump - Dump interchromosomal contacts. (Can be skipped if output exists)                       |"
	  echo "| 3) Calder Saddle - Compute compartment segregation scores (Can be skipped if output exists)       |"
	  echo "| 4) Insulation - Compute insulation scores. (Can be skipped if output exists)                      |"
	  echo "| 5) LCS analysis - Requires output of all the above steps.                                         |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| USAGE: ./Run_LCS_pipeline.sh -i DATA directory with hic pipeline output. It serves as output      |"
	  echo "|                      path for Calder, Dump, Saddle and Insulation scores.                         |"
	  echo "|                      -i INPUT_DIR                                                                 |"
	  echo "|                      -c COMPARISON_TABLE [tsv, csv, xlsx] desired comparisons.                    |"
	  echo "|                      -o OUTPUT directory for plots                                                |"
	  echo "|                      -g GENOME -b BIN_SIZE -p THREADS                        # CALDER flags       |"
	  echo "|---------------------------------------------------------------------------------------------------|"
          echo "|---------------------------- CALDER MODULE --------------------------------------------------------|"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| 1) Calder parameters: 									    |"
	  echo "|    -b HARDCODED: 50000.  BIN_SIZE - Resolution at which to run Calder.                            |"
	  echo "|    -g GENOME  [hg19 or hg38]									    |"
	  echo "|    -c HARDCODED: List of chromosomes to include. All autosomes + X.				    |"
	  echo "|    -p NPROC - Number of threads to use. Default 4.						    |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "|--------------------------- DUMP MODULE -----------------------------------------------------------|"
	  echo "|---------------------------------------------------------------------------------------------------|"	  
	  echo "| 2) Run Dump Module parameters:								    |"
	  echo "|    Dump interchromosomal contacts at HARDCODED 1 000 000 000 resolution.			    |"
	  echo "|    - GENOME [hg19 or hg38]									    |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "|---------------------------- SADDLE MODULE --------------------------------------------------------|"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| 3) Run compartment segregation computation. This is the largest module.			    |"
	  echo "|    Requires Calder output (subcompartments bed file) and the .mcool file.		            |"
	  echo "|    Compute compartment segregation scores. Outputs cscores tsv file.                              |"
	  echo "|    -b HARDCODED: 50000.  BIN_SIZE - Resolution of Calder.			                    |"
	  echo "|    -r Hic resoultion to use for superpixelling: HARDCODED 1 000 000.			            |"
	  echo "|    - Number of tiles for compartment segregation is set to 16.                                    |"
	  echo "|      Modify N_TILE parameter to change.                                                           |"
	  echo "|    - Excluded chromosomes are read from ref/chroms_to_exclude.txt file.                           |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "|----------------------------- INSULATION MODULE ---------------------------------------------------|"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| 4) Run Insulation Module:                                                                         |"
	  echo "|    -w Window size for insulation scores. Default 1000000.                                         |"
	  echo "|    -n Normalisation for insulation scores [VC_SQRT or KR].                                        |"
	  echo "|    -b BIN_SIZE - Resolution of Calder. Default 50000. Minimum is 10000. 			    |"
	  echo "|    Fixed params for insulation analysis (Hardcoded in the script):                                |"
	  echo "|    - Output format is bigwig.  DELTA set to 4,  TRIM_MEAN = 0                                     |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| Sample dir structure required for LCS analysis:                                                   |"
	  echo "|            ./{hic_sample}/1_subcomp_50000 						            |"
	  echo "|                           2_trans                                                                 |"
	  echo "|                           3_saddle 				     	                            |"
	  echo "|                           4_insulation 						    	    |"	
	  echo "|                           aligned (hic_pipeline output) 					    |"
	  echo "|                                        						            |"
	  echo "|------------------------------ LCS MAIN ANALYSIS --------------------------------------------------|"
	  echo "| 5) Run LCS analysis: -i Input data                                                                |"
	  echo "|                      -o Output directory for plots.                                               |"
	  echo "|                      -c Comparison table.                                                         |"
	  echo "|                      -g Genome [hg19 or hg38].                                                    |"
	  echo "|                      -r Hic resolution.                                                           |"
	  echo "|                      -w Window size of insulation scores.                                         |"
	  echo "|                      --intra_chr_res Intra chromosomal resolution. [-d]                           |"
	  echo "|---------------------------------------------------------------------------------------------------|"
}

while getopts ":c:g:i:o:p:h" opt; do
    case ${opt} in 
	c ) COMPARISON_TABLE=$OPTARG;;                              # Comparison table
	g )
	    # Genome
	    if [ $OPTARG == "hg19" ]; then
		GEN_PATH=$SCRIPT_DIR/ref/hg19.main.chrom.sizes
		GENOME="hg19"
	    elif [ $OPTARG == "hg38" ]; then
		GEN_PATH=$SCRIPT_DIR/ref/hg38.main.chrom.sizes
		GENOME="hg38"
	    elif [ $OPTARG == "mm10" ]; then
		GEN_PATH=$SCRIPT_DIR/ref/mm10.main.chrom.sizes
		GENOME="mm10"
	    else
		echo "Invalid genome. Supported genomes are hg19 and hg38"
		exit 1
	    fi
	    ;;
	i ) INPUT_DIR=$OPTARG
	    if [ ! -d $INPUT_DIR ]; then
		echo "Input directory does not exist"
		exit 1
	    fi
	    ;;
	p ) NPROC=$OPTARG;;
	o ) OUT_PATH=$OPTARG
	    if [ ! -d $OUT_PATH ]; then
		mkdir -p $OUT_PATH
	    fi
	    ;;
	: )
	   echo "Option -$OPTARG requires an argument." >&2
	   exit 1
	   ;;
	h ) usage
	    exit 1 ;;
	? )
 	   echo "Invalid option passed: $OPTARG" 1>&2
	   exit 1 ;;
    esac
done

# HARDCODED DEFAULTS
CHROMS=$(cut -f1 $SCRIPT_DIR/ref/hg19.main.chrom.sizes |\
	     grep -v chrY |\
	     sort | uniq |\
	     tr '\n' , | sed 's/chr//g' ) # Chromosomes to include - used for Calder, pass explicitly to avoid issues with non-standard chromosomes.

CHROMS=${CHROMS%,}
echo "Using chromosomes: $CHROMS"

EXCLUDED_CHROMOSOMES=$(cat $SCRIPT_DIR/ref/chroms_to_exclude.txt )                                                      
EXCLUDED_CHROMOSOMES=$(echo "$EXCLUDED_CHROMOSOMES" | sed 's/^chr//' )
EXCLUDED_CHROMOSOMES=${EXCLUDED_CHROMOSOMES%,}

CALDER_RES=50000
HIC_RES=1000000

# TRANS CONTACTS ANALYSIS
INTER_CHR_RES=1000000000

# Segregation analysis 
N_TILES=16

## INSULATION - same as original WGD paper
INS_OUTPUT_FORMAT="bigwig"
INS_NORMALISATION='KR'
INS_TRIM_MEAN=0
BOUNDARY_DELTA_WINDOW=4
INS_W_SIZE=1000000
INTRA_CHR_RES=1000000

ask_user () {
        echo "$1 output found."
	echo "Do you want run $1 module? Type 'yes' or 'y' to proceed..."
	read in_p
	echo $in_p
	if [ $in_p == 'yes' ] || [ $in_p == 'y' ] ; then
	    return 0
	else
	    return 1
	fi
}

# Screen samples for required input. Collect for further processing.
SAMPLES=()
RUN_CALDER=()
RUN_DUMP=()
RUN_SADDLE=()
RUN_INSULATION=()
  
while IFS= read -r line; do
    # check if sample exists in input dir
    if [ ! -f "$INPUT_DIR/$line/aligned/inter_30.hic" ]; then
	echo "Error: Sample $line does not appear to have hic file!"
	echo "Please check your input directory and comparison table!"
	exit 1
    else
	SAMPLES+=("$line")
    fi
done < <(cut -f1-2 $COMPARISON_TABLE | tr '\n' '\n' | tr '\t' '\n' | sort | uniq)

# Step 2: Check for missing required inputs
for sample in "${SAMPLES[@]}"; do
    echo "Checking sample $sample for missing inputs..."

    # Check if Calder output is missing
    if ! find -L "$INPUT_DIR" -type d -path "*/$sample/1_subcomp_${CALDER_RES}/sub_compartments" | grep -q .; then
        RUN_CALDER+=("$sample")
    fi

    # Check if binsize dump is missing
    if ! find -L "$INPUT_DIR" -type f -path "*/$sample/*binsize_${INTER_CHR_RES}.tsv" | grep -q .; then
        RUN_DUMP+=("$sample")
    fi

    # Check if Saddle output is missing
    if ! find -L "$INPUT_DIR" -type f -path "*/$sample/3_saddle/compartment_scores.tsv" | grep -q .; then
        RUN_SADDLE+=("$sample")
    fi

    # Check if Insulation output is missing
    if ! find -L "$INPUT_DIR" -type d -path "*/$sample/4_insulation" | grep -q .; then
        RUN_INSULATION+=("$sample")
    fi
done

run_module_calder() {
    calder --input $1 \
	   --bin_size $2 \
	   --nproc $3 \
	   -c $CHROMS \
	   --genome $GENOME \
	   --outpath $4 \
	   --adaptive
    if [ ! -f "$4/sub_compartments/all_sub_compartments.bed" ]; then
     	echo "Error: Calder subcompartments bed file not found!"
	cp "$4/sub_compartments/all_sub_compartments.bed" "$4/sub_compartments/all_sub_compartments.bed.bak"
	# Adjust 1-indexed fields: columns 2 and 7 (convert to 0-based for BED format); See comments section.
	awk '{print $1"\t"($2-1)"\t"$3"\t"$4"\t"$5"\t"$6"\t"($7-1)"\t"$8"\t"$9}' \
	    "$4/sub_compartments/all_sub_compartments.bed.bak" > "$4/sub_compartments/tmp.bed"
	sort -k1,1 -k2,2n "$4/sub_compartments/tmp.bed" > "$4/sub_compartments/all_sub_compartments.bed"
	rm "$4/sub_compartments/tmp.bed"
    fi
}

run_module_dump() {
    python src/DUMP_COOL.py --resolution ${1} \
	   --genome ${GENOME} \
	   "${2}::/resolutions/10000000" \
	   ${3}
}

run_module_saddle() {
	python src/CALDER_SADDLE_FAST.py ${1} ${2} ${3} \
	   --hic_resolution ${HIC_RES} \
	   --compartments_resolution ${CALDER_RES} \
	   --n_tiles ${N_TILES} \
	   --excludeChroms ${EXCLUDED_CHROMOSOMES}
}

run_module_insulation() {
    # Compute insulation scores
    fanc insulation $1@${CALDER_RES}@${INS_NORMALISATION} \
	 $3 \
	 -o "bigwig" \
	 --window-sizes ${INS_W_SIZE} \
	 --offset 0 \
	 --geom-mean \
	 --trim-mean ${INS_TRIM_MEAN}
    for f in `find -L ${2} -name "$(basename ${3})_*.bigwig"`
    do
	mv ${f} ${3}.bigWig
    done
    
    bigWigToBedGraph ${3}.bigWig ${3}.bedGraph.tmp
    
    cat ${3}.bedGraph.tmp | \
	awk 'BEGIN {F = "\t"; OFS = "\t"} {if($4 != "nan"){gsub(/^chr/, "", $1); print "chr"$1,$2,$3,$4;}}'\
	    > ${3}.bedGraph
    rm ${3}.bedGraph.tmp
	
    # Compute insualtion boundaries
    fanc boundaries ${3}.bigWig \
	 ${3}_boundaries.bed \
	 --window-sizes ${INS_W_SIZE} \
	 --delta ${BOUNDARY_DELTA_WINDOW}
}

display_msg_run () {
	echo "--------------------------------------------------------------------------------"
	echo " Run module $1 on $2 with $3 bp resolution"
	echo "--------------------------------------------------------------------------------"
}

display_msg_skip () {
	echo "--------------------------------------------------------------------------------"
	echo "Skipping module $1"
	echo "--------------------------------------------------------------------------------"
}
display_msg_present () {
	echo "--------------------------------------------------------------------------------"
	echo "$1 output already present for all samples. Skipping $1 module"
	echo "--------------------------------------------------------------------------------"
}

echo "---------------------------------------------------------------------------------------------------"
echo " ------------------ Running modules 1-4: Calder, Dump, Saddle, Insulation ---------------------------"
echo "---------------------------------------------------------------------------------------------------"

# check every RUN_ array and in case it is not empty ask user whether to run modules for given samples.
if [ ${#RUN_CALDER[@]} -gt 0 ]; then
    echo "Detected missing subcompartment calls for samples: ${RUN_CALDER[@]}"
    ask_user "Calder"
    if [ $? -eq 0 ]; then
	for sample in "${RUN_CALDER[@]}"; do
	    display_msg_run "Calder" $sample $CALDER_RES
	    HIC_PATH=$( realpath $(find -L $INPUT_DIR/$sample -name "*30.hic") ) # Get hic file
	    CALDER_OUT=${INPUT_DIR}/${sample}/1_subcomp_${CALDER_RES}
	    if [ ! -d $CALDER_OUT ]; then mkdir -p $CALDER_OUT; fi
	    run_module_calder ${HIC_PATH} ${CALDER_RES} ${NPROC} ${CALDER_OUT}
	done
    else
	display_msg_skip "Calder"
    fi
else
    display_msg_present "Calder"
fi

if [ ${#RUN_DUMP[@]} -gt 0 ]; then
	echo "Detected missing inter-chromosomal dumps for samples: ${RUN_DUMP[@]}"
	ask_user "Dump"
	if [ $? -eq 0 ]; then
	    for sample in "${RUN_DUMP[@]}"; do
		display_msg_run "Dump" $sample $INTER_CHR_RES
		COOL_PATH=$( realpath $(find -L $INPUT_DIR/$sample -name "*30.mcool") ) # Get cool file
		DUMP_OUT="${INPUT_DIR}/${sample}/2_trans/"
		if [ ! -d $DUMP_OUT ]; then mkdir -p $DUMP_OUT; fi
		run_module_dump ${INTER_CHR_RES} ${COOL_PATH} "${DUMP_OUT}/binsize_${INTER_CHR_RES}.tsv"
	    done
	else
	    display_msg_skip "Dump"
	fi
else
    display_msg_present "Dump"
fi

if [ ${#RUN_SADDLE[@]} -gt 0 ]; then
	echo "Detected missing Calder Saddle output for samples: ${RUN_SADDLE[@]}"
	ask_user "Saddle"
	if [ $? -eq 0 ]; then
	for sample in "${RUN_SADDLE[@]}"; do
	    display_msg_run "Saddle" $sample $HIC_RES
	    COOL_PATH=$( realpath $(find -L $INPUT_DIR/$sample -name "*30.mcool") ) # Get cool file
	    CALDER_OUT=${INPUT_DIR}/${sample}/1_subcomp_${CALDER_RES}
	    COMP_PATH=${CALDER_OUT}/sub_compartments/all_sub_compartments.bed
	    SADDLE_OUT=${INPUT_DIR}/${sample}/3_saddle
	    if [ ! -d $SADDLE_OUT ]; then mkdir -p $SADDLE_OUT; fi
	    run_module_saddle ${COOL_PATH} ${COMP_PATH} ${SADDLE_OUT}
	done
	else
	    display_msg_skip "Saddle"
	fi
else
    display_msg_present "Saddle"
fi

if [ ${#RUN_INSULATION[@]} -gt 0 ]; then
	echo "Detected missing insulation scores for samples: ${RUN_INSULATION[@]}"
	ask_user "Insulation"
	if [ $? -eq 0 ]; then
	for sample in "${RUN_INSULATION[@]}"; do
	    display_msg_run "Insulation" $sample $CALDER_RES
	    HIC_PATH=$( realpath $(find -L $INPUT_DIR/$sample -name "*30.hic") ) # Get hic file
	    INS_PATH=${INPUT_DIR}/${sample}/4_insulation
	    if [ ! -d $INS_PATH ]; then mkdir -p $INS_PATH; fi
	    OUT_F=${INS_PATH}/insulation_scores_r${CALDER_RES}_${INS_NORMALISATION}_w${INS_W_SIZE}
	    run_module_insulation ${HIC_PATH} ${INS_PATH} ${OUT_F}
	done
	else
	    display_msg_skip "Insulation"
	fi
else
    display_msg_present "Insulation"
fi

echo "---------------------------------------------------------------------------------------------------"
echo " ------------------ All required input for LCS analysis is present --------------------------------"
echo "---------------------------------------------------------------------------------------------------"
echo " ------------------ Run Genome-wide LCS analysis -------------------------------------------------"
echo "---------------------------------------------------------------------------------------------------"

# Main analysis starts here
python ${SCRIPT_DIR}/src/lcs_main_dev.py \
        --INPUT_DIR ${INPUT_DIR} \
        --comparison ${COMPARISON_TABLE} \
        --intra_chr_res ${INTRA_CHR_RES} \
	--INTER_CHR_RES ${INTER_CHR_RES} \
	-g $(realpath $GEN_PATH) \
	-r $HIC_RES \
	-w $INS_W_SIZE \
        -o $OUT_PATH

echo "---------------------------------------------------------------------------------------------------"
echo " ------------------ Run Chromosome-level compartment segregation analysis ---------------------------"
echo "---------------------------------------------------------------------------------------------------"

CHROM_LCS="${OUT_PATH}/Chrom_LCS"
if [ ! -d $CHROM_LCS ]; then mkdir -p $CHROM_LCS; fi
python ${SCRIPT_DIR}/src/LCS_subcomp_chrom.py \
		--data ${INPUT_DIR} \
		--comparisons ${COMPARISON_TABLE} \
		--output ${CHROM_LCS}

echo "---------------------------------------------------------------------------------------------------"
echo " ------------------ LCS PIPELINE COMPLETED -------------------------------------------------------"
echo "---------------------------------------------------------------------------------------------------"
