#!/bin/bash

# Wrapper around Luca's LCS pipeline.
# Compute calder compartments, dump interchromosomal contacts, compute cscores and insulation boundaries/scores.
# AUTHOR: Karol Piera and Luca Nanni
# Email: karol.piera@unil.ch

SCRIPT_DIR=$(dirname $0)

# Parse arguments to variables
# Add most important parameters to getopts.
# Optionally one can trun fixed params such as genome into OPTARG as well.

usage () {
    	  echo "|---------------------------------------------------------------------------------------------------|"  
	  echo "| Run Loss of chromatin segregation analysis.                                                       |"
	  echo "| This script is a wrapper around Luca's LCS pipeline. It is meant to be run directly on            |"
	  echo "| hic pipeline output directory. See test comparison.tsv for comparison.                            |"
	  echo "| The main analysis is found in src/LCS_main.py.                                                    |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| It requires output of 4 prior modules: Calder compartment calls, DUMP interchromosomal contacts,  |"
	  echo "| compartment segregation scores from Calder Saddle and insulation (scores and boundires), 	    |"
	  echo "| from insulation module.                                                                           |"
	  echo "| Module 1 and module 4 are  implemented in fanc. Other modules are implemented in python.          |"
	  echo "| Scripts can be found in src/ directory.                                                           |"
	  echo "| IMPORTANT: To save time, if output exists, user is asked if they want to re-run modules 1-4!      |"
	  echo "|---------------------------------------------------------------------------------------------------|"    
	  echo "| USAGE: ./Run_LCS_pipeline.sh -i DATA -o PLOTS_PATH -c COMPARISON_TABLE       # In/Out, Comparison |" 
	  echo "|                      -g GENOME -b BIN_SIZE -p THREADS                        # CALDER flags       |"
	  echo "|                      -D DUMP_RES                                             # Dump flag          |"
	  echo "|                      -r HIC_RES -b BIN_SIZE                                  # Calder Saddle flag |"
	  echo "|                      -w INS_W_SIZE -n INS_NORMALISATION                      # Insulation flag    |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| See verbose descritpion below:                                                                    |"
	  echo "| -i Path to the directory with hic pipeline results.                                               |"
	  echo "|    IMPORTANT: this path serves as output path for Calder, Dump, Saddle and Insulation scores.     |"
	  echo "| -c TSV table with desired comparison. Sample names HAVE TO match hic pipeline output results!!!!  |"
          echo "|---------------------------------------------------------------------------------------------------|"
          echo "| 1) Run Calder Module:                                                                             |"
	  echo "|    -b BIN_SIZE - Resolution for Calder. Default 50000. Minimum is 10000.  			    |"
	  echo "|    (Used also for insulation analysis and calder saddle script.)                                  |"
	  echo "|    -g Genome [hg19 or hg38]. Default hg19. See ./ref for files.                                   |"
	  echo "|    -p Thread control. Default 1.								    |"
          echo "|    - Adaptive is set to true.                                                                     |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| 2) Run Dump Module: 	                                                    			    |"
	  echo "|   -D Resolution at which to dump interchromosomal contacts. Defualt is 1000000000.                |"
	  echo "|   -g GENOME                                                                                       |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| 3) Run Calder_Saddle. This module requires output from Calder.				    |"
	  echo "|    Compute compartment segregation scores. Outputs cscores tsv file.                              |"
	  echo "|    -b BIN_SIZE - Calder resoultion                                                                |"
	  echo "|    -r Hic resoultion to use. Default 1 000 000.                                                   |" 
	  echo "|    - Number of tiles for compartment segregation is set to 16.                                    |"
	  echo "|      Modify N_TILE parameter to change.                                                           |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| 4) Run Insulation Module:                                                                         |"
	  echo "|    -w Window size for insulation scores. Default 1000000.                                         |"
	  echo "|    -n Normalisation for insulation scores [VC_SQRT or KR].                                        |"
	  echo "|    -b BIN_SIZE - Resolution of Calder. Default 50000. Minimum is 10000. 			    |"
	  echo "|    Fixed params for insulation analysis (Modify in the script):                                   |"
	  echo "|    - Output format is bigwig.  DELTA set to 4,  TRIM_MEAN = 0                                     |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| One can skip above steps if needed (comment out core functions) and run LCS analysis directly.    |"
	  echo "| Sample dir structure required for LCS analysis:                                                   |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "|            ./{hic_sample}/1_CalderOutput_50000 						    |"
	  echo "|                           2_Dumps                                                                 |"
	  echo "|                           3_CalderSaddleOutput_r1000000 				     	    |"
	  echo "|                           4_InsulationOutput 						    	    |"	
	  echo "|                           aligned (hic_pipeline output) 					    |"
	  echo "|                                        						            |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| 5) Run LCS analysis: -c Comparisons table in tsv format.                                          |"
	  echo "|                      -o PLOTS PATH                                                                |"
	  echo "|                      -i Input data as above                                                       |"
	  echo "|---------------------------------------------------------------------------------------------------|"
}

declare -a DUMP_RES 
while getopts ":b:c:d:D:g:i:o:p:r:w:n:h" opt; do
    case ${opt} in
	b ) CALDER_RES=$OPTARG
	    if [ $CALDER_RES -lt 10000 ]; then
		echo "Calder resolution must be at least 10000"
		exit 1
	    fi
	    ;;
	c )
	    COMPARISON_TABLE=$OPTARG
	    echo "Extract comparisons from ${COMPARISON_TABLE}"
	    # check file passed and is a tsv
	    if [ $(echo $COMPARISON_TABLE | awk -F. '{print $NF}') == "tsv" ]; then
		if [ ! -f $COMPARISON_TABLE ]; then
		    echo "Comparison table does not exist"
		    exit 1
		fi
	    else
		echo "Invalid file format. Please provide a tsv file"
		exit 1
	    fi
		    
	    ;;
	d )
	    # Dump resolution for inter-chromosomal contacts
	    DUMP_RES+=($OPTARG)
	    ;;
	D )
	    # Dump resolution for intra-chromosomal contacts
	    DUMP_RES+=($OPTARG)
	    ;;
	g )
	    # Genome
	    if [ $OPTARG == "hg19" ]; then
		GEN_PATH=$SCRIPT_DIR/ref/hg19.main.chrom.sizes
		GENOME="hg19"
	    elif [ $OPTARG == "hg38" ]; then
		GEN_PATH=$SCRIPT_DIR/ref/hg38.main.chrom.sizes
		GENOME="hg38"
	    else
		echo "Invalid genome. Supported genomes are hg19 and hg38"
		exit 1
	    fi
	    ;;
	i )
	    # INPUT. Results directory
	    INPUT_DIR=$OPTARG
	    if [ ! -d $INPUT_DIR ]; then
		echo "Input directory does not exist"
		exit 1
	    fi
	    ;;
   
	p ) NPROC=$OPTARG;;
	r ) HIC_RES=$OPTARG;;
	w ) INS_W_SIZE=$OPTARG;;
	n)  INS_NORMALISATION=$OPTARG
	    if [ $INS_NORMALISATION != "KR" ] && [ $INS_NORMALISATION != "VC_SQRT" ]; then
		echo "Invalid normalisation. Choose either KR or VC_SQRT"
		exit 1
	    fi  
	    ;;
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
	    exit 1
	    ;;
	? )
 	   echo "Invalid option passed: $OPTARG" 1>&2
	   exit 1
	   ;;
    esac
done

#if [ "$INPUT_DIR == "" " ] || [ "$OUT_PATH" == "" ] || [ "$COMPARISON_TABLE" == "" ]; then
#    echo "Missing required arguments"
#    usage
#-   exit 1
#fi

# Fixed parameters are those that do not change across analysis very often
# If neeeded, they can be added to getpos as well to customize the input better

# Fixexd parameters - DUMP
CHROMS=$(cut -f1 $SCRIPT_DIR/ref/hg19.main.chrom.sizes | grep -v chrY | sort | uniq | tr '\n' , ) # Chromosomes to include 
EXCLUDED_CHROMOSOMES="chrY,chrM,chrMT,Y,M,MT"                                                     # Chromosomes to exclude in DUMP and SADDLE

# Fixed parameters - SADDLE
N_TILES=16
# Fixed parameters - INSUALTION
INS_OUTPUT_FORMAT="bigwig"
INS_TRIM_MEAN=0
BOUNDARY_DELTA_WINDOW=4

# Env
#source ~/mambaforge/etc/profile.d/conda.sh
#source ~/mambaforge/etc/profile.d/mamba.sh
mamba activate lcs_pipeline

# Define functions
# Ask user if they want to proceed - helper function
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

# Since Calder, Dump, Saddle and Insulation are essential for the LCS analysis but costly on time,
# we let user to decide if they want to re-run them or skip if the output is already present.
if [ $(find $INPUT_DIR -name *30.hic | wc -l) -eq 0 ]; then
	echo "No hic files found in $INPUT_DIR"
	exit 1
fi


if [ $(find $INPUT_DIR -name *sub*.bed | wc -l) -gt 0 ]; then
    ask_user "Calder"
    if [ $? -eq 0 ]; then
	Calder="Calder"
    else
	Calder=""
    fi
else
    Calder="Calder"
fi

if [ $(find $INPUT_DIR -type f -name binsize_${DUMP_RES[0]}.tsv | wc -l) -gt 0 ] &&
   [ $(find $INPUT_DIR -type f -name binsize_${DUMP_RES[1]}.tsv | wc -l) -gt 0 ]; then
    ask_user "Dump"
    if [ $? -eq 0 ]; then
	Dump="Dump"
    else
	Dump=""
    fi
else
    Dump="Dump"
fi

if [ $(find $INPUT_DIR -type d -name 3_CalderSaddleOutput* | wc -l) -gt 0 ]; then
    ask_user "Saddle and Insulation"
    if [ $? -eq 0 ]; then
	Saddle="Saddle"
    else
	Saddle=""
    fi
else
    Saddle="Saddle"
fi

if [ $(find $INPUT_DIR -type d -name 4_InsulationOutput | wc -l) -gt 0 ]; then
	ask_user "Insulation"
	if [ $? -eq 0 ]; then
	    Insulation="Insulation"
	else
	    Insulation="Stop"
	fi
else
    Insulation="Insulation"
fi


# Define each module. LCS scritpt is executed only if all the modules were run or data strcutre with corresponding output is present.
# Define each module. LCS_main.py scritpt is executed in either of the two cases:
# 1) all the modules were run in consectutive manner
# 2) module output was found in the input directory and user decided to skip the module.

run_module_calder() {
    calder --input $1 \
	   --bin_size $2 \
	   --nproc $3 \
	   -c $CHROMS \
	   --genome $GENOME \
	   --outpath $4 \
	   --adaptive
}

run_module_dump() {
    python src/DUMP_COOL.py --resolution ${1} \
	   --genome ${GENOME} \
	   "${2}::/resolutions/10000000" \
	   ${3}
}

run_module_saddle() {
    python src/CALDER_SADDLE.py ${1} ${2} ${3} \
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
    for f in `find ${2} -name "$(basename ${3})_*.bigwig"`
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

 
# Loop over hic pipeline results (input directory)
for dir in $INPUT_DIR/*; 
do
    if [ ! -d $dir ]; then
	continue
    else
	# Get files
	sample=$(basename $dir)
	HIC_PATH=$( realpath $(find $dir -name "*30.hic") ) # Get hic file
	COOL_PATH=$( realpath $(find $dir -name "*30.mcool") ) # Get cool file
	
	
	CALDER_OUT=${dir}/1_CalderOutput_${CALDER_RES}
	
	if [ ! -d $CALDER_OUT ]; then
	    mkdir -p $CALDER_OUT
	fi

	if [ $Calder == "Calder" ]; then
	    echo "---------------------------------------------"
	    echo "Run calder on $sample"
	    echo "---------------------------------------------"
	    run_module_calder ${HIC_PATH} ${CALDER_RES} ${NPROC} ${CALDER_OUT}
	else
	    echo "Skipping Calder for $sample sample"
	fi
		    	
	if [ "$Dump" == "Dump" ]; then
	    echo "---------------------------------------------"
	    echo "Dump interchromosomal contacts at resolutions ${DUMP_RES[0]} and ${DUMP_RES[1]}"
	    echo "---------------------------------------------"
	    DUMP_OUT="${dir}/2_Dumps/"
	    if [ ! -d $DUMP_OUT ]; then
		mkdir -p $DUMP_OUT
	    fi
	    for res in ${DUMP_RES[@]}; do
		echo "Dumping interchromosomal contacts at resolution $res"
		run_module_dump $res ${COOL_PATH} "${DUMP_OUT}/binsize_${res}.tsv"
	    done
	else
	    echo "Skipping Dump for $sample sample"
	fi
	
	if [ -f ${CALDER_OUT}/sub_compartments/all_sub_compartments.bed ]; then 
	    # Double-check if compartment file exsists (neccessary for Calder Saddle)
	    COMP_PATH=${CALDER_OUT}/sub_compartments/all_sub_compartments.bed
	    SADDLE_OUT=${dir}/3_CalderSaddleOutput_r${HIC_RES}
	    
	    if [ "$Saddle" == "Saddle" ]; then
		echo "---------------------------------------------"
		echo "Run calder saddle on $sample using resolution ${HIC_RES} hic resolution"
		echo "---------------------------------------------"
		if [ ! -d $SADDLE_OUT ]; then
		    mkdir -p $SADDLE_OUT
		fi
 		run_module_saddle ${COOL_PATH} ${COMP_PATH} ${SADDLE_OUT} 
	    else
		echo "Skipping Saddle for $sample"
	    fi
	else
	    echo "No compartment file found for $sample sample. Skipping Saddle"
	fi

	if [ $Insulation == "Insulation" ]; then
	    INS_PATH=${dir}/4_InsulationOutput
	
	    if [ ! -d $INS_PATH ]; then
		mkdir -p $INS_PATH
	    fi
	    OUT_F=${INS_PATH}/insulation_scores_r${CALDER_RES}_${INS_NORMALISATION}_w${INS_W_SIZE}
	
	    echo "---------------------------------------------"
	    echo "Run insulation module on $sample"
	    echo "---------------------------------------------"
	    run_module_insulation ${HIC_PATH} ${INS_PATH} ${OUT_F}
	else
	    echo "Skipping Insulation for $sample"
	fi    
    fi
done
echo "Done"
echo "---------------------------------------------------------------------------------------------------"
echo "---> Running LCS_main.py analysis for following comparisons:  					 |"
echo "---------------------------------------------------------------------------------------------------"

cat $COMPARISON_TABLE

if [ ! -d $OUT_PATH ]; then
	mkdir -p $OUT_PATH
fi

# Run LCS_main.py main analysis
python ${SCRIPT_DIR}/src/LCS_main.py \
        --INPUT_DIR ${INPUT_DIR} \
        --comparison ${COMPARISON_TABLE} \
        --intra_chr_res "${DUMP_RES[0]}" \
	--INTER_CHR_RES "${DUMP_RES[1]}" \
	-g $GENOME \
	-r $HIC_RES \
	-w $INS_W_SIZE \
        -o $OUT_PATH

# Merge plots into one report for convenience 
for res_dir in $OUT_PATH/*; do
    if [[  "$res_dir" == "Interchromosomal_OE" ]]; then
	# Skip
	echo "Skipping $res_dir"
    else
	echo "Merging reports for $res_dir"
	gs -o ${res_dir}/$(basename ${res_dir})_MergedReport.pdf \
	   -dNOPAUSE -sDEVICE=pdfwrite -dAutoRotatePages=/None -sPAPERSIZE=a4 \
	   -dPDFSETTINGS=/prepress \
	   -c "<</Orientation 3>> setpagedevice" \
	   -f ${res_dir}/*pdf
    fi
done
mamba deactivate
