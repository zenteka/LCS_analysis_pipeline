#!/bin/bash

# Wrapper around Luca's LCS pipeline
# AUTHOR: Karol Piera and Luca Nanni
# Work inprogress
# Email: karol.piera@unil.ch

#set -e
SCRIPT_DIR=$(dirname $0)

# Parse arguments to variables
# Add most important parameters to getopts.
# Optionally one can trun fixed params such as genome into OPTARG as well.
# TODO:
## 1) Add more fixed parameters to getopts -d to control dump resolution, -g to control genome


usage () {
    	  echo "|---------------------------------------------------------------------------------------------------|"  
	  echo "| Run LCS PIPELINE.                                                                                 |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| Usage: ./cmd -i DATA -c COMPARISONS -o PLOTS_PATH -b BIN_SIZE -p THREADS -n NORM -w WINDOW_SIZE   |"
	  echo "| See verbose descritpion below:                                                                    |"
	  echo "| -i Path to the directory with hic pipeline results.                                               |"
	  echo "|    IMPORTANT: this path serves as output path for Calder, Dump, Saddle and Insulation scores.     |"
	  echo "| -c TSV table with desired comparison. Sample names HAVE TO match hic pipeline output results!!!!  |"
          echo "|---------------------------------------------------------------------------------------------------|"
          echo "| 1) Run Calder: -b BIN SIZE (also used to compute insulation scores)                               |"
	  echo "|                -p No. of threads to use. Genome is set to hg19 (modify in the script).            |"
          echo "|                -Adaptive is set to true.                                                          |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| 2) Dump interchromosomal contacts: Check DUMP_RES varaible to modify.                             |"
	  echo "|    -d Dump resolution (minimal resolution). Default 10000                                         |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| 3) Run Saddle:                                                                                    |"
	  echo "|    Run Saddle script for each resolution that we was used to dump matrices.                       |"
	  echo "|    N. Tiles is set to 16. Modify in the script.                                                   |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| 4) Run Insulation scores: -w Insulation window size.                                              |"
	  echo "|                           -n Normalisation for insulation scores. [KR or VC_SQRT]                 |"
	  echo "|                           -b resolution (reused BIN SIZE for Calder)                              |"
	  echo "|    Fixed params for insulation analysis (Modify in the script):                                   |"
	  echo "|    - Output format is bigwig.  DELTA set to 4,  TRIM_MEAN = 0                                     |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| One can skip above steps if needed (comment out core functions) and run LCS analysis directly.    |"
	  echo "| Sample dir structure required for LCS analysis:                                                   |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "|            ./{hic_sample}/1_CalderOutput_50000 						    |"
	  echo "|                           2_Dumps                                                                 |"
	  echo "|                           3_CalderSaddleOutput_r100000 					    |"
	  echo "|                           3_CalderSaddleOutput_r1000000 				     	    |"
	  echo "|                           3_CalderSaddleOutput_r10000000 				   	    |"
	  echo "|                           4_InsulationOutput 						    	    |"	
	  echo "|                           aligned 						             	    |"
	  echo "|                                        						            |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| 5) Run LCS analysis: -c Comparisons                                                               |"
	  echo "|                      -o PLOTS PATH                                                                |"
	  echo "|                      -i Input data as above                                                       |"
	  echo "|---------------------------------------------------------------------------------------------------|"
}

while getopts ":b:i:c:p:w:n:o:h" opt; do
    case ${opt} in
	c )
	    # LCS comparisons
	    LCS_COMPARISON=$OPTARG
	    ;;
	b )
	    # Calder resolution 
	    CALDER_RES=$OPTARG
	    if [ $CALDER_RES -lt 10000 ]; then
		echo "Calder resolution must be at least 10000"
		exit 1
	    fi
	    ;;
	g)
	    # Genome
	    if [ $OPTARG == "hg19" ]; then
		GEN_PATH=$SCRIPT_DIR/ref/hg19.main.chrom.sizes
		GENOME="hg19"
	    elif [ $OPTARG == "hg38" ]; then
		GEN_PATH=$SCRIPT_DIR/ref/hg38.main.chrom.sizes
		GENOME="hg38"
	    else
		echo "Invalid genome. Provide either hg19 or hg38"
		exit 1
	    fi
	    ;;
	i )
	    # Results directory
	    INPUT_DIR=$OPTARG
	    ;;
   
	p)
	  # Thread control
	    NPROC=$OPTARG	  
	    ;;
	w)
	    # Insulation window size
	    INSULATION_WINDOW_SIZE=$OPTARG
	    ;;
	n)
	    # Normalisation for insulation scores
	    INS_NORMALISATION=$OPTARG
            
	    # check if normalisation is valid
	    if [ $INS_NORMALISATION != "KR" ] && [ $INS_NORMALISATION != "VC_SQRT" ]; then
		echo "Invalid normalisation. Choose either KR or VC_SQRT"
		exit 1
	    fi  
	    ;;
	o)
	    # Output directory
	    OUT_PATH=$OPTARG
	    ;;
       * | h )
	   usage
	    exit 1
	    ;;
       \? )
 	    echo "Invalid option: $OPTARG" 1>&2
	    exit 1
	    ;;
  esac
done


# Fixed parameters are those that do not change across analysis very often
# If neeeded, they can be added to getpos as well to customize the input better

# Fixed parameters CALDER
GENOME="hg38"

# Fixexd parameters - DUMP
DUMP_RES=("100000" "1000000" "10000000" ) # Dump resolutions
EXCLUDED_CHROMOSOMES="chrM,chrY" # Chromosomes to exclude

# Fixed parameters - SADDLE
N_TILES=16

# Fixed parameters - INSUALTION
INS_OUTPUT_FORMAT="bigwig"
INS_TRIM_MEAN=0
BOUNDARY_DELTA_WINDOW=4

# Env
source ~/mambaforge/etc/profile.d/conda.sh
source ~/mambaforge/etc/profile.d/mamba.sh
mamba activate lcs_pipeline

# Define functions
# Ask user if they want to proceed
ask_user () {
        echo "$1 output found."
	echo "Do you want run $1 module? Type 'yes' or 'y' to proceed..."
	read in_p
	echo $in_p
	if [ $in_p == 'yes' ] || [ $in_p == 'y' ]; then
	    return 0
	else
	    return 1
	fi
}

run_module_calder() {
    calder --input $1 \
	   --bin_size $2 \
	   --nproc $3 \
	   --genome $GENOME \
	   --outpath $4 \
	   --adaptive
}

run_module_dump() {
	python 2.0_DUMP_COOL.py --resolution ${1} \
	   --excludeChroms ${EXCLUDED_CHROMOSOMES} \
	   "${2}::/resolutions/${DUMP_RES[0]}" \
	   ${3}
}

run_saddle() {
	python 2.1_CALDER_SADDLE.py ${1} ${2} ${3} \
	   --hic_resolution ${4} \
	   --compartments_resolution ${5} \
	   --n_tiles ${6} \
	   --excludeChroms ${EXCLUDED_CHROMOSOMES}
}

run_insulation() {
    # Compute insulation scores
    fanc insulation $1@${CALDER_RES}@${INS_NORMALISATION} \
	 $3 \
	 -o "bigwig" \
	 --window-sizes ${INSULATION_WINDOW_SIZE} \
	 --offset 0 \
	 --geom-mean \
	 --trim-mean ${INS_TRIM_MEAN}
    
    for f in `find ${2} -name "$(basename ${3})_*.bigwig"`
    do
	mv ${f} ${3}.bigWig
    done
    
    bigWigToBedGraph ${3}.bigWig ${3}.bedGraph.tmp
    
    cat ${3}.bedGraph.tmp | \
	awk 'BEGIN {F = "\t"; OFS = "\t"} {if($4 != "nan"){gsub(/^chr/, "", $1); print "chr"$1,$2,$3,$4;}}' > ${3}.bedGraph
    rm ${3}.bedGraph.tmp
	
    # Compute insualtion boundaries
    fanc boundaries ${3}.bigWig \
	 ${3}_boundaries.bed \
	 --window-sizes ${INSULATION_WINDOW_SIZE} \
	 --delta ${BOUNDARY_DELTA_WINDOW}
}

# Make some check before running modules
# Since Calder, Dump, Saddle and Insulation are essential for the LCS analysis but costly on time,
# we let user to decide if they want to re-run them or skip if the output is already present.

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

if [ $(find $INPUT_DIR -type d -name 2_Dumps | wc -l) -gt 0 ]; then
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

# Loop over hic pipeline results:
for dir in $INPUT_DIR/*; 
do
    if [ ! -d $dir ]; then
	continue
    else
	# Get files
	sample=$(basename $dir)
	HIC_PATH=$( realpath $(find $dir -name "*30.hic") ) # Get hic file
	COOL_PATH=$( realpath $(find $dir -name "*30.mcool") ) # Get cool file

	# Store 
	SAM_ARRAY+=($sample)
	HIC_ARRAY+=($HIC_PATH)
	COOL_ARRAY+=($COOL_PATH)
	
	# Run Calder
	# Activate conda environment
	       
	echo "---------------------------------------------"
	echo "Run calder on $sample"
	echo "---------------------------------------------"
	
	CALDER_OUT=${dir}/1_CalderOutput_${CALDER_RES}
	
	if [ ! -d $CALDER_OUT ]; then
	    mkdir -p $CALDER_OUT
	fi

	if [ $Calder == "Calder" ]; then
	    run_module_calder ${HIC_PATH} ${CALDER_RES} ${NPROC} ${CALDER_OUT}
	else
	    echo "Skipping Calder for $sample"
	fi
	
	# Dump interchromosomal contacts for given resolutions
	for res in ${DUMP_RES[@]};
	do
	    echo "---------------------------------------------"
	    echo "Dump interchromosomal contacts at resolution ${res}"
	    echo "---------------------------------------------"
	    
            DUMP_OUT="${dir}/2_Dumps/"
	    if [ ! -d $DUMP_OUT ]; then
		mkdir -p $DUMP_OUT
	    fi
	    if [ $Dump == "Dump" ]; then
		run_module_dump $res ${COOL_PATH} "${DUMP_OUT}/binsize_${res}.tsv"
	    else
		echo "Skipping Dump for $sample"
	    fi
            
	    # Run CALDER saddle for each his resolution
	    echo "---------------------------------------------"
	    echo "Run calder saddle on $sample at resolution ${res}"
	    echo "---------------------------------------------"
	    
	    # Check to bed on the save side
	    if [ -f ${CALDER_OUT}/sub_compartments/all_sub_compartments.bed ]; then

		# Get compartemt calls
		COMP_PATH=${CALDER_OUT}/sub_compartments/all_sub_compartments.bed

		SADDLE_OUT=${dir}/3_CalderSaddleOutput_r${res}
		
		if [ ! -d $SADDLE_OUT ]; then
		    mkdir -p $SADDLE_OUT
		fi
		if [ $Saddle == "Saddle" ]; then
 		    run_saddle ${COOL_PATH} ${COMP_PATH} ${SADDLE_OUT} ${res} ${CALDER_RES} ${N_TILES}
		else
		    echo "Skipping Saddle for $sample"
		fi
		
	    fi
	done
	
	echo "Completed CALDER, Dump and Saddle for $sample"
	echo "---------------------------------------------"
	echo "Compute insulation scores and insualtion boundries"
	echo "---------------------------------------------"
	
	INS_PATH=${dir}/4_InsulationOutput
	
	if [ ! -d $INS_PATH ]; then
	    mkdir -p $INS_PATH
	fi
	
	OUT_F=${INS_PATH}/insulation_scores_r${CALDER_RES}_${INS_NORMALISATION}_w${INSULATION_WINDOW_SIZE}
	if [ $Saddle == "Saddle" ]; then
	    run_insulation ${HIC_PATH} ${INS_PATH} ${OUT_F}
	else
	    echo "Skipping Insulation for $sample"
	fi
    fi    
done

echo "Done"
echo "---------------------------------------------------------------------------------------------------"
echo "---> Run LCS analysis using python 3.0_LCS.py with the following arguments:                        |"
echo "---> -c $LCS_COMPARISON -i $INPUT_DIR -o $OUT_PATH                                                 |"
echo "---------------------------------------------------------------------------------------------------"

if [ ! -d $OUT_PATH ]; then
	mkdir -p $OUT_PATH
fi

# Run main LCS analysis
python 3.0_RUN_LCS.py -c $LCS_COMPARISON \
        -g $GENOME \
        -d $INPUT_DIR \
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
