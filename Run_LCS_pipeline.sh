#!/bin/bash
# Wrapper around Luca's LCS pipeline that is part of the WGD project.
# AUTHOR: Karol Piera and Luca Nanni
# Email: karol.piera@unil.ch

SCRIPT_DIR=$(dirname $0)

# LCS is meant to be run on hic pipeline output. Desired comparison should be provided in a tsv, csv or xlsx format.
# LCS pipeline is a module-based analysis. Ouptut of module 1 to 4 is required to run the main LCS analysis.

usage () {
    	  echo "|---------------------------------------------------------------------------------------------------|"  
	  echo "| Run Loss of chromatin segregation analysis.                                                       |"
	  echo "| This script is a wrapper around WGD LCS pipeline. It is meant to be run directly on               |"
	  echo "| hic pipeline output. See test comparison.tsv for a sample file.                                   |"
	  echo "| The main analysis is found in src/LCS_main.py.                                                    |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| The pipeline consists of 5 steps. Output of first 4 steps is REQUIRED to run the main analysis.   |"
	  echo "| The steps are:                                                                                    |"
	  echo "| 1) Calder - Compute compartments and subcompartments. (Can be skipped if output exists)           |"
	  echo "| 2) Dump - Dump interchromosomal contacts. (Can be skipped if output exists)                       |"
	  echo "| 3) Calder Saddle - Compute compartment segregation scores (Can be skipped if output exists)       |"
	  echo "| 4) Insulation - Compute insulation scores. (Can be skipped if output exists)                      |"
	  echo "| 5) LCS analysis - Reuires output of all the above steps.                                          |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| Each module is described below.                                                                   |"
	  echo "| IMPORTANT: To save time, if output exists, user is asked if they want to re-run modules 1-4!      |"
	  echo "|---------------------------------------------------------------------------------------------------|"    
	  echo "| USAGE: ./Run_LCS_pipeline.sh -i DATA -o PLOTS_PATH -c COMPARISON_TABLE       # In/Out, Comparison |" 
	  echo "|                      -g GENOME -b BIN_SIZE -p THREADS                        # CALDER flags       |"
	  echo "|                      -D DUMP_RES                                             # Dump flag          |"
	  echo "|                      -r HIC_RES -b BIN_SIZE                                  # Calder Saddle flags|"
	  echo "|                      -w INS_W_SIZE -n INS_NORMALISATION                      # Insulation flags   |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| See verbose descritpion below:                                                                    |"
	  echo "| -i Path to the directory with hic pipeline results.                                               |"
	  echo "|    IMPORTANT: this path serves as output path for Calder, Dump, Saddle and Insulation scores.     |"
	  echo "| -c TSV table with desired comparison. Sample names HAVE TO match hic pipeline output results!!!!  |"
          echo "|---------------------------- CALDER MODULE --------------------------------------------------------|"
          echo "| 1) Run Calder Module:                                                                             |"
	  echo "|    -b BIN_SIZE - Resolution for Calder. Default 50000. Minimum is 10000.  			    |"
	  echo "|    Chromosomes will be set according to the selected genome and file provided in the              |"
	  echo "|    ${SCRIPT_DIR}/ref directory.                                                                   |"
	  echo "|    (Used also for insulation analysis and calder saddle script.)                                  |"
	  echo "|    -g Genome [hg19 or hg38]. Default hg19. See ./ref for files.                                   |"
	  echo "|    -p Thread control. Default 1.								    |"
          echo "|    - Adaptive is set to true.                                                                     |"
	  echo "|--------------------------- DUMP MODULE -----------------------------------------------------------|"
	  echo "| 2) Run Dump Module: 	                                                    			    |"
	  echo "|   -D Resolution at which to dump interchromosomal contacts. Defualt is 1 000 000 000.             |"
	  echo "|   -g GENOME                                                                                       |"
	  echo "|---------------------------- SADDLE MODULE --------------------------------------------------------|"
	  echo "| 3) Run Calder_Saddle. This module requires output from Calder.				    |"
	  echo "|    Compute compartment segregation scores. Outputs cscores tsv file.                              |"
	  echo "|    -b BIN_SIZE - Calder resoultion                                                                |"
	  echo "|    -r Hic resoultion to use. Default 1 000 000.                                                   |" 
	  echo "|    - Number of tiles for compartment segregation is set to 16.                                    |"
	  echo "|      Modify N_TILE parameter to change.                                                           |"
	  echo "|----------------------------- INSULATION MODULE ---------------------------------------------------|"
	  echo "| 4) Run Insulation Module:                                                                         |"
	  echo "|    -w Window size for insulation scores. Default 1000000.                                         |"
	  echo "|    -n Normalisation for insulation scores [VC_SQRT or KR].                                        |"
	  echo "|    -b BIN_SIZE - Resolution of Calder. Default 50000. Minimum is 10000. 			    |"
	  echo "|    Fixed params for insulation analysis (Modify in the script):                                   |"
	  echo "|    - Output format is bigwig.  DELTA set to 4,  TRIM_MEAN = 0                                     |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| One can skip above steps if needed (comment out core functions) and run LCS analysis directly.    |"
	  echo "| Sample dir structure required for LCS analysis:                                                   |"
	  echo "|                                                                                                   |"
	  echo "|            ./{hic_sample}/1_CalderOutput_50000 						    |"
	  echo "|                           2_Dumps                                                                 |"
	  echo "|                           3_CalderSaddleOutput_r1000000 				     	    |"
	  echo "|                           4_InsulationOutput 						    	    |"	
	  echo "|                           aligned (hic_pipeline output) 					    |"
	  echo "|                                        						            |"
	  echo "|------------------------------ LCS MAIN ANALYSIS --------------------------------------------------|"
	  echo "| 5) Run LCS analysis: -i Input data as above.                                                      |"
	  echo "|                      -o Output directory for plots.                                               |"
	  echo "|                      -c Comparison table.                                                         |"
	  echo "|                      -g Genome [hg19 or hg38].                                                    |"
	  echo "|                      -r Hic resolution.                                                           |"
	  echo "|                      -w Window size for insulation scores.                                        |"
	  echo "|                      --intra_chr_res Intra chromosomal resolution.                                |"
	  echo "|                      --INTER_CHR_RES Inter chromosomal resolution.                                |"
	  echo "|---------------------------------------------------------------------------------------------------|"

}

while getopts ":b:c:d:D:g:i:o:p:r:w:n:h" opt; do
    case ${opt} in
	b ) CALDER_RES=$OPTARG
	    if [ $CALDER_RES -lt 10000 ]; then
		echo "Calder resolution must be at least 10000"
		exit 1
	    fi
	    ;;
	c ) COMPARISON_TABLE=$OPTARG;;                              # Comparison table
       	d ) INTRA_CHR_RES=$OPTARG;;                                 # Intra chromosomal resolution
	D ) INTER_CHR_RES=$OPTARG;; 				    # Inter chromosomal resolution
	d )
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

# Fixed parameters:
CHROMS=$(cut -f1 $SCRIPT_DIR/ref/hg19.main.chrom.sizes | grep -v chrY | sort | uniq | tr '\n' , ) # Chromosomes to include - used for Calder
EXCLUDED_CHROMOSOMES="chrY,chrM,chrMT,Y,M,MT"                                                     # Chromosomes to exclude in DUMP and SADDLE
## SADDLE
N_TILES=16
## INSULATION
INS_OUTPUT_FORMAT="bigwig"
INS_TRIM_MEAN=0
BOUNDARY_DELTA_WINDOW=4

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

if [ $(find $INPUT_DIR -type f -name binsize_${INTER_CHR_RES}.tsv | wc -l) -gt 0 ]; then
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

display_msg_run () {
	echo "--------------------------------------------------------------------------------"
	echo " Run module $1 on $2 with $3 bp resolution"
	echo "--------------------------------------------------------------------------------"
}

display_msg_skip () {
	echo "--------------------------------------------------------------------------------"
	echo "Skipping $1 for $2"
	echo "--------------------------------------------------------------------------------"
}
 
# Run main modules: Calder, Dump, Saddle, Insulation
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
	if [ "$Calder" == "Calder" ]; then
	    display_msg_run $Calder $sample $CALDER_RES
	    run_module_calder ${HIC_PATH} ${CALDER_RES} ${NPROC} ${CALDER_OUT}
	else
	    display_msg_skip "Calder" $sample 
	fi
		    	
	if [ "$Dump" == "Dump" ]; then
            
	    echo "---------------------------------------------"
	    echo "Dump interchromosomal contacts at resolutions ${DUMP_RES[0]} and ${DUMP_RES[1]}"
	    echo "---------------------------------------------"
            
	    DUMP_OUT="${dir}/2_Dumps/"
	    if [ ! -d $DUMP_OUT ]; then
		mkdir -p $DUMP_OUT
	    fi
	    display_msg_run $Dump $sample $INTER_CHR_RES
	    run_module_dump ${INTER_CHR_RES} ${COOL_PATH} "${DUMP_OUT}/binsize_${INTER_CHR_RES}.tsv"
	    for res in ${DUMP_RES[@]}; do
		echo "Dumping interchromosomal contacts at resolution $res"
		run_module_dump $res ${COOL_PATH} "${DUMP_OUT}/binsize_${res}.tsv"
	    done
	else
	    display_msg_skip "Dump" $sample
	fi
	
	if [ -f ${CALDER_OUT}/sub_compartments/all_sub_compartments.bed ]; then 
	    # Double-check if compartment file exsists (neccessary for Calder Saddle)
	    COMP_PATH=${CALDER_OUT}/sub_compartments/all_sub_compartments.bed
	    SADDLE_OUT=${dir}/3_CalderSaddleOutput_r${HIC_RES}
	    
	    if [ "$Saddle" == "Saddle" ]; then
		if [ ! -d $SADDLE_OUT ]; then
		    mkdir -p $SADDLE_OUT
		fi
		display_msg_run $Saddle $sample $HIC_RES
 		run_module_saddle ${COOL_PATH} ${COMP_PATH} ${SADDLE_OUT} 
	    else
		display_msg_skip "Saddle" $sample
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
	    display_msg_run $Insulation $sample $CALDER_RES
	    run_module_insulation ${HIC_PATH} ${INS_PATH} ${OUT_F}
	else
	    display_msg_skip "Insulation" $sample
	fi    
    fi
done

echo "---------------------------------------------------------------------------------------------------"
echo " ------------------ All modules run successfully --------------------------------------------------"
echo " ------------------ Run LCS_main.py main analysis -------------------------------------------------"
echo "---------------------------------------------------------------------------------------------------"

# Run LCS_main.py main analysis
python ${SCRIPT_DIR}/src/LCS_main.py \
        --INPUT_DIR ${INPUT_DIR} \
        --comparison ${COMPARISON_TABLE} \
        --intra_chr_res ${INTRA_CHR_RES} \
	--INTER_CHR_RES ${INTER_CHR_RES} \
	-g $(realpath $GEN_PATH) \
	-r $HIC_RES \
	-w $INS_W_SIZE \
        -o $OUT_PATH
if [ $? -ne 0 ]; then
	echo "Error in LCS_main.py"
	exit 1
fi
    
# If on UNIX levrage ghostscript to merge reports
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
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
fi
