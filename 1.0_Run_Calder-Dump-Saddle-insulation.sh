#!/bin/bash

# Wrapper around Luca's LCS pipeline
# AUTHOR: Karol Piera and Luca Nanni
# DATE: 2024-03-24
# WOrk-inporgress

set -e

# Parse arguments to variables

while getopts ":a:i:c:p:w:n:o:h" opt; do
    case ${opt} in
	a )
	# LCS comparison
	LCS_COMPARISON=$OPTARG
	;;
	#b )
	## Control comparison
	#CONTROL_COMPARISON=$OPTARG
	#;;
	c )
	  # Calder resolution 
	  CALDER_RES=$OPTARG
	  if [ $CALDER_RES -lt 10000 ]; then
	      echo "Calder resolution must be at least 10000"
	      exit 1
	  fi
	  ;;
	i )
	  # Results directory
	  INPUT_DIR=$OPTARG
	  ;;
   
	p)
	  # Number of processees 
	  NPROC=$OPTARG
	  
	  ;;
	w)
	  # Insulation window size parameter
	  INSULATION_WINDOW_SIZE=$OPTARG
	  ;;
	n)
	  # Choose normalisation for insulation scores
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
	  echo "|---------------------------------------------------------------------------------------------------|"  
	  echo "| Scritpts runs full LCS analysis on the output of the hic pipeline.                                |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  echo "| Usage: ./cmd -i INPUTT DIRECTORY; output of the hic pipeline.                                     |"
	  echo "| Calder flags: -c Calder resolution; -p Number of processes. Adaptive is set to True.              |"
	  echo "| Insulation scores flags: -w Insulation window size; -n Insulation normalisation. 		    |"
	  echo "| -c Calder resolution is also passed as hic resolution to compute insulation scores.		    |"
	  echo "|---------------------------------------------------------------------------------------------------|"
	  exit 0
	  ;;
      \? )
	  echo "Usage: cmd -i INPUT_DIR"
	  exit 1
	  ;;
  esac
done


# Fixed parameters are those that do not change across analysis very often
# If neeeded, they can be added to getpos as well to customize the input better

# Fixed parameters CALDER
GENOME="hg19"

# Fixexd parameters - DUMP
DUMP_RES=("100000" "1000000" "10000000") # Dump resolutions
EXCLUDED_CHROMOSOMES="chrM,chrY" # Chromosomes to exclude

# Fixed parameters - SADDLE
N_TILES=16

# Fixed parameters - INSUALTION
INS_OUTPUT_FORMAT="bigwig"
INS_TRIM_MEAN=0
BOUNDARY_DELTA_WINDOW=4

# Initalize arrays 
SAM_ARRAY=()
HIC_ARRAY=()
COOL_ARRAY=()

source ~/mambaforge/etc/profile.d/conda.sh
source ~/mambaforge/etc/profile.d/mamba.sh
mamba activate calder

# Define functions
run_calder() {
    calder --input $1 \
	   --bin_size $2 \
	   --nproc $3 \
	   --genome $GENOME \
	   --outpath $4 \
	   --adaptive
}

run_dump() {
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


# Loop over results
for dir in $INPUT_DIR/*; 
do
    if [ ! -d $dir ]; then
	continue
    else
	# Collect variables
	sample=$(basename $dir)
	HIC_PATH=$( realpath $(find $dir -name "*30.hic") ) # Get hic file
	COOL_PATH=$( realpath $(find $dir -name "*30.mcool") ) # Get cool file

	# Store
	SAM_ARRAY+=($sample)
	HIC_ARRAY+=($HIC_PATH)z``
	COOL_ARRAY+=($COOL_PATH)

	# Get file path
	
	# Run Calder
	# Activate conda environment
	
       
	echo "---------------------------------------------"
	echo "Running calder on $sample"
	echo "---------------------------------------------"
	
	CALDER_OUTPATH=${dir}/1_CalderOutput_${CALDER_RES}
	
	if [ ! -d $CALDER_OUTPATH ]; then
	    mkdir -p $CALDER_OUTPATH
	fi
		
	run_calder ${HIC_PATH} ${CALDER_RES} ${NPROC} ${CALDER_OUTPATH}
	
	# Dump interchromosomal contacts at given resolutions
	mamba activate hic_pipeline	

	for res in ${DUMP_RES[@]};
	do
	    echo "---------------------------------------------"
	    echo "Dumping interchromosomal contacts at resolution ${res}"
	    echo "---------------------------------------------"
	    
	    echo "Dumping interchromosomal contacts at resolution ${res}"
	    
	    if [ ! -d ${dir}/dumps ]; then
		mkdir -p ${dir}/dumps
	    fi
	    
	    outPath="${dir}/2_Dumps/binsize_${res}.tsv"
	    mkdir -p $(dirname ${outPath})
	    
	    # Run CALDER saddle for each his resolution
	    echo "---------------------------------------------"
	    echo "Running calder saddle on $sample at resolution ${res}"
	    echo "---------------------------------------------"
	    
	    # Check to bed on the save side
	    if [ -f ${CALDER_OUTPATH}/sub_compartments/all_sub_compartments.bed ]; then

		# Get compartemt calls
		COMP_PATH=${CALDER_OUTPATH}/sub_compartments/all_sub_compartments.bed

		SADDLE_OUT=${dir}/3_CalderSaddleOutput_r${res}
		
		if [ ! -d $SADDLE_OUT ]; then
		    mkdir -p $SADDLE_OUT
		fi	
 		run_saddle ${COOL_PATH} ${COMP_PATH} ${SADDLE_OUT} ${res} ${CALDER_RES} ${N_TILES} 
	    fi
	done
	
	echo "Completed CALDER, Dump and Saddle for $sample"
	echo "---------------------------------------------"
	echo "Run insulation scores"
	
	INS_PATH=${dir}/4_InsulationOutput
	
	if [ ! -d $INS_PATH ]; then
	    mkdir -p $INS_PATH
	fi
	
	OUT_F=${INS_PATH}/insulation_scores_r${CALDER_RES}_${INS_NORMALISATION}_w${INSULATION_WINDOW_SIZE}
	
	fanc insulation ${HIC_PATH}@${CALDER_RES}@${INS_NORMALISATION} \
	      ${OUT_F} \
	      -o "bigwig" \
	      --window-sizes ${INSULATION_WINDOW_SIZE} \
	      --offset 0 \
	      --geom-mean \
	      --trim-mean ${INS_TRIM_MEAN}
	for f in `find ${INS_PATH} -name "$(basename ${OUT_F})_*.bigwig"`
	do
	    mv ${f} ${OUT_F}.bigWig
	done
	
	bigWigToBedGraph ${OUT_F}.bigWig ${OUT_F}.bedGraph.tmp
	cat ${OUT_F}.bedGraph.tmp | \
	    awk 'BEGIN {F = "\t"; OFS = "\t"} {if($4 != "nan"){gsub(/^chr/, "", $1); print "chr"$1,$2,$3,$4;}}' > ${OUT_F}.bedGraph
	rm ${OUT_F}.bedGraph.tmp
	
	# Compute insualtion boundaries
	fanc boundaries ${OUT_F}.bigWig \
		 ${OUT_F}_boundaries.bed \
		 --window-sizes ${INSULATION_WINDOW_SIZE} \
		 --delta ${BOUNDARY_DELTA_WINDOW}  

    fi
done

# Generate tsv for LCS analysis
for i in "${!SAM_ARRAY[@]}"; do
    echo "Index: $i; Sample: ${SAM_ARRAY[$i]}; Path: ${PATH_ARRAY[$i]}"
    echo -e "${PATH_ARRAY[$i]}\t${SAM_ARRAY[$i]}\tNONE" >> analyzed_samples_table.tsv
done

echo "Done"
echo "---------------------------------------------------------------------------------------------------"
echo "---> Run LCS analysis using python 3.0_LCS.py with the following arguments:                        |"
echo "---> -i $INPUT_DIR -t samples_table.tsv -o LCS_output                                              |"
echo "---------------------------------------------------------------------------------------------------"

# Make results section in the current directory
if [ ! -d $OUT_PATH ]; then
	mkdir -p $OUT_PATH
fi

# Run main LCS analysis
python 3.0_RUN_LCS.py -a $LCS_COMPARISON \
       -d $INPUT_DIR \
       -o $OUT_PATH
## Merge reports
#for res in $OUT_PATH/*; do
#    gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=${res}_Merged_Report.pdf ${res}/*.pdf
#    mv to repeorts
#done
mamba deactivate
