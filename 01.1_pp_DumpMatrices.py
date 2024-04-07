dimport pandas as pd
import os

import subprocess
from cooler import Cooler

# 
starting_binsize = 100000
resolutions = [1000000000]
output_path = 'data/processed/hic_dumped_matrices'
hic_samples = pd.read_excel("hic_samples_karol.xlsx").set_index('hic_sample_name').dropna(subset=['cool_path'])

os.makedirs(output_path, exist_ok=True)

failed_samples = []

chromosomes_to_exclude = {
	"Lozano2021_Watermelon_Diploid": ["CG_Chr00"],
	"Lozano2021_Watermelon_Tetraploid": ["CG_Chr00"]
}

for sample in hic_samples.index:
	print("-"*20)
	print(sample)
	cool_path = hic_samples.loc[sample, 'cool_path'] 
	if not os.path.isfile(cool_path):
		print("* MCOOL file not available")
		continue
	cool_path += f"::/resolutions/{starting_binsize}"
	genome = hic_samples.loc[sample, 'genome']
	spath = os.path.join(output_path, sample)
	os.makedirs(spath, exist_ok=True)

	for res in resolutions:
		print(f"* Resolution: {res} bp")
		outpath = os.path.join(spath, f"dump_{res}.tsv")
		if os.path.isfile(outpath):
			print("** Already present")
			continue
                
		CMD = f"python scripts_luca/scripts/dump_cool.py --resolution {res}"

		if sample in chromosomes_to_exclude.keys():
			excChroms = ",".join(chromosomes_to_exclude[sample])
			CMD += f" --excludeChroms {excChroms}"
		else:
			CMD += f" --excludeChroms chrY,chrMT,chrM,Y,M,MT"

		CMD += f" {cool_path} {outpath}"
		try:
			subprocess.run(CMD, shell=True, check=True)
		except subprocess.CalledProcessError:
			failed_samples.append( (sample, res) )
	print("-"*20)

