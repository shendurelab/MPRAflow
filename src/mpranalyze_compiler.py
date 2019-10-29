#contributed by Tal Ashuach

import os
import re
import sys
from collections import defaultdict

annot_pattern = re.compile("^([DR]NA).*\(condition (.*), batch (.*)\)$")
def get_annot(head):
	m = annot_pattern.match(head)
	if m is not None:
		return m.group(1,2,3)
		
def parse(in_path, out_directory):
	with open(in_path, 'r') as infile:
		header = next(infile).strip().split('\t')
		## parse header
		header_annot = [get_annot(h) for h in header[3:]]
		dna_annot = [h for h in header_annot if h[0] == 'DNA']
		rna_annot = [h for h in header_annot if h[0] == 'RNA']
		n_dna_obs = len(dna_annot)
		n_rna_obs = len(rna_annot)

		rna_dict = defaultdict(list)
		dna_dict = defaultdict(list)
		## 
		for l in infile:
			l = l.strip().split('\t')
			seq_id = l[0]
			if seq_id == 'no_BC':
				continue
			
			## add DNA counts to dictionary
			dna_dict[seq_id].append(l[3:(3 + n_dna_obs)])

			## add RNA counts to dictionary
			rna_dict[seq_id].append(l[(3 + n_dna_obs):]) 

	n_bc = max([len(x) for x in rna_dict.values()])

	## write output DNA annotations
	dna_colnames = []
	with open(os.path.join(out_directory, 'dna_annot.tsv'), 'w') as ofile:
		ofile.write('\t'.join(["sample", "type", "condition", "batch", "barcode"]) + '\n')
		for i in range(1, n_bc + 1):
			for x in dna_annot:
				sample_name = '_'.join(list(x) + [str(i)])
				dna_colnames.append(sample_name)
				ofile.write('\t'.join([sample_name] + list(x) + [str(i)]) + '\n')

	## write output DNA counts
	with open(os.path.join(out_directory, 'dna_counts.tsv'), 'w') as ofile:
		ofile.write('\t'.join(['seq_id'] + dna_colnames) + '\n')
		for seq_id in dna_dict:
			ofile.write('\t'.join(
				[seq_id] + # sequence ID
				[count for bc in dna_dict[seq_id] for count in bc] + # flattened list of observations
				(['0'] * n_dna_obs * (n_bc - len(dna_dict[seq_id]))) # zero padding
				) + '\n')

	## write output RNA annotations
	rna_colnames = []
	with open(os.path.join(out_directory, 'rna_annot.tsv'), 'w') as ofile:
		ofile.write('\t'.join(["sample", "type", "condition", "batch", "barcode"]) + '\n')
		for i in range(1, n_bc + 1):
			for x in rna_annot:
				sample_name = '_'.join(list(x) + [str(i)])
				rna_colnames.append(sample_name)
				ofile.write('\t'.join([sample_name] + list(x) + [str(i)]) + '\n')

	## write output RNA counts
	with open(os.path.join(out_directory, 'rna_counts.tsv'), 'w') as ofile:
		ofile.write('\t'.join(['seq_id'] + rna_colnames) + '\n')
		for seq_id in rna_dict:
			ofile.write('\t'.join(
				[seq_id] + # sequence ID
				[count for bc in rna_dict[seq_id] for count in bc] + # flattened list of observations
				(['0'] * n_rna_obs * (n_bc - len(rna_dict[seq_id]))) # zero padding
				) + '\n')

if __name__ == '__main__':
	parse(sys.argv[1], sys.argv[2])


