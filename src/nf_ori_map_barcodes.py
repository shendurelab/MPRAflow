

import itertools
import os
import pandas as pd
import pickle
import pysam
import sys

from collections import Counter, defaultdict
from tqdm import tqdm

##CMD python ori_map_barcodes.py data_assoc/  

#verify trailing arguments
project_dir=sys.argv[1]
fastq_in=sys.argv[2]
n_fastq_f=sys.argv[3]
bamfile=sys.argv[4]
n_bam_f=sys.argv[5]
prefix=sys.argv[6]
mapq=sys.argv[7]
baseq=sys.argv[8]
cigar=sys.argv[9]

print(project_dir)
print(fastq_in)
#print(n_fastq)
print(bamfile)
#print(n_bam)

if cigar=='n':
    cigar=''


print(n_fastq_f)
print(n_bam_f)
with open(n_fastq_f, 'r') as file:
    n_fastq = file.read().replace('\n','')

#f=open(n_fastq_f, 'r')
#n_fastq=f.read().replace('\n','')

with open(n_bam_f, 'r') as file:
    n_bam = file.read().replace('\n','')

print('counts')
print(n_bam)
print(n_fastq)
#map coords to BCs and filter based on given read/map quality and an exact cigar string if provided
def get_coords_to_barcodes(fastq_in, n_fastq,bamfile,n_bam,mapq=30,baseq=30,cigar=''):
    coords_to_barcodes_fn = f'{prefix}_coords_to_barcodes.pickle'

    if not os.path.exists(coords_to_barcodes_fn):
        #fastq = pysam.FastxFile(f'{project_dir}/R-PL-ass_S2_R2_001.fastq.gz')
        #fastq = pysam.FastxFile(f'{project_dir}/Gracie-assoc_S4_R2_001.fastq')
        fastq = pysam.FastxFile(fastq_in)
        #n_fastq_records = 1205322 # (wc -l) / 4
        print(n_fastq)
        n_fastq_records = float(n_fastq)/4 # (wc -l) / 4
        
        #bam = pysam.AlignmentFile(f'{project_dir}/library.bam', 'rb')
        #bam = pysam.AlignmentFile(f'{project_dir}/mem_Gracie-assoc_S4_merged.bam', 'rb')
        #n_bam_records = 712525

        bam = pysam.AlignmentFile(bamfile, 'rb')
        n_bam_records = float(n_bam)
        #print((bam))   
 
        query_to_coords = {}
        bad_pairs = 0
        poor_quality = 0
        print('start')
        #get names of the aligned reads that havae a high enough quality and match sequence
        for i, read in tqdm(enumerate(bam), 'paired-end reads', total = n_bam_records):
            #print(read.cigarstring)
            #print(read)
            #print(read.reference_name)
            #print(isinstance(read.reference_name, str))
            
            #if not read.is_read1:
            #    continue
            
            #if not read.is_proper_pair:
            #    bad_pairs += 1
            #    continue
            
            if read.mapping_quality < int(mapq):
                poor_quality += 1
                continue
            #only save exact matches and high read quality 		
            #if read.cigarstring == '146M':
            #if read.cigarstring == '179M' or read.cigarstring == '178M':
            
            ## filter reads with too low of map quality and exact cigar match if provided
            if read.mapping_quality > int(mapq):
                if cigar == "":
                    if isinstance(read.reference_name, str):
                        query_to_coords[read.query_name] = read.reference_name
                else:
                    if read.cigarstring==cigar:
                        if isinstance(read.reference_name, str):
                            query_to_coords[read.query_name] = read.reference_name           

        print(f'bad pairs: {bad_pairs} poor quality: {poor_quality}')
        #print(query_to_coords)
        #print(read.reference_name)
        #print(fastq)
        #print(bam)
        print("start")
        coords_to_barcodes = defaultdict(list)
        ###get get barcodes for the aligned reads that passed QC
        for i, barcode in tqdm(enumerate(fastq), 'barcodes', total = n_fastq_records):
            #print(i)
            #print(barcode.quality)
            #print(barcode.name)
            min_base_quality = min(map(ord, barcode.quality))
            
            #filter BCs with poor base quality
            if barcode.name in query_to_coords and min_base_quality >= int(baseq):
                coord = query_to_coords[barcode.name]
                coords_to_barcodes[coord].append(barcode.sequence)
                #print(coord)
        #get unique set of BCs
        #for k, v in coords_to_barcodes.items():   
            #print('set')
            #print(len(set(v)))
            #print(len(v))
        pickle.dump(coords_to_barcodes, open(coords_to_barcodes_fn, 'wb'))
    return pickle.load(open(coords_to_barcodes_fn, 'rb'))

def save_barcodes_per_candidate(d, fn):
    (
        pd.Series(d, name = 'n_barcodes')
        .str.len()
        .rename_axis('coord')
        .reset_index()
        .to_feather(fn)
    )


def save_non_repeated_barcodes_per_candidate(d, fn):
    save_barcodes_per_candidate(
        pd.Series(d, name = 'n_barcodes').apply(set),
        fn
    )





#project_dir = '/ye/yelabstore3/gracieg/scMPRA/pilot_mixing/slimsdata.genomecenter.ucdavis.edu/Data/96qdysy9rg/Unaligned/Project_NAFI_FI48'

coords_to_barcodes = get_coords_to_barcodes(fastq_in, n_fastq,bamfile,n_bam)


save_barcodes_per_candidate(coords_to_barcodes, f'{prefix}_barcodes_per_candidate.feather')
save_non_repeated_barcodes_per_candidate(coords_to_barcodes, f'{prefix}_barcodes_per_candidate-no_repeats.feather')


counter = Counter(itertools.chain.from_iterable(coords_to_barcodes.values()))
pickle.dump(counter, open(f'{prefix}_barcode_counts.pickle', 'wb'))

jackpotted_barcodes = {barcode for barcode, count in counter.most_common(50)}
filtered_coords_to_barcodes = {}
for coord, barcodes in tqdm(coords_to_barcodes.items(), 'filtered barcodes'):
    filtered_coords_to_barcodes[coord] = list(filter(lambda x: x not in jackpotted_barcodes, barcodes))

save_barcodes_per_candidate(filtered_coords_to_barcodes, f'{prefix}_barcodes_per_candidate-no_jackpots.feather')
save_non_repeated_barcodes_per_candidate(filtered_coords_to_barcodes, f'{prefix}_barcodes_per_candidate-no_repeats-no_jackpots.feather')


