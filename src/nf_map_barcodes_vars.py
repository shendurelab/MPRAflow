#From Sean Whalen with adaptation by Gracie Gordon
import pandas as pd
import pysam
import sys

from itertools import chain


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


#get commandline arguments
variant_mode=True

project_dir=sys.argv[1]
fastq_in=sys.argv[2]
n_fastq_f=sys.argv[3]
library_bam_fn=sys.argv[4]
variants_fn=sys.argv[5]
pref=sys.argv[6]
mapq=sys.argv[7]
baseq=sys.argv[8]

query_to_coords = {} # holds barcode associations
#library_bam_fn = f'{config["index_dir"]}/library.bam' # filename of aligned reads (paired ends first merged to single end)

# cigar operators
BAM_CMATCH     = 0 # M
BAM_CINS       = 1 # I
BAM_CDEL       = 2 # D
BAM_CREF_SKIP  = 3 # N
BAM_CSOFT_CLIP = 4 # S
BAM_CHARD_CLIP = 5 # H
BAM_CPAD       = 6 # P
BAM_CEQUAL     = 7 # =
BAM_CDIFF      = 8 # X
BAM_CBACK      = 9 # B

# load simple text file of variants
# expected columns:
# reference_name    variant_positions   ref_bases   alt_bases
# reference name is the geneomic coordinate of each library insert
# variant_positions is the comma separated list of indices where we are testing SNVs
# ref_bases is a comma separated list of reference bases at these positions
# alt_bases is a comma separated list of alternate allele bases at these positions
# if there is only a single SNV per insert, then just have the position with no comma
variants = pd.read_csv(
    variants_fn,
    sep = '\t',
    index_col = 'reference_name'
)

# create a simple mapping from reference name to a python list of integer variant positions
reference_name_to_variant_positions = (
    variants
    ['variant_positions']
    .apply(
        lambda positions:
            list(map(int, positions.split(',')))
            if isinstance(positions, str) else [positions]
    )
    .to_dict()
)

# convert comma separated lists of ref and alt bases into python lists
variants['ref_bases'] = variants['ref_bases'].str.split(',')
variants['alt_bases'] = variants['alt_bases'].str.split(',')

# loop over each aligned read, check QC, check cigar scores at optional variant positions, and add association from sequencer id to mapped insert coordinate if everything looks good
for read in pysam.AlignmentFile(library_bam_fn, 'rb'):
    if read.is_unmapped or read.mapping_quality < config['min_mapping_quality']:
        continue

    reference_name = read.reference_name

    if variant_mode:
        # get position of variants for current insert
        variant_positions = reference_name_to_variant_positions[reference_name]

        # start building an expanded list of cigar operator codes that we can index into using variant positions.
        # if read is clipped and doesn't start at position 0, add a buffer of any character at the start so our indices line up
        cigar = ['X'] * read.reference_start

        # adjust the cigar string for inserts, deletions, and other cigar operators so that the original variant indices line up.
        # this was simplified from my code and could need tweaking
        for cigar_operator, cigar_length in read.cigartuples:
            if cigar_operator == BAM_CINS:
                pass
            elif cigar_operator in {BAM_CDEL, BAM_CREF_SKIP, BAM_CSOFT_CLIP, BAM_CEQUAL, BAM_CDIFF}:
                cigar.extend([cigar_operator] * cigar_length)
            else:
                raise Exception(f'unsupported cigar operator: {cigar_operator} in {read.cigartuples}')

        # if any variant position extends past read length, skip it. can happen if alignment is messy, but should be rare
        if max(variant_positions) >= len(cigar):
            continue

        # get cigar scores only at position of variants: doing this with a list/itemgetter is much faster than creating numpy arrays
        positions_getter = itemgetter(*variant_positions)
        variants_cigar = positions_getter(cigar)
        if isinstance(variants_cigar, int):
                variants_cigar = [variants_cigar]

        # we expect a match at all variant positions
        expected_variants_cigar = [BAM_CEQUAL] * len(variant_positions)

        # if any of the observed cigar scores is not a match, skip the read
        if variants_cigar != expected_variants_cigar:
            continue

    # if we passed QC, and the sequence matched at all variant positions in variant mode,
    # then associate the read id from sequencer to the mapped coordinate of the insert
    query_to_coords[read.query_name] = reference_name


coords_to_barcodes = defaultdict(list)
###get get barcodes for the aligned reads that passed QC
for i, barcode in tqdm(enumerate(fastq), 'barcodes', total = n_fastq_records):

    min_base_quality = min(map(ord, barcode.quality))

    #filter BCs with poor base quality
    if barcode.name in query_to_coords and min_base_quality >= int(baseq):
        coord = query_to_coords[barcode.name]
        coords_to_barcodes[coord].append(barcode.sequence)



pickle.dump(coords_to_barcodes, open(coords_to_barcodes_fn, 'wb'))

counter = Counter(itertools.chain.from_iterable(coords_to_barcodes.values()))
pickle.dump(counter, open(f'{prefix}_barcode_counts.pickle', 'wb'))
jackpotted_barcodes = {barcode for barcode, count in counter.most_common(50)}
filtered_coords_to_barcodes = {}
for coord, barcodes in tqdm(coords_to_barcodes.items(), 'filtered barcodes'):
    filtered_coords_to_barcodes[coord] = list(filter(lambda x: x not in jackpotted_barcodes, barcodes))

save_barcodes_per_candidate(filtered_coords_to_barcodes, f'{prefix}_barcodes_per_candidate-no_jackpots.feather')
save_non_repeated_barcodes_per_candidate(filtered_coords_to_barcodes, f'{prefix}_barcodes_per_candidate-no_repeats-no_jackpots.feather')


