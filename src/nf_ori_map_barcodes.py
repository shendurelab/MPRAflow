#!/usr/bin/env python

"""
Barcode association code by Sean Whalen (sean.whalen@gladstone.ucsf.edu), Pollard Lab, 2019-2020
Variant processing concept by Katie Pollard
"""

import argparse
import os
import pandas as pd
import pickle
import pysam

from operator import itemgetter
from tqdm import tqdm


# removes barcodes with low coverage per candidate
def filter_low_coverage_barcodes(barcodes_coords):
    return (
        barcodes_coords
        .groupby(['coord', 'barcode'])
        .size()
        .rename('barcode_coverage')
        .reset_index()
        .query('barcode_coverage >= @args.min_barcode_coverage')
    )


# keep only the most frequent candidate mapping per barcode
def filter_minority_associations(barcodes_coords):
    percent_association = (
        barcodes_coords
        .groupby('barcode')
        ['coord']
        .value_counts(normalize = True)
    )
    return (
        percent_association
        [percent_association >= args.min_barcode_mapping_percent]
        .rename('barcode_mapping_percent')
        .reset_index()
    )


parser = argparse.ArgumentParser()
parser.add_argument('--bam_fn', required = True)
parser.add_argument('--barcode_fn', required = True)
parser.add_argument('--variant_fn', default = None)
parser.add_argument('--output_fn', required = True)
parser.add_argument('--min_mapping_quality', default = 30, type = int)
parser.add_argument('--min_barcode_quality', default = 30, type = int)
parser.add_argument('--min_barcode_entropy', default = 0.75, type = float)
parser.add_argument('--min_barcode_coverage', default = 3, type = int)
parser.add_argument('--min_barcode_mapping_percent', default = 0.75, type = float)
parser.add_argument('--max_reference_start', default = 2, type = int)
parser.add_argument('--n_bam_records', default = None, type = int)
parser.add_argument('--n_barcode_records', default = None, type = int)
args = parser.parse_args()

assert 0.51 <= args.min_barcode_mapping_percent <= 1

if args.variant_fn != None:
    variant_mode_enabled = os.path.exists(args.variant_fn)
    print(f'variant mode: {variant_mode_enabled}')
else:
    variant_mode_enabled = False 

# cigar operators for variant mode
# https://samtools.github.io/hts-specs/SAMv1.pdf
BAM_CMATCH     = 0 # M (alignment match)
BAM_CINS       = 1 # I (insertion)
BAM_CDEL       = 2 # D (deletion)
BAM_CREF_SKIP  = 3 # N (skipped reference)
BAM_CSOFT_CLIP = 4 # S (soft clipping)
BAM_CHARD_CLIP = 5 # H (hard clipping)
BAM_CPAD       = 6 # P (padding)
BAM_CEQUAL     = 7 # = (sequence match)
BAM_CDIFF      = 8 # X (sequence difference)
consumes_query_and_reference_ops = {BAM_CMATCH, BAM_CEQUAL, BAM_CDIFF}
consumes_query_ops = {BAM_CINS, BAM_CSOFT_CLIP}

if variant_mode_enabled:
    print('setting mapping quality to 0 so alleles with low scores are not filtered')
    args.min_mapping_quality = 0

    # load simple text file of variants
    # expected columns:
    # reference_name    variant_positions   ref_bases   alt_bases
    # reference name is the geneomic coordinate of each library insert
    # variant_positions is the comma separated list of indices where we are testing SNVs
    # ref_bases is a comma separated list of reference bases at these positions
    # alt_bases is a comma separated list of alternate allele bases at these positions
    # if there is only a single SNV per insert, then just have the position with no comma
    variants = pd.read_csv(
        args.variant_fn,
        sep = '\t',
        index_col = 'reference_name'
    )
    variant_reference_names = set(variants.index)

    # create a simple mapping from reference name to a python list of integer variant positions:
    # splitting the position string by comma, map each element to an integer, and convert back to a list.
    # if the variant column is for single variants, pandas will correctly infer dtype=int,
    # so create a single element list for consistent downstream logic.
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

    # convert comma separated lists of ref and alt bases into python objects
    variants['ref_bases'] = (
        variants
        ['ref_bases']
        .str.split(',')
        .apply(tuple)
    )
    variants['alt_bases'] = (
        variants
        ['alt_bases']
        .str.split(',')
        .apply(tuple)
    )

# holds association from query name (sequencer read id) to mapped genomic coordinate.
# query name is later used to associate coordinates with barcodes.
query_to_coord = {}

# associate read ids from sequencer with mapped genomic coordinates
for read in tqdm(pysam.AlignmentFile(args.bam_fn, 'rb'), total = args.n_bam_records):
    # ensure read is mapped and mapping quality is above threshold
    if read.is_unmapped or read.mapping_quality < args.min_mapping_quality or read.reference_start > args.max_reference_start:
        continue

    # if in variant mode and current read is a variant, check for sequence matches at variant positions
    if variant_mode_enabled and read.reference_name in variant_reference_names:
        # variant positions are given with respect to the reference,
        # but cigar string is given with respect to the query.
        # need to adjust query string so it lines up with reference
        # to check for sequence matches at variant positions.
        reference = read.get_reference_sequence().upper()
        query = ''
        offset = 0
        for cigar_operator, cigar_length in read.cigartuples:
            if cigar_operator in consumes_query_and_reference_ops:
                query += reference[offset:(offset + cigar_length)]
                offset += cigar_length
            elif cigar_operator in consumes_query_ops:
                offset += cigar_length

        variant_positions = reference_name_to_variant_positions[read.reference_name]
        if read.reference_start > 0:
            variant_positions = [_ - read.reference_start for _ in variant_positions]

        if max(variant_positions) >= min(len(reference), len(query)):
            continue

        reference_variant_bases = itemgetter(*variant_positions)(reference)
        if isinstance(reference_variant_bases, str):
            reference_variant_bases = (reference_variant_bases,)

        query_variant_bases = itemgetter(*variant_positions)(query)
        if isinstance(query_variant_bases, str):
            query_variant_bases = (query_variant_bases,)

        expected_reference_variant_bases = variants.loc[read.reference_name, 'ref_bases']
        expected_alt_variant_bases = variants.loc[read.reference_name, 'alt_bases']

        assert reference_variant_bases in {expected_reference_variant_bases, expected_alt_variant_bases}

        # discard if query has a sequence mismatch at any variant position
        if query_variant_bases != reference_variant_bases:
            continue

    # perform association: read.query_name is sequencer read id, read.reference_name is mapped genomic coordinate
    query_to_coord[read.query_name] = read.reference_name

# associate barcodes to mapped genomic coordinates (using the read id -> coordinate associations we just created)
barcodes_coords = [
    (barcode.sequence, query_to_coord[barcode.name])
    for barcode in tqdm(pysam.FastxFile(args.barcode_fn), total = args.n_barcode_records)
    if barcode.name in query_to_coord and min(barcode.get_quality_array()) >= args.min_barcode_quality
]
barcodes_coords = pd.DataFrame(barcodes_coords, columns = ['barcode', 'coord'])
del query_to_coord

print('filtering')
print(barcodes_coords)
barcodes_coords = filter_low_coverage_barcodes(barcodes_coords)
barcodes_coords = filter_minority_associations(barcodes_coords)

print(barcodes_coords)

#filtered_barcode_to_coord = (
#    barcodes_coords
#    [['barcode', 'coord']]
#    .set_index('barcode')
#    .squeeze()
#    .to_dict()
#)

filtered_barcode_to_coord = (
    barcodes_coords
    [['coord', 'barcide']]
    .set_index('coord')
    .squeeze()
    .to_dict()
)

print(filtered_barcode_to_coord)
pickle.dump(
    filtered_barcode_to_coord,
    open(args.output_fn, 'wb')
)

#add in plotting script here. 



