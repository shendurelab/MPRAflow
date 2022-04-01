import argparse
import csv
import itertools
import typing
import pandas as pd
import Bio.SeqIO as bsio
import Bio.Seq as bseq


class Namespace(argparse.Namespace):
    in_design: typing.TextIO
    in_labels: typing.TextIO
    variants: typing.TextIO
    out_design: typing.TextIO
    out_labels: typing.TextIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('in_design', type=argparse.FileType('r'))
    parser.add_argument('in_labels', type=argparse.FileType('r'))
    parser.add_argument('variants', type=argparse.FileType('r'))
    parser.add_argument('out_design', type=argparse.FileType('w'))
    parser.add_argument('out_labels', type=argparse.FileType('w'))
    args = parser.parse_args(namespace=Namespace())

    # Read the design
    design = pd.DataFrame([[fa.name, fa.seq] for fa in bsio.parse(args.in_design, 'fasta')], columns=['name', 'sequence'])

    # Read the labels
    labels = pd.read_table(args.in_labels, names=['name', 'label'])

    # Read the mutations
    mutations = pd.read_table(args.variants, names=['name', 'pos', 'reference', 'alternate'], header=0)

    # Merge them all together
    merged_df = design.merge(labels, how='outer', on='name').merge(mutations, how='outer', on='name').fillna(value={'pos': '', 'reference': '', 'alternate': '', 'label': 'unknown'})

    # Split mutations
    out_labs = []
    out_fa = []

    def split_rows(row):
        if row['pos']:
            variants = [(int(x), y, z) for x, y, z in zip(row['pos'].split(','), row['reference'].split(','), row['alternate'].split(','))]
            for is_mut in itertools.product(range(2), repeat=len(variants)):
                seq = row['sequence']
                name = row['name']
                for flag, (pos, ref, alt) in zip(is_mut, variants):
                    seq = seq[:pos - 1] + [ref, alt][flag] + seq[pos:]
                    name += f':{pos}={[ref, alt][flag]}'
                out_labs.append([name, row['label']])
                out_fa.append(bsio.SeqRecord(seq, name, name, row['label']))
        else:
            seq = row['sequence']
            out_labs.append([row['name'], row['label']])
            out_fa.append(bsio.SeqRecord(seq, row['name'], row['name'], row['label']))

    merged_df.apply(split_rows, axis=1)

    # Write output
    bsio.write(out_fa, args.out_design, 'fasta')
    pd.DataFrame(out_labs).to_csv(args.out_labels, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE)


if __name__ == '__main__':
    main()
