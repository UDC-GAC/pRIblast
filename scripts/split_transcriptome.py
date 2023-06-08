#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse

def read_fasta_file(file_name):
    seqs = {}
    with open(file_name) as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].rstrip('\n')
                seqs[name] = ''
            else:
                seqs[name] += line.rstrip('\n')
    return seqs

def write_fasta_file(seqs, file_name):
    with open(file_name, 'w') as f:
        for name, seq in seqs.items():
            pretty_seq = '\n'.join(
                [seq[i:i + 60] for i in range(0, len(seq), 60)]
            )
            f.write(f'>{name}\n')
            f.write(f'{pretty_seq}\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Split a fasta file into two files, one with lncRNA and
one with RNA. The input file must have the name of the lncRNA sequences
containing the word "lncRNA".'''
    )
    parser.add_argument(
        'fasta_file',
        help='Fasta file to be split.',
    )
    args = parser.parse_args()

    if not os.path.exists(args.fasta_file):
        raise Exception(f'{args.fasta_file} does not exist.')
    if not os.path.isfile(args.fasta_file):
        raise Exception(f'{args.fasta_file} is not a file.')

    RNA_seqs = read_fasta_file(args.fasta_file)
    lncRNA_seqs = {
        name: seq
        for name, seq in RNA_seqs.items()
        if "lncRNA" in name
    }

    input_dir = os.path.dirname(args.fasta_file)
    input_file_name = os.path.basename(args.fasta_file)
    input_file_extension = os.path.splitext(input_file_name)[1]
    input_file_name = os.path.splitext(input_file_name)[0]

    output_lncRNA_file_name = f'{input_file_name}_lncRNA{input_file_extension}'
    output_RNA_file_name = f'{input_file_name}_RNA{input_file_extension}'

    output_lncRNA_file_path = os.path.join(input_dir, output_lncRNA_file_name)
    output_RNA_file_path = os.path.join(input_dir, output_RNA_file_name)

    write_fasta_file(RNA_seqs, output_RNA_file_path)
    write_fasta_file(lncRNA_seqs, output_lncRNA_file_path)
