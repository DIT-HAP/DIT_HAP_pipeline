
from pathlib import Path
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
from multiprocessing import Pool
from utils.genome_utils import parse_gff_file, extract_coding_gene, extract_CDS_intron, cal_accumlated_CDS_bases

def main(args):

    # read gff file using pandas
    gff_df = parse_gff_file(args.gff)

    args.output.parent.mkdir(parents=True, exist_ok=True)

    # extract coding gene
    region_of_coding_genes = extract_coding_gene(gff_df)

    # extract CDS and intron
    representative_CDS_and_intron = extract_CDS_intron(gff_df)

    # chunk the representative_CDS_and_intron
    chunk_size = int(len(representative_CDS_and_intron)/args.cores)
    chunks = np.array_split(representative_CDS_and_intron, len(representative_CDS_and_intron) // chunk_size + 1)

    # expand the region to each nucleotide
    pool = Pool(args.cores)
    chunked_nucleotide_df = pool.map(generate_nucleotide_df, chunks)
    pool.close()

    all_nucleotide_df = pd.concat(chunked_nucleotide_df, ignore_index=True)

    gene_lists = list(all_nucleotide_df.groupby("Systematic ID"))
    gene_df_lists = [gene_df for gene, gene_df in gene_lists]

    pool = Pool(args.cores)
    accumulated_CDS_base = pool.map(cal_accumlated_CDS_bases, gene_df_lists)
    pool.close()

    all_nucleotide_with_accumulation = pd.concat(accumulated_CDS_base, ignore_index=True)
   
    all_nucleotide_with_accumulation.to_csv(args.output, sep="\t", index=False, header=True)

def expand_region_to_each_nucleotide(row):
    
    start_coordinates = range(row['Start'], row['End'])
    end_coordinates = range(row['Start']+1, row['End']+1)

    nucleotide_df = pd.DataFrame()
    for col in row.index:
        nucleotide_df[col] = [row.loc[col]] * len(start_coordinates)
    nucleotide_df['Start'] = start_coordinates
    nucleotide_df['End'] = end_coordinates
    nucleotide_df['Length'] = 1

    return nucleotide_df

def generate_nucleotide_df(chunk):

    sub_nucleotide_df = pd.DataFrame(columns=chunk.columns)
    for index, row in chunk.iterrows():
        i_nucleotide_df = expand_region_to_each_nucleotide(row)
        sub_nucleotide_df = pd.concat([sub_nucleotide_df, i_nucleotide_df], ignore_index=True)

    return sub_nucleotide_df

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Create bed file for each nucleotide of the coding region in the genome'
    )
    # parser.add_argument(
    #     '-f',
    #     '--fasta',
    #     required=True,
    #     type=Path,
    #     help='Genome file in fasta format'
    # )
    parser.add_argument(
        '-g',
        '--gff',
        required=True,
        type=Path,
        help='GFF file for the genome'
    )
    parser.add_argument(
        "-c",
        "--cores",
        required=True,
        type=int,
        help="The number of cores to use."
    )
    parser.add_argument(
        '-o',
        '--output',
        required=True,
        type=Path,
        help='Output bed file for each nucleotide of the coding region in the genome'
    )

    args = parser.parse_args()

    main(args)