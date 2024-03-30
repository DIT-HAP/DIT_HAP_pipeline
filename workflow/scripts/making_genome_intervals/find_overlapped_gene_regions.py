import argparse
import pybedtools
from pybedtools import BedTool
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge overlapping genes in a BED file."
    )
    parser.add_argument("-i", "--input", help="Path to the input BED file")
    parser.add_argument("-o", "--output", help="Path to the output BED file")
    return parser.parse_args()


def main():
    args = parse_args()
    merge_overlapping_genes(args.input, args.output)


def merge_overlapping_genes(input_file, output_file):
    # Load the bed file
    coding_gene = BedTool(input_file)

    # Find overlapping intervals
    overlaps = coding_gene.intersect(
        coding_gene, wa=True, wb=True, header=True)

    # Process overlaps
    results = overlaps.each(find_overlapping_regions)

    # drop duplicate lines
    results = (
        results.saveas()
        .to_dataframe(names=["#Chr", "Start", "End", "Systematic ID", "Type", "Strand"])
        .drop_duplicates(subset=["#Chr", "Start", "End"], keep="first")
    )

    # overlapped genes
    overlapped_region = results[results["Type"] == "Overlapping genes"].copy()
    overlapped_region["Length"] = (
        overlapped_region["End"] - overlapped_region["Start"]
    )

    # Save the results
    overlapped_region.to_csv(output_file, sep="\t", index=False, header=True)


def find_overlapping_regions(feature):
    chr_a, chr_b = feature[0], feature[7]
    start_a, start_b = int(feature[1]), int(feature[8])
    end_a, end_b = int(feature[2]), int(feature[9])
    name_a, name_b = feature[3], feature[10]
    score_a, score_b = feature[4], feature[11]
    strand_a, strand_b = feature[5], feature[12]

    chrom = chr_a
    start = max(start_a, start_b)
    end = min(end_a, end_b)

    if name_a != name_b:
        name = name_a + "," + name_b
        score = "Overlapping genes"
        strand = strand_a + "," + strand_b
        return pybedtools.create_interval_from_list(
            [chrom, str(start), str(end), name, score, strand]
        )
    else:
        return pybedtools.create_interval_from_list(
            [chrom, str(start), str(end), name_a, score_a, strand_a]
        )


if __name__ == "__main__":
    main()
