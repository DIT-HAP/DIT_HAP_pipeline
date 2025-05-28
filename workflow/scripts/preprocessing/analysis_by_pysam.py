import pysam
import argparse
from collections import defaultdict
import concurrent.futures # Added for multi-threading

def get_read_info(read, unique_tags_list):
    """Extracts relevant information from a pysam AlignedSegment object, including tags, strand, NCIGAR, chrom, pos, and flag."""
    if read is None:
        info = {
            "MAPQ": "N/A",
            "LEN": "N/A",
            "CIGAR": "N/A",
            "STRAND": "N/A",
            "NCIGAR": "N/A",
            "CHROM": "N/A",
            "POS": "N/A",
            "FLAG": "N/A",
        }
        for tag_name in unique_tags_list:
            info[tag_name] = "N/A"
        return info
    
    read_specific_tags = dict(read.get_tags())
    strand = "-" if read.is_reverse else "+"
    
    ncigar = "N/A"
    if read.cigartuples:
        ncigar = len(read.cigartuples)
    elif read.cigarstring: # if cigartuples is None but cigarstring exists (e.g., "*")
        ncigar = 0 # Or treat as N/A if "*" means no alignment for CIGAR ops count
    # If cigarstring is also None or empty, ncigar remains "N/A"

    chrom = read.reference_name if read.reference_name is not None else "N/A"
    # pysam uses 0-based coordinates for reference_start. SAM format is 1-based.
    pos = str(read.reference_start + 1) if read.reference_start is not None and read.reference_start != -1 else "N/A"
    flag = str(read.flag) if read.flag is not None else "N/A"

    info = {
        "MAPQ": read.mapping_quality,
        "LEN": read.query_alignment_length,
        "CIGAR": read.cigarstring if read.cigarstring else "N/A",
        "STRAND": strand,
        "NCIGAR": str(ncigar), # Ensure NCIGAR is a string for consistent output
        "CHROM": chrom,
        "POS": pos,
        "FLAG": flag,
    }
    
    for tag_name in unique_tags_list:
        value = read_specific_tags.get(tag_name, "N/A")
        # Ensure consistent string representation, especially for booleans or None
        if value is True:
            info[tag_name] = "True"
        elif value is False:
            info[tag_name] = "False"
        elif value is None: # Should be caught by .get(tag_name, "N/A") but as an extra check
            info[tag_name] = "N/A"
        else:
            info[tag_name] = str(value) # Convert other types to string
            
    return info

def process_qname_pair(qname, read1, read2, unique_tags_list):
    """Processes a single qname pair and returns the formatted output line data, including all details."""
    is_proper_pair = "N/A"
    if read1 and read1.is_paired:
        is_proper_pair = "Yes" if read1.is_proper_pair else "No"
    elif read2 and read2.is_paired:
        is_proper_pair = "Yes" if read2.is_proper_pair else "No"
    elif (read1 and not read1.is_paired) or (read2 and not read2.is_paired):
        is_proper_pair = "Single_End_Or_Flag_Issue"

    r1_info = get_read_info(read1, unique_tags_list)
    r2_info = get_read_info(read2, unique_tags_list)

    output_line = [qname]
    output_line.extend([
        str(r1_info["MAPQ"]), 
        str(r1_info["LEN"]), 
        r1_info["CIGAR"],
        r1_info["STRAND"],
        r1_info["NCIGAR"],
        r1_info["CHROM"],
        r1_info["POS"],
        r1_info["FLAG"],
    ])
    for tag_name in unique_tags_list:
        output_line.append(r1_info[tag_name]) # Already stringified in get_read_info

    output_line.extend([
        str(r2_info["MAPQ"]), 
        str(r2_info["LEN"]), 
        r2_info["CIGAR"],
        r2_info["STRAND"],
        r2_info["NCIGAR"],
        r2_info["CHROM"],
        r2_info["POS"],
        r2_info["FLAG"],
    ])
    for tag_name in unique_tags_list:
        output_line.append(r2_info[tag_name]) # Already stringified in get_read_info
    
    output_line.append(is_proper_pair)
    return output_line

def main():
    parser = argparse.ArgumentParser(description="Extract comprehensive summary for read pairs from a BAM/SAM file, including tags, strand, NCIGAR, chrom, pos, and flag.")
    parser.add_argument("-i", "--input_bam", required=True, help="Path to the input BAM/SAM file.")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output tab-delimited file.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use for processing. (Default: 4)")
    args = parser.parse_args()

    print(f"Using {args.threads} threads for processing.")
    print("Reading BAM/SAM file, grouping reads by QNAME, and collecting all unique tags...")

    read_pairs = defaultdict(lambda: [None, None])
    all_tags_set = set()
    
    try:
        samfile = pysam.AlignmentFile(args.input_bam, "rb" if args.input_bam.endswith(".bam") else "r", threads=args.threads)
    except ValueError as e:
        print(f"Error opening SAM/BAM file: {e}")
        print("Please ensure the file is a valid SAM/BAM file and the index is present if it's a BAM file and you are seeking.")
        return
    except FileNotFoundError:
        print(f"Error: Input file not found at {args.input_bam}")
        return

    read_count = 0
    for read in samfile:
        read_count += 1
        if read_count % 1000000 == 0:
            print(f"  Processed {read_count // 1000000}M alignments for grouping and tag collection...")

        if not read.is_unmapped: # Collect tags only from mapped reads that might be part of a pair
            current_read_tags = read.get_tags()
            for tag_name, _ in current_read_tags:
                all_tags_set.add(tag_name)

        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue

        qname = read.query_name
        if read.is_read1:
            if read_pairs[qname][0] is None:
                 read_pairs[qname][0] = read
        elif read.is_read2:
            if read_pairs[qname][1] is None:
                read_pairs[qname][1] = read
    
    samfile.close()
    
    sorted_unique_tags = sorted(list(all_tags_set))
    if not sorted_unique_tags:
        print("Warning: No optional tags found in the BAM/SAM file.")
    else:
        print(f"Found unique tags: {sorted_unique_tags}")

    print(f"Finished grouping reads. Total alignments scanned: {read_count}.")
    print(f"Found {len(read_pairs)} unique query names to process.")
    print("Processing pairs and writing to output...")

    # Dynamically build header
    header_fields = ["QueryName", 
                     "R1_MAPQ", "R1_LEN", "R1_CIGAR", "R1_Strand", "R1_NCIGAR", 
                     "R1_Chrom", "R1_Pos", "R1_Flag"]
    for tag_name in sorted_unique_tags:
        header_fields.append(f"R1_{tag_name}")
    
    header_fields.extend(["R2_MAPQ", "R2_LEN", "R2_CIGAR", "R2_Strand", "R2_NCIGAR",
                          "R2_Chrom", "R2_Pos", "R2_Flag"])
    for tag_name in sorted_unique_tags:
        header_fields.append(f"R2_{tag_name}")
    
    header_fields.append("Is_Proper_Pair")

    with open(args.output_file, "w") as outfile:
        outfile.write("\t".join(header_fields) + "\n")

        sorted_qnames = sorted(read_pairs.keys())
        results_to_write = [None] * len(sorted_qnames)
        qname_to_index = {qname: i for i, qname in enumerate(sorted_qnames)}

        with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
            future_to_qname = {
                executor.submit(process_qname_pair, qname, read_pairs[qname][0], read_pairs[qname][1], sorted_unique_tags): qname 
                for qname in sorted_qnames
            }
            
            processed_count = 0
            for future in concurrent.futures.as_completed(future_to_qname):
                qname = future_to_qname[future]
                try:
                    result_line_data = future.result()
                    results_to_write[qname_to_index[qname]] = result_line_data
                    processed_count +=1
                    if processed_count % 100000 == 0:
                        print(f"  Processed {processed_count} QNAMEs for output generation...")
                except Exception as exc:
                    print(f'{qname} generated an exception: {exc}')
                    # Ensure the error line has the correct number of columns
                    results_to_write[qname_to_index[qname]] = [qname] + ["ERROR"] * (len(header_fields) - 1)
        
        for result_line_data in results_to_write:
            if result_line_data:
                outfile.write("\t".join(map(str, result_line_data)) + "\n") # Ensure all elements are strings

    print(f"Processed {len(sorted_qnames)} unique query names.")
    print(f"Output written to {args.output_file}")
    print("Memory Usage Note: For very large BAM/SAM files, consider pre-sorting by QNAME (e.g., 'samtools sort -n') for optimal memory performance, as this script currently loads all read pairs into memory before processing.")

if __name__ == "__main__":
    main()
