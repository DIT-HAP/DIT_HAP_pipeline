import pysam
import argparse
from collections import defaultdict
# import concurrent.futures # No longer needed for streaming

def get_read_info(read, unique_tags_list):
    """Extracts relevant information from a pysam AlignedSegment object, including tags, strand, NCIGAR, chrom, pos (0-based), ref_start (0-based), ref_end (0-based), and flag."""
    if read is None:
        info = {
            "MAPQ": "N/A",
            "LEN": "N/A",
            "CIGAR": "N/A",
            "STRAND": "N/A",
            "NCIGAR": "N/A",
            "CHROM": "N/A",
            "POS": "N/A",         # 0-based
            "Ref_Start": "N/A", # 0-based
            "Ref_End": "N/A",   # 0-based (exclusive end)
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
    elif read.cigarstring: 
        ncigar = 0 

    chrom = read.reference_name if read.reference_name is not None else "N/A"
    # Output 0-based POS (leftmost alignment position)
    pos = str(read.reference_start) if read.reference_start is not None and read.reference_start != -1 else "N/A"
    flag = str(read.flag) if read.flag is not None else "N/A"
    # Output 0-based reference_start
    ref_start_val = str(read.reference_start) if read.reference_start is not None else "N/A"
    # Output 0-based reference_end (exclusive end of alignment)
    ref_end_val = str(read.reference_end) if read.reference_end is not None else "N/A"

    info = {
        "MAPQ": read.mapping_quality,
        "LEN": read.query_alignment_length,
        "CIGAR": read.cigarstring if read.cigarstring else "N/A",
        "STRAND": strand,
        "NCIGAR": str(ncigar),
        "CHROM": chrom,
        "POS": pos, # Now 0-based
        "Ref_Start": ref_start_val, # 0-based
        "Ref_End": ref_end_val,     # 0-based (exclusive end)
        "FLAG": flag,
    }
    
    for tag_name in unique_tags_list:
        value = read_specific_tags.get(tag_name, "N/A")
        if value is True:
            info[tag_name] = "True"
        elif value is False:
            info[tag_name] = "False"
        elif value is None:
            info[tag_name] = "N/A"
        else:
            info[tag_name] = str(value)
            
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
        r1_info["POS"],         # 0-based
        r1_info["Ref_Start"], # 0-based
        r1_info["Ref_End"],   # 0-based (exclusive end)
        r1_info["FLAG"],
    ])
    for tag_name in unique_tags_list:
        output_line.append(r1_info[tag_name])

    output_line.extend([
        str(r2_info["MAPQ"]), 
        str(r2_info["LEN"]), 
        r2_info["CIGAR"],
        r2_info["STRAND"],
        r2_info["NCIGAR"],
        r2_info["CHROM"],
        r2_info["POS"],         # 0-based
        r2_info["Ref_Start"], # 0-based
        r2_info["Ref_End"],   # 0-based (exclusive end)
        r2_info["FLAG"],
    ])
    for tag_name in unique_tags_list:
        output_line.append(r2_info[tag_name])
    
    output_line.append(is_proper_pair)
    return output_line

def main():
    parser = argparse.ArgumentParser(description="Extract comprehensive summary for read pairs from a BAM/SAM file, including tags, strand, NCIGAR, chrom, pos (0-based), ref_start (0-based), ref_end (0-based), and flag. Assumes qname-sorted input for memory efficiency.")
    parser.add_argument("-i", "--input_bam", required=True, help="Path to the input BAM/SAM file (must be qname-sorted).")
    parser.add_argument("-o", "--output_file", required=True, help="Path to the output tab-delimited file.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads for BAM decompression. (Default: 4)")
    args = parser.parse_args()

    print(f"Using {args.threads} threads for BAM decompression.")
    print("Processing BAM/SAM file in a memory-friendly streaming mode (requires qname-sorted input).")
    print("Using a fixed list of tags: ['AS', 'MC', 'MD', 'MQ', 'NM', 'SA', 'XA', 'XS']")

    # Predefined list of tags
    sorted_unique_tags = ['AS', 'MC', 'MD', 'MQ', 'NM', 'SA', 'XA', 'XS']
    # Ensure it's sorted for consistent header order, though the provided list is already sorted.
    sorted_unique_tags.sort()

    # The first pass for tag collection has been removed.
    # Directly build header using the fixed list of tags.

    header_fields = ["QueryName", 
                     "R1_MAPQ", "R1_LEN", "R1_CIGAR", "R1_Strand", "R1_NCIGAR", 
                     "R1_Chrom", "R1_Pos", "R1_Ref_Start", "R1_Ref_End", "R1_Flag"]
    for tag_name in sorted_unique_tags:
        header_fields.append(f"R1_{tag_name}")
    
    header_fields.extend(["R2_MAPQ", "R2_LEN", "R2_CIGAR", "R2_Strand", "R2_NCIGAR",
                          "R2_Chrom", "R2_Pos", "R2_Ref_Start", "R2_Ref_End", "R2_FLAG"])
    for tag_name in sorted_unique_tags:
        header_fields.append(f"R2_{tag_name}")
    
    header_fields.append("Is_Proper_Pair")

    print("Processing reads and writing to output...") # Changed from "Second pass:..."
    with open(args.output_file, "w") as outfile:
        outfile.write("\t".join(header_fields) + "\n")

        current_qname = None
        current_r1 = None
        current_r2 = None
        processed_qname_count = 0
        read_count_pass = 0 # Renamed from read_count_pass2

        try:
            # This is now the only pass over the BAM file
            samfile_process = pysam.AlignmentFile(args.input_bam, "rb" if args.input_bam.endswith(".bam") else "r", threads=args.threads)
            
            for read in samfile_process:
                read_count_pass +=1
                if read_count_pass % 2000000 == 0:
                    print(f"  Scanned {read_count_pass // 1000000}M alignments...")

                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                qname = read.query_name

                if qname != current_qname:
                    if current_qname is not None:
                        output_line_data = process_qname_pair(current_qname, current_r1, current_r2, sorted_unique_tags)
                        outfile.write("\t".join(map(str, output_line_data)) + "\n")
                        processed_qname_count += 1
                        if processed_qname_count % 500000 == 0:
                             print(f"  Processed and wrote {processed_qname_count} QNAME groups...")
                    
                    current_qname = qname
                    current_r1 = None
                    current_r2 = None
                
                if read.is_read1:
                    if current_r1 is None: 
                        current_r1 = read
                elif read.is_read2:
                    if current_r2 is None:
                        current_r2 = read
            
            if current_qname is not None:
                output_line_data = process_qname_pair(current_qname, current_r1, current_r2, sorted_unique_tags)
                outfile.write("\t".join(map(str, output_line_data)) + "\n")
                processed_qname_count += 1
            
            samfile_process.close()

        except ValueError as e:
            print(f"Error opening or processing SAM/BAM file: {e}") # Simplified error message
            return 
        except FileNotFoundError:
            print(f"Error: Input file not found at {args.input_bam}")
            return
            
    print(f"Finished processing. Total QNAME groups written: {processed_qname_count}. Total alignments scanned: {read_count_pass}.") # Updated variable name
    print(f"Output written to {args.output_file}")
    print("Memory Usage Note: This script processes BAM/SAM files in a memory-efficient, streaming fashion (assuming qname-sorted input) using a fixed list of tags.")

if __name__ == "__main__":
    main()
