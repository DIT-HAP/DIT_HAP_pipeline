import pandas as pd
import argparse
from collections import defaultdict

def process_insertion_counts(input_tsv, output_tsv, chunk_size=500000):
    """
    Process the TSV output from parse_bam_to_tsv.py to generate insertion count summary using chunking.
    
    For each read pair:
    - If R1_Strand is '+', use R1_Ref_Start as the insertion coordinate
    - If R1_Strand is '-', use R1_Ref_End as the insertion coordinate
    
    Output: Chr, Coordinate, +, -
    """
    
    print(f"Processing TSV file in chunks of {chunk_size:,} rows...")
    
    # Initialize counters and data structures
    insertion_counts = defaultdict(lambda: {'+': 0, '-': 0})
    total_rows = 0
    valid_rows = 0
    invalid_rows = 0
    chunk_count = 0
    
    # Required columns
    required_columns = ['R1_Strand', 'R1_Chrom', 'R1_Ref_Start', 'R1_Ref_End']
    
    print("Starting chunked processing...")
    
    try:
        # Create chunk iterator
        chunk_iterator = pd.read_csv(
            input_tsv, 
            sep='\t',
            chunksize=chunk_size,
            na_values=['N/A', 'NA', '']
        )
        
        # Process each chunk
        for chunk_df in chunk_iterator:
            chunk_df = chunk_df[required_columns].copy()
            chunk_count += 1
            chunk_rows = len(chunk_df)
            total_rows += chunk_rows
            
            if chunk_count % 10 == 0:
                print(f"  Processing chunk {chunk_count}, total rows processed: {total_rows:,}")
            
            # Validate required columns on first chunk
            if chunk_count == 1:
                missing_columns = [col for col in required_columns if col not in chunk_df.columns]
                if missing_columns:
                    raise ValueError(f"Missing required columns: {missing_columns}")
                
                print(f"Found all required columns: {required_columns}")
                print(f"First chunk size: {chunk_rows:,} rows")
            
            # Filter out rows where essential data is missing
            valid_chunk = chunk_df[
                (chunk_df['R1_Strand'].notna()) & 
                (chunk_df['R1_Chrom'].notna()) & 
                (chunk_df['R1_Ref_Start'].notna()) & 
                (chunk_df['R1_Ref_End'].notna()) &
                (chunk_df['R1_Strand'].isin(['+', '-']))
            ].copy()
            
            chunk_valid_rows = len(valid_chunk)
            chunk_invalid_rows = chunk_rows - chunk_valid_rows
            valid_rows += chunk_valid_rows
            invalid_rows += chunk_invalid_rows
            
            if chunk_count == 1 or chunk_count % 10 == 0:
                print(f"    Chunk {chunk_count}: {chunk_valid_rows:,}/{chunk_rows:,} valid rows ({chunk_valid_rows/chunk_rows*100:.1f}%)")
            
            # Process valid rows in this chunk
            if chunk_valid_rows > 0:
                # Determine insertion coordinate based on strand
                def get_insertion_coordinate(row):
                    try:
                        if row['R1_Strand'] == '+':
                            # for +, TTAA[Genome] use the last A position
                            return int(row['R1_Ref_Start'] + 4)
                        elif row['R1_Strand'] == '-':
                            # for -, [Genome]AATT use the T position, is already the last A position
                            return int(row['R1_Ref_End'])
                        else:
                            return None
                    except (ValueError, TypeError):
                        return None
                
                valid_chunk['Insertion_Coordinate'] = valid_chunk.apply(get_insertion_coordinate, axis=1)
                
                # Remove rows where coordinate determination failed
                valid_chunk = valid_chunk[valid_chunk['Insertion_Coordinate'].notna()]
                
                # Count insertions for this chunk
                for _, row in valid_chunk.iterrows():
                    try:
                        chrom = row['R1_Chrom']
                        coord = int(row['Insertion_Coordinate'])
                        strand = row['R1_Strand']
                        
                        key = (chrom, coord)
                        insertion_counts[key][strand] += 1
                    except (ValueError, TypeError, KeyError):
                        # Skip rows with invalid data
                        continue
        
        print(f"\nCompleted processing {chunk_count} chunks")
        print(f"Total rows processed: {total_rows:,}")
        print(f"Valid rows: {valid_rows:,} ({valid_rows/total_rows*100:.1f}%)")
        print(f"Invalid rows: {invalid_rows:,} ({invalid_rows/total_rows*100:.1f}%)")
        print(f"Found insertions at {len(insertion_counts):,} unique chromosome-coordinate combinations")
        
    except Exception as e:
        print(f"Error during chunked processing: {e}")
        raise
    
    # Convert accumulated counts to output format
    print("Preparing output table...")
    output_data = []
    
    for (chrom, coord), strand_counts in insertion_counts.items():
        output_data.append({
            'Chr': chrom,
            'Coordinate': coord,
            '+': strand_counts['+'],
            '-': strand_counts['-']
        })
    
    if not output_data:
        print("Warning: No insertion sites found!")
        # Create empty output file with headers
        empty_df = pd.DataFrame(columns=['Chr', 'Coordinate', '+', '-'])
        empty_df.to_csv(output_tsv, sep='\t', index=False)
        return empty_df
    
    # Sort by chromosome and coordinate
    output_df = pd.DataFrame(output_data)
    
    output_df = output_df.sort_values(['Chr', 'Coordinate'])
    
    # Write output
    print(f"Writing {len(output_df):,} rows to {output_tsv}")
    output_df.to_csv(output_tsv, sep='\t', index=False)
    
    # Print summary statistics
    total_plus = output_df['+'].sum()
    total_minus = output_df['-'].sum()
    total_insertions = total_plus + total_minus
    
    print(f"\nSummary:")
    print(f"Total unique insertion sites: {len(output_df):,}")
    print(f"Total insertions on + strand: {total_plus:,}")
    print(f"Total insertions on - strand: {total_minus:,}")
    print(f"Total insertions: {total_insertions:,}")
    
    # Additional statistics
    if len(output_df) > 0:
        sites_with_both_strands = len(output_df[(output_df['+'] > 0) & (output_df['-'] > 0)])
        sites_plus_only = len(output_df[(output_df['+'] > 0) & (output_df['-'] == 0)])
        sites_minus_only = len(output_df[(output_df['+'] == 0) & (output_df['-'] > 0)])
        
        print(f"\nStrand distribution:")
        print(f"Sites with both strands: {sites_with_both_strands:,}")
        print(f"Sites with + strand only: {sites_plus_only:,}")
        print(f"Sites with - strand only: {sites_minus_only:,}")
    
    return output_df

def main():
    parser = argparse.ArgumentParser(description="Process TSV output from parse_bam_to_tsv.py to generate insertion count summary by chromosome, coordinate, and strand using memory-efficient chunking.")
    parser.add_argument("-i", "--input_tsv", required=True, help="Path to the input TSV file from parse_bam_to_tsv.py")
    parser.add_argument("-o", "--output_tsv", required=True, help="Path to the output TSV file with insertion counts")
    parser.add_argument("-c", "--chunk_size", type=int, default=500000, help="Number of rows to process per chunk (default: 500,000)")
    
    args = parser.parse_args()
    
    print("Using pandas chunking for memory-efficient processing")
    print(f"Pandas version: {pd.__version__}")
    print(f"Chunk size: {args.chunk_size:,} rows")
    
    try:
        result_df = process_insertion_counts(args.input_tsv, args.output_tsv, args.chunk_size)
        print(f"\nProcessing completed successfully!")
        print(f"Output written to: {args.output_tsv}")
        
    except Exception as e:
        print(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main()) 