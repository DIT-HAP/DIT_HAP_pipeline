# Set pandas options to display all columns
import pandas as pd
import numpy as np
import argparse
import warnings
warnings.filterwarnings('ignore')

def filter_read_pairs(input_filepath, output_filepath, filter_config):
    """
    Loads read pair data using pandas chunking, filters R1 and R2 independently based on separate configurable quality criteria,
    and saves the filtered data with descriptive statistics.
    
    Args:
        input_filepath (str): Path to input TSV file
        output_filepath (str): Path to output TSV file
        filter_config (dict): Dictionary with separate R1 and R2 filter parameters
    """
    print(f"INFO: Loading data from: {input_filepath}")
    
    # Print filter configuration
    print("\n" + "="*60)
    print("FILTER CONFIGURATION")
    print("="*60)
    print("R1 Filters:")
    print(f"  - MAPQ threshold: {filter_config['r1_mapq_threshold']}")
    print(f"  - NCIGAR value: {filter_config['r1_ncigar_value']}")
    print(f"  - NM threshold: {filter_config['r1_nm_threshold']}")
    print(f"  - Require no supplementary (SA): {filter_config['r1_no_sa']}")
    print(f"  - Require no secondary (XA): {filter_config['r1_no_xa']}")
    print("R2 Filters:")
    print(f"  - MAPQ threshold: {filter_config['r2_mapq_threshold']}")
    print(f"  - NCIGAR value: {filter_config['r2_ncigar_value']}")
    print(f"  - NM threshold: {filter_config['r2_nm_threshold']}")
    print(f"  - Require no supplementary (SA): {filter_config['r2_no_sa']}")
    print(f"  - Require no secondary (XA): {filter_config['r2_no_xa']}")
    print("Pair-level Filters:")
    print(f"  - Require proper pairs: {filter_config['require_proper_pair']}")
    
    try:
        
        # Get chunk size from config
        chunk_size = filter_config.get('chunk_size', 50000)
        print(f"\nINFO: Processing file in chunks of {chunk_size:,} rows...")
        
        # Initialize counters
        total_rows = 0
        filtered_rows = 0
        chunk_count = 0
        first_chunk = True
        
        # Process data in chunks
        print(f"\nINFO: Starting chunked processing...")
        chunk_iterator = pd.read_csv(
            input_filepath, 
            sep='\t',
            na_values=['N/A', 'NA', ''],
            chunksize=chunk_size
        )
        
        # Process each chunk
        for chunk_df in chunk_iterator:
            chunk_count += 1
            chunk_rows_before = len(chunk_df)
            total_rows += chunk_rows_before
            
            if chunk_count % 10 == 0:
                print(f"  Processing chunk {chunk_count}, total rows processed: {total_rows:,}")
            
            # Display original data info for first chunk
            if first_chunk:
                print("\n" + "="*60)
                print("ORIGINAL DATA INFORMATION")
                print("="*60)
                print(f"Columns: {len(chunk_df.columns)}")
                print(f"First chunk size: {len(chunk_df):,} rows")
                
                print("\nColumn Data Types:")
                for col, dtype in chunk_df.dtypes.items():
                    print(f"  {col}: {dtype}")
                
                print("\nSample of first few rows:")
                print(chunk_df.head(3).to_string())
                
                first_chunk = False
            
            # Build filter conditions for this chunk
            filter_mask = pd.Series([True] * len(chunk_df), index=chunk_df.index)
            
            # R1 filters
            if filter_config['r1_mapq_threshold'] is not None:
                filter_mask = filter_mask & (chunk_df['R1_MAPQ'] >= filter_config['r1_mapq_threshold'])
            if filter_config['r1_ncigar_value'] is not None:
                filter_mask = filter_mask & (chunk_df['R1_NCIGAR'] == filter_config['r1_ncigar_value'])
            if filter_config['r1_nm_threshold'] is not None:
                filter_mask = filter_mask & (chunk_df['R1_NM'] <= filter_config['r1_nm_threshold'])
            if filter_config['r1_no_sa']:
                filter_mask = filter_mask & (chunk_df['R1_SA'].isna() | (chunk_df['R1_SA'] == 'N/A'))
            if filter_config['r1_no_xa']:
                filter_mask = filter_mask & (chunk_df['R1_XA'].isna() | (chunk_df['R1_XA'] == 'N/A'))
            
            # R2 filters
            if filter_config['r2_mapq_threshold'] is not None:
                filter_mask = filter_mask & (chunk_df['R2_MAPQ'] >= filter_config['r2_mapq_threshold'])
            if filter_config['r2_ncigar_value'] is not None:
                filter_mask = filter_mask & (chunk_df['R2_NCIGAR'] == filter_config['r2_ncigar_value'])
            if filter_config['r2_nm_threshold'] is not None:
                filter_mask = filter_mask & (chunk_df['R2_NM'] <= filter_config['r2_nm_threshold'])
            if filter_config['r2_no_sa']:
                filter_mask = filter_mask & (chunk_df['R2_SA'].isna() | (chunk_df['R2_SA'] == 'N/A'))
            if filter_config['r2_no_xa']:
                filter_mask = filter_mask & (chunk_df['R2_XA'].isna() | (chunk_df['R2_XA'] == 'N/A'))
            
            # Proper pair filter
            if filter_config['require_proper_pair']:
                filter_mask = filter_mask & (chunk_df['Is_Proper_Pair'].str.capitalize() == 'Yes')
            
            # Apply filters to chunk
            filtered_chunk = chunk_df[filter_mask]
            chunk_filtered_rows = len(filtered_chunk)
            filtered_rows += chunk_filtered_rows
            
            # Write filtered chunk to output (append mode after first chunk)
            if chunk_count == 1:
                # First chunk: write with header
                filtered_chunk.to_csv(output_filepath, sep='\t', index=False, mode='w')
                print(f"\nINFO: Created output file: {output_filepath}")
            else:
                # Subsequent chunks: append without header
                filtered_chunk.to_csv(output_filepath, sep='\t', index=False, mode='a', header=False)
            
            # Print progress for this chunk
            if chunk_count == 1 or chunk_count % 10 == 0:
                retention_rate = chunk_filtered_rows / chunk_rows_before * 100 if chunk_rows_before > 0 else 0
                print(f"    Chunk {chunk_count}: {chunk_filtered_rows:,}/{chunk_rows_before:,} rows retained ({retention_rate:.1f}%)")
        
        print(f"\nINFO: Completed processing {chunk_count} chunks")
        
        # Final statistics
        print("\n" + "="*60)
        print("FILTERING SUMMARY")
        print("="*60)
        print(f"Total chunks processed: {chunk_count}")
        print(f"Original read pairs: {total_rows:,}")
        print(f"Filtered read pairs: {filtered_rows:,}")
        print(f"Removed read pairs: {total_rows - filtered_rows:,}")
        
        if total_rows > 0:
            retention_rate = filtered_rows / total_rows * 100
            print(f"Overall retention rate: {retention_rate:.2f}%")
        
        print(f"Output written to: {output_filepath}")
        
        # Display sample of filtered data
        if filtered_rows > 0:
            print("\n" + "="*60)
            print("FILTERED DATA SAMPLE")
            print("="*60)
            try:
                sample_df = pd.read_csv(output_filepath, sep='\t', nrows=5)
                print("First 5 rows of filtered data:")
                print(sample_df.to_string())
            except Exception as e:
                print(f"Could not read sample of filtered data: {e}")
        
    except FileNotFoundError:
        print(f"ERROR: File not found: {input_filepath}")
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()

def main():
    parser = argparse.ArgumentParser(
        description="Filter read pairs independently for R1 and R2 with separate configurable criteria using pandas chunking for memory efficiency.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("-i", "--input-file", required=True, help="Input TSV file with read pair data")
    parser.add_argument("-o", "--output-file", required=True, help="Output TSV file for filtered data")
    
    # Chunking configuration
    parser.add_argument("-c", "--chunk-size", type=int, default=50000, help="Number of rows to process per chunk (default: 50,000)")
    
    # R1 filter parameters
    parser.add_argument("--r1-mapq-threshold", type=float, default=10,
                       help="Minimum mapping quality score for R1 reads")
    parser.add_argument("--r1-ncigar-value", type=int, default=1,
                       help="Required NCIGAR value for R1 reads (1 = single alignment)")
    parser.add_argument("--r1-nm-threshold", type=int, default=3,
                       help="Maximum number of mismatches allowed for R1 reads")
    parser.add_argument("--r1-require-no-supplementary", action="store_true",
                       help="Require R1 reads to have no supplementary alignments (SA tag)")
    parser.add_argument("--r1-require-no-secondary", action="store_true",
                       help="Require R1 reads to have no secondary alignments (XA tag)")
    
    # R2 filter parameters
    parser.add_argument("--r2-mapq-threshold", type=float, default=10,
                       help="Minimum mapping quality score for R2 reads")
    parser.add_argument("--r2-ncigar-value", type=int, default=5,
                       help="Required NCIGAR value for R2 reads")
    parser.add_argument("--r2-nm-threshold", type=int, default=10,
                       help="Maximum number of mismatches allowed for R2 reads")
    parser.add_argument("--r2-require-no-supplementary", action="store_true",
                       help="Require R2 reads to have no supplementary alignments (SA tag)")
    parser.add_argument("--r2-require-no-secondary", action="store_true",
                       help="Require R2 reads to have no secondary alignments (XA tag)")
    
    # Pair-level filters
    parser.add_argument("--require-proper-pair", action="store_true",
                       help="Require proper pairs")
    
    # Disable specific filters for R1
    parser.add_argument("--r1-disable-mapq", action="store_true",
                       help="Disable MAPQ filtering for R1")
    parser.add_argument("--r1-disable-ncigar", action="store_true",
                       help="Disable NCIGAR filtering for R1")
    parser.add_argument("--r1-disable-nm", action="store_true",
                       help="Disable NM (mismatch) filtering for R1")
    
    # Disable specific filters for R2
    parser.add_argument("--r2-disable-mapq", action="store_true",
                       help="Disable MAPQ filtering for R2")
    parser.add_argument("--r2-disable-ncigar", action="store_true",
                       help="Disable NCIGAR filtering for R2")
    parser.add_argument("--r2-disable-nm", action="store_true",
                       help="Disable NM (mismatch) filtering for R2")
    
    args = parser.parse_args()
    
    print("Using pandas with chunking for stable data processing")
    print(f"Pandas version: {pd.__version__}")
    print(f"Chunk size: {args.chunk_size:,} rows")
    
    # Build filter configuration with separate R1 and R2 parameters
    filter_config = {
        # R1 parameters
        'r1_mapq_threshold': None if args.r1_disable_mapq else args.r1_mapq_threshold,
        'r1_ncigar_value': None if args.r1_disable_ncigar else args.r1_ncigar_value,
        'r1_nm_threshold': None if args.r1_disable_nm else args.r1_nm_threshold,
        'r1_no_sa': args.r1_require_no_supplementary,
        'r1_no_xa': args.r1_require_no_secondary,
        # R2 parameters
        'r2_mapq_threshold': None if args.r2_disable_mapq else args.r2_mapq_threshold,
        'r2_ncigar_value': None if args.r2_disable_ncigar else args.r2_ncigar_value,
        'r2_nm_threshold': None if args.r2_disable_nm else args.r2_nm_threshold,
        'r2_no_sa': args.r2_require_no_supplementary,
        'r2_no_xa': args.r2_require_no_secondary,
        # Pair-level parameters
        'require_proper_pair': args.require_proper_pair,
        # Chunking parameters
        'chunk_size': args.chunk_size
    }
    
    filter_read_pairs(args.input_file, args.output_file, filter_config)

if __name__ == "__main__":
    main()