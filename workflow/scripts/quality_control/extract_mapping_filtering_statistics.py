
# %%
"""
Extract filtering statistics from read pair filtering log files.
This script parses log files to extract PBL and PBR filtering statistics
including original counts, filtered counts, and retention rates.
"""

from pathlib import Path
import sys
import pandas as pd
import re
import argparse

def extract_filtering_statistics(log_files, output_file):
    """
    Extract filtering statistics from all log files in the specified folder.
    
    Args:
        log_folder (str): Path to the folder containing log files
        output_file (str): Path to save the output statistics CSV
    """
    print(f"*** Found {len(log_files)} log files with filtering statistics")
    
    filtering_statistics = {}
    
    for log_file in log_files:
        log_file_stem = Path(log_file).stem
        print(f"**** Processing: {log_file_stem}")
        
        # Initialize dictionary for this log file
        filtering_statistics[log_file_stem] = {}
        
        try:
            with open(log_file, "r") as f:
                content = f.read()
            
            # Find all filtering summary sections
            # Pattern to match the entire filtering summary block
            summary_pattern = re.compile(
                r"============================================================\s*\n"
                r"FILTERING SUMMARY\s*\n"
                r"============================================================\s*\n"
                r"Total chunks processed: (\d+)\s*\n"
                r"Original read pairs: ([\d,]+)\s*\n"
                r"Filtered read pairs: ([\d,]+)\s*\n"
                r"Removed read pairs: ([\d,]+)\s*\n"
                r"Overall retention rate: ([\d.]+)%\s*\n"
                r"Output written to: (.+?\.(?:PBL|PBR)\.filtered\.tsv)"
            )
            
            matches = summary_pattern.findall(content)
            
            if matches:
                print(f"*** Found {len(matches)} filtering summary sections")
                
                for match in matches:
                    chunks_processed = int(match[0])
                    original_pairs = int(match[1].replace(',', ''))
                    filtered_pairs = int(match[2].replace(',', ''))
                    removed_pairs = int(match[3].replace(',', ''))
                    retention_rate = float(match[4])/100
                    output_path = match[5]
                    
                    # Determine if this is PBL or PBR based on output path
                    if ".PBL.filtered.tsv" in output_path:
                        suffix = "PBL"
                    elif ".PBR.filtered.tsv" in output_path:
                        suffix = "PBR"
                    else:
                        print(f"*** Warning: Could not determine PBL/PBR from output path: {output_path}")
                        continue
                    
                    # Store statistics
                    filtering_statistics[log_file_stem].update({
                        f"chunks_processed_{suffix}": chunks_processed,
                        f"original_read_pairs_{suffix}": original_pairs,
                        f"filtered_read_pairs_{suffix}": filtered_pairs,
                        f"removed_read_pairs_{suffix}": removed_pairs,
                        f"retention_rate_{suffix}": retention_rate,
                        f"output_file_{suffix}": output_path
                    })
                    
                    print(f"***   {suffix}: {original_pairs:,} -> {filtered_pairs:,} ({retention_rate*100:.2f}% retained)")
                
            else:
                print(f"*** No filtering summary sections found in: {log_file}")
                
        except Exception as e:
            print(f"*** Error processing {log_file}: {str(e)}")
            continue
    
    # Convert to DataFrame
    if filtering_statistics:
        filtering_statistics_df = pd.DataFrame(filtering_statistics).T
        
        # Sort columns for better readability
        pbl_cols = [col for col in filtering_statistics_df.columns if col.endswith('_PBL')]
        pbr_cols = [col for col in filtering_statistics_df.columns if col.endswith('_PBR')]
        all_cols = sorted(pbl_cols) + sorted(pbr_cols)
        
        filtering_statistics_df = filtering_statistics_df.reindex(columns=all_cols)
        
        # Calculate totals
        if 'original_read_pairs_PBL' in filtering_statistics_df.columns and 'original_read_pairs_PBR' in filtering_statistics_df.columns:
            filtering_statistics_df['total_original_pairs'] = (
                filtering_statistics_df['original_read_pairs_PBL'] + 
                filtering_statistics_df['original_read_pairs_PBR']
            )
        
        if 'filtered_read_pairs_PBL' in filtering_statistics_df.columns and 'filtered_read_pairs_PBR' in filtering_statistics_df.columns:
            filtering_statistics_df['total_filtered_pairs'] = (
                filtering_statistics_df['filtered_read_pairs_PBL'] + 
                filtering_statistics_df['filtered_read_pairs_PBR']
            )
            
        # Calculate overall retention rate
        filtering_statistics_df['overall_retention_rate'] = (
            filtering_statistics_df['total_filtered_pairs'] / 
            filtering_statistics_df['total_original_pairs']
        ).round(4)
        
        # Save to file
        filtering_statistics_df = filtering_statistics_df.rename_axis("Sample", axis=0).sort_index()
        filtering_statistics_df.to_csv(output_file, sep="\t", index=True, float_format="%.2f")
        print(f"\n*** Statistics extracted and saved to: {output_file}")
        print(f"*** Processed {len(filtering_statistics_df)} samples")
        
        # Display summary
        print("\n*** SUMMARY STATISTICS ***")
        print(filtering_statistics_df[['total_original_pairs', 'total_filtered_pairs', 'overall_retention_rate']].describe())
        
        return filtering_statistics_df
    else:
        print("*** No statistics extracted from any log files")
        return None

def main():
    parser = argparse.ArgumentParser(description="Extract filtering statistics from read pair filtering log files")
    parser.add_argument("-i", "--input", required=True, type=Path, nargs="+", help="Path to the folder containing log files")
    parser.add_argument("-o", "--output", required=True, type=Path, help="Path to save the output statistics CSV")
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    input_files = [str(input_file) for input_file in args.input]
    
    # Extract statistics
    stats_df = extract_filtering_statistics(input_files, args.output)
    
    if stats_df is not None:
        print(f"\n*** Extraction completed successfully!")
        print(f"*** Output saved to: {args.output}")
    else:
        print(f"\n*** Extraction failed - no data found")
        sys.exit(1)

if __name__ == "__main__":
    main()

# %%
