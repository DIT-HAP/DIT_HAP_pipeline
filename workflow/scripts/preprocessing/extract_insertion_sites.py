"""
Enhanced insertion site extraction with Pydantic validation and Loguru logging.

Processes TSV output from BAM parsing to identify and count transposon insertion
sites based on read alignment coordinates and strand orientation.
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from loguru import logger
from pydantic import BaseModel, Field, field_validator


# ======================== Configuration & Models ========================

class InsertionSite(BaseModel):
    """Model for a single insertion site."""
    
    chromosome: str = Field(..., min_length=1, description="Chromosome name")
    coordinate: int = Field(..., ge=0, description="0-based coordinate")
    plus_count: int = Field(0, ge=0, description="Count on plus strand")
    minus_count: int = Field(0, ge=0, description="Count on minus strand")
    
    class Config:
        frozen = True
    
    @property
    def total_count(self) -> int:
        """Total insertions at this site."""
        return self.plus_count + self.minus_count
    
    @property
    def has_both_strands(self) -> bool:
        """Whether site has insertions on both strands."""
        return self.plus_count > 0 and self.minus_count > 0


class ProcessingConfig(BaseModel):
    """Configuration for insertion extraction."""
    
    input_file: Path = Field(..., description="Input TSV file path")
    output_file: Path = Field(..., description="Output TSV file path")
    chunk_size: int = Field(500000, ge=10000, le=5000000, description="Rows per chunk")
        
    @field_validator('input_file')
    def validate_input_exists(cls, v):
        if not v.exists():
            raise ValueError(f"Input file not found: {v}")
        return v
    
    @field_validator('output_file')
    def validate_output_dir(cls, v):
        output_dir = v.parent
        if not output_dir.exists():
            logger.info(f"Creating output directory: {output_dir}")
            output_dir.mkdir(parents=True, exist_ok=True)
        return v
    
    class Config:
        frozen = True


class ExtractionStats(BaseModel):
    """Statistics from insertion extraction."""
    
    total_rows: int = Field(..., ge=0, description="Total input rows")
    valid_rows: int = Field(..., ge=0, description="Valid rows processed")
    invalid_rows: int = Field(..., ge=0, description="Invalid/skipped rows")
    unique_sites: int = Field(..., ge=0, description="Unique insertion sites")
    total_plus_insertions: int = Field(..., ge=0, description="Total + strand insertions")
    total_minus_insertions: int = Field(..., ge=0, description="Total - strand insertions")
    sites_both_strands: int = Field(..., ge=0, description="Sites with both strands")
    sites_plus_only: int = Field(..., ge=0, description="Sites with + strand only")
    sites_minus_only: int = Field(..., ge=0, description="Sites with - strand only")
    chunks_processed: int = Field(..., ge=0, description="Number of chunks processed")
    
    @property
    def total_insertions(self) -> int:
        """Total number of insertions."""
        return self.total_plus_insertions + self.total_minus_insertions
    
    @property
    def validity_rate(self) -> float:
        """Percentage of valid rows."""
        if self.total_rows == 0:
            return 0.0
        return (self.valid_rows / self.total_rows) * 100
    
    class Config:
        frozen = True


# ======================== Logging Setup ========================

def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for insertion extraction."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )


# ======================== Core Processing Functions ========================

@logger.catch
def calculate_insertion_coordinate(row: pd.Series) -> Optional[int]:
    """
    Calculate insertion coordinate based on strand orientation.
    
    For + strand: TTAA[Genome] - use position after TTAA (ref_start + 4)
    For - strand: [Genome]TTAA - use position at end (ref_end)
    
    Args:
        row: Pandas series with read alignment data
        
    Returns:
        Insertion coordinate or None if invalid
    """
    try:
        strand = row['R1_Strand']
        
        if strand == '+':
            # For + strand, TTAA[Genome] use the position after TTAA
            return int(row['R1_Ref_Start']) + 4
        elif strand == '-':
            # For - strand, [Genome]TTAA use the end position
            return int(row['R1_Ref_End'])
        else:
            return None
            
    except (ValueError, TypeError, KeyError):
        return None


@logger.catch
def process_chunk(
    chunk: pd.DataFrame,
    chunk_num: int,
    required_columns: List[str]
) -> Tuple[Dict[Tuple[str, int], Dict[str, int]], int, int]:
    """
    Process a single chunk of data to extract insertion sites.
    
    Args:
        chunk: DataFrame chunk to process
        chunk_num: Chunk number for logging
        required_columns: Required column names
        
    Returns:
        Tuple of (insertion_counts, valid_rows, invalid_rows)
    """
    chunk_rows = len(chunk)
    
    # Validate columns on first chunk
    if chunk_num == 1:
        missing_columns = [col for col in required_columns if col not in chunk.columns]
        if missing_columns:
            raise ValueError(f"Missing required columns: {missing_columns}")
        logger.debug(f"Found all required columns: {required_columns}")
    
    # Select only required columns
    chunk = chunk[required_columns].copy()
    
    # Filter valid rows
    valid_mask = (
        chunk['R1_Strand'].notna() & 
        chunk['R1_Chrom'].notna() & 
        chunk['R1_Ref_Start'].notna() & 
        chunk['R1_Ref_End'].notna() &
        chunk['R1_Strand'].isin(['+', '-'])
    )
    
    valid_chunk = chunk[valid_mask].copy()
    valid_rows = len(valid_chunk)
    invalid_rows = chunk_rows - valid_rows
    
    if chunk_num == 1 or chunk_num % 10 == 0:
        retention_rate = (valid_rows / chunk_rows * 100) if chunk_rows > 0 else 0
        logger.info(
            f"Chunk {chunk_num}: {valid_rows:,}/{chunk_rows:,} valid rows "
            f"({retention_rate:.1f}%)"
        )
    
    # Calculate insertion coordinates
    insertion_counts = defaultdict(lambda: {'+': 0, '-': 0})
    
    if valid_rows > 0:
        valid_chunk['Insertion_Coordinate'] = valid_chunk.apply(
            calculate_insertion_coordinate, axis=1
        )
        
        # Remove rows where coordinate calculation failed
        valid_chunk = valid_chunk[valid_chunk['Insertion_Coordinate'].notna()]
        
        # Count insertions
        for _, row in valid_chunk.iterrows():
            try:
                chrom = row['R1_Chrom']
                coord = int(row['Insertion_Coordinate'])
                strand = row['R1_Strand']
                
                key = (chrom, coord)
                insertion_counts[key][strand] += 1
                
            except (ValueError, TypeError, KeyError):
                continue
    
    return dict(insertion_counts), valid_rows, invalid_rows


@logger.catch
def extract_insertion_sites(config: ProcessingConfig) -> ExtractionStats:
    """
    Main function to extract insertion sites from aligned reads.
    
    Args:
        config: Validated processing configuration
        
    Returns:
        ExtractionStats object with processing statistics
    """
    logger.info(f"Processing TSV file: {config.input_file}")
    logger.info(f"Chunk size: {config.chunk_size:,} rows")
    
    # Initialize counters
    insertion_counts = defaultdict(lambda: {'+': 0, '-': 0})
    total_rows = 0
    total_valid_rows = 0
    total_invalid_rows = 0
    chunk_count = 0
    
    # Required columns for processing
    required_columns = ['R1_Strand', 'R1_Chrom', 'R1_Ref_Start', 'R1_Ref_End']
    
    logger.info("Starting chunked processing...")
    
    try:
        # Create chunk iterator
        chunk_iterator = pd.read_csv(
            config.input_file,
            sep='\t',
            chunksize=config.chunk_size,
            na_values=['N/A', 'NA', '']
        )
        
        # Process each chunk
        for chunk_df in chunk_iterator:
            chunk_count += 1
            total_rows += len(chunk_df)
            
            if chunk_count % 10 == 0:
                logger.info(f"Processing chunk {chunk_count}, total rows: {total_rows:,}")
            
            # Process chunk
            chunk_counts, valid_rows, invalid_rows = process_chunk(
                chunk_df, chunk_count, required_columns
            )
            
            # Accumulate counts
            for key, strand_counts in chunk_counts.items():
                insertion_counts[key]['+'] += strand_counts['+']
                insertion_counts[key]['-'] += strand_counts['-']
            
            total_valid_rows += valid_rows
            total_invalid_rows += invalid_rows
        
        logger.success(f"Completed processing {chunk_count} chunks")
        
    except Exception as e:
        logger.error(f"Error during chunked processing: {e}")
        raise
    
    # Convert to output format
    logger.info("Preparing output table...")
    
    if not insertion_counts:
        logger.warning("No insertion sites found!")
        # Create empty output with headers
        empty_df = pd.DataFrame(columns=['Chr', 'Coordinate', '+', '-'])
        empty_df.to_csv(config.output_file, sep='\t', index=False)
        
        return ExtractionStats(
            total_rows=total_rows,
            valid_rows=total_valid_rows,
            invalid_rows=total_invalid_rows,
            unique_sites=0,
            total_plus_insertions=0,
            total_minus_insertions=0,
            sites_both_strands=0,
            sites_plus_only=0,
            sites_minus_only=0,
            chunks_processed=chunk_count
        )
    
    # Create insertion site objects
    insertion_sites = []
    for (chrom, coord), strand_counts in insertion_counts.items():
        site = InsertionSite(
            chromosome=chrom,
            coordinate=coord,
            plus_count=strand_counts['+'],
            minus_count=strand_counts['-']
        )
        insertion_sites.append(site)
    
    # Convert to DataFrame
    output_data = [
        {
            'Chr': site.chromosome,
            'Coordinate': site.coordinate,
            '+': site.plus_count,
            '-': site.minus_count
        }
        for site in insertion_sites
    ]
    
    output_df = pd.DataFrame(output_data)
    output_df = output_df.sort_values(['Chr', 'Coordinate'])
    
    # Write output
    logger.info(f"Writing {len(output_df):,} insertion sites to {config.output_file}")
    output_df.to_csv(config.output_file, sep='\t', index=False)
    
    # Calculate statistics
    total_plus = sum(site.plus_count for site in insertion_sites)
    total_minus = sum(site.minus_count for site in insertion_sites)
    sites_both = sum(1 for site in insertion_sites if site.has_both_strands)
    sites_plus_only = sum(1 for site in insertion_sites 
                          if site.plus_count > 0 and site.minus_count == 0)
    sites_minus_only = sum(1 for site in insertion_sites 
                           if site.minus_count > 0 and site.plus_count == 0)
    
    stats = ExtractionStats(
        total_rows=total_rows,
        valid_rows=total_valid_rows,
        invalid_rows=total_invalid_rows,
        unique_sites=len(insertion_sites),
        total_plus_insertions=total_plus,
        total_minus_insertions=total_minus,
        sites_both_strands=sites_both,
        sites_plus_only=sites_plus_only,
        sites_minus_only=sites_minus_only,
        chunks_processed=chunk_count
    )
    
    # Display summary
    logger.info("=" * 60)
    logger.info("EXTRACTION SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total rows processed: {stats.total_rows:,}")
    logger.info(f"Valid rows: {stats.valid_rows:,} ({stats.validity_rate:.1f}%)")
    logger.info(f"Invalid rows: {stats.invalid_rows:,}")
    logger.success(f"Unique insertion sites: {stats.unique_sites:,}")
    
    logger.info("\nInsertion counts:")
    logger.info(f"  Total insertions: {stats.total_insertions:,}")
    logger.info(f"  Plus strand: {stats.total_plus_insertions:,}")
    logger.info(f"  Minus strand: {stats.total_minus_insertions:,}")
    
    logger.info("\nStrand distribution:")
    logger.info(f"  Sites with both strands: {stats.sites_both_strands:,}")
    logger.info(f"  Sites with + only: {stats.sites_plus_only:,}")
    logger.info(f"  Sites with - only: {stats.sites_minus_only:,}")
    
    logger.success(f"Output saved to: {config.output_file}")
    
    return stats


# ======================== Main Entry Point ========================

def main():
    """Main entry point for the script."""
    setup_logging()
    
    parser = argparse.ArgumentParser(
        description="Extract insertion sites from aligned read TSV files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "-i", "--input_tsv",
        required=True,
        help="Input TSV file from BAM parsing"
    )
    parser.add_argument(
        "-o", "--output_tsv",
        required=True,
        help="Output TSV file with insertion counts"
    )
    parser.add_argument(
        "-c", "--chunk_size",
        type=int,
        default=500000,
        help="Number of rows to process per chunk"
    )
    
    args = parser.parse_args()
    
    logger.info(f"Pandas version: {pd.__version__}")
    
    try:
        # Create and validate configuration
        config = ProcessingConfig(
            input_file=Path(args.input_tsv),
            output_file=Path(args.output_tsv),
            chunk_size=args.chunk_size
        )
        
        # Extract insertion sites
        stats = extract_insertion_sites(config)
        
        logger.success("Processing completed successfully!")
        return 0
        
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())