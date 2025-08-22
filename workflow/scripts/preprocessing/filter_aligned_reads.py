"""
Enhanced read pair filtering with Pydantic validation and Loguru logging.

Filters aligned read pairs from BAM-derived TSV files based on configurable
quality criteria for R1 and R2 reads independently, using chunked processing
for memory efficiency.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
import numpy as np
from loguru import logger
from pydantic import BaseModel, Field, field_validator


# ======================== Configuration & Models ========================

class FilterThresholds(BaseModel):
    """Validation model for read filtering thresholds."""
    
    mapq_threshold: Optional[float] = Field(None, ge=0, le=255, description="Minimum MAPQ score")
    ncigar_value: Optional[int] = Field(None, ge=0, description="Required NCIGAR value")
    nm_threshold: Optional[int] = Field(None, ge=0, description="Maximum mismatches allowed")
    no_sa: bool = Field(False, description="Require no supplementary alignments")
    no_xa: bool = Field(False, description="Require no secondary alignments")
    
    class Config:
        frozen = True


class FilterConfig(BaseModel):
    """Complete filtering configuration with validation."""
    
    input_file: Path = Field(..., description="Input TSV file path")
    output_file: Path = Field(..., description="Output TSV file path")
    chunk_size: int = Field(50000, ge=1000, le=10000000, description="Rows per chunk")
    r1_filters: FilterThresholds = Field(default_factory=FilterThresholds)
    r2_filters: FilterThresholds = Field(default_factory=FilterThresholds)
    require_proper_pair: bool = Field(False, description="Require proper pairs")
    
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


class FilteringStats(BaseModel):
    """Statistics from filtering operation."""
    
    total_rows: int = Field(..., ge=0)
    filtered_rows: int = Field(..., ge=0)
    removed_rows: int = Field(..., ge=0)
    retention_rate: float = Field(..., ge=0, le=100)
    chunks_processed: int = Field(..., ge=0)
    
    class Config:
        frozen = True


# ======================== Logging Setup ========================

def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for read filtering."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )


# ======================== Core Filtering Functions ========================

@logger.catch
def build_filter_mask(
    chunk: pd.DataFrame,
    r1_filters: FilterThresholds,
    r2_filters: FilterThresholds,
    require_proper_pair: bool
) -> pd.Series:
    """
    Build boolean mask for filtering read pairs.
    
    Args:
        chunk: DataFrame chunk to filter
        r1_filters: R1 filtering thresholds
        r2_filters: R2 filtering thresholds
        require_proper_pair: Whether to require proper pairs
        
    Returns:
        Boolean series indicating which rows pass filters
    """
    # Initialize mask with all True
    filter_mask = pd.Series([True] * len(chunk), index=chunk.index)
    
    # Apply R1 filters
    if r1_filters.mapq_threshold is not None:
        filter_mask &= (chunk['R1_MAPQ'] >= r1_filters.mapq_threshold)
        
    if r1_filters.ncigar_value is not None:
        filter_mask &= (chunk['R1_NCIGAR'] == r1_filters.ncigar_value)
        
    if r1_filters.nm_threshold is not None:
        filter_mask &= (chunk['R1_NM'] <= r1_filters.nm_threshold)
        
    if r1_filters.no_sa:
        filter_mask &= (chunk['R1_SA'].isna() | (chunk['R1_SA'] == 'N/A'))
        
    if r1_filters.no_xa:
        filter_mask &= (chunk['R1_XA'].isna() | (chunk['R1_XA'] == 'N/A'))
    
    # Apply R2 filters
    if r2_filters.mapq_threshold is not None:
        filter_mask &= (chunk['R2_MAPQ'] >= r2_filters.mapq_threshold)
        
    if r2_filters.ncigar_value is not None:
        filter_mask &= (chunk['R2_NCIGAR'] == r2_filters.ncigar_value)
        
    if r2_filters.nm_threshold is not None:
        filter_mask &= (chunk['R2_NM'] <= r2_filters.nm_threshold)
        
    if r2_filters.no_sa:
        filter_mask &= (chunk['R2_SA'].isna() | (chunk['R2_SA'] == 'N/A'))
        
    if r2_filters.no_xa:
        filter_mask &= (chunk['R2_XA'].isna() | (chunk['R2_XA'] == 'N/A'))
    
    # Apply proper pair filter
    if require_proper_pair:
        filter_mask &= (chunk['Is_Proper_Pair'].str.capitalize() == 'Yes')
    
    return filter_mask


def display_filter_configuration(config: FilterConfig) -> None:
    """Display filtering configuration in formatted output."""
    logger.info("=" * 60)
    logger.info("FILTER CONFIGURATION")
    logger.info("=" * 60)
    
    logger.info("R1 Filters:")
    if config.r1_filters.mapq_threshold is not None:
        logger.info(f"  - MAPQ threshold: {config.r1_filters.mapq_threshold}")
    if config.r1_filters.ncigar_value is not None:
        logger.info(f"  - NCIGAR value: {config.r1_filters.ncigar_value}")
    if config.r1_filters.nm_threshold is not None:
        logger.info(f"  - NM threshold: {config.r1_filters.nm_threshold}")
    logger.info(f"  - Require no supplementary (SA): {config.r1_filters.no_sa}")
    logger.info(f"  - Require no secondary (XA): {config.r1_filters.no_xa}")
    
    logger.info("R2 Filters:")
    if config.r2_filters.mapq_threshold is not None:
        logger.info(f"  - MAPQ threshold: {config.r2_filters.mapq_threshold}")
    if config.r2_filters.ncigar_value is not None:
        logger.info(f"  - NCIGAR value: {config.r2_filters.ncigar_value}")
    if config.r2_filters.nm_threshold is not None:
        logger.info(f"  - NM threshold: {config.r2_filters.nm_threshold}")
    logger.info(f"  - Require no supplementary (SA): {config.r2_filters.no_sa}")
    logger.info(f"  - Require no secondary (XA): {config.r2_filters.no_xa}")
    
    logger.info("Pair-level Filters:")
    logger.info(f"  - Require proper pairs: {config.require_proper_pair}")
    logger.info(f"  - Chunk size: {config.chunk_size:,} rows")


@logger.catch
def process_chunk(
    chunk: pd.DataFrame,
    chunk_num: int,
    config: FilterConfig,
    first_chunk: bool
) -> Tuple[pd.DataFrame, bool]:
    """
    Process a single chunk of data.
    
    Args:
        chunk: DataFrame chunk to process
        chunk_num: Chunk number for logging
        config: Filtering configuration
        first_chunk: Whether this is the first chunk
        
    Returns:
        Tuple of (filtered chunk, is_first_chunk_after_processing)
    """
    chunk_rows_before = len(chunk)
    
    # Display info for first chunk
    if first_chunk:
        logger.info("=" * 60)
        logger.info("ORIGINAL DATA INFORMATION")
        logger.info("=" * 60)
        logger.info(f"Columns: {len(chunk.columns)}")
        logger.info(f"First chunk size: {chunk_rows_before:,} rows")
        
        logger.debug("Column Data Types:")
        for col, dtype in chunk.dtypes.items():
            logger.debug(f"  {col}: {dtype}")
    
    # Build and apply filter mask
    filter_mask = build_filter_mask(
        chunk,
        config.r1_filters,
        config.r2_filters,
        config.require_proper_pair
    )
    
    filtered_chunk = chunk[filter_mask]
    chunk_filtered_rows = len(filtered_chunk)
    
    # Log progress
    if chunk_num == 1 or chunk_num % 10 == 0:
        retention_rate = (chunk_filtered_rows / chunk_rows_before * 100 
                         if chunk_rows_before > 0 else 0)
        logger.info(
            f"Chunk {chunk_num}: {chunk_filtered_rows:,}/{chunk_rows_before:,} "
            f"rows retained ({retention_rate:.1f}%)"
        )
    
    return filtered_chunk, False


@logger.catch
def filter_read_pairs(config: FilterConfig) -> FilteringStats:
    """
    Main function to filter read pairs using chunked processing.
    
    Args:
        config: Validated filtering configuration
        
    Returns:
        FilteringStats object with processing statistics
    """
    logger.info(f"Loading data from: {config.input_file}")
    display_filter_configuration(config)
    
    # Initialize counters
    total_rows = 0
    filtered_rows = 0
    chunk_count = 0
    first_chunk = True
    
    logger.info(f"Processing file in chunks of {config.chunk_size:,} rows...")
    
    try:
        # Create chunk iterator
        chunk_iterator = pd.read_csv(
            config.input_file,
            sep='\t',
            na_values=['N/A', 'NA', ''],
            chunksize=config.chunk_size
        )
        
        # Process each chunk
        for chunk_df in chunk_iterator:
            chunk_count += 1
            total_rows += len(chunk_df)
            
            if chunk_count % 10 == 0:
                logger.info(f"Processing chunk {chunk_count}, total rows: {total_rows:,}")
            
            # Process chunk
            filtered_chunk, first_chunk = process_chunk(
                chunk_df, chunk_count, config, first_chunk
            )
            filtered_rows += len(filtered_chunk)
            
            # Write filtered chunk
            if chunk_count == 1:
                filtered_chunk.to_csv(
                    config.output_file, sep='\t', index=False, mode='w'
                )
                logger.info(f"Created output file: {config.output_file}")
            else:
                filtered_chunk.to_csv(
                    config.output_file, sep='\t', index=False, 
                    mode='a', header=False
                )
        
        logger.info(f"Completed processing {chunk_count} chunks")
        
        # Calculate statistics
        removed_rows = total_rows - filtered_rows
        retention_rate = filtered_rows / total_rows * 100 if total_rows > 0 else 0
        
        stats = FilteringStats(
            total_rows=total_rows,
            filtered_rows=filtered_rows,
            removed_rows=removed_rows,
            retention_rate=retention_rate,
            chunks_processed=chunk_count
        )
        
        # Display summary
        logger.info("=" * 60)
        logger.info("FILTERING SUMMARY")
        logger.info("=" * 60)
        logger.success(f"Total chunks processed: {stats.chunks_processed}")
        logger.info(f"Original read pairs: {stats.total_rows:,}")
        logger.info(f"Filtered read pairs: {stats.filtered_rows:,}")
        logger.info(f"Removed read pairs: {stats.removed_rows:,}")
        logger.success(f"Overall retention rate: {stats.retention_rate:.2f}%")
        logger.info(f"Output written to: {config.output_file}")
        
        # Display sample of filtered data
        if filtered_rows > 0:
            try:
                sample_df = pd.read_csv(config.output_file, sep='\t', nrows=5)
                logger.debug(f"Sample of filtered data (first 5 rows):")
                logger.debug(f"Shape: {sample_df.shape}")
            except Exception as e:
                logger.warning(f"Could not read sample of filtered data: {e}")
        
        return stats
        
    except Exception as e:
        logger.error(f"Error processing file: {e}")
        raise


# ======================== Main Entry Point ========================

def main():
    """Main entry point for the script."""
    setup_logging()
    
    parser = argparse.ArgumentParser(
        description="Filter aligned read pairs with configurable quality criteria",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument(
        "-i", "--input-file", 
        required=True, 
        help="Input TSV file with read pair data"
    )
    parser.add_argument(
        "-o", "--output-file", 
        required=True, 
        help="Output TSV file for filtered data"
    )
    
    # Chunking configuration
    parser.add_argument(
        "-c", "--chunk-size", 
        type=int, 
        default=50000, 
        help="Number of rows to process per chunk"
    )
    
    # R1 filter parameters
    parser.add_argument(
        "--r1-mapq-threshold", 
        type=float, 
        default=10,
        help="Minimum MAPQ score for R1 reads"
    )
    parser.add_argument(
        "--r1-ncigar-value", 
        type=int, 
        default=1,
        help="Required NCIGAR value for R1 reads"
    )
    parser.add_argument(
        "--r1-nm-threshold", 
        type=int, 
        default=3,
        help="Maximum mismatches for R1 reads"
    )
    parser.add_argument(
        "--r1-require-no-supplementary", 
        action="store_true",
        help="Require R1 reads to have no supplementary alignments"
    )
    parser.add_argument(
        "--r1-require-no-secondary", 
        action="store_true",
        help="Require R1 reads to have no secondary alignments"
    )
    
    # R2 filter parameters
    parser.add_argument(
        "--r2-mapq-threshold", 
        type=float, 
        default=10,
        help="Minimum MAPQ score for R2 reads"
    )
    parser.add_argument(
        "--r2-ncigar-value", 
        type=int, 
        default=5,
        help="Required NCIGAR value for R2 reads"
    )
    parser.add_argument(
        "--r2-nm-threshold", 
        type=int, 
        default=10,
        help="Maximum mismatches for R2 reads"
    )
    parser.add_argument(
        "--r2-require-no-supplementary", 
        action="store_true",
        help="Require R2 reads to have no supplementary alignments"
    )
    parser.add_argument(
        "--r2-require-no-secondary", 
        action="store_true",
        help="Require R2 reads to have no secondary alignments"
    )
    
    # Pair-level filters
    parser.add_argument(
        "--require-proper-pair", 
        action="store_true",
        help="Require proper pairs"
    )
    
    # Disable specific filters
    parser.add_argument(
        "--r1-disable-mapq", 
        action="store_true",
        help="Disable MAPQ filtering for R1"
    )
    parser.add_argument(
        "--r1-disable-ncigar", 
        action="store_true",
        help="Disable NCIGAR filtering for R1"
    )
    parser.add_argument(
        "--r1-disable-nm", 
        action="store_true",
        help="Disable NM filtering for R1"
    )
    parser.add_argument(
        "--r2-disable-mapq", 
        action="store_true",
        help="Disable MAPQ filtering for R2"
    )
    parser.add_argument(
        "--r2-disable-ncigar", 
        action="store_true",
        help="Disable NCIGAR filtering for R2"
    )
    parser.add_argument(
        "--r2-disable-nm", 
        action="store_true",
        help="Disable NM filtering for R2"
    )
    
    args = parser.parse_args()
    
    logger.info(f"Pandas version: {pd.__version__}")
    logger.info(f"NumPy version: {np.__version__}")
    
    try:
        # Build filter configuration
        r1_filters = FilterThresholds(
            mapq_threshold=None if args.r1_disable_mapq else args.r1_mapq_threshold,
            ncigar_value=None if args.r1_disable_ncigar else args.r1_ncigar_value,
            nm_threshold=None if args.r1_disable_nm else args.r1_nm_threshold,
            no_sa=args.r1_require_no_supplementary,
            no_xa=args.r1_require_no_secondary
        )
        
        r2_filters = FilterThresholds(
            mapq_threshold=None if args.r2_disable_mapq else args.r2_mapq_threshold,
            ncigar_value=None if args.r2_disable_ncigar else args.r2_ncigar_value,
            nm_threshold=None if args.r2_disable_nm else args.r2_nm_threshold,
            no_sa=args.r2_require_no_supplementary,
            no_xa=args.r2_require_no_secondary
        )
        
        config = FilterConfig(
            input_file=Path(args.input_file),
            output_file=Path(args.output_file),
            chunk_size=args.chunk_size,
            r1_filters=r1_filters,
            r2_filters=r2_filters,
            require_proper_pair=args.require_proper_pair
        )
        
        # Process the file
        stats = filter_read_pairs(config)
        
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()