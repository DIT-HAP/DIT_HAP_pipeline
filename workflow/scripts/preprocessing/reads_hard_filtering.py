"""
Enhanced hard filtering for insertion reads with Pydantic validation and Loguru logging.

Filters insertion reads based on minimum read count thresholds at initial timepoints,
supporting reproducible scientific analysis with comprehensive logging and validation.
"""

import argparse
import sys
from pathlib import Path
from typing import Tuple

import pandas as pd
from loguru import logger
from pydantic import BaseModel, Field, field_validator


# ======================== Data Models ========================

class FilteringConfig(BaseModel):
    """Configuration for hard filtering operation."""
    
    input_file: Path = Field(..., description="Input TSV file with insertion reads")
    output_file: Path = Field(..., description="Output TSV file for filtered reads")
    initial_timepoint: str = Field(..., description="Initial timepoint column name for filtering")
    cutoff_threshold: int = Field(..., ge=0, description="Minimum read count threshold")
    
    @field_validator('input_file')
    def validate_input_file(cls, v: Path) -> Path:
        """Validate that input file exists and is readable."""
        if not v.exists():
            raise ValueError(f"Input file not found: {v}")
        if not v.suffix.lower() in ['.tsv', '.txt']:
            logger.warning(f"Input file may not be TSV format: {v.suffix}")
        return v
    
    @field_validator('output_file')
    def validate_output_path(cls, v: Path) -> Path:
        """Validate output directory exists or create it."""
        output_dir = v.parent
        if not output_dir.exists():
            logger.info(f"Creating output directory: {output_dir}")
            output_dir.mkdir(parents=True, exist_ok=True)
        return v
    
    @field_validator('initial_timepoint')
    def validate_timepoint(cls, v: str) -> str:
        """Validate timepoint is not empty."""
        if not v or not v.strip():
            raise ValueError("Initial timepoint cannot be empty")
        return v.strip()
    
    class Config:
        frozen = True


class FilteringStats(BaseModel):
    """Statistics from the filtering operation."""
    
    total_insertions: int = Field(..., ge=0, description="Total insertions before filtering")
    retained_insertions: int = Field(..., ge=0, description="Insertions retained after filtering")
    removed_insertions: int = Field(..., ge=0, description="Insertions removed by filtering")
    retention_rate: float = Field(..., ge=0, le=100, description="Percentage of insertions retained")
    samples_processed: int = Field(..., ge=0, description="Number of samples processed")
    
    class Config:
        frozen = True


class SampleFilteringResult(BaseModel):
    """Result for filtering a single sample."""
    
    sample_name: str = Field(..., description="Sample identifier")
    total_insertions: int = Field(..., ge=0)
    retained_insertions: int = Field(..., ge=0)
    retention_rate: float = Field(..., ge=0, le=100)


# ======================== Logging Configuration ========================

def setup_logging(log_level: str = "INFO") -> None:
    """Configure Loguru logging for the filtering script."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )


# ======================== Core Filtering Functions ========================

@logger.catch
def load_insertion_data(input_file: Path) -> pd.DataFrame:
    """
    Load insertion reads from TSV file.
    
    Args:
        input_file: Path to input TSV file
        
    Returns:
        DataFrame with multi-index (Chr, Coordinate, Strand, Target)
    """
    logger.info(f"Loading insertion data from: {input_file}")
    
    try:
        df = pd.read_csv(
            input_file, 
            sep="\t", 
            index_col=[0, 1, 2, 3], 
            header=[0, 1]
        )
        logger.success(f"Loaded {df.shape[0]:,} insertions with {df.shape[1]} columns")
        
        # Validate structure
        if not isinstance(df.columns, pd.MultiIndex):
            raise ValueError("Expected MultiIndex columns with (Sample, Timepoint) structure")
        
        if len(df.columns.levels) != 2:
            raise ValueError("Expected columns to have 2 levels: Sample and Timepoint")
            
        return df
        
    except Exception as e:
        logger.error(f"Failed to load data: {e}")
        raise


@logger.catch
def validate_timepoint_exists(df: pd.DataFrame, timepoint: str) -> None:
    """
    Validate that the specified timepoint exists in the data.
    
    Args:
        df: DataFrame with MultiIndex columns
        timepoint: Timepoint to validate
    """
    available_timepoints = df.columns.get_level_values(1).unique()
    if timepoint not in available_timepoints:
        logger.error(f"Timepoint '{timepoint}' not found in data")
        logger.error(f"Available timepoints: {list(available_timepoints)}")
        raise ValueError(f"Invalid timepoint: {timepoint}")
    
    logger.debug(f"Validated timepoint '{timepoint}' exists in data")


@logger.catch
def filter_sample_insertions(
    sample_data: pd.DataFrame, 
    initial_timepoint: str, 
    cutoff: int
) -> SampleFilteringResult:
    """
    Filter insertions for a single sample based on read count threshold.
    
    Args:
        sample_data: DataFrame for a single sample
        initial_timepoint: Column name for initial timepoint
        cutoff: Minimum read count threshold
        
    Returns:
        SampleFilteringResult with filtering statistics
    """
    sample_name = sample_data.columns.get_level_values(0)[0]
    total_insertions = len(sample_data)
    
    # Apply filtering
    mask = sample_data[(sample_name, initial_timepoint)] >= cutoff
    filtered_data = sample_data[mask]
    retained_insertions = len(filtered_data)
    
    retention_rate = (retained_insertions / total_insertions * 100 
                     if total_insertions > 0 else 0)
    
    logger.debug(
        f"Sample {sample_name}: {retained_insertions:,}/{total_insertions:,} "
        f"insertions retained ({retention_rate:.2f}%)"
    )
    
    return SampleFilteringResult(
        sample_name=sample_name,
        total_insertions=total_insertions,
        retained_insertions=retained_insertions,
        retention_rate=retention_rate
    )


@logger.catch
def apply_hard_filtering(
    df: pd.DataFrame, 
    config: FilteringConfig
) -> Tuple[pd.DataFrame, FilteringStats]:
    """
    Apply hard filtering across all samples.
    
    Args:
        df: Input DataFrame with insertion data
        config: Filtering configuration
        
    Returns:
        Tuple of (filtered DataFrame, FilteringStats)
    """
    logger.info("Starting hard filtering process...")
    
    # Validate timepoint exists
    validate_timepoint_exists(df, config.initial_timepoint)
    
    # Display initial data info
    logger.info("=" * 60)
    logger.info("INITIAL DATA SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total insertions: {df.shape[0]:,}")
    logger.info(f"Total samples: {len(df.columns.get_level_values(0).unique())}")
    logger.info(f"Timepoints: {list(df.columns.get_level_values(1).unique())}")
    logger.info(f"Initial timepoint: {config.initial_timepoint}")
    logger.info(f"Cutoff threshold: {config.cutoff_threshold}")
    
    # Process each sample
    filtered_samples = {}
    sample_results = []
    
    for sample_name, sample_data in df.groupby(level="Sample", axis=1):
        logger.debug(f"Processing sample: {sample_name}")
        
        result = filter_sample_insertions(
            sample_data, 
            config.initial_timepoint, 
            config.cutoff_threshold
        )
        sample_results.append(result)
        
        if result.retained_insertions > 0:
            filtered_samples[sample_name] = sample_data[
                sample_data[(sample_name, config.initial_timepoint)] >= config.cutoff_threshold
            ]
    
    # Combine filtered samples
    if not filtered_samples:
        logger.warning("No samples retained any insertions after filtering")
        filtered_df = pd.DataFrame(columns=df.columns)
    else:
        filtered_df = pd.concat(filtered_samples.values(), axis=1)
    
    # Calculate overall statistics
    total_insertions = sum(r.total_insertions for r in sample_results)
    retained_insertions = sum(r.retained_insertions for r in sample_results)
    removed_insertions = total_insertions - retained_insertions
    retention_rate = (retained_insertions / total_insertions * 100 
                     if total_insertions > 0 else 0)
    
    stats = FilteringStats(
        total_insertions=total_insertions,
        retained_insertions=retained_insertions,
        removed_insertions=removed_insertions,
        retention_rate=retention_rate,
        samples_processed=len(sample_results)
    )
    
    return filtered_df, stats


# ======================== Main Entry Point ========================

def create_argument_parser() -> argparse.ArgumentParser:
    """Create argument parser with detailed help."""
    parser = argparse.ArgumentParser(
        description="Filter insertion reads by hard filtering based on read count thresholds",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python reads_hard_filtering.py -i raw_reads.tsv -o filtered_reads.tsv -itp 0h -c 5
  python reads_hard_filtering.py -i counts.tsv -o filtered.tsv --init-timepoint YES0 --cutoff 10
        """
    )
    
    parser.add_argument(
        "-i", "--input",
        type=Path,
        required=True,
        help="Input TSV file with insertion reads (multi-index: Chr, Coordinate, Strand, Target)")
    
    parser.add_argument(
        "-o", "--output",
        type=Path,
        required=True,
        help="Output TSV file for filtered insertion reads")
    
    parser.add_argument(
        "-itp", "--init-timepoint",
        type=str,
        required=True,
        help="Initial timepoint column name for filtering (e.g., '0h', 'YES0')")
    
    parser.add_argument(
        "-c", "--cutoff",
        type=int,
        required=True,
        help="Minimum read count threshold at initial timepoint")
    
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level (default: INFO)")
    
    return parser


@logger.catch
def main():
    """Main entry point for the hard filtering script."""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.log_level)
    
    logger.info("=" * 70)
    logger.info("HARD FILTERING FOR INSERTION READS")
    logger.info("=" * 70)
    
    try:
        # Create configuration
        config = FilteringConfig(
            input_file=args.input,
            output_file=args.output,
            initial_timepoint=args.init_timepoint,
            cutoff_threshold=args.cutoff
        )
        
        # Display configuration
        logger.info("Configuration:")
        logger.info(f"  Input file: {config.input_file}")
        logger.info(f"  Output file: {config.output_file}")
        logger.info(f"  Initial timepoint: {config.initial_timepoint}")
        logger.info(f"  Cutoff threshold: {config.cutoff_threshold}")
        
        # Load data
        df = load_insertion_data(config.input_file)
        
        # Apply filtering
        filtered_df, stats = apply_hard_filtering(df, config)
        
        # Save results
        logger.info("Saving filtered results...")
        filtered_df.to_csv(config.output_file, sep="\t", header=True, index=True)
        logger.success(f"Results saved to: {config.output_file}")
        
        # Display summary
        logger.info("=" * 70)
        logger.info("FILTERING SUMMARY")
        logger.info("=" * 70)
        logger.info(f"Total insertions: {stats.total_insertions:,}")
        logger.info(f"Retained insertions: {stats.retained_insertions:,}")
        logger.info(f"Removed insertions: {stats.removed_insertions:,}")
        logger.success(f"Retention rate: {stats.retention_rate:.2f}%")
        logger.info(f"Samples processed: {stats.samples_processed}")
        
        # Log completion
        logger.success("Hard filtering completed successfully!")
        
    except Exception as e:
        logger.error(f"Hard filtering failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()