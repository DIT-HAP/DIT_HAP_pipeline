"""
Enhanced strand-specific insertion merging with Pydantic validation and Loguru logging.

Merges insertion counts from PBL (left primer) and PBR (right primer) reads,
organizing them by chromosome, coordinate, and strand for downstream analysis.
"""

import argparse
import sys
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from loguru import logger
from pydantic import BaseModel, Field, field_validator


# ======================== Configuration & Models ========================

class InsertionData(BaseModel):
    """Model for insertion data from a single source (PBL or PBR)."""
    
    dataframe: pd.DataFrame
    source: str = Field(..., pattern="^(PBL|PBR)$", description="Data source")
    site_count: int = Field(..., ge=0, description="Number of insertion sites")
    
    class Config:
        arbitrary_types_allowed = True
        frozen = True
    
    @field_validator('dataframe')
    def validate_columns(cls, v):
        """Ensure dataframe has required columns."""
        required_cols = ["Chr", "Coordinate", "+", "-"]
        missing_cols = [col for col in required_cols if col not in v.columns]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        return v


class MergeConfig(BaseModel):
    """Configuration for merging insertion data."""
    
    pbl_file: Path = Field(..., description="PBL insertion file path")
    pbr_file: Path = Field(..., description="PBR insertion file path")
    output_file: Path = Field(..., description="Output file path")
    
    @field_validator('pbl_file', 'pbr_file')
    def validate_input_exists(cls, v, field):
        if not v.exists():
            logger.warning(f"{field.name} not found: {v}")
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


class MergeStats(BaseModel):
    """Statistics from merging operation."""
    
    pbl_sites: int = Field(..., ge=0, description="PBL insertion sites")
    pbr_sites: int = Field(..., ge=0, description="PBR insertion sites")
    merged_sites: int = Field(..., ge=0, description="Total merged sites")
    plus_strand_sites: int = Field(..., ge=0, description="Plus strand sites")
    minus_strand_sites: int = Field(..., ge=0, description="Minus strand sites")
    shared_sites: int = Field(..., ge=0, description="Sites in both PBL and PBR")
    pbl_only_sites: int = Field(..., ge=0, description="PBL-only sites")
    pbr_only_sites: int = Field(..., ge=0, description="PBR-only sites")
    
    class Config:
        frozen = True
    
    @property
    def overlap_percentage(self) -> float:
        """Percentage of sites shared between PBL and PBR."""
        total_unique = self.pbl_sites + self.pbr_sites - self.shared_sites
        if total_unique == 0:
            return 0.0
        return (self.shared_sites / total_unique) * 100


# ======================== Logging Setup ========================

def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for insertion merging."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )


# ======================== Core Functions ========================

@logger.catch
def load_insertion_data(file_path: Path, source: str) -> InsertionData:
    """
    Load insertion data from TSV file.
    
    Args:
        file_path: Path to insertion TSV file
        source: Data source identifier (PBL or PBR)
        
    Returns:
        InsertionData object with validated dataframe
    """
    if not file_path.exists():
        logger.warning(f"{source} file not found: {file_path}")
        # Return empty dataframe with proper columns
        empty_df = pd.DataFrame(columns=["Chr", "Coordinate", "+", "-"])
        return InsertionData(
            dataframe=empty_df,
            source=source,
            site_count=0
        )
    
    logger.info(f"Loading {source} file: {file_path}")
    
    try:
        df = pd.read_csv(file_path, sep="\t", header=0)
        site_count = len(df)
        
        logger.success(f"{source} data loaded: {site_count:,} insertion sites")
        
        # Validate and create InsertionData object
        insertion_data = InsertionData(
            dataframe=df,
            source=source,
            site_count=site_count
        )
        
        # Log basic statistics
        if site_count > 0:
            plus_total = df["+"].sum()
            minus_total = df["-"].sum()
            logger.debug(f"{source} strand counts: + = {plus_total:,}, - = {minus_total:,}")
        
        return insertion_data
        
    except Exception as e:
        logger.error(f"Error loading {source} file: {e}")
        raise


@logger.catch
def merge_insertion_data(
    pbl_data: InsertionData,
    pbr_data: InsertionData
) -> Tuple[pd.DataFrame, MergeStats]:
    """
    Merge PBL and PBR insertion data by coordinate and reorganize by strand.
    
    The merging logic:
    - For + strand insertions: PBL(-) reads map to left, PBR(+) reads map to right
    - For - strand insertions: PBL(+) reads map to left, PBR(-) reads map to right
    
    Args:
        pbl_data: PBL insertion data
        pbr_data: PBR insertion data
        
    Returns:
        Tuple of (merged dataframe, merge statistics)
    """
    logger.info("Merging PBL and PBR insertion data...")
    
    pbl_df = pbl_data.dataframe
    pbr_df = pbr_data.dataframe
    
    # Perform outer merge on Chr and Coordinate
    merged_df = pd.merge(
        pbl_df,
        pbr_df,
        how="outer",
        on=["Chr", "Coordinate"],
        suffixes=("_PBL", "_PBR"),
    )
    
    # Fill NaN with 0 for count columns
    merged_df.fillna(0, inplace=True)
    
    # Calculate statistics before reorganization
    shared_mask = (merged_df["+_PBL"] + merged_df["-_PBL"] > 0) & \
                  (merged_df["+_PBR"] + merged_df["-_PBR"] > 0)
    shared_sites = shared_mask.sum()
    
    pbl_only_mask = (merged_df["+_PBL"] + merged_df["-_PBL"] > 0) & \
                    (merged_df["+_PBR"] + merged_df["-_PBR"] == 0)
    pbl_only_sites = pbl_only_mask.sum()
    
    pbr_only_mask = (merged_df["+_PBL"] + merged_df["-_PBL"] == 0) & \
                    (merged_df["+_PBR"] + merged_df["-_PBR"] > 0)
    pbr_only_sites = pbr_only_mask.sum()
    
    logger.debug(f"Merge results: {len(merged_df):,} unique coordinate pairs")
    logger.debug(f"Shared sites: {shared_sites:,}")
    logger.debug(f"PBL-only sites: {pbl_only_sites:,}")
    logger.debug(f"PBR-only sites: {pbr_only_sites:,}")
    
    # Reorganize by strand
    # Plus strand: PBL(-) and PBR(+)
    plus_insertion = merged_df[["Chr", "Coordinate", "-_PBL", "+_PBR"]].copy()
    plus_insertion["Strand"] = "+"
    plus_insertion.rename(
        columns={"-_PBL": "PBL", "+_PBR": "PBR"}, inplace=True
    )
    
    # Minus strand: PBL(+) and PBR(-)
    minus_insertion = merged_df[["Chr", "Coordinate", "+_PBL", "-_PBR"]].copy()
    minus_insertion["Strand"] = "-"
    minus_insertion.rename(
        columns={"+_PBL": "PBL", "-_PBR": "PBR"}, inplace=True
    )
    
    # Combine and format final dataframe
    final_df = pd.concat([plus_insertion, minus_insertion], axis=0)
    final_df = final_df.set_index(["Chr", "Coordinate", "Strand"])
    
    # Convert to integer and sort
    final_df = final_df.astype(int)
    final_df = final_df.sort_index()
    
    # Add total read count
    final_df["Reads"] = final_df["PBL"] + final_df["PBR"]
    
    # Filter out rows with no reads (optional, keeps all coordinates)
    # final_df = final_df[final_df["Reads"] > 0]
    
    # Calculate final statistics
    plus_sites = (final_df.loc[final_df.index.get_level_values("Strand") == "+", "Reads"] > 0).sum()
    minus_sites = (final_df.loc[final_df.index.get_level_values("Strand") == "-", "Reads"] > 0).sum()
    
    stats = MergeStats(
        pbl_sites=pbl_data.site_count,
        pbr_sites=pbr_data.site_count,
        merged_sites=len(final_df[final_df["Reads"] > 0]),
        plus_strand_sites=plus_sites,
        minus_strand_sites=minus_sites,
        shared_sites=shared_sites,
        pbl_only_sites=pbl_only_sites,
        pbr_only_sites=pbr_only_sites
    )
    
    return final_df, stats


@logger.catch
def save_merged_data(
    merged_df: pd.DataFrame,
    output_path: Path,
    stats: MergeStats
) -> None:
    """
    Save merged data to TSV file and display statistics.
    
    Args:
        merged_df: Merged dataframe
        output_path: Output file path
        stats: Merge statistics
    """
    # Save to file
    logger.info(f"Writing merged data to: {output_path}")
    merged_df.to_csv(output_path, sep='\t', index=True, header=True)
    
    # Display comprehensive statistics
    logger.info("=" * 60)
    logger.info("MERGE SUMMARY")
    logger.info("=" * 60)
    
    logger.info("Input files:")
    logger.info(f"  PBL sites: {stats.pbl_sites:,}")
    logger.info(f"  PBR sites: {stats.pbr_sites:,}")
    
    logger.info("\nMerge results:")
    logger.success(f"  Total merged sites: {stats.merged_sites:,}")
    logger.info(f"  Plus strand sites: {stats.plus_strand_sites:,}")
    logger.info(f"  Minus strand sites: {stats.minus_strand_sites:,}")
    
    logger.info("\nSite overlap:")
    logger.info(f"  Shared sites: {stats.shared_sites:,}")
    logger.info(f"  PBL-only sites: {stats.pbl_only_sites:,}")
    logger.info(f"  PBR-only sites: {stats.pbr_only_sites:,}")
    logger.info(f"  Overlap percentage: {stats.overlap_percentage:.1f}%")
    
    # Additional statistics from the dataframe
    if len(merged_df) > 0:
        total_reads = merged_df["Reads"].sum()
        pbl_reads = merged_df["PBL"].sum()
        pbr_reads = merged_df["PBR"].sum()
        
        logger.info("\nRead counts:")
        logger.info(f"  Total reads: {total_reads:,}")
        logger.info(f"  PBL reads: {pbl_reads:,}")
        logger.info(f"  PBR reads: {pbr_reads:,}")
        
        if total_reads > 0:
            logger.info(f"  PBL contribution: {pbl_reads/total_reads*100:.1f}%")
            logger.info(f"  PBR contribution: {pbr_reads/total_reads*100:.1f}%")
    
    logger.success(f"Output saved to: {output_path}")


@logger.catch
def merge_insertions(config: MergeConfig) -> MergeStats:
    """
    Main function to merge PBL and PBR insertions.
    
    Args:
        config: Validated merge configuration
        
    Returns:
        MergeStats object with merge statistics
    """
    # Load PBL and PBR data
    pbl_data = load_insertion_data(config.pbl_file, "PBL")
    pbr_data = load_insertion_data(config.pbr_file, "PBR")
    
    # Merge the data
    merged_df, stats = merge_insertion_data(pbl_data, pbr_data)
    
    # Save results
    save_merged_data(merged_df, config.output_file, stats)
    
    return stats


# ======================== Main Entry Point ========================

def main():
    """Main entry point for the script."""
    setup_logging()
    
    parser = argparse.ArgumentParser(
        description="Merge insertion counts from PBL and PBR reads",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "-il", "--inputPBL",
        type=Path,
        required=True,
        help="Input file with PBL insertion counts"
    )
    parser.add_argument(
        "-ir", "--inputPBR",
        type=Path,
        required=True,
        help="Input file with PBR insertion counts"
    )
    parser.add_argument(
        "-o", "--output",
        type=Path,
        required=True,
        help="Output file for merged insertions"
    )
    
    args = parser.parse_args()
    
    logger.info(f"Pandas version: {pd.__version__}")
    logger.info(f"NumPy version: {np.__version__}")
    
    try:
        # Create and validate configuration
        config = MergeConfig(
            pbl_file=args.inputPBL,
            pbr_file=args.inputPBR,
            output_file=args.output
        )
        
        # Perform merging
        stats = merge_insertions(config)
        
        logger.success("Merging completed successfully!")
        return 0
        
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())