"""
Enhanced timepoint concatenation with Pydantic validation and Loguru logging.

Concatenates insertion count data across multiple timepoints for the same sample,
adding target sequence information from the reference genome.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from Bio import SeqIO
from loguru import logger
from pydantic import BaseModel, Field, field_validator


# ======================== Configuration & Models ========================

class TimepointConfig(BaseModel):
    """Configuration for timepoint concatenation."""
    
    sample_name: str = Field(..., min_length=1, description="Sample identifier")
    input_files: List[Path] = Field(..., min_items=1, description="Input insertion files")
    timepoints: List[str] = Field(..., min_items=1, description="Timepoint names")
    genome_file: Path = Field(..., description="Reference genome FASTA file")
    output_pbl: Path = Field(..., description="Output PBL file path")
    output_pbr: Path = Field(..., description="Output PBR file path")
    output_reads: Path = Field(..., description="Output reads file path")
    
    @field_validator('input_files')
    @classmethod
    def validate_input_files(cls, v):
        missing = [f for f in v if not f.exists()]
        if missing:
            raise ValueError(f"Input files not found: {missing}")
        return v
    
    @field_validator('timepoints')
    @classmethod
    def validate_timepoint_count(cls, v, info):
        # In Pydantic v2, we need to access data through info.data
        if hasattr(info, 'data') and 'input_files' in info.data:
            input_files = info.data['input_files']
            if len(v) != len(input_files):
                raise ValueError(
                    f"Number of timepoints ({len(v)}) must match number of input files "
                    f"({len(input_files)})"
                )
        return v
    
    @field_validator('genome_file')
    @classmethod
    def validate_genome_exists(cls, v):
        if not v.exists():
            raise ValueError(f"Genome file not found: {v}")
        return v
    
    @field_validator('output_pbl', 'output_pbr', 'output_reads')
    @classmethod
    def validate_output_dir(cls, v):
        output_dir = v.parent
        if not output_dir.exists():
            logger.info(f"Creating output directory: {output_dir}")
            output_dir.mkdir(parents=True, exist_ok=True)
        return v
    
    class Config:
        frozen = True


class ConcatenationStats(BaseModel):
    """Statistics from concatenation operation."""
    
    num_timepoints: int = Field(..., ge=1, description="Number of timepoints")
    num_insertions: int = Field(..., ge=0, description="Total unique insertions")
    num_chromosomes: int = Field(..., ge=0, description="Number of chromosomes")
    total_pbl_reads: int = Field(..., ge=0, description="Total PBL reads")
    total_pbr_reads: int = Field(..., ge=0, description="Total PBR reads")
    total_reads: int = Field(..., ge=0, description="Total combined reads")
    
    class Config:
        frozen = True


# ======================== Logging Setup ========================

def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for timepoint concatenation."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )


# ======================== Core Functions ========================

@logger.catch
def load_reference_genome(genome_path: Path) -> Dict:
    """
    Load reference genome sequences.
    
    Args:
        genome_path: Path to genome FASTA file
        
    Returns:
        Dictionary of sequence records
    """
    logger.info(f"Loading reference genome: {genome_path}")
    
    try:
        ref_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
        logger.success(f"Loaded {len(ref_dict)} sequences from genome")
        
        # Log chromosome names
        chroms = list(ref_dict.keys())[:5]
        if len(ref_dict) > 5:
            logger.debug(f"Chromosomes: {chroms} ... and {len(ref_dict) - 5} more")
        else:
            logger.debug(f"Chromosomes: {chroms}")
        
        return ref_dict
        
    except Exception as e:
        logger.error(f"Error loading genome: {e}")
        raise


@logger.catch
def extract_target_sequence(
    chrom: str,
    coordinate: int,
    ref_dict: Dict
) -> str:
    """
    Extract TTAA target sequence from reference.
    
    Args:
        chrom: Chromosome name
        coordinate: Insertion coordinate
        ref_dict: Reference genome dictionary
        
    Returns:
        4bp target sequence (TTAA or variant)
    """
    try:
        if chrom not in ref_dict:
            logger.warning(f"Chromosome {chrom} not found in reference")
            return "NNNN"
        
        # Extract 4bp target sequence (TTAA)
        # Coordinate is 1-based, convert to 0-based
        start = coordinate - 4
        end = coordinate
        
        if start < 0:
            logger.warning(f"Coordinate {coordinate} too close to chromosome start")
            return "NNNN"
        
        seq = str(ref_dict[chrom].seq[start:end])
        
        if len(seq) != 4:
            logger.warning(f"Could not extract 4bp at {chrom}:{coordinate}")
            return "NNNN"
        
        return seq.upper()
        
    except Exception as e:
        logger.warning(f"Error extracting target at {chrom}:{coordinate}: {e}")
        return "NNNN"


@logger.catch
def concatenate_timepoints(
    config: TimepointConfig,
    ref_dict: Dict
) -> Tuple[pd.DataFrame, ConcatenationStats]:
    """
    Concatenate insertion data across timepoints.
    
    Args:
        config: Validated configuration
        ref_dict: Reference genome dictionary
        
    Returns:
        Tuple of (concatenated dataframe, statistics)
    """
    logger.info(f"Concatenating {len(config.timepoints)} timepoints")
    
    # Load all timepoint files
    dfs = []
    for i, (file, tp) in enumerate(zip(config.input_files, config.timepoints)):
        logger.debug(f"Loading timepoint {tp} from {file}")
        
        try:
            df = pd.read_csv(file, header=0, index_col=[0, 1, 2], sep="\t")
            logger.debug(f"  Loaded {len(df)} insertions for {tp}")
            dfs.append(df)
        except Exception as e:
            logger.error(f"Error loading {file}: {e}")
            raise
    
    # Concatenate all timepoints
    logger.info("Concatenating dataframes...")
    concatenated = pd.concat(
        dfs,
        axis=1,
        keys=config.timepoints,
        join="outer"
    )
    
    # Sort by timepoint names and coordinates
    concatenated = concatenated.sort_index(
        level=0, axis=1, key=lambda x: x.str.lower()
    ).sort_index(axis=0)
    
    logger.success(f"Concatenated {len(concatenated)} unique insertion sites")
    
    # Add target sequence information
    logger.info("Adding target sequences...")
    target_sequences = []
    
    for idx in concatenated.index:
        chrom = idx[0]
        coordinate = idx[1]
        target = extract_target_sequence(chrom, coordinate, ref_dict)
        target_sequences.append(target)
    
    # Add target as new index level
    concatenated = concatenated.set_index(
        pd.Series(target_sequences, name="Target", index=concatenated.index),
        append=True
    )
    
    # Count unique targets
    unique_targets = concatenated.index.get_level_values("Target").unique()
    logger.info(f"Found {len(unique_targets)} unique target sequences")
    
    # Log target distribution if interesting
    target_counts = concatenated.index.get_level_values("Target").value_counts()
    if "TTAA" in target_counts.index:
        ttaa_fraction = target_counts["TTAA"] / len(concatenated) * 100
        logger.info(f"TTAA targets: {target_counts['TTAA']} ({ttaa_fraction:.1f}%)")
    
    # Calculate read totals before creating frozen stats object
    total_pbl_reads = 0
    total_pbr_reads = 0
    total_reads = 0
    
    if "PBL" in concatenated.columns.get_level_values(1):
        pbl_data = concatenated.xs("PBL", level=1, axis=1)
        total_pbl_reads = int(pbl_data.sum().sum())
    
    if "PBR" in concatenated.columns.get_level_values(1):
        pbr_data = concatenated.xs("PBR", level=1, axis=1)
        total_pbr_reads = int(pbr_data.sum().sum())
    
    if "Reads" in concatenated.columns.get_level_values(1):
        reads_data = concatenated.xs("Reads", level=1, axis=1)
        total_reads = int(reads_data.sum().sum())
    
    # Create statistics object with calculated values (can't modify after creation due to frozen=True)
    stats = ConcatenationStats(
        num_timepoints=len(config.timepoints),
        num_insertions=len(concatenated),
        num_chromosomes=concatenated.index.get_level_values("Chr").nunique(),
        total_pbl_reads=total_pbl_reads,
        total_pbr_reads=total_pbr_reads,
        total_reads=total_reads
    )
    
    return concatenated, stats


@logger.catch
def save_concatenated_data(
    concatenated: pd.DataFrame,
    config: TimepointConfig,
    stats: ConcatenationStats
) -> None:
    """
    Save concatenated data to output files.
    
    Args:
        concatenated: Concatenated dataframe
        config: Configuration
        stats: Concatenation statistics
    """
    logger.info("Saving concatenated data...")
    
    # Save PBL data
    if "PBL" in concatenated.columns.get_level_values(1):
        pbl_data = concatenated.xs("PBL", level=1, axis=1)
        pbl_data.fillna(0).astype(int).to_csv(
            config.output_pbl, index=True, sep="\t"
        )
        logger.success(f"Saved PBL data to {config.output_pbl}")
    else:
        logger.warning("No PBL data found in concatenated results")
    
    # Save PBR data
    if "PBR" in concatenated.columns.get_level_values(1):
        pbr_data = concatenated.xs("PBR", level=1, axis=1)
        pbr_data.fillna(0).astype(int).to_csv(
            config.output_pbr, index=True, sep="\t"
        )
        logger.success(f"Saved PBR data to {config.output_pbr}")
    else:
        logger.warning("No PBR data found in concatenated results")
    
    # Save Reads data
    if "Reads" in concatenated.columns.get_level_values(1):
        reads_data = concatenated.xs("Reads", level=1, axis=1)
        reads_data.fillna(0).astype(int).to_csv(
            config.output_reads, index=True, sep="\t"
        )
        logger.success(f"Saved Reads data to {config.output_reads}")
    else:
        logger.warning("No Reads data found in concatenated results")
    
    # Display summary statistics
    logger.info("=" * 60)
    logger.info("CONCATENATION SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Sample: {config.sample_name}")
    logger.info(f"Timepoints: {', '.join(config.timepoints)}")
    logger.success(f"Unique insertion sites: {stats.num_insertions:,}")
    logger.info(f"Chromosomes: {stats.num_chromosomes}")
    
    if stats.total_reads > 0:
        logger.info(f"\nRead counts:")
        logger.info(f"  PBL reads: {stats.total_pbl_reads:,}")
        logger.info(f"  PBR reads: {stats.total_pbr_reads:,}")
        logger.info(f"  Total reads: {stats.total_reads:,}")


# ======================== Main Entry Point ========================

def main():
    """Main entry point for the script."""
    setup_logging()
    
    parser = argparse.ArgumentParser(
        description="Concatenate insertion data across timepoints",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("-s", "--sample", required=True, help="Sample name")
    parser.add_argument(
        "-i", "--input-files",
        dest="input_files",
        nargs="+",
        type=Path,
        required=True,
        help="Input insertion count files"
    )
    parser.add_argument(
        "-tp", "--timepoints",
        dest="timepoints",
        nargs="+",
        type=str,
        required=True,
        help="Timepoint names"
    )
    parser.add_argument(
        "-g", "--genome",
        type=Path,
        required=True,
        help="Reference genome FASTA file"
    )
    parser.add_argument(
        "-ol", "--outputPBL",
        type=Path,
        required=True,
        help="Output PBL file"
    )
    parser.add_argument(
        "-or", "--outputPBR",
        type=Path,
        required=True,
        help="Output PBR file"
    )
    parser.add_argument(
        "-o", "--outputReads",
        type=Path,
        required=True,
        help="Output Reads file"
    )
    
    args = parser.parse_args()
    
    logger.info(f"Pandas version: {pd.__version__}")
    
    try:
        # Create and validate configuration
        config = TimepointConfig(
            sample_name=args.sample,
            input_files=args.input_files,
            timepoints=args.timepoints,
            genome_file=args.genome,
            output_pbl=args.outputPBL,
            output_pbr=args.outputPBR,
            output_reads=args.outputReads
        )
        
        # Load reference genome
        ref_dict = load_reference_genome(config.genome_file)
        
        # Concatenate timepoints
        concatenated, stats = concatenate_timepoints(config, ref_dict)
        
        # Save results
        save_concatenated_data(concatenated, config, stats)
        
        logger.success("Concatenation completed successfully!")
        return 0
        
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())