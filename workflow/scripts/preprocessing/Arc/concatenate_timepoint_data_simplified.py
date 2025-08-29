"""
Simplified timepoint concatenation with cleaner logic.

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
        if hasattr(info, 'data') and 'input_files' in info.data:
            input_files = info.data['input_files']
            if len(v) != len(input_files):
                raise ValueError(
                    f"Timepoints ({len(v)}) must match input files ({len(input_files)})"
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
        v.parent.mkdir(parents=True, exist_ok=True)
        return v


class ConcatenationStats(BaseModel):
    """Statistics from concatenation operation."""
    
    num_timepoints: int = Field(..., ge=1)
    num_insertions: int = Field(..., ge=0)
    num_chromosomes: int = Field(..., ge=0)
    total_pbl_reads: int = Field(..., ge=0)
    total_pbr_reads: int = Field(..., ge=0)
    total_reads: int = Field(..., ge=0)


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

def load_reference_genome(genome_path: Path) -> Dict:
    """Load reference genome sequences."""
    logger.info(f"Loading reference genome: {genome_path}")
    
    try:
        ref_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
        logger.success(f"Loaded {len(ref_dict)} sequences from genome")
        
        # Log first few chromosome names
        chroms = list(ref_dict.keys())[:5]
        if len(ref_dict) > 5:
            logger.debug(f"Chromosomes: {chroms} ... and {len(ref_dict) - 5} more")
        else:
            logger.debug(f"Chromosomes: {chroms}")
        
        return ref_dict
        
    except Exception as e:
        logger.error(f"Error loading genome: {e}")
        raise


def extract_target_sequence(chrom: str, coordinate: int, ref_dict: Dict) -> str:
    """Extract TTAA target sequence from reference."""
    try:
        if chrom not in ref_dict:
            logger.warning(f"Chromosome {chrom} not found in reference")
            return "NNNN"
        
        # Extract 4bp target sequence (TTAA)
        start = coordinate - 4  # Convert 1-based to 0-based
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


def load_timepoint_files(input_files: List[Path], timepoints: List[str]) -> List[pd.DataFrame]:
    """Load all timepoint files into dataframes."""
    dfs = []
    for file, tp in zip(input_files, timepoints):
        logger.debug(f"Loading timepoint {tp} from {file}")
        
        try:
            df = pd.read_csv(file, header=0, index_col=[0, 1, 2], sep="\t")
            logger.debug(f"  Loaded {len(df)} insertions for {tp}")
            dfs.append(df)
        except Exception as e:
            logger.error(f"Error loading {file}: {e}")
            raise
    
    return dfs


def calculate_read_totals(concatenated: pd.DataFrame) -> Tuple[int, int, int]:
    """Calculate total read counts from concatenated data."""
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
    
    return total_pbl_reads, total_pbr_reads, total_reads


def concatenate_timepoints(config: TimepointConfig, ref_dict: Dict) -> Tuple[pd.DataFrame, ConcatenationStats]:
    """Concatenate insertion data across timepoints."""
    logger.info(f"Concatenating {len(config.timepoints)} timepoints")
    
    # Load all timepoint files
    dfs = load_timepoint_files(config.input_files, config.timepoints)
    
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
    target_sequences = [
        extract_target_sequence(idx[0], idx[1], ref_dict)
        for idx in concatenated.index
    ]
    
    # Add target as new index level
    concatenated = concatenated.set_index(
        pd.Series(target_sequences, name="Target", index=concatenated.index),
        append=True
    )
    
    # Log target distribution
    unique_targets = concatenated.index.get_level_values("Target").unique()
    logger.info(f"Found {len(unique_targets)} unique target sequences")
    
    target_counts = concatenated.index.get_level_values("Target").value_counts()
    if "TTAA" in target_counts.index:
        ttaa_fraction = target_counts["TTAA"] / len(concatenated) * 100
        logger.info(f"TTAA targets: {target_counts['TTAA']} ({ttaa_fraction:.1f}%)")
    
    # Calculate read totals
    total_pbl_reads, total_pbr_reads, total_reads = calculate_read_totals(concatenated)
    
    # Create statistics
    stats = ConcatenationStats(
        num_timepoints=len(config.timepoints),
        num_insertions=len(concatenated),
        num_chromosomes=concatenated.index.get_level_values("Chr").nunique(),
        total_pbl_reads=total_pbl_reads,
        total_pbr_reads=total_pbr_reads,
        total_reads=total_reads
    )
    
    return concatenated, stats


def save_data_by_type(data: pd.DataFrame, output_path: Path, data_type: str) -> None:
    """Save data of specific type to output file."""
    if data_type in data.columns.get_level_values(1):
        filtered_data = data.xs(data_type, level=1, axis=1)
        filtered_data.fillna(0).astype(int).to_csv(output_path, index=True, sep="\t")
        logger.success(f"Saved {data_type} data to {output_path}")
    else:
        logger.warning(f"No {data_type} data found in concatenated results")


def save_concatenated_data(concatenated: pd.DataFrame, config: TimepointConfig, stats: ConcatenationStats) -> None:
    """Save concatenated data to output files."""
    logger.info("Saving concatenated data...")
    
    # Save each data type
    save_data_by_type(concatenated, config.output_pbl, "PBL")
    save_data_by_type(concatenated, config.output_pbr, "PBR")
    save_data_by_type(concatenated, config.output_reads, "Reads")
    
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