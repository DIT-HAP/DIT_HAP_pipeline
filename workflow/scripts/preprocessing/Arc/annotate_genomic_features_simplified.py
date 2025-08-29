"""
Simplified genomic feature annotation for transposon insertion sites.

Annotates insertions with genomic features including genes, intergenic regions,
and coding sequences. Calculates distances to start/stop codons and determines
affected amino acid residues.
"""

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from loguru import logger
from pybedtools import BedTool
from pydantic import BaseModel, Field, field_validator


# ======================== Configuration & Models ========================

class GenomicFeature(BaseModel):
    """Model for a genomic feature annotation."""
    
    chromosome: str = Field(..., min_length=1, description="Chromosome")
    coordinate: int = Field(..., ge=0, description="Insertion coordinate")
    strand: str = Field(..., pattern="^[+-]$", description="Strand")
    target: str = Field(..., description="Target sequence")
    feature_type: str = Field(..., description="Feature type")
    gene_name: Optional[str] = Field(None, description="Gene name if applicable")
    distance_to_start: Optional[float] = Field(None, description="Distance to start codon")
    distance_to_stop: Optional[float] = Field(None, description="Distance to stop codon")
    fraction_to_start: Optional[float] = Field(None, ge=0, le=1, description="Fraction to start")
    fraction_to_stop: Optional[float] = Field(None, ge=0, le=1, description="Fraction to stop")
    residue_affected: Optional[int] = Field(None, ge=1, description="Affected amino acid")
    residue_frame: Optional[int] = Field(None, ge=0, le=2, description="Reading frame")
    insertion_direction: Optional[str] = Field(None, description="Forward/Reverse")
    
    class Config:
        frozen = True


class AnnotationConfig(BaseModel):
    """Configuration for genomic annotation."""
    
    input_file: Path = Field(..., description="Input insertion file")
    genome_region_file: Path = Field(..., description="Genome region BED file")
    output_file: Path = Field(..., description="Output annotation file")
    
    @field_validator('input_file', 'genome_region_file')
    def validate_input_exists(cls, v, field):
        if not v.exists():
            raise ValueError(f"{field.name} not found: {v}")
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


class AnnotationStats(BaseModel):
    """Statistics from annotation process."""
    
    total_insertions: int = Field(..., ge=0, description="Total insertion sites")
    annotated_insertions: int = Field(..., ge=0, description="Successfully annotated")
    coding_insertions: int = Field(..., ge=0, description="In coding regions")
    intergenic_insertions: int = Field(..., ge=0, description="In intergenic regions")
    unique_genes: int = Field(..., ge=0, description="Unique genes affected")
    forward_insertions: int = Field(..., ge=0, description="Forward direction")
    reverse_insertions: int = Field(..., ge=0, description="Reverse direction")
    
    class Config:
        frozen = True
    
    @property
    def coding_percentage(self) -> float:
        """Percentage of insertions in coding regions."""
        if self.total_insertions == 0:
            return 0.0
        return (self.coding_insertions / self.total_insertions) * 100


# ======================== Logging Setup ========================

def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for genomic annotation."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False
    )


# ======================== Core Annotation Functions ========================

@logger.catch
def load_insertion_data(input_path: Path) -> pd.DataFrame:
    """Load insertion data from TSV or CSV file."""
    logger.info(f"Loading insertion data from {input_path}")
    
    file_ext = input_path.suffix.lower()
    
    try:
        if file_ext == ".tsv":
            df = pd.read_csv(input_path, header=0, usecols=[0, 1, 2, 3], sep="\t")
        elif file_ext == ".csv":
            df = pd.read_csv(input_path, header=0, usecols=[0, 1, 2, 3])
        else:
            raise ValueError(f"Unsupported file format: {file_ext}")
        
        # Rename coordinate column for BED processing
        df = df.rename(columns={"Coordinate": "End"})
        df.insert(1, "Start", df["End"])
        
        logger.success(f"Loaded {len(df)} insertion sites")
        return df
        
    except Exception as e:
        logger.error(f"Error loading insertion data: {e}")
        raise


@logger.catch
def load_genome_regions(region_path: Path) -> pd.DataFrame:
    """Load genome region annotations."""
    logger.info(f"Loading genome regions from {region_path}")
    
    try:
        df = pd.read_csv(region_path, sep="\t", header=0)
        df = df.rename(columns={"#Chr": "Chr"})
        
        logger.success(f"Loaded {len(df)} genome regions")
        return df
        
    except Exception as e:
        logger.error(f"Error loading genome regions: {e}")
        raise


@logger.catch
def calculate_codon_distances(row: pd.Series) -> pd.Series:
    """Calculate distances to start and stop codons."""
    if row["Type"] == "Intergenic region":
        return pd.Series([np.nan] * 4, index=[
            "Distance_to_start_codon", "Distance_to_stop_codon",
            "Fraction_to_start_codon", "Fraction_to_stop_codon"
        ])
    
    if row["Strand_Interval"] == "+":
        return pd.Series([
            row["Distance_to_region_start"],
            row["Distance_to_region_end"],
            row["Fraction_to_region_start"],
            row["Fraction_to_region_end"]
        ], index=[
            "Distance_to_start_codon", "Distance_to_stop_codon",
            "Fraction_to_start_codon", "Fraction_to_stop_codon"
        ])
    else:
        return pd.Series([
            row["Distance_to_region_end"],
            row["Distance_to_region_start"],
            row["Fraction_to_region_end"],
            row["Fraction_to_region_start"]
        ], index=[
            "Distance_to_start_codon", "Distance_to_stop_codon",
            "Fraction_to_start_codon", "Fraction_to_stop_codon"
        ])


@logger.catch
def calculate_affected_residue(row: pd.Series) -> pd.Series:
    """Calculate affected amino acid residue and reading frame."""
    if row["Type"] in ["Intergenic region", "Non-coding gene"]:
        return pd.Series([np.nan, np.nan], index=["Residue_affected", "Residue_frame"])
    
    cds_base = float(row.get("Accumulated_CDS_bases", 0))
    
    if row["Feature"] == "CDS":
        if row["Strand_Interval"] == "+":
            cds_base += int(row["Coordinate"]) - int(row["Start_Interval"])
        else:
            cds_base += int(row["End_Interval"] - row["Coordinate"])
    
    residue_affected = cds_base // 3 + 1
    residue_frame = cds_base % 3
    
    return pd.Series([residue_affected, residue_frame], index=["Residue_affected", "Residue_frame"])


@logger.catch
def assign_insertion_direction(row: pd.Series) -> str:
    """Determine insertion direction relative to gene."""
    if row["Type"] == "Intergenic region":
        return np.nan
    return "Forward" if row["Strand"] == row["Strand_Interval"] else "Reverse"


@logger.catch
def drop_boundary_duplicates(sub_df: pd.DataFrame) -> pd.DataFrame:
    """Remove duplicated insertions at gene boundaries."""
    # Remove insertions exactly at boundaries
    boundary_mask = (
        (sub_df["Distance_to_start_codon"] == 0) |
        (sub_df["Distance_to_stop_codon"] == 0)
    )
    sub_df = sub_df[~boundary_mask]
    
    # If only one type remains, return it
    if sub_df["Type"].nunique() == 1:
        return sub_df
    
    # Otherwise, prefer coding regions over intergenic
    coding_mask = sub_df["Type"] != "Intergenic region"
    if coding_mask.any():
        return sub_df[coding_mask]
    
    return sub_df


@logger.catch
def annotate_insertions(
    insertions_df: pd.DataFrame,
    regions_df: pd.DataFrame
) -> Tuple[pd.DataFrame, AnnotationStats]:
    """Annotate insertions with genomic features."""
    logger.info("Annotating insertions with genomic features...")
    
    # Convert to BedTool objects and intersect
    insertions_bed = BedTool.from_dataframe(insertions_df)
    regions_bed = BedTool.from_dataframe(regions_df)
    
    # Get column names for merging
    insertion_cols = insertions_df.columns.tolist()
    region_cols = [
        f"{col}_Interval" if col in insertion_cols else col
        for col in regions_df.columns
    ]
    
    # Intersect insertions with regions
    intersected = insertions_bed.intersect(regions_bed, wa=True, wb=True)
    annotated_df = intersected.to_dataframe(names=insertion_cols + region_cols)
    
    # Clean up DataFrame
    annotated_df.drop(columns=["Start"], inplace=True)
    annotated_df.rename(columns={"End": "Coordinate"}, inplace=True)
    annotated_df.replace(r"^\.$", np.nan, inplace=True, regex=True)
    
    # Calculate distances to region boundaries
    annotated_df["Distance_to_region_start"] = annotated_df["Coordinate"] - annotated_df["ParentalRegion_start"]
    annotated_df["Distance_to_region_end"] = annotated_df["ParentalRegion_end"] - annotated_df["Coordinate"]
    annotated_df["Fraction_to_region_start"] = annotated_df["Distance_to_region_start"] / annotated_df["ParentalRegion_length"]
    annotated_df["Fraction_to_region_end"] = annotated_df["Distance_to_region_end"] / annotated_df["ParentalRegion_length"]
    
    # Apply all calculations at once
    calculations = [
        annotated_df.apply(calculate_codon_distances, axis=1),
        annotated_df.apply(calculate_affected_residue, axis=1),
        pd.DataFrame({"Insertion_direction": annotated_df.apply(assign_insertion_direction, axis=1)})
    ]
    
    # Concatenate all results
    for calc_result in calculations:
        annotated_df = pd.concat([annotated_df, calc_result], axis=1)
    
    # Remove boundary duplicates
    annotated_df = (
        annotated_df.groupby(["Chr", "Coordinate", "Strand"])
        .apply(drop_boundary_duplicates)
        .reset_index(drop=True)
    )
    
    # Calculate statistics
    stats = AnnotationStats(
        total_insertions=len(insertions_df),
        annotated_insertions=len(annotated_df),
        coding_insertions=len(annotated_df[annotated_df["Type"] != "Intergenic region"]),
        intergenic_insertions=len(annotated_df[annotated_df["Type"] == "Intergenic region"]),
        unique_genes=annotated_df["GeneName"].nunique() if "GeneName" in annotated_df.columns else 0,
        forward_insertions=len(annotated_df[annotated_df["Insertion_direction"] == "Forward"]),
        reverse_insertions=len(annotated_df[annotated_df["Insertion_direction"] == "Reverse"])
    )
    
    return annotated_df, stats


@logger.catch
def save_annotations(
    annotated_df: pd.DataFrame,
    output_path: Path,
    stats: AnnotationStats
) -> None:
    """Save annotated insertions to file."""
    logger.info(f"Saving annotations to {output_path}")
    
    annotated_df.to_csv(
        output_path,
        index=False,
        header=True,
        float_format="%.3f",
        sep="\t"
    )
    
    # Display statistics
    logger.info("=" * 60)
    logger.info("ANNOTATION SUMMARY")
    logger.info("=" * 60)
    logger.info(f"Total insertions: {stats.total_insertions:,}")
    logger.success(f"Annotated insertions: {stats.annotated_insertions:,}")
    
    logger.info("\nFeature distribution:")
    logger.info(f"  Coding regions: {stats.coding_insertions:,} ({stats.coding_percentage:.1f}%)")
    logger.info(f"  Intergenic regions: {stats.intergenic_insertions:,}")
    
    if stats.unique_genes > 0:
        logger.info(f"\nGene impact:")
        logger.info(f"  Unique genes affected: {stats.unique_genes:,}")
        logger.info(f"  Forward insertions: {stats.forward_insertions:,}")
        logger.info(f"  Reverse insertions: {stats.reverse_insertions:,}")
    
    logger.success(f"Annotations saved to {output_path}")


# ======================== Main Entry Point ========================

def main():
    """Main entry point for the script."""
    setup_logging()
    
    parser = argparse.ArgumentParser(
        description="Annotate insertion sites with genomic features",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "-i", "--input",
        dest="input",
        type=Path,
        required=True,
        help="Input insertion file (TSV or CSV)"
    )
    parser.add_argument(
        "-g", "--genome-region",
        dest="genome_region",
        type=Path,
        required=True,
        help="Genome region BED file"
    )
    parser.add_argument(
        "-o", "--output",
        dest="output",
        type=Path,
        required=True,
        help="Output annotation file"
    )
    
    args = parser.parse_args()
    
    try:
        # Create and validate configuration
        config = AnnotationConfig(
            input_file=args.input,
            genome_region_file=args.genome_region,
            output_file=args.output
        )
        
        # Load data
        insertions_df = load_insertion_data(config.input_file)
        regions_df = load_genome_regions(config.genome_region_file)
        
        # Annotate insertions
        annotated_df, stats = annotate_insertions(insertions_df, regions_df)
        
        # Save results
        save_annotations(annotated_df, config.output_file, stats)
        
        logger.success("Annotation completed successfully!")
        return 0
        
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())