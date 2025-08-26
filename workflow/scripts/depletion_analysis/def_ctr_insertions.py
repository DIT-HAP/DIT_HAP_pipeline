"""
Use selected insertions as control insertions for depletion analysis.

Combine the insertion tables and annotation tables, and select the control insertions.
The control insertions are defined as:
- In the intergenic regions
- Distance to the nearest gene is greater than 500 bp
"""

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from loguru import logger
from pydantic import BaseModel, Field, field_validator


# ======================== Configuration & Models ========================


class ControlInsertionConfig(BaseModel):
    """Configuration for control insertion selection."""

    insertion_table: Path = Field(..., description="Input insertion table")
    annotation_table: Path = Field(..., description="Input annotation table")
    output_file: Path = Field(..., description="Output file")

    @field_validator("insertion_table", "annotation_table")
    def validate_input_exists(cls, v, field):
        """Validate that input files exist."""
        if not v.exists():
            raise ValueError(f"Input file {v} does not exist")
        return v

    @field_validator("output_file")
    def validate_output_path(cls, v: Path) -> Path:
        """Validate output directory exists or create it."""
        output_dir = v.parent
        if not output_dir.exists():
            logger.info(f"Creating output directory: {output_dir}")
            output_dir.mkdir(parents=True, exist_ok=True)
        return v

    class Config:
        frozen = True


# ======================== Logging Setup ========================


def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for control insertion selection."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="{time:HH:mm:ss} | {level: <8} | {message}",
        level=log_level,
        colorize=False,
    )


# ======================== Core Functions ========================
@logger.catch
def load_and_preprocess_data(
    counts_file: Path, annotations_file: Path
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load and preprocess the insertion and annotation tables.

    Args:
        counts_file (Path): Path to the insertion table
        annotations_file (Path): Path to the annotation table

    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]:
            - counts_df: Insertion table with counts
            - insertion_annotations: Annotation table
    """
    counts_df = pd.read_csv(
        counts_file, index_col=[0, 1, 2, 3], header=[0, 1], sep="\t"
    )

    # Remove the rows with any NA value
    counts_df = counts_df.dropna(axis=0, how="any").copy()

    # Load and process annotations
    insertion_annotations = pd.read_csv(
        annotations_file, index_col=[0, 1, 2, 3], sep="\t"
    )

    return counts_df, insertion_annotations


def get_control_insertions(
    counts_df: pd.DataFrame, insertion_annotations: pd.DataFrame
) -> pd.DataFrame:
    ctr_insertions = insertion_annotations.query(
        "Type == 'Intergenic region' and Distance_to_region_start > 500 and Distance_to_region_end > 500"
    )

    ctr_insertions = ctr_insertions[ctr_insertions.index.isin(counts_df.index)].drop_duplicates(keep="first")

    return ctr_insertions


def main():
    parser = argparse.ArgumentParser(
        description="Select control insertions for depletion analysis"
    )
    parser.add_argument(
        "-i",
        "--insertion_table",
        type=Path,
        required=True,
        help="Path to the insertion table",
    )
    parser.add_argument(
        "-a",
        "--annotation_table",
        type=Path,
        required=True,
        help="Path to the annotation table",
    )
    parser.add_argument(
        "-o", "--output_file", type=Path, required=True, help="Path to the output file"
    )
    args = parser.parse_args()

    setup_logging()

    try:
        config = ControlInsertionConfig(
            insertion_table=args.insertion_table,
            annotation_table=args.annotation_table,
            output_file=args.output_file,
        )

        counts_df, insertion_annotations = load_and_preprocess_data(
            config.insertion_table, config.annotation_table
        )

        ctr_insertions = get_control_insertions(counts_df, insertion_annotations)

        ctr_insertions.to_csv(config.output_file, sep="\t", index=True)

        logger.success("Control insertion selection completed successfully!")
        return 0
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        return 1


if __name__ == "__main__":
    main()
