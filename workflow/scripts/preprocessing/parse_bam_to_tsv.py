"""
Enhanced BAM to TSV parser with Pydantic validation and Loguru logging.

Extracts comprehensive summary for read pairs from BAM/SAM files,
including tags, strand, NCIGAR, chrom, pos (0-based), ref_start (0-based),
ref_end (0-based), and flag information.
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pysam
from loguru import logger
from pydantic import BaseModel, Field, validator


# ======================== Configuration & Models ========================

class ReadInfo(BaseModel):
    """Validated read information from pysam AlignedSegment."""
    
    mapq: int = Field(default=0, ge=0, le=255, description="Mapping quality")
    length: Optional[int] = Field(default=None, ge=0, description="Query alignment length")
    cigar: str = Field(default="N/A", description="CIGAR string")
    strand: str = Field(default="N/A", regex="^[+-]$|^N/A$", description="Read strand")
    n_cigar: int = Field(default=0, ge=0, description="Number of CIGAR operations")
    chrom: str = Field(default="N/A", description="Reference chromosome")
    pos: str = Field(default="N/A", description="0-based position")
    ref_start: str = Field(default="N/A", description="0-based reference start")
    ref_end: str = Field(default="N/A", description="0-based reference end (exclusive)")
    flag: int = Field(default=0, ge=0, description="SAM flag")
    tags: Dict[str, str] = Field(default_factory=dict, description="Read tags")
    
    class Config:
        frozen = True
        
    @validator('strand')
    def validate_strand(cls, v):
        if v not in ['+', '-', 'N/A']:
            logger.warning(f"Unexpected strand value: {v}")
        return v


class ReadPairInfo(BaseModel):
    """Validated information for a read pair."""
    
    qname: str = Field(..., min_length=1, description="Query name")
    read1: Optional[ReadInfo] = None
    read2: Optional[ReadInfo] = None
    is_proper_pair: str = Field(default="N/A", description="Proper pair status")
    
    class Config:
        frozen = True
    
    @validator('is_proper_pair')
    def validate_proper_pair(cls, v):
        valid_values = ["Yes", "No", "Single_End_Or_Flag_Issue", "N/A"]
        if v not in valid_values:
            logger.warning(f"Unexpected proper pair value: {v}")
        return v


class ProcessingConfig(BaseModel):
    """Configuration for BAM processing."""
    
    input_bam: Path = Field(..., description="Input BAM/SAM file path")
    output_file: Path = Field(..., description="Output TSV file path")
    threads: int = Field(default=4, ge=1, le=32, description="Number of threads")
    tag_list: List[str] = Field(
        default=['AS', 'MC', 'MD', 'MQ', 'NM', 'SA', 'XA', 'XS'],
        description="Tags to extract"
    )
    
    @validator('input_bam')
    def validate_input_exists(cls, v):
        if not v.exists():
            raise ValueError(f"Input file not found: {v}")
        return v
    
    @validator('output_file')
    def validate_output_dir(cls, v):
        output_dir = v.parent
        if not output_dir.exists():
            logger.info(f"Creating output directory: {output_dir}")
            output_dir.mkdir(parents=True, exist_ok=True)
        return v
    
    class Config:
        frozen = True


# ======================== Core Functions ========================

def setup_logging(log_level: str = "INFO") -> None:
    """Configure loguru for BAM processing."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | {message}",
        level=log_level,
        colorize=True
    )


@logger.catch
def extract_read_info(read: Optional[pysam.AlignedSegment], tag_list: List[str]) -> ReadInfo:
    """
    Extract and validate information from a pysam AlignedSegment.
    
    Args:
        read: AlignedSegment object or None
        tag_list: List of tags to extract
        
    Returns:
        Validated ReadInfo object
    """
    if read is None:
        return ReadInfo(tags={tag: "N/A" for tag in tag_list})
    
    # Extract basic information
    strand = "-" if read.is_reverse else "+"
    
    # Calculate NCIGAR
    n_cigar = 0
    if read.cigartuples:
        n_cigar = len(read.cigartuples)
    elif read.cigarstring:
        n_cigar = 0
    
    # Extract tags
    read_tags = dict(read.get_tags())
    formatted_tags = {}
    
    for tag_name in tag_list:
        value = read_tags.get(tag_name, "N/A")
        if value is True:
            formatted_tags[tag_name] = "True"
        elif value is False:
            formatted_tags[tag_name] = "False"
        elif value is None:
            formatted_tags[tag_name] = "N/A"
        else:
            formatted_tags[tag_name] = str(value)
    
    # Create validated ReadInfo
    return ReadInfo(
        mapq=read.mapping_quality if read.mapping_quality is not None else 0,
        length=read.query_alignment_length,
        cigar=read.cigarstring if read.cigarstring else "N/A",
        strand=strand,
        n_cigar=n_cigar,
        chrom=read.reference_name if read.reference_name is not None else "N/A",
        pos=str(read.reference_start) if read.reference_start is not None and read.reference_start != -1 else "N/A",
        ref_start=str(read.reference_start) if read.reference_start is not None else "N/A",
        ref_end=str(read.reference_end) if read.reference_end is not None else "N/A",
        flag=read.flag if read.flag is not None else 0,
        tags=formatted_tags
    )


@logger.catch
def determine_proper_pair_status(
    read1: Optional[pysam.AlignedSegment],
    read2: Optional[pysam.AlignedSegment]
) -> str:
    """
    Determine if reads form a proper pair.
    
    Args:
        read1: First read of the pair
        read2: Second read of the pair
        
    Returns:
        Proper pair status string
    """
    if read1 and read1.is_paired:
        return "Yes" if read1.is_proper_pair else "No"
    elif read2 and read2.is_paired:
        return "Yes" if read2.is_proper_pair else "No"
    elif (read1 and not read1.is_paired) or (read2 and not read2.is_paired):
        return "Single_End_Or_Flag_Issue"
    return "N/A"


@logger.catch
def process_read_pair(
    qname: str,
    read1: Optional[pysam.AlignedSegment],
    read2: Optional[pysam.AlignedSegment],
    tag_list: List[str]
) -> ReadPairInfo:
    """
    Process a read pair and return validated information.
    
    Args:
        qname: Query name
        read1: First read
        read2: Second read
        tag_list: Tags to extract
        
    Returns:
        Validated ReadPairInfo object
    """
    is_proper_pair = determine_proper_pair_status(read1, read2)
    
    r1_info = extract_read_info(read1, tag_list)
    r2_info = extract_read_info(read2, tag_list)
    
    return ReadPairInfo(
        qname=qname,
        read1=r1_info,
        read2=r2_info,
        is_proper_pair=is_proper_pair
    )


def format_output_line(pair_info: ReadPairInfo, tag_list: List[str]) -> List[str]:
    """
    Format ReadPairInfo into output line.
    
    Args:
        pair_info: Validated read pair information
        tag_list: Tags to include in output
        
    Returns:
        List of strings for TSV output
    """
    output_line = [pair_info.qname]
    
    # Read 1 information
    r1 = pair_info.read1 or ReadInfo(tags={tag: "N/A" for tag in tag_list})
    output_line.extend([
        str(r1.mapq),
        str(r1.length) if r1.length is not None else "N/A",
        r1.cigar,
        r1.strand,
        str(r1.n_cigar),
        r1.chrom,
        r1.pos,
        r1.ref_start,
        r1.ref_end,
        str(r1.flag)
    ])
    for tag in tag_list:
        output_line.append(r1.tags.get(tag, "N/A"))
    
    # Read 2 information
    r2 = pair_info.read2 or ReadInfo(tags={tag: "N/A" for tag in tag_list})
    output_line.extend([
        str(r2.mapq),
        str(r2.length) if r2.length is not None else "N/A",
        r2.cigar,
        r2.strand,
        str(r2.n_cigar),
        r2.chrom,
        r2.pos,
        r2.ref_start,
        r2.ref_end,
        str(r2.flag)
    ])
    for tag in tag_list:
        output_line.append(r2.tags.get(tag, "N/A"))
    
    output_line.append(pair_info.is_proper_pair)
    return output_line


def build_header(tag_list: List[str]) -> List[str]:
    """
    Build the header for the output TSV file.
    
    Args:
        tag_list: Tags to include in header
        
    Returns:
        List of header field names
    """
    header_fields = [
        "QueryName",
        "R1_MAPQ", "R1_LEN", "R1_CIGAR", "R1_Strand", "R1_NCIGAR",
        "R1_Chrom", "R1_Pos", "R1_Ref_Start", "R1_Ref_End", "R1_Flag"
    ]
    
    for tag_name in tag_list:
        header_fields.append(f"R1_{tag_name}")
    
    header_fields.extend([
        "R2_MAPQ", "R2_LEN", "R2_CIGAR", "R2_Strand", "R2_NCIGAR",
        "R2_Chrom", "R2_Pos", "R2_Ref_Start", "R2_Ref_End", "R2_FLAG"
    ])
    
    for tag_name in tag_list:
        header_fields.append(f"R2_{tag_name}")
    
    header_fields.append("Is_Proper_Pair")
    return header_fields


# ======================== Main Processing ========================

@logger.catch
def process_bam_file(config: ProcessingConfig) -> None:
    """
    Main function to process BAM/SAM file and write TSV output.
    
    Args:
        config: Validated processing configuration
    """
    logger.info(f"Starting BAM processing with {config.threads} threads")
    logger.info(f"Input: {config.input_bam}")
    logger.info(f"Output: {config.output_file}")
    logger.info(f"Extracting tags: {', '.join(config.tag_list)}")
    
    # Ensure tags are sorted for consistent output
    sorted_tags = sorted(config.tag_list)
    
    # Build header
    header_fields = build_header(sorted_tags)
    
    # Process BAM file
    with open(config.output_file, "w") as outfile:
        # Write header
        outfile.write("\t".join(header_fields) + "\n")
        
        # Initialize tracking variables
        current_qname = None
        current_r1 = None
        current_r2 = None
        processed_qname_count = 0
        read_count = 0
        
        try:
            # Open BAM/SAM file
            mode = "rb" if str(config.input_bam).endswith(".bam") else "r"
            samfile = pysam.AlignmentFile(
                str(config.input_bam),
                mode,
                threads=config.threads
            )
            
            logger.info("Processing reads in streaming mode (requires qname-sorted input)")
            
            for read in samfile:
                read_count += 1
                
                if read_count % 2000000 == 0:
                    logger.info(f"Processed {read_count // 1000000}M alignments")
                
                # Skip unmapped, secondary, and supplementary alignments
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                
                qname = read.query_name
                
                # Process completed pair when encountering new qname
                if qname != current_qname:
                    if current_qname is not None:
                        pair_info = process_read_pair(
                            current_qname, current_r1, current_r2, sorted_tags
                        )
                        output_line = format_output_line(pair_info, sorted_tags)
                        outfile.write("\t".join(output_line) + "\n")
                        
                        processed_qname_count += 1
                        if processed_qname_count % 500000 == 0:
                            logger.info(f"Written {processed_qname_count} read pairs")
                    
                    # Reset for new qname
                    current_qname = qname
                    current_r1 = None
                    current_r2 = None
                
                # Store reads
                if read.is_read1:
                    if current_r1 is None:
                        current_r1 = read
                elif read.is_read2:
                    if current_r2 is None:
                        current_r2 = read
            
            # Process final pair
            if current_qname is not None:
                pair_info = process_read_pair(
                    current_qname, current_r1, current_r2, sorted_tags
                )
                output_line = format_output_line(pair_info, sorted_tags)
                outfile.write("\t".join(output_line) + "\n")
                processed_qname_count += 1
            
            samfile.close()
            
        except Exception as e:
            logger.error(f"Error processing BAM file: {e}")
            raise
    
    logger.success(f"Processing complete!")
    logger.info(f"Total alignments scanned: {read_count:,}")
    logger.info(f"Total read pairs written: {processed_qname_count:,}")
    logger.info(f"Output saved to: {config.output_file}")


def main():
    """Main entry point for the script."""
    setup_logging()
    
    parser = argparse.ArgumentParser(
        description="Extract comprehensive summary for read pairs from BAM/SAM files"
    )
    parser.add_argument(
        "-i", "--input_bam",
        required=True,
        help="Path to the input BAM/SAM file (must be qname-sorted)"
    )
    parser.add_argument(
        "-o", "--output_file",
        required=True,
        help="Path to the output tab-delimited file"
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=4,
        help="Number of threads for BAM decompression (default: 4)"
    )
    
    args = parser.parse_args()
    
    try:
        # Create and validate configuration
        config = ProcessingConfig(
            input_bam=Path(args.input_bam),
            output_file=Path(args.output_file),
            threads=args.threads
        )
        
        # Process BAM file
        process_bam_file(config)
        
    except ValueError as e:
        logger.error(f"Configuration error: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()