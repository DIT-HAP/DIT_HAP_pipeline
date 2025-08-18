"""
GFF3 file parsing utilities for genomic feature extraction.
"""

import re
import pandas as pd
from typing import Optional, Dict, Any
from pathlib import Path


def extract_transcript_id(
    row: pd.Series,
    id_pattern: re.Pattern = re.compile(r"ID=([^:;]+)"),
    parent_pattern: re.Pattern = re.compile(r"Parent=([^:;]+)"),
    feature_types: tuple = ("CDS", "five_prime_UTR", "three_prime_UTR", "intron")
) -> Optional[str]:
    """
    Extract transcript ID from GFF attributes.
    
    For CDS and UTR features, extracts the Parent attribute as transcript ID.
    For other features, extracts the ID attribute directly.
    
    Parameters
    ----------
    row : pd.Series
        GFF row with 'Type' and 'Attributes' columns
    id_pattern : re.Pattern
        Compiled regex for ID extraction
    parent_pattern : re.Pattern  
        Compiled regex for Parent extraction
    feature_types : tuple
        Feature types that use Parent as transcript ID
        
    Returns
    -------
    Optional[str]
        Transcript ID or None if not found
        
    Examples
    --------
    >>> row = pd.Series({'Type': 'CDS', 'Attributes': 'Parent=SPAC1002.01.1'})
    >>> extract_transcript_id(row)
    'SPAC1002.01.1'
    """
    attributes = row.get('Attributes', '')
    feature_type = row.get('Type', '')
    
    if feature_type in feature_types:
        # For CDS, UTR, intron: use Parent
        match = parent_pattern.search(attributes)
    else:
        # For genes, mRNA: use ID
        match = id_pattern.search(attributes)
    
    return match.group(1) if match else None


class GFFParser:
    """Parser for GFF3 files with genomic feature extraction."""
    
    def __init__(self):
        self.id_pattern = re.compile(r"ID=([^:;]+)")
        self.parent_pattern = re.compile(r"Parent=([^:;]+)")
        
    def parse_gff_file(self, gff_file: Path) -> pd.DataFrame:
        """
        Parse GFF3 file into pandas DataFrame.
        
        Parameters
        ----------
        gff_file : Path
            Path to GFF3 file
            
        Returns
        -------
        pd.DataFrame
            Parsed GFF data with standard columns
        """
        if not gff_file.exists():
            raise FileNotFoundError(f"GFF file not found: {gff_file}")
        
        # Standard GFF3 column names
        columns = [
            'Chromosome', 'Source', 'Type', 'Start', 'End', 
            'Score', 'Strand', 'Phase', 'Attributes'
        ]
        
        try:
            # Read GFF, skip comments and headers
            gff_df = pd.read_csv(
                gff_file,
                sep='\t',
                comment='#',
                names=columns,
                na_values='.',
                dtype={'Chromosome': str, 'Start': int, 'End': int}
            )
            
            print(f"✅ Loaded {len(gff_df)} GFF records from {gff_file.name}")
            return gff_df
            
        except Exception as e:
            print(f"❌ Failed to parse GFF file: {e}")
            return pd.DataFrame()
    
    def extract_features_with_transcripts(self, gff_df: pd.DataFrame) -> pd.DataFrame:
        """
        Extract transcript IDs and create enhanced GFF DataFrame.
        
        Parameters
        ----------
        gff_df : pd.DataFrame
            Raw GFF data
            
        Returns
        -------
        pd.DataFrame
            GFF data with Transcript_ID column
        """
        # Add transcript ID column
        gff_df['Transcript_ID'] = gff_df.apply(
            lambda row: extract_transcript_id(
                row, self.id_pattern, self.parent_pattern
            ),
            axis=1
        )
        
        # Filter out rows without transcript IDs
        valid_features = gff_df.dropna(subset=['Transcript_ID'])
        
        print(f"✅ Extracted {len(valid_features)} features with transcript IDs")
        return valid_features
    
    def get_gene_categories(self, gff_df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        """
        Categorize genes by type (coding vs non-coding).
        
        Parameters
        ----------
        gff_df : pd.DataFrame
            GFF data with features
            
        Returns
        -------
        Dict[str, pd.DataFrame]
            Dictionary with 'coding' and 'noncoding' DataFrames
        """
        coding_types = ['CDS', 'mRNA']
        noncoding_types = ['tRNA', 'rRNA', 'snoRNA', 'snRNA', 'lncRNA']
        
        # Get unique genes by transcript ID and type
        unique_genes = gff_df.drop_duplicates(subset=['Transcript_ID', 'Type'])
        
        # Categorize genes
        coding_genes = unique_genes[unique_genes['Type'].isin(coding_types)]
        noncoding_genes = unique_genes[unique_genes['Type'].isin(noncoding_types)]
        
        print(f"✅ Found {len(coding_genes)} coding and {len(noncoding_genes)} non-coding genes")
        
        return {
            'coding': coding_genes,
            'noncoding': noncoding_genes
        }