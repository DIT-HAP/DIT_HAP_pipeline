"""
Feature processing utilities for genomic analysis.
"""

import pandas as pd
import numpy as np
from typing import Dict, Optional, Tuple


def calculate_accumulated_cds_bases(transcript_features_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate accumulated CDS bases for features within a transcript.
    
    Sorts features by start position. For negative strand genes, reverses 
    the sorted order to process features in transcription direction (3' to 5').
    
    Parameters
    ----------
    transcript_features_df : pd.DataFrame
        Features for a single transcript with Start, End, Strand columns
        
    Returns
    -------
    pd.DataFrame
        Features with CDS_accumulated_bases column
        
    Notes
    -----
    Accumulated bases represent the cumulative CDS length from the
    transcription start, accounting for strand direction.
    """
    if transcript_features_df.empty:
        return transcript_features_df
        
    # Sort by genomic start position
    sorted_features = transcript_features_df.sort_values('Start').copy()
    
    # For negative strand, reverse order for transcription direction
    strand = sorted_features['Strand'].iloc[0]
    if strand == '-':
        sorted_features = sorted_features.iloc[::-1].copy()
    
    # Calculate feature lengths and accumulated bases
    sorted_features['Feature_length'] = (
        sorted_features['End'] - sorted_features['Start'] + 1
    )
    
    # Cumulative sum for CDS bases only
    cds_mask = sorted_features['Type'] == 'CDS'
    cds_lengths = sorted_features['Feature_length'].where(cds_mask, 0)
    sorted_features['CDS_accumulated_bases'] = cds_lengths.cumsum()
    
    return sorted_features


class FeatureProcessor:
    """Processor for converting genomic features to analysis-ready formats."""
    
    def __init__(self):
        pass
    
    def features_to_bed_format(
        self, 
        transcript_features_df: pd.DataFrame,
        gene_type_label: str,
        peptide_length_map: Dict[str, int]
    ) -> pd.DataFrame:
        """
        Convert GFF features for a transcript to BED-like format.
        
        Parameters
        ----------
        transcript_features_df : pd.DataFrame
            Features for a single transcript (CDS, intron, etc.)
        gene_type_label : str
            Type label for the gene (e.g., 'protein_coding')
        peptide_length_map : Dict[str, int]
            Mapping of transcript IDs to peptide lengths
            
        Returns
        -------
        pd.DataFrame
            BED-formatted features with genomic annotations
            
        Notes
        -----
        Creates comprehensive annotations including:
        - Genomic coordinates (BED format)
        - Feature type and properties
        - Transcript and gene information
        - Accumulated CDS positions
        """
        if transcript_features_df.empty:
            return pd.DataFrame()
        
        # Calculate accumulated CDS bases
        features_with_cds = calculate_accumulated_cds_bases(transcript_features_df)
        
        # Extract basic information
        transcript_id = features_with_cds['Transcript_ID'].iloc[0]
        chromosome = features_with_cds['Chromosome'].iloc[0]
        strand = features_with_cds['Strand'].iloc[0]
        
        # Get peptide length
        peptide_length = peptide_length_map.get(transcript_id, 0)
        
        bed_records = []
        
        for _, feature in features_with_cds.iterrows():
            # BED coordinates (0-based start)
            bed_start = feature['Start'] - 1
            bed_end = feature['End']
            
            # Create feature annotation
            record = {
                '#Chr': chromosome,
                'Start': bed_start,
                'End': bed_end,
                'Systematic_ID': transcript_id,
                'Feature_type': feature['Type'],
                'Strand': strand,
                'Gene_type': gene_type_label,
                'Feature_length': feature['Feature_length'],
                'CDS_accumulated_bases': feature['CDS_accumulated_bases'],
                'Peptide_length': peptide_length,
                'GFF_Start': feature['Start'],
                'GFF_End': feature['End']
            }
            
            bed_records.append(record)
        
        return pd.DataFrame(bed_records)
    
    def identify_primary_transcripts(
        self, 
        coding_features_df: pd.DataFrame,
        length_threshold: float = 0.8
    ) -> pd.DataFrame:
        """
        Identify primary transcripts based on peptide length.
        
        Parameters
        ----------
        coding_features_df : pd.DataFrame
            All coding features with peptide length information
        length_threshold : float
            Minimum fraction of max length to be considered primary
            
        Returns
        -------
        pd.DataFrame
            Primary transcript features only
            
    Notes
        -----
        Primary transcripts are identified as the longest isoform
        per gene, or isoforms within the length threshold of the longest.
        """
        if coding_features_df.empty:
            return pd.DataFrame()
        
        # Group by systematic ID to find max peptide length per gene
        gene_stats = coding_features_df.groupby('Systematic_ID').agg({
            'Peptide_length': 'max'
        }).reset_index()
        gene_stats.columns = ['Systematic_ID', 'Max_peptide_length']
        
        # Merge back with features
        features_with_max = coding_features_df.merge(
            gene_stats, on='Systematic_ID', how='left'
        )
        
        # Calculate length ratio and identify primary transcripts
        features_with_max['Length_ratio'] = (
            features_with_max['Peptide_length'] / 
            features_with_max['Max_peptide_length']
        )
        
        # Flag primary transcripts
        primary_mask = features_with_max['Length_ratio'] >= length_threshold
        features_with_max['Primary_transcript_flag'] = primary_mask.map({True: 'Yes', False: 'No'})
        
        # Return only primary transcripts
        primary_transcripts = features_with_max[primary_mask].copy()
        
        print(f"✅ Identified {len(primary_transcripts)} primary transcript features")
        return primary_transcripts
    
    def process_noncoding_genes(
        self, 
        noncoding_genes_df: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Process non-coding genes to BED format.
        
        Parameters
        ----------
        noncoding_genes_df : pd.DataFrame
            Non-coding gene features from GFF
            
        Returns
        -------
        pd.DataFrame
            BED-formatted non-coding gene features
        """
        if noncoding_genes_df.empty:
            return pd.DataFrame()
        
        bed_records = []
        
        for _, gene in noncoding_genes_df.iterrows():
            record = {
                '#Chr': gene['Chromosome'],
                'Start': gene['Start'] - 1,  # BED 0-based
                'End': gene['End'],
                'Systematic_ID': gene['Transcript_ID'],
                'Feature_type': gene['Type'],
                'Strand': gene['Strand'],
                'Gene_type': 'non_coding',
                'Feature_length': gene['End'] - gene['Start'] + 1,
                'GFF_Start': gene['Start'],
                'GFF_End': gene['End']
            }
            bed_records.append(record)
        
        print(f"✅ Processed {len(bed_records)} non-coding genes")
        return pd.DataFrame(bed_records)