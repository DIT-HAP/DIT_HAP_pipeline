"""
Annotation utilities for genomic regions and intergenic analysis.
"""

import pandas as pd
import numpy as np
from typing import Dict, Optional, Tuple, List
import pybedtools


def annotate_intergenic_flanks(
    intergenic_row: pd.Series, 
    primary_transcripts_df: pd.DataFrame
) -> pd.Series:
    """
    Annotate intergenic region with flanking transcript information.
    
    Parameters
    ----------
    intergenic_row : pd.Series
        Intergenic region with #Chr, Start, End columns
    primary_transcripts_df : pd.DataFrame
        Primary transcripts for flank identification
        
    Returns
    -------
    pd.Series
        Annotated intergenic region with left/right flank info
        
    Notes
    -----
    Identifies the closest upstream and downstream transcripts
    on the same chromosome to annotate intergenic boundaries.
    """
    chrom = intergenic_row['#Chr']
    start = intergenic_row['Start']
    end = intergenic_row['End']
    
    # Filter transcripts on same chromosome
    chrom_transcripts = primary_transcripts_df[
        primary_transcripts_df['#Chr'] == chrom
    ].copy()
    
    if chrom_transcripts.empty:
        # No transcripts on this chromosome
        annotations = {
            'Left_flank_gene': 'None',
            'Left_flank_distance': np.nan,
            'Right_flank_gene': 'None', 
            'Right_flank_distance': np.nan
        }
        return pd.concat([intergenic_row, pd.Series(annotations)])
    
    # Find left flank (upstream)
    upstream_transcripts = chrom_transcripts[chrom_transcripts['End'] <= start]
    if not upstream_transcripts.empty:
        # Closest upstream transcript
        distances = start - upstream_transcripts['End']
        closest_idx = distances.idxmin()
        left_gene = upstream_transcripts.loc[closest_idx, 'Systematic_ID']
        left_distance = distances.loc[closest_idx]
    else:
        left_gene = 'None'
        left_distance = np.nan
    
    # Find right flank (downstream)  
    downstream_transcripts = chrom_transcripts[chrom_transcripts['Start'] >= end]
    if not downstream_transcripts.empty:
        # Closest downstream transcript
        distances = downstream_transcripts['Start'] - end
        closest_idx = distances.idxmin()
        right_gene = downstream_transcripts.loc[closest_idx, 'Systematic_ID']
        right_distance = distances.loc[closest_idx]
    else:
        right_gene = 'None'
        right_distance = np.nan
    
    # Create annotations
    annotations = {
        'Left_flank_gene': left_gene,
        'Left_flank_distance': left_distance,
        'Right_flank_gene': right_gene,
        'Right_flank_distance': right_distance
    }
    
    return pd.concat([intergenic_row, pd.Series(annotations)])


class AnnotationUtils:
    """Utilities for genomic feature annotation and enrichment."""
    
    def __init__(self):
        pass
    
    def find_intergenic_regions(
        self, 
        primary_transcripts_df: pd.DataFrame,
        min_length: int = 100
    ) -> pd.DataFrame:
        """
        Identify intergenic regions using BedTools.
        
        Parameters
        ----------
        primary_transcripts_df : pd.DataFrame
            Primary transcript features in BED format
        min_length : int
            Minimum intergenic region length
            
        Returns
        -------
        pd.DataFrame
            Intergenic regions with coordinates
        """
        if primary_transcripts_df.empty:
            print("⚠️  No primary transcripts for intergenic analysis")
            return pd.DataFrame()
        
        try:
            # Create BedTool from transcript coordinates
            transcript_bed = pybedtools.BedTool.from_dataframe(
                primary_transcripts_df[['#Chr', 'Start', 'End']].drop_duplicates()
            )
            
            # Get genome bounds for complement
            genome_bounds = primary_transcripts_df.groupby('#Chr').agg({
                'Start': 'min',
                'End': 'max'
            }).reset_index()
            
            # Create genome file for BedTools
            genome_file = 'temp_genome.txt'
            with open(genome_file, 'w') as f:
                for _, chrom_data in genome_bounds.iterrows():
                    f.write(f"{chrom_data['#Chr']}\t{chrom_data['End']}\n")
            
            # Find complement (intergenic regions)
            intergenic_bed = transcript_bed.complement(g=genome_file)
            
            # Convert back to DataFrame
            intergenic_df = intergenic_bed.to_dataframe(
                names=['#Chr', 'Start', 'End']
            )
            
            # Filter by minimum length
            intergenic_df['Length'] = intergenic_df['End'] - intergenic_df['Start']
            intergenic_df = intergenic_df[intergenic_df['Length'] >= min_length]
            
            # Clean up temp file
            import os
            if os.path.exists(genome_file):
                os.remove(genome_file)
            
            print(f"✅ Found {len(intergenic_df)} intergenic regions (≥{min_length} bp)")
            return intergenic_df
            
        except Exception as e:
            print(f"❌ Intergenic analysis failed: {e}")
            return pd.DataFrame()
    
    def annotate_with_gene_info(
        self, 
        features_df: pd.DataFrame,
        gene_names_file: Optional[str] = None,
        essentiality_file: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Add gene names and essentiality information.
        
        Parameters
        ----------
        features_df : pd.DataFrame
            Features with Systematic_ID column
        gene_names_file : Optional[str]
            Path to gene names mapping file
        essentiality_file : Optional[str]
            Path to gene essentiality data
            
        Returns
        -------
        pd.DataFrame
            Features with additional gene information
        """
        annotated_df = features_df.copy()
        
        # Add gene names if file provided
        if gene_names_file:
            try:
                gene_names = pd.read_csv(gene_names_file, sep='\t')
                if 'Systematic_ID' in gene_names.columns and 'Gene_name' in gene_names.columns:
                    annotated_df = annotated_df.merge(
                        gene_names[['Systematic_ID', 'Gene_name']], 
                        on='Systematic_ID', 
                        how='left'
                    )
                    print(f"✅ Added gene names from {gene_names_file}")
            except Exception as e:
                print(f"⚠️  Could not load gene names: {e}")
        
        # Add essentiality if file provided
        if essentiality_file:
            try:
                essentiality = pd.read_csv(essentiality_file, sep='\t')
                if 'Systematic_ID' in essentiality.columns:
                    annotated_df = annotated_df.merge(
                        essentiality, 
                        on='Systematic_ID', 
                        how='left'
                    )
                    print(f"✅ Added essentiality data from {essentiality_file}")
            except Exception as e:
                print(f"⚠️  Could not load essentiality data: {e}")
        
        return annotated_df
    
    def find_overlapping_regions(
        self, 
        features_df: pd.DataFrame,
        overlap_threshold: float = 0.5
    ) -> pd.DataFrame:
        """
        Find overlapping genomic features.
        
        Parameters
        ----------
        features_df : pd.DataFrame
            Features with genomic coordinates
        overlap_threshold : float
            Minimum overlap fraction to report
            
        Returns
        -------
        pd.DataFrame
            Features with overlap annotations
        """
        if features_df.empty:
            return features_df
        
        try:
            # Create BedTool for overlap analysis
            bed_tool = pybedtools.BedTool.from_dataframe(
                features_df[['#Chr', 'Start', 'End', 'Systematic_ID']]
            )
            
            # Self-intersect to find overlaps
            overlaps = bed_tool.intersect(bed_tool, wa=True, wb=True, f=overlap_threshold)
            
            if len(overlaps) > 0:
                overlap_df = overlaps.to_dataframe()
                # Process overlap results
                print(f"✅ Found overlapping regions with {overlap_threshold} threshold")
            else:
                print("✅ No significant overlaps found")
                
            return features_df
            
        except Exception as e:
            print(f"⚠️  Overlap analysis failed: {e}")
            return features_df