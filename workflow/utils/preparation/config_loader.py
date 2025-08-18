"""
Configuration management for genomic feature analysis.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Any


@dataclass
class AnalysisConfig:
    """Configuration parameters for genomic feature analysis."""
    
    # File paths
    gff_file: Path
    peptide_stats_file: Path
    gene_names_file: Path
    gene_essentiality_file: Path
    output_dir: Path
    
    # Analysis parameters  
    min_cds_length: int = 50
    intergenic_min_length: int = 100
    primary_transcript_threshold: float = 0.8
    
    # Gene categories
    coding_types: tuple = ('CDS', 'mRNA', 'gene')
    noncoding_types: tuple = ('tRNA', 'rRNA', 'snoRNA', 'snRNA', 'lncRNA')
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        # Ensure output directory exists
        self.output_dir.mkdir(exist_ok=True, parents=True)
        
        # Validate required files exist
        for file_path in [self.gff_file, self.peptide_stats_file]:
            if not file_path.exists():
                raise FileNotFoundError(f"Required file not found: {file_path}")
    
    @classmethod
    def from_dict(cls, config_dict: Dict[str, Any]) -> 'AnalysisConfig':
        """Create config from dictionary with Path conversion."""
        # Convert string paths to Path objects
        path_fields = ['gff_file', 'peptide_stats_file', 'gene_names_file', 
                      'gene_essentiality_file', 'output_dir']
        
        for field in path_fields:
            if field in config_dict and isinstance(config_dict[field], str):
                config_dict[field] = Path(config_dict[field])
        
        return cls(**config_dict)


def load_default_config() -> AnalysisConfig:
    """Load default configuration for the analysis."""
    return AnalysisConfig(
        gff_file=Path("../data/schizosaccharomyces_pombe.chr.gff3"),
        peptide_stats_file=Path("../data/peptide_stats.tsv"),
        gene_names_file=Path("../data/gene_IDs_names.tsv"),
        gene_essentiality_file=Path("../data/gene_essentiality.tsv"),
        output_dir=Path("../results/genomic_features")
    )