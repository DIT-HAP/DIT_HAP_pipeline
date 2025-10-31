# %% ============================================= Import Libraries =============================================
import sys
from pathlib import Path
from dataclasses import dataclass
from loguru import logger
import gffutils
import gffutils.biopython_integration as bp
from Bio.Seq import Seq, translate
# %% =========================================== Constants ============================================
POMBASE_VERSION = "2025-10-01"
# %% ========================================== Configuration ============================================
@dataclass
class config:
    """Configuration parameters for PomBase feature functions."""

    PomBase_resource_dir: Path = Path(__file__).parent.parent.parent / "resources" / "pombase_data" / POMBASE_VERSION

    def __post_init__(self):
        """Ensure that the resource directory exists."""
        if not self.PomBase_resource_dir.exists():
            raise FileNotFoundError(f"PomBase resource directory not found: {self.PomBase_resource_dir}")
        
        # sequence and annotation files
        self.fasta_file: str = str(self.PomBase_resource_dir / "genome_sequence_and_features" / "Schizosaccharomyces_pombe_all_chromosomes.fa")
        self.fai_file: str = str(self.PomBase_resource_dir / "genome_sequence_and_features" / "Schizosaccharomyces_pombe_all_chromosomes.fa.fai")
        self.gff3_file: str = str(self.PomBase_resource_dir / "genome_sequence_and_features" / "Schizosaccharomyces_pombe_all_chromosomes.gff3")
        self.database_file: str = str(self.PomBase_resource_dir / "genome_sequence_and_features" / "Schizosaccharomyces_pombe_all_chromosomes.db")


# %% ============================================ Logger Setup ============================================
def setup_logger():
    """Set up the logger for the module."""
    logger.remove()
    logger.add(
        sys.stderr,
        format="{time:YYYY-MM-DD HH:mm:ss} | {level:<8} | {message}",
        level="DEBUG"
    )

# %% ============================================ Functions ============================================
cfg = config()
# %%

 

db = gffutils.create_db(cfg.gff3_file, cfg.database_file, force=True)
db = gffutils.FeatureDB(cfg.database_file)

# %%
mRNAs = list(db.features_of_type('mRNA'))

coding_genes = set([
    gene.attributes.get("Parent")[0] for gene in mRNAs
])
# %%
