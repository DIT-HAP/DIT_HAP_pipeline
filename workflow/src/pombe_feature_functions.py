#============================================= Import Libraries =============================================
from __future__ import annotations
import re
from pathlib import Path
from dataclasses import dataclass
from typing import Literal
from loguru import logger
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq, translate
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import gc_fraction, GC123, seq3
import gffutils
from codonbias.scores import EffectiveNumberOfCodons

#============================================== Constants ==============================================
POMBASE_VERSION = "2025-10-01"

CHR1_LEFT_TELEMERE_END = 29663
CHR1_RIGHT_TELEMERE_START = 5554844
CHR2_LEFT_TELEMERE_END = 39186
CHR2_RIGHT_TELEMERE_START = 4500619
CHR3_LEFT_RIBOSOMAL_DNA_END = 23130
CHR3_RIGHT_RIBOSOMAL_DNA_START = 2440994
CHROMOSOME_END = {
    "left": {
        "I": CHR1_LEFT_TELEMERE_END,
        "II": CHR2_LEFT_TELEMERE_END,
        "III": CHR3_LEFT_RIBOSOMAL_DNA_END
    },
    "right": {
        "I": CHR1_RIGHT_TELEMERE_START,
        "II": CHR2_RIGHT_TELEMERE_START,
        "III": CHR3_RIGHT_RIBOSOMAL_DNA_START
    }
}

CENTROMERE_POSITIONS = {
    "I": (3753687, 3789421),
    "II": (1602418, 1644747),
    "III": (1070904, 1137003)
}

#========================================== Configuration ============================================
@dataclass
class config:
    """Configuration parameters for PomBase feature functions."""

    PomBase_resource_dir: Path = Path(__file__).parent.parent.parent / "resources" / "pombase_data" / POMBASE_VERSION

    def __post_init__(self):
        """Ensure that the resource directory exists."""
        if not self.PomBase_resource_dir.exists():
            raise FileNotFoundError(f"PomBase resource directory not found: {self.PomBase_resource_dir}")
        
        # sequence and annotation files
        self.gene_meta_file: Path = self.PomBase_resource_dir / "Gene_metadata" / "gene_IDs_names_products.tsv"
        self.protein_meta: pd.DataFrame = pd.read_csv(str(self.PomBase_resource_dir / "Protein_features" / "peptide_stats.tsv"), sep="\t", index_col=[0])
        self.fasta_file: str = str(self.PomBase_resource_dir / "genome_sequence_and_features" / "Schizosaccharomyces_pombe_all_chromosomes.fa")
        self.genome_sequence_dict: dict[str, Seq] = SeqIO.to_dict(SeqIO.parse(self.fasta_file, "fasta"))
        self.genome_length_dict: dict[str, int] = {k: len(v) for k, v in self.genome_sequence_dict.items()}
        self.fai_file: str = str(self.PomBase_resource_dir / "genome_sequence_and_features" / "Schizosaccharomyces_pombe_all_chromosomes.fa.fai")
        self.gff3_file: str = str(self.PomBase_resource_dir / "genome_sequence_and_features" / "Schizosaccharomyces_pombe_all_chromosomes.gff3")
        self.database_file: str = str(self.PomBase_resource_dir / "genome_sequence_and_features" / "Schizosaccharomyces_pombe_all_chromosomes.db")
        self.primary_peptide_length: dict[str, int] = self.protein_meta["Residues"].to_dict()
#=========================================== Helper Functions ============================================
@logger.catch
def determine_primary_candidate(gene_id: str, mRNA_id: str, peptide_length: int, primary_peptide_length: int) -> bool:
    """Determine if the mRNA is a primary transcript based on peptide length."""
    # Check if the gene ID is one of the specific cases (same coding region but different UTRs)
    if gene_id in ["SPBC119.04", "SPBC17A3.07"]:
        if mRNA_id.endswith(".1"):
            return True
        else:
            return False
    # Specific cases where wrong annotation or NNNs in sequence
    elif gene_id in ["SPAC212.11", "SPAC2E12.05", "SPAC977.01", "SPMIT.03", "SPMIT.06", "SPMIT.08"]:
        if mRNA_id.endswith(".1"):
            return True
        else:
            return False
    else:
        if peptide_length == primary_peptide_length:
            return True
        else:
            return False

@logger.catch
def calculate_aliphatic_index_biopython(protein_sequence: str) -> float:
    """Calculate the aliphatic index of a protein sequence using Biopython's ProteinAnalysis."""
    analysis = ProteinAnalysis(protein_sequence)
    aa_percent = analysis.amino_acids_percent
    X_ala = aa_percent.get('A', 0) * 100  # Alanine
    X_val = aa_percent.get('V', 0) * 100  # Valine
    X_leu = aa_percent.get('L', 0) * 100  # Leucine
    X_ile = aa_percent.get('I', 0) * 100  # Isoleucine
    
    # Calculate aliphatic index
    aliphatic_index = X_ala + 2.9 * X_val + 3.9 * (X_leu + X_ile)
    
    return round(aliphatic_index, 3)

#============================================ Dataclass ============================================
@dataclass
class DNA_level_features:
    """Dataclass to hold DNA-level features."""
    Gene_id: str
    mRNA_id: str
    Chromosome: Literal["chr_II_telomeric_gap", "I", "II", "III", "mating_type_region", "mitochondrial"]
    Start: int
    End: int
    Strand: Literal["+", "-"]
    Abs_distance_from_telomere: float
    Relative_distance_from_telomere: float
    Abs_distance_from_centromere: float
    Relative_distance_from_centromere: float
    Gene_length: int
    GC_content_of_gene: float
    CDS_number: int
    GC_content_of_CDS: float
    Fraction_of_CDS: float
    GC3: float
    Containing_intron: bool
    Intron_number: int
    GC_content_of_intron: float
    Total_intron_length: int
    Average_intron_length: float
    Length_of_first_intron: int
    GC_contents_of_first_intron: float
    ENC: float
    Peptide_length: int
    Primary_peptide_length: int
    Primary_candidate: bool

    @classmethod
    def from_gffutils_feature(cls, mRNA: gffutils.Feature, db: gffutils.FeatureDB, cfg: config) -> DNA_level_features:
        """Create a DNA_level_features instance from a gffutils Feature."""

        # Extract the CDS features for the mRNA
        if mRNA.strand == "+":
            CDSs = list(db.children(mRNA, featuretype='CDS', order_by='start'))
            start = getattr(CDSs[0], 'start', 0)
            end = getattr(CDSs[-1], 'end', 0)
        else:
            CDSs = list(db.children(mRNA, featuretype='CDS', reverse=True, order_by='start'))
            start = getattr(CDSs[0], 'end', 0)
            end = getattr(CDSs[-1], 'start', 0)

        # basic features
        gene_id = mRNA.attributes.get("Parent")[0]
        mRNA_id = mRNA.id
        chrom = mRNA.chrom
        strand = mRNA.strand
        midpoint = (start + end) // 2

        # telemere and centromere related features
        abs_distance_from_telomere = min(abs(midpoint - CHROMOSOME_END["left"].get(chrom, np.nan)), abs(midpoint - CHROMOSOME_END["right"].get(chrom, np.nan)))
        relative_distance_from_telomere = round(abs_distance_from_telomere / cfg.genome_length_dict[chrom], 3)
        abs_distance_from_centromere = abs(midpoint - np.mean(CENTROMERE_POSITIONS.get(chrom, (np.nan, np.nan))))
        relative_distance_from_centromere = round(abs_distance_from_centromere / cfg.genome_length_dict[chrom], 3)

        # feature of the sequence from start codon to stop codon
        gene_length = abs(end - start) + 1
        GC_content_of_gene = round(gc_fraction(cfg.genome_sequence_dict[chrom][min(start, end):max(start, end)]), 3)

        # CDS features
        CDS_number = len(CDSs)
        CDS_sequence = ""
        for cds in CDSs:
            CDS_sequence += cds.sequence(cfg.fasta_file)
        GC_content_of_CDS = round(gc_fraction(CDS_sequence), 3)
        Fraction_of_CDS = round(len(CDS_sequence) / gene_length, 3)
        GC3 = round(GC123(Seq(CDS_sequence))[-1], 3)

        # Intron features
        introns = list(db.children(mRNA, featuretype='intron', order_by='start'))
        Containing_intron = len(introns) > 0
        Intron_number = len(introns)
        intron_sequences = [intron.sequence(cfg.fasta_file) for intron in introns]
        intron_sequence = "".join(intron_sequences)
        GC_content_of_intron = round(gc_fraction(intron_sequence), 3)
        Total_intron_length = len(intron_sequence)
        Average_intron_length = Total_intron_length / Intron_number if Intron_number > 0 else 0
        Length_of_first_intron = len(intron_sequences[0]) if Intron_number > 0 else 0
        GC_contents_of_first_intron = round(gc_fraction(intron_sequences[0]), 3) if Intron_number > 0 else 0.0

        # ENC for CUB
        ENC_model = EffectiveNumberOfCodons(mean="unweighted")
        ENC = np.round(ENC_model.get_score(CDS_sequence), 2)

        # Peptide lengths
        Peptide_length = len(Seq(CDS_sequence).translate(to_stop=True))
        primary_peptide_length = cfg.primary_peptide_length[gene_id]
        primary_candidate = determine_primary_candidate(gene_id, mRNA_id, Peptide_length, primary_peptide_length)

        return cls(
            Gene_id=gene_id,
            mRNA_id=mRNA_id,
            Chromosome=chrom,
            Start=start,
            End=end,
            Strand=strand,
            Abs_distance_from_telomere=abs_distance_from_telomere,
            Relative_distance_from_telomere=relative_distance_from_telomere,
            Abs_distance_from_centromere=abs_distance_from_centromere,
            Relative_distance_from_centromere=relative_distance_from_centromere,
            Gene_length=gene_length,
            GC_content_of_gene=GC_content_of_gene,
            CDS_number=CDS_number,
            GC_content_of_CDS=GC_content_of_CDS,
            Fraction_of_CDS=Fraction_of_CDS,
            GC3=GC3,
            Containing_intron=Containing_intron,
            Intron_number=Intron_number,
            GC_content_of_intron=GC_content_of_intron,
            Total_intron_length=Total_intron_length,
            Average_intron_length=Average_intron_length,
            Length_of_first_intron=Length_of_first_intron,
            GC_contents_of_first_intron=GC_contents_of_first_intron,
            ENC=ENC,
            Peptide_length=Peptide_length,
            Primary_peptide_length=primary_peptide_length,
            Primary_candidate=primary_candidate
        )

#============================================ Functions ============================================
@logger.catch
def extract_protein_features_from_peptide_sequence(peptide_fasta_file: Path, return_redundant_meta: bool = False) -> pd.DataFrame:
    """Extract protein-level features from a peptide sequence FASTA file."""
    records = []
    aa_content = {}
    aa_percent = {}
    for record in SeqIO.parse(peptide_fasta_file, "fasta"):
        gene_id = re.search(r"(\S+)\.\d:pep$", record.id).groups()[0]
        sequence = str(record.seq).rstrip("*")
        if "*" in sequence:
            print(f"Warning: Stop codon found in sequence of {gene_id}. Truncating at first stop codon.")
            sequence = sequence.split("*")[0]
        analysis = ProteinAnalysis(sequence)
        aa_content.update({gene_id: analysis.count_amino_acids()})
        aa_percent.update({gene_id: analysis.amino_acids_percent})
        protein_features = {
            "Gene_id": gene_id,
            "aromaticity": analysis.aromaticity(),
            "aliphatic_index": calculate_aliphatic_index_biopython(sequence),
            "gravy": analysis.gravy(),
            "flexibility": np.mean(analysis.flexibility()),
            "instability_index": analysis.instability_index(),
            "monoisotopic": analysis.monoisotopic
        }
        protein_features.update(
            dict(zip(("molar_extinction_reduced", "molar_extinction_cystines"), analysis.molar_extinction_coefficient()))
        )
        protein_features.update(
            dict(zip(("Helix_fraction", "Turn_fraction", "Sheet_fraction"), analysis.secondary_structure_fraction()))
        )
        if return_redundant_meta:
            protein_features.update({
                "charge_at_pH": analysis.charge_at_pH(7.0),
                "isoelectric_point": analysis.isoelectric_point(),
                "length": len(sequence),
                "molecular_weight(kDa)": analysis.molecular_weight() / 1000,
            })
        records.append(protein_features)

    aa_content_df = pd.DataFrame.from_dict(aa_content, orient='index')
    aa_percent_df = pd.DataFrame.from_dict(aa_percent, orient='index')

    aa_content_df.columns = [f"aa_count_{seq3(col)}" for col in aa_content_df.columns]
    aa_percent_df.columns = [f"aa_percent_{seq3(col)}" for col in aa_percent_df.columns]

    records_df = pd.DataFrame(records).set_index("Gene_id")
    records_df = records_df.join(aa_content_df).join(aa_percent_df).reset_index()
    return records_df