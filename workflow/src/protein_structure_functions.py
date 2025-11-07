# ================================ Imports =================================
import gzip
from pathlib import Path
from loguru import logger
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

# ================================ Constants =================================


# =============================== Protein Structure Functions =================================
@logger.catch
def extract_pLDDT(structure_file: Path | str) -> list[float]:
    """Extracts the per-residue pLDDT scores from a PDB or mmCIF file."""

    if isinstance(structure_file, str):
        structure_file = Path(structure_file)

    # compressed or not
    if structure_file.name.endswith(".gz"):
        f = gzip.open(structure_file, "rt")
    else:
        f = open(structure_file, "r")

    # PDB or mmCIF
    if structure_file.name.rstrip(".gz").lower().endswith(".pdb"):
        parser = PDBParser()
        structure = parser.get_structure(structure_file.stem, f)
    elif structure_file.name.rstrip(".gz").lower().endswith(".cif"):
        parser = MMCIFParser()
        structure = parser.get_structure(structure_file.stem, f)
    else:
        raise ValueError(f"Unknown file format: {structure_file.name}")

    # extract pLDDT
    pLDDT = []
    for residue in structure.get_residues():
        if residue.has_id("CA"):
            pLDDT.append(residue["CA"].bfactor)

    f.close()

    return pLDDT

@logger.catch
def extract_pLDDT_pdb_gz(structure_file: Path | str) -> list[float]:
    """Extracts the per-residue pLDDT scores from a pdb.gz file."""
    if isinstance(structure_file, str):
        structure_file = Path(structure_file)

    f = gzip.open(structure_file, "rt")
    parser = PDBParser()
    structure = parser.get_structure(structure_file.stem, f)

    # extract pLDDT
    pLDDT = []
    for residue in structure.get_residues():
        pLDDT.append(residue["CA"].bfactor)

    f.close()

    return pLDDT

@logger.catch
def extract_pLDDT_pdb(structure_file: Path | str) -> list[float]:
    """Extracts the per-residue pLDDT scores from a PDB or mmCIF file."""

    if isinstance(structure_file, str):
        structure_file = Path(structure_file)

    parser = PDBParser()
    structure = parser.get_structure(structure_file.stem, structure_file)

    # extract pLDDT
    pLDDT = []
    for residue in structure.get_residues():
        pLDDT.append(residue["CA"].bfactor)

    return pLDDT


def extract_protein_seq_pdb_gz(structure_file: Path | str) -> str:
    """Extracts the residue sequence from a pdb.gz file."""
    if isinstance(structure_file, str):
        structure_file = Path(structure_file)
    
    f = gzip.open(structure_file, "rt")
    parser = PDBParser()
    structure = parser.get_structure(structure_file.stem, f)

    ppb = PPBuilder()
    seq = ppb.build_peptides(structure)[0].get_sequence()

    f.close()

    return seq