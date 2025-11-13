# ================================ Imports =================================
import gzip
from pathlib import Path
from loguru import logger
from typing import Literal
from tqdm import tqdm
import numpy as np
import pandas as pd
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

@logger.catch
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

@logger.catch
def pLDDT_statistics_report(
    structure_dir: Path,
    structure_format: Literal["pdb", "pdb.gz", "cif", "cif.gz", "mixed"] = "pdb.gz"
) -> pd.DataFrame:
    """Generates a report of pLDDT statistics for all protein structures in a directory."""
    all_pdb_files = list((structure_dir).glob(f"*.{structure_format}"))
    tqdm_iterator = tqdm(all_pdb_files)
    pLDDT_records = []

    for pdb_file in tqdm_iterator:
        uniprot_id = pdb_file.name.split("-F1-")[0].split("AF-")[1]
        match structure_format:
            case "pdb":
                pLDDT = np.array(extract_pLDDT_pdb(pdb_file))
            case "pdb.gz":
                pLDDT = np.array(extract_pLDDT_pdb_gz(pdb_file))
            case "cif" | "cif.gz" | "mixed":
                pLDDT = np.array(extract_pLDDT(pdb_file))
            case _:
                raise ValueError(f"Unsupported structure format: {structure_format}")
        length_protein = len(pLDDT)
        mean_pLDDT = np.mean(pLDDT)
        std_pLDDT = np.std(pLDDT)
        cv_pLDDT = std_pLDDT / mean_pLDDT if mean_pLDDT != 0 else np.nan
        disorder_fraction = np.sum(pLDDT < 50) / length_protein
        pLDDT_records.append({
            "uniprot_id": uniprot_id,
            "protein_length": length_protein,
            "pLDDT": ",".join(pLDDT.astype(str)),
            "mean_pLDDT": round(mean_pLDDT, 3),
            "std_pLDDT": round(std_pLDDT, 3),
            "cv_pLDDT": round(cv_pLDDT, 3),
            "disorder_fraction": round(disorder_fraction, 3),
        })

    pLDDTs = pd.DataFrame(pLDDT_records)
    return pLDDTs