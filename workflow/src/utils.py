# ================================ Imports =================================
from pathlib import Path
import pandas as pd

# ================================= Utility Functions =================================
def read_file(
    file: Path,
    **kwargs,
) -> pd.DataFrame:
    """Read a file into a pandas DataFrame based on file extension."""
    if "tsv" in file.name:
        return pd.read_csv(file, sep="\t", **kwargs)
    elif "bed" in file.name:
        return pd.read_csv(file, sep="\t", **kwargs)
    elif "csv" in file.name:
        return pd.read_csv(file, sep=",", **kwargs)
    elif "xlsx" in file.name:
        return pd.read_excel(file, **kwargs)
    else:
        raise ValueError(f"Unsupported file type: {file.name}")
