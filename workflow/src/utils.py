# ================================ Imports =================================
from pathlib import Path
from typing import List, Optional
import pandas as pd

# ================================= Utility Functions =================================
def read_file(
    file: Path,
    index_col: Optional[List[int]] = None,
    header: Optional[List[int] | str] = "infer",
    **kwargs,
) -> pd.DataFrame:
    """Read a file into a pandas DataFrame based on file extension."""
    if "tsv" in file.name:
        return pd.read_csv(file, index_col=index_col, header=header, sep="\t", **kwargs)
    elif "bed" in file.name:
        return pd.read_csv(file, index_col=index_col, header=header, sep="\t", **kwargs)
    elif "csv" in file.name:
        return pd.read_csv(file, index_col=index_col, header=header, sep=",", **kwargs)
    elif "xlsx" in file.name:
        if header == "infer":
            header = 0
        return pd.read_excel(file, index_col=index_col, header=header, **kwargs)
    else:
        raise ValueError(f"Unsupported file type: {file.name}")