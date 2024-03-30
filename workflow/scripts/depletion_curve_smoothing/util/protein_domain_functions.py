import numpy as np
import pandas as pd

def assign_protein_domain(row, domain):
    Gene = row["Systematic ID"]
    if Gene in domain["Systematic ID"].values:
        residue = int(row["Residue_affected"])
        domains = domain[domain["Systematic ID"] == Gene].copy()
        for idx, row in domains.iterrows():
            for iDomain in row.loc["domain_residues"].split(","):
                domain_start = int(iDomain.split("-")[0])
                domain_end = int(iDomain.split("-")[1])
                domain_id = row.loc["domain_id"]
                if residue >= domain_start and residue <= domain_end:
                    return domain_id, row.loc["domain_residues"]
        return np.nan, np.nan
    else:
        return np.nan, np.nan