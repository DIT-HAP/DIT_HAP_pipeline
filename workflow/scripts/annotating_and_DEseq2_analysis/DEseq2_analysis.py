# %%
from pathlib import Path
import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

def load_and_preprocess_data(counts_file, annotations_file):
    # Load counts data
    counts_df = pd.read_csv(counts_file, index_col=[0, 1, 2, 3], header=[0, 1])
    counts_df.columns = ["_".join(col) for col in counts_df.columns]
    counts_df.index = ["=".join(map(str, index)) for index in counts_df.index]
    counts_df = counts_df.astype(int).T

    # Create metadata
    metadata = pd.DataFrame()
    metadata["sample"] = counts_df.index
    metadata["condition"] = [idx.split("_")[1] for idx in counts_df.index]
    metadata["group"] = [idx.split("_")[0] for idx in counts_df.index]
    metadata.set_index("sample", inplace=True)

    # Remove NA values
    counts_df = counts_df.loc[:, ~counts_df.isna().any(axis=0)].copy()

    # Load and process annotations
    insertion_annotations = pd.read_csv(annotations_file, index_col=[0, 1, 2, 3])
    insertion_annotations.index = ["=".join(map(str, index)) for index in insertion_annotations.index]

    return counts_df, metadata, insertion_annotations

def merge_counts_and_annotations(counts_df, insertion_annotations):
    counts_annotations = counts_df.T.merge(
        insertion_annotations, left_index=True, right_index=True, how="left")
    counts_annotations = counts_annotations[~counts_annotations.index.duplicated(keep='first')]
    return counts_annotations

def get_control_insertions(counts_annotations):
    return counts_annotations[
        (counts_annotations["Type"] == "Intergenic region") & 
        (counts_annotations["Distance_to_region_start"] > 500) & 
        (counts_annotations["Distance_to_region_end"] > 500)
    ].index

def create_deseq_dataset(counts_df, metadata, control_insertions):
    inference = DefaultInference(n_cpus=36)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors="condition",
        refit_cooks=True,
        inference=inference,
        ref_level=["condition", "0h"],
        min_replicates=3,
    )
    
    dds.fit_size_factors(control_genes=control_insertions)
    dds.fit_genewise_dispersions()
    dds.fit_dispersion_trend()
    dds.fit_dispersion_prior()
    dds.fit_MAP_dispersions()
    dds.fit_LFC()
    dds.calculate_cooks()
    if dds.refit_cooks:
        dds.refit()
    
    return dds

def perform_differential_analysis(dds, timepoints):
    stat_res = {}
    inference = DefaultInference(n_cpus=36)
    
    for tp in timepoints:
        stat_res[tp] = DeseqStats(
            dds, contrast=["condition", "0h", tp], inference=inference, 
            cooks_filter=True, independent_filter=True, quiet=True
        )
        stat_res[tp].summary()
        # Uncomment the following line if you want to perform LFC shrinkage
        # stat_res[tp].lfc_shrink(coeff=f"condition_{tp}_vs_0h")
    
    return stat_res

def plot_ma(stat_res, output_dir):
    for tp, res in stat_res.items():
        res.plot_MA(save_path=Path(output_dir) / f"MA_{tp}.png")

def concatenate_results(stat_res, timepoints):
    result_df = {tp: stat_res[tp].results_df for tp in timepoints}
    concated_results = pd.concat(result_df, axis=1)
    concated_results.index = pd.MultiIndex.from_tuples(
        concated_results.index.str.split("=").tolist())
    return concated_results

def parse_args():
    parser = argparse.ArgumentParser(description="Perform differential expression analysis on insertion counts.")
    parser.add_argument("-i", "--counts_file", type=Path, required=True, help="Path to the counts file.")
    parser.add_argument("-a", "--annotations_file", type=Path, required=True, help="Path to the annotations file.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output file.")
    parser.add_argument("-t", "--initial_timepoint", type=str, required=True, help="Initial timepoint to analyze.")
    return parser.parse_args()

def main():
    args = parse_args()

    counts_df, metadata, insertion_annotations = load_and_preprocess_data(args.counts_file, args.annotations_file)
    counts_annotations = merge_counts_and_annotations(counts_df, insertion_annotations)
    control_insertions = get_control_insertions(counts_annotations)

    timepoints = metadata["condition"].unique().tolist()
    timepoints.remove(args.initial_timepoint)
    
    dds = create_deseq_dataset(counts_df, metadata, control_insertions)
    dds.plot_dispersions()
    
    stat_res = perform_differential_analysis(dds, timepoints)
    plot_ma(stat_res, args.output.parent)
    
    concated_results = concatenate_results(stat_res, timepoints)
    concated_results.to_csv(Path(args.output))

if __name__ == "__main__":
    main()

# %%
