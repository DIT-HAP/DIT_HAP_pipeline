# %%
from pathlib import Path
import numpy as np
import pandas as pd
import argparse
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

def load_and_preprocess_data(counts_file, annotations_file):
    # Load counts data
    counts_df = pd.read_csv(counts_file, index_col=[0, 1, 2, 3], header=[0, 1], sep="\t")
    counts_df_index_names = counts_df.index.names
    counts_df_columns_names = counts_df.columns.names

    counts_df.columns = ["#".join(col) for col in counts_df.columns]
    counts_df.index = ["=".join(map(str, index)) for index in counts_df.index]
    counts_df = counts_df.astype(int).T

    # Create metadata
    metadata = pd.DataFrame()
    metadata["sample"] = counts_df.index
    metadata["condition"] = [idx.split("#")[1] for idx in counts_df.index]
    metadata["group"] = [idx.split("#")[0] for idx in counts_df.index]
    metadata.set_index("sample", inplace=True)

    # Remove NA values
    counts_df = counts_df.loc[:, ~counts_df.isna().any(axis=0)].copy()

    # Load and process annotations
    insertion_annotations = pd.read_csv(annotations_file, index_col=[0, 1, 2, 3], sep="\t")
    insertion_annotations_index_names = insertion_annotations.index.names
    insertion_annotations.index = ["=".join(map(str, index)) for index in insertion_annotations.index]

    return counts_df, metadata, insertion_annotations, counts_df_index_names, counts_df_columns_names, insertion_annotations_index_names

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

def create_deseq_dataset(counts_df, metadata, control_insertions, initial_timepoint="0h"):
    inference = DefaultInference(n_cpus=36)
    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design="~condition",
        refit_cooks=True,
        inference=inference,
        min_replicates=7,
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

def perform_differential_analysis(dds, timepoints, initial_timepoint="0h"):
    stat_res = {}
    inference = DefaultInference(n_cpus=36)
    
    for tp in timepoints:
        stat_res[tp] = DeseqStats(
            dds, contrast=["condition", tp, initial_timepoint], inference=inference, 
            cooks_filter=True, independent_filter=True, quiet=True
        )
        stat_res[tp].summary()
        # Uncomment the following line if you want to perform LFC shrinkage
        # stat_res[tp].lfc_shrink(coeff=f"condition[T.{tp}]")
    
    return stat_res

def plot_ma(stat_res, output_dir):
    for tp, res in stat_res.items():
        res.plot_MA(save_path=Path(output_dir) / f"MA_{tp}.png")

def concatenate_results(stat_res, timepoints):
    result_df = {}
    for tp in timepoints:
        result_df[tp] = stat_res[tp].results_df
        result_df[tp]["log2FoldChange"] = -result_df[tp]["log2FoldChange"]
        result_df[tp]["stat"] = -result_df[tp]["stat"]
    concated_results = pd.concat(result_df, axis=1)
    concated_results.index = pd.MultiIndex.from_tuples(
        concated_results.index.str.split("=").tolist())
    # Convert string format numbers to numeric values in the MultiIndex
    new_index = []
    for idx in concated_results.index:
        chr_name = idx[0]
        # Convert coordinate from string to integer
        coordinate = int(idx[1]) if idx[1].isdigit() else idx[1]
        strand = idx[2]
        target = idx[3]
        new_index.append((chr_name, coordinate, strand, target))
    
    # Create a new MultiIndex with the converted values
    concated_results.index = pd.MultiIndex.from_tuples(
        new_index, names=concated_results.index.names)
    return concated_results

def transform_index_to_multiindex(dds, layer_name):
    df = pd.DataFrame(dds.layers[layer_name], index=dds.obs.index.tolist(), columns=dds.var.index.tolist()).T
    df.index = pd.MultiIndex.from_tuples(df.index.str.split("=").tolist())
    new_index = []
    for idx in df.index:
        chr_name = idx[0]
        # Convert coordinate from string to integer
        coordinate = int(idx[1]) if idx[1].isdigit() else idx[1]
        strand = idx[2]
        target = idx[3]
        new_index.append((chr_name, coordinate, strand, target))
    
    # Create a new MultiIndex with the converted values
    df.index = pd.MultiIndex.from_tuples(
        new_index, names=df.index.names)
    df.columns = pd.MultiIndex.from_tuples(df.columns.str.split("#").tolist())

    return df

def parse_args():
    parser = argparse.ArgumentParser(description="Perform differential expression analysis on insertion counts.")
    parser.add_argument("-i", "--counts_file", type=Path, required=True, help="Path to the counts file.")
    parser.add_argument("-t", "--initial_timepoint", type=str, required=True, help="Initial timepoint to analyze.")
    parser.add_argument("-a", "--annotations_file", type=Path, required=True, help="Path to the annotations file.")
    parser.add_argument("-o", "--output", type=Path, required=True, help="Path to the output file.")
    return parser.parse_args()

def main():

    args = parse_args()

    print("*** Loading and preprocessing data...")
    counts_df, metadata, insertion_annotations, counts_df_index_names, counts_df_columns_names, insertion_annotations_index_names = load_and_preprocess_data(args.counts_file, args.annotations_file)

    print("The metadata for analysis:")
    print(metadata)

    print("*** Merging counts and annotations...")
    counts_annotations = merge_counts_and_annotations(counts_df, insertion_annotations)

    print("*** Getting control insertions...")
    control_insertions = get_control_insertions(counts_annotations)
    print("*** The number of control insertions:", control_insertions.shape[0])


    timepoints = metadata["condition"].unique().tolist()
    timepoints.remove(args.initial_timepoint)
    print("*** The control timepoint:", args.initial_timepoint)
    print("*** The timepoints for analysis:", timepoints)
    
    print("*** Creating DESeq2 dataset...")
    dds = create_deseq_dataset(counts_df, metadata, control_insertions, args.initial_timepoint)
    print("*** Plotting dispersions...")
    dds.plot_dispersions()

    print("*** Transforming index to multiindex...")
    normalized_counts = transform_index_to_multiindex(dds, "normed_counts")
    normalized_counts = normalized_counts.rename_axis(counts_df_index_names, axis=0).rename_axis(counts_df_columns_names, axis=1)
    normalized_counts.to_csv(Path(args.output.parent) / "normed_counts.tsv", index=True, float_format="%.3f", sep="\t")

    print("*** Transforming count_X to multiindex...")
    count_X = pd.DataFrame(dds.X, index=dds.obs.index.tolist(), columns=dds.var.index.tolist()).T
    count_X.index = pd.MultiIndex.from_tuples(count_X.index.str.split("=").tolist())
    count_X.columns = pd.MultiIndex.from_tuples(count_X.columns.str.split("#").tolist())
    count_X = count_X.rename_axis(counts_df_index_names, axis=0).rename_axis(counts_df_columns_names, axis=1)
    count_X.to_csv(Path(args.output.parent) / "count_X.tsv", index=True, float_format="%.3f", sep="\t")

    # refit_count_norm = transform_index_to_multiindex(dds.filtered_genes, "normed_counts")
    # dds.filtered_genes.to_csv(Path(args.output.parent) / "refit_count_norm.csv", index=True, float_format="%.3f")

    cooks_df = transform_index_to_multiindex(dds, "cooks")
    cooks_df = cooks_df.rename_axis(counts_df_index_names, axis=0).rename_axis(counts_df_columns_names, axis=1)
    cooks_df.to_csv(Path(args.output.parent) / "cooks.tsv", index=True, float_format="%.3f", sep="\t")

    print("*** Performing differential analysis...")
    stat_res = perform_differential_analysis(dds, timepoints, args.initial_timepoint)
    plot_ma(stat_res, args.output.parent)

    print("*** Concatenating results...")
    concated_results = concatenate_results(stat_res, timepoints)
    concated_results = concated_results.rename_axis(insertion_annotations_index_names, axis=0).rename_axis(["Timepoint", "Statistic"], axis=1)

    # add the metrics (baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) for the initial timepoint
    print("*** Adding the metrics for the initial timepoint...")
    print("All timepoints share the same baseMean")
    baseMean_initial = concated_results.xs("baseMean", axis=1, level="Statistic").iloc[:,0]
    print("Set the initial log2FoldChange to 0")
    log2FoldChange_initial = 0
    print("Set the initial lfcSE to NaN")
    lfcSE_initial = np.nan
    print("Set the initial stat to NaN")
    stat_initial = np.nan
    print("Set the initial pvalue to 1")
    pvalue_initial = 1
    print("Set the initial padj to 1")
    padj_initial = 1
    print("Insert the initial timepoint metrics to the first column...")
    concated_results.insert(0, (args.initial_timepoint,"padj"), padj_initial)
    concated_results.insert(0, (args.initial_timepoint,"pvalue"), pvalue_initial)
    concated_results.insert(0, (args.initial_timepoint,"stat"), stat_initial)
    concated_results.insert(0, (args.initial_timepoint,"lfcSE"), lfcSE_initial)
    concated_results.insert(0, (args.initial_timepoint,"log2FoldChange"), log2FoldChange_initial)
    concated_results.insert(0, (args.initial_timepoint,"baseMean"), baseMean_initial)

    concated_results = concated_results.round(3)
    print("*** Saving results...")
    concated_results.to_csv(Path(args.output), sep="\t")

    # baseMean,log2FoldChange,lfcSE,stat,pvalue,padj

    baseMean_df = concated_results.xs("baseMean", axis=1, level="Statistic")
    baseMean_df.to_csv(Path(args.output.parent) / "baseMean.tsv", index=True, sep="\t")

    LFC_df = concated_results.xs("log2FoldChange", axis=1, level="Statistic")
    LFC_df.to_csv(Path(args.output.parent) / "LFC.tsv", index=True, sep="\t")

    lfcSE_df = concated_results.xs("lfcSE", axis=1, level="Statistic")
    lfcSE_df.to_csv(Path(args.output.parent) / "lfcSE.tsv", index=True, sep="\t")

    stat_df = concated_results.xs("stat", axis=1, level="Statistic")
    stat_df.to_csv(Path(args.output.parent) / "stat.tsv", index=True, sep="\t")

    pvalue_df = concated_results.xs("pvalue", axis=1, level="Statistic")
    pvalue_df.to_csv(Path(args.output.parent) / "pvalue.tsv", index=True, sep="\t")

    padj_df = concated_results.xs("padj", axis=1, level="Statistic")
    padj_df.to_csv(Path(args.output.parent) / "padj.tsv", index=True, sep="\t")

if __name__ == "__main__":
    main()
