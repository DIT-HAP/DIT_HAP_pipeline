import sys
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def main(args):
    # get arguments
    cutoff = float(args.cutoff)
    level_type = args.type
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    # read files
    if level_type == "Site":
        use_index_col = [0, 1, 2, 3]
    elif level_type == "Bin":
        use_index_col = [0, 1, 2, 3, 4]
    elif level_type == "Domain":
        use_index_col = [0, 1, 2]
    elif level_type == "Gene":
        use_index_col = [0]

    weighted_M = pd.read_csv(args.weighted_M, index_col=use_index_col, header=0)
    slope = pd.read_csv(args.input_slope, index_col=use_index_col, header=0)
    lag = pd.read_csv(args.input_lag, index_col=use_index_col, header=0)
    stat = pd.read_csv(args.input_stat, index_col=use_index_col, header=0)

    timepoints = slope.columns.str.split("_")
    timepoints = timepoints[0]

    DG_raw = weighted_M[timepoints[-1]].apply(
        lambda row: find_interval_across_cutoff(
            row, timepoints, cutoff
        )
    )

    DG_DR_DL = DG_raw.apply(
        lambda row: assign_DR_DL_func(row, slope, lag, stat), axis=1
    )

    DG_DR_DL = DG_DR_DL.astype(
        {
            "Individual insertions": int,
            "Confidence score": float,
            "DG": str,
            "DR": float,
            "DL": float,
        }
    )

    DG_DR_DL.to_csv(output, index=True, header=True, float_format="%.3f")


def find_interval_across_cutoff(row, timepoints, cutoff):

    has_lag_or_not = "NO"
    DG = "_".join(timepoints)
    if row >= cutoff:
        has_lag_or_not = "YES"
    return pd.Series(
        [has_lag_or_not, DG], index=["has_lag_or_not", "DG"]
    )


def assign_DR_DL_func(row, slope, lag, stat):

    row_id = row.name
    istat = stat.loc[row_id]
    has_DL_or_not = row["has_lag_or_not"]
    DG = row["DG"]
    DR = slope.loc[row_id, DG]
    if has_DL_or_not == "NO":
        DL = np.nan
    else:
        DL = lag.loc[row_id, DG]

    DR_DL = pd.Series([DG, DR, DL], index=["DG", "DR", "DL"])

    return pd.concat([istat, DR_DL])


if __name__ == "__main__":
    # add arguments
    parser = argparse.ArgumentParser(description="Assign DR and DL")
    parser.add_argument(
        "-wm",
        "--weighted-M",
        dest="weighted_M",
        type=Path,
        required=True,
        help="Weighted M file",
    )
    parser.add_argument(
        "-is",
        "--input-slope",
        dest="input_slope",
        type=Path,
        required=True,
        help="Input slope file",
    )
    parser.add_argument(
        "-il",
        "--input-lag",
        dest="input_lag",
        type=Path,
        required=True,
        help="Input lag file",
    )
    parser.add_argument(
        "-istat",
        "--input-stat",
        dest="input_stat",
        type=Path,
        required=True,
        help="Input stat file",
    )
    parser.add_argument(
        "-c", "--cutoff", dest="cutoff", required=True, help="Cutoff value"
    )
    parser.add_argument(
        "-t", "--type", dest="type", required=True, help="Type of level"
    )
    parser.add_argument(
        "-o", "--output", dest="output", type=Path, required=True, help="Output file"
    )

    # parse arguments
    args = parser.parse_args()

    main(args)
