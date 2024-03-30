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
    early_timepoints = [i[0] for i in timepoints]
    late_timepoints = [i[1] for i in timepoints]
    timepoints = early_timepoints + [late_timepoints[-1]]

    DG_raw = weighted_M[late_timepoints].apply(
        lambda row: find_interval_across_cutoff(
            row, early_timepoints, late_timepoints, cutoff
        ),
        axis=1,
    )

    DG_DR_DL = DG_raw.apply(
        lambda row: uses_DG_or_nextDG(row, slope, lag, stat), axis=1
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


def find_interval_across_cutoff(row, early_timepoints, late_timepoints, cutoff):
    loc_gt_cutoff = np.where(row > cutoff)[0]

    has_lag_or_not = "NO"
    if loc_gt_cutoff.size == 0:
        etp = late_timepoints[-2]
        ltp = late_timepoints[-1]
        DG = etp + "_" + ltp
        next_DG = DG
    elif loc_gt_cutoff[0] == row.size - 1:
        has_lag_or_not = "YES"
        etp = late_timepoints[-2]
        ltp = late_timepoints[-1]
        DG = etp + "_" + ltp
        next_DG = DG
    else:
        has_lag_or_not = "YES"
        etp = early_timepoints[loc_gt_cutoff[0]]
        ltp = late_timepoints[loc_gt_cutoff[0]]
        lltp = late_timepoints[loc_gt_cutoff[0] + 1]
        DG = etp + "_" + ltp
        next_DG = ltp + "_" + lltp
    return pd.Series(
        [has_lag_or_not, DG, next_DG], index=["has_lag_or_not", "DG", "next_DG"]
    )


def uses_DG_or_nextDG(row, slope, lag, stat):
    row_id = row.name
    istat = stat.loc[row_id]
    has_DL_or_not = row["has_lag_or_not"]
    slope_DG = slope.loc[row_id, row["DG"]]
    slope_next_DG = slope.loc[row_id, row["next_DG"]]
    lag_DG = lag.loc[row_id, row["DG"]]
    lag_next_DG = lag.loc[row_id, row["next_DG"]]
    if slope_next_DG > slope_DG:
        DG = row["next_DG"]
        DR = slope_next_DG
        if has_DL_or_not == "YES":
            DL = lag_next_DG
        else:
            DL = np.nan
    else:
        DG = row["DG"]
        DR = slope_DG
        if has_DL_or_not == "YES":
            DL = lag_DG
        else:
            DL = np.nan

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
        "-o", "--output", dest="output", required=True, help="Output file"
    )

    # parse arguments
    args = parser.parse_args()

    main(args)
