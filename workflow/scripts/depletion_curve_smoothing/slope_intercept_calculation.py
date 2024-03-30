"""

"""

import sys
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
from numpy import trapz


def main(args):
    generation_file = Path(args.generation_file)
    M_value_file = Path(args.M_value_file)
    level_type = args.level_type
    start_timepoint = args.start_timepoint
    output_slope = Path(args.output_slope)
    output_intercept = Path(args.output_intercept)
    output_lag = Path(args.output_lag)
    output_area = Path(args.output_area)
    output_slope.parent.mkdir(parents=True, exist_ok=True)
    output_intercept.parent.mkdir(parents=True, exist_ok=True)
    output_lag.parent.mkdir(parents=True, exist_ok=True)

    generation = (
        pd.read_csv(generation_file, index_col=[0], header=[0]).mean(axis=1).round(3)
    )
    if level_type == "Site":
        M_values = pd.read_csv(M_value_file, index_col=[0, 1, 2, 3], header=[0])
    elif level_type == "Bin":
        M_values = pd.read_csv(M_value_file, index_col=[0, 1, 2, 3, 4], header=[0])
    elif level_type == "Domain":
        M_values = pd.read_csv(M_value_file, index_col=[0, 1, 2], header=[0])
    elif level_type == "Gene":
        M_values = pd.read_csv(M_value_file, index_col=[0], header=[0])

    area_and_last_M = M_values.apply(
        calculate_area,
        axis=1,
        generation=generation,
        start_timepoint=start_timepoint,
    )

    slope_intercept_lag = M_values.apply(
        calculate_slope_and_intercept,
        axis=1,
        generation=generation,
        start_timepoint=start_timepoint,
    )

    area_and_last_M.to_csv(
        output_area, header=True, index=True, float_format="%.3f"
    )
    slope_intercept_lag.xs("slope", axis=1).to_csv(
        output_slope, header=True, index=True, float_format="%.3f"
    )
    slope_intercept_lag.xs("intercept", axis=1).to_csv(
        output_intercept, header=True, index=True, float_format="%.3f"
    )
    slope_intercept_lag.xs("lag", axis=1).to_csv(
        output_lag, header=True, index=True, float_format="%.3f"
    )



def calculate_slope_and_intercept(row, generation, start_timepoint):
    GMs = (
        pd.concat([generation, row], axis=1, keys=["G", "M"])
        .sort_values("G", ascending=True)
        .loc[start_timepoint:]
    )

    GM_early_timepoint = GMs.index[:-1]
    GM_late_timepoint = GMs.index[1:]

    slope_intercept_lag = {
        "slope": pd.Series(),
        "intercept": pd.Series(),
        "lag": pd.Series(),
    }

    for etp, ltp in zip(GM_early_timepoint, GM_late_timepoint):
        tp_pair = etp + "_" + ltp
        slope = (GMs.loc[ltp, "M"] - GMs.loc[etp, "M"]) / (
            GMs.loc[ltp, "G"] - GMs.loc[etp, "G"]
        )
        intercept = GMs.loc[etp, "M"] - slope * GMs.loc[etp, "G"]
        if slope == 0:
            lag = np.nan
        else:
            lag = -intercept / slope

        slope_intercept_lag["slope"].loc[tp_pair] = slope
        slope_intercept_lag["intercept"].loc[tp_pair] = intercept
        slope_intercept_lag["lag"].loc[tp_pair] = lag

    return pd.concat(slope_intercept_lag, axis=0)

def calculate_area(row, generation, start_timepoint):

    GMs = (
        pd.concat([generation, row], axis=1, keys=["G", "M"])
        .sort_values("G", ascending=True)
        .loc[start_timepoint:]
    )
    
    x = GMs["G"]
    y = GMs["M"]

    last_M = y.iloc[-1]

    area = trapz(y, x)

    return pd.Series([area, last_M], index=["area", "last_M"])


if __name__ == "__main__":
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Calculate the slope and intercept of the linear regression of M values and generation."
    )
    parser.add_argument(
        "-G",
        "--generation",
        dest="generation_file",
        required=True,
        type=Path,
        help="File of generation",
    )
    parser.add_argument(
        "-M",
        "--M-values",
        dest="M_value_file",
        required=True,
        type=Path,
        help="File of M values",
    )
    parser.add_argument(
        "-lt",
        "--level-type",
        dest="level_type",
        required=True,
        type=str,
        help="Level type",
    )
    parser.add_argument(
        "-stp",
        "--start-timepoint",
        dest="start_timepoint",
        required=True,
        type=str,
        help="Start timepoint",
    )
    parser.add_argument(
        "-os",
        "--output-slope",
        dest="output_slope",
        required=True,
        type=Path,
        help="Output file for slope",
    )
    parser.add_argument(
        "-oi",
        "--output-intercept",
        dest="output_intercept",
        required=True,
        type=Path,
        help="Output file for intercept",
    )
    parser.add_argument(
        "-ol",
        "--output-lag",
        dest="output_lag",
        required=True,
        type=Path,
        help="Output file for lag",
    )
    parser.add_argument(
        "-oa",
        "--output-area",
        dest="output_area",
        required=True,
        type=Path,
        help="Output file for area",
    )

    args = parser.parse_args()

    main(args)
