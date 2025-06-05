import pandas as pd
import numpy as np


def calculate_weight_averaged_M(
    index, M_values, transformaed_CS_values, initial_timepoint
):
    if isinstance(index, tuple):
        index = [index]
    else:
        index = index.tolist()

    sub_Ms = M_values.loc[index]
    sub_Ms = sub_Ms.stack(level="Sample").stack(level="Timepoint").to_frame("M")
    CS_index = sub_Ms.index.droplevel(level="Timepoint")
    sub_CS = transformaed_CS_values.loc[CS_index, "Confidence_score"].tolist()

    sub_Ms["Confidence_score"] = sub_CS

    weighted_Ms = sub_Ms.groupby("Timepoint").apply(
        lambda sub_df: np.average(sub_df["M"], weights=sub_df["Confidence_score"])
    )

    confidence_score = sub_Ms.xs(initial_timepoint, level="Timepoint")[
        "Confidence_score"
    ]

    statistics = pd.Series()
    individual_insertions = confidence_score.shape[0]
    sum_confidence_score = confidence_score.sum()
    statistics.loc["Individual insertions"] = individual_insertions
    statistics.loc["Confidence score"] = sum_confidence_score

    merge_data = pd.concat(
        [weighted_Ms, statistics], axis=0, keys=["Weighted M", "Statistic"]
    )

    return merge_data
