
import warnings
from util.calculate_weight_averaged_M import calculate_weight_averaged_M
from scipy.optimize import curve_fit, minimize
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import statsmodels.api as sm
lowess = sm.nonparametric.lowess
# warnings.filterwarnings("error", category=RuntimeWarning)

# def the function for calculating the slope and the x-intercept for given two points


def slope_intercept(x1, y1, x2, y2):
    slope = (y2 - y1) / (x2 - x1)
    if slope == 0:
        x_intercept = 0
    else:
        x_intercept = x1 - y1 / slope
    return round(slope, 3), round(x_intercept, 3)


def sigmoid_GM_curve(G, DL, DR, LIM):

    try:
        M = LIM / (1 + np.exp(4*DR / LIM*(DL-G)+2))
    except Warning:
        M = np.array([0]*len(G))

    return M


def cauchy_loss_customed(params, G, Y_true):

    LIM, DR, DL = params
    y_pred = sigmoid_GM_curve(G, LIM, DR, DL)
    residuals = Y_true - y_pred

    return np.sum(np.log1p(residuals**2))


def curve_fitting(
    index,
    GM_df,
    fatol,
    useWeightedM=True,
):
    insertions_df = GM_df[GM_df.index.isin(index)].copy()
    insertions_df_no_0h = insertions_df[insertions_df["G"] > 0].copy()
    smoothed = lowess(insertions_df_no_0h["M"], insertions_df_no_0h["G"],
                      frac=0.3, return_sorted=False, it=10, is_sorted=False)
    insertions_df_no_0h["Smoothed_M"] = smoothed

    if useWeightedM:
        WG = GM_df.reset_index(drop=True)[["Sample", "Timepoint", "G"]].drop_duplicates().groupby("Timepoint").apply(
            lambda x: np.average(x["G"])).sort_values()
        WM = insertions_df_no_0h.groupby("Timepoint").apply(
            lambda x: np.average(
                x["Smoothed_M"], weights=x["CS"])
        )
        WG = WG[1:]

    try:
        init_guess = [0, 0, 0.1]
        xdata = np.array(WG)
        ydata = np.array(WM)
        result = minimize(cauchy_loss_customed, init_guess, args=(xdata, ydata), options={
                          "maxiter": 1000000, 'disp': False, 'fatol': fatol}, method="Nelder-Mead")

        sigmoid_DL, sigmoid_DR, sigmoid_LIM = result.x
        ypred = sigmoid_GM_curve(xdata, sigmoid_DL, sigmoid_DR, sigmoid_LIM)

        two_point_paris = zip(xdata[:-1], ypred[:-1], xdata[1:], ypred[1:])
        DR_DLs = [slope_intercept(x1, y1, x2, y2)
                  for x1, y1, x2, y2 in two_point_paris]
        # define the DR as the largest value of the slope and the DL as the same index of the largest value of the slope
        DR, DL = max(DR_DLs, key=lambda x: x[0])
        LIM = sigmoid_LIM

    except RuntimeError:
        sigmoid_DL, sigmoid_DR, sigmoid_LIM = np.nan, np.nan, np.nan
        DR, DL = np.nan, np.nan
        LIM = np.nan

    try:
        confidence_score = insertions_df.set_index(["Sample", "Timepoint"], append=True).unstack(
            "Sample")["CS"].reset_index(drop=True).drop_duplicates().sum(axis=1)

        individual_insertions = confidence_score.shape[0]
        sum_confidence_score = confidence_score.sum()
    except KeyError:
        confidence_score = pd.Series()
        individual_insertions = 0
        sum_confidence_score = 0

    statistics = pd.Series()
    statistics.loc["Individual insertions"] = individual_insertions
    statistics.loc["Confidence score"] = sum_confidence_score
    statistics.loc["DL"] = DL
    statistics.loc["DR"] = DR
    statistics.loc["LIM"] = LIM
    statistics.loc["sigmoid_DL"] = sigmoid_DL
    statistics.loc["sigmoid_DR"] = sigmoid_DR
    statistics.loc["sigmoid_LIM"] = sigmoid_LIM
    statistics.loc["WM_YES0"] = WM.iloc[0]
    statistics.loc["WM_YES1"] = WM.iloc[1]
    statistics.loc["WM_YES2"] = WM.iloc[2]
    statistics.loc["WM_YES3"] = WM.iloc[3]
    statistics.loc["WM_YES4"] = WM.iloc[4]

    return statistics
