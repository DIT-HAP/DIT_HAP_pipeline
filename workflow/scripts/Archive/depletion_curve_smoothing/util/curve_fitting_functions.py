
import matplotlib.pyplot as plt
import warnings
from util.calculate_weight_averaged_M import calculate_weight_averaged_M
from scipy.optimize import curve_fit, minimize
import pandas as pd
from pathlib import Path
import numpy as np
from scipy.interpolate import interp1d
import statsmodels.api as sm
lowess = sm.nonparametric.lowess
warnings.filterwarnings("error", category=RuntimeWarning)

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


def smoothing_comparison_plot(insertions_df_no_0h, smoothed, level_id, generation, fitting, prediction):

    insertions_df_no_0h["Smoothed_M"] = smoothed
    if type(level_id) == str:
        level_id = [level_id]
    else:
        level_id = list(level_id)
    id = "_".join([str(i) for i in level_id])
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.scatter(insertions_df_no_0h["G"], insertions_df_no_0h["M"],
               edgecolor="blue", s=insertions_df_no_0h["CS"], alpha=0.5, facecolor="none", label="Before LOWESS")
    ax.scatter(insertions_df_no_0h["G"],
               smoothed, edgecolor="red", s=insertions_df_no_0h["CS"], alpha=0.5, facecolor="none", label="After LOWESS")
    ax.set_ylim(-3, 10)
    # create a random string for the name of the plot
    ax.legend()
    ax.set_xlabel("Generation")
    ax.set_ylabel("M")

    Gs = generation[1:]
    Ms = insertions_df_no_0h.groupby("Timepoint").apply(
        lambda x: np.average(
            x["M"], weights=x["CS"])
    )
    smoothed_Ms = insertions_df_no_0h.groupby("Timepoint").apply(
        lambda x: np.average(
            x["Smoothed_M"], weights=x["CS"])
    )
    ax.plot(Gs, Ms, label="Weighted M Before LOWESS", color="blue", alpha=0.3)
    ax.plot(Gs, smoothed_Ms, label="Weighted M After LOWESS",
            color="red", alpha=0.3)
    output_dir = Path(f"tmp/{fitting}_{prediction}")
    output_dir.mkdir(exist_ok=True, parents=True)
    fig.savefig(output_dir/f"{fitting}_{prediction}_{id}.png")
    plt.close()


def curve_fitting(
    level_id,
    generation,
    insertions_in_current_level,
    fatol,
    useWeightedM=True,
    useLOWESS=True,
    fitting="points5",
    prediction="points5"
):

    insertions_df_no_0h = insertions_in_current_level[insertions_in_current_level["G"] > 0].copy(
    )

    if useLOWESS:
        smoothed = lowess(insertions_df_no_0h["M"], insertions_df_no_0h["G"],
                          frac=0.3, return_sorted=False, it=10, is_sorted=False)
        smoothing_comparison_plot(
            insertions_df_no_0h, smoothed, level_id, generation, fitting, prediction)
        insertions_df_no_0h["M"] = smoothed

    if useWeightedM:
        Gs = generation[1:]
        Ms = insertions_df_no_0h.groupby("Timepoint").apply(
            lambda x: np.average(
                x["M"], weights=x["CS"])
        )
    else:
        Gs = insertions_df_no_0h["G"]
        Ms = insertions_df_no_0h["M"]
    try:
        init_guess = [0, 0, 0.1]
        xdata = np.array(Gs)
        ydata = np.array(Ms)

        f = interp1d(xdata, ydata)
        inter_x = (xdata[:-1] + xdata[1:])/2

        if fitting == "points5":
            xnew = xdata
            ynew = ydata
        elif fitting == "points7":
            xnew = np.append(xdata, [inter_x[0], inter_x[-1]])
            xnew = np.sort(xnew)
            ynew = f(xnew)

        result = minimize(cauchy_loss_customed, init_guess, args=(xnew, ynew), options={
                          "maxiter": 1000000, 'disp': False, 'fatol': fatol}, method="Nelder-Mead")

        sigmoid_DL, sigmoid_DR, sigmoid_LIM = result.x

        if prediction == "points5":
            xpred = xdata
        elif prediction == "points9":
            xpred = np.append(xdata, inter_x)
            xpred = np.sort(xpred)
        elif prediction == "points9g":
            xpred = np.arange(1, 14, 2)
        elif prediction == "points13g":
            xpred = np.arange(1, 14, 1)

        ypred = sigmoid_GM_curve(
            xpred, sigmoid_DL, sigmoid_DR, sigmoid_LIM)

        two_point_paris = zip(
            xpred[:-1], ypred[:-1], xpred[1:], ypred[1:])
        DR_DLs = [slope_intercept(x1, y1, x2, y2)
                  for x1, y1, x2, y2 in two_point_paris]
        # define the DR as the largest value of the slope and the DL as the same index of the largest value of the slope
        DR, DL = max(DR_DLs, key=lambda x: np.abs(x[0]))
        LIM = sigmoid_GM_curve(
            np.array([13.5]), sigmoid_DL, sigmoid_DR, sigmoid_LIM)[0]

    except RuntimeError:
        sigmoid_DL, sigmoid_DR, sigmoid_LIM = np.nan, np.nan, np.nan
        DR, DL = np.nan, np.nan
        LIM = np.nan

    try:
        confidence_score = insertions_in_current_level.set_index(["Sample", "Timepoint"], append=True).unstack(
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
    statistics.loc["WM_YES0"] = Ms.iloc[0]
    statistics.loc["WM_YES1"] = Ms.iloc[1]
    statistics.loc["WM_YES2"] = Ms.iloc[2]
    statistics.loc["WM_YES3"] = Ms.iloc[3]
    statistics.loc["WM_YES4"] = Ms.iloc[4]

    return statistics
