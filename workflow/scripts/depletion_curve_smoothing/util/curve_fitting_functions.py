
from util.calculate_weight_averaged_M import calculate_weight_averaged_M
from scipy.optimize import curve_fit, minimize
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("error", category=RuntimeWarning)


def function_before_DL(G, DL, DR):
    return 0


def function_after_DL_V1(G, DL, DR):
    numerator = 2 ** G
    denominator = 2 ** (DL + (G - DL) * (1 - DR)) * 63 / 64 + 2 ** G / 64
    M = np.log2(numerator / denominator)
    return M


def function_after_DL_V2(G, DL, DR, LIM):
    numerator = 2 ** G
    denominator = 2 ** (DL + (G - DL) * (1 - DR)) * \
        (LIM-1) / LIM + 2 ** G / LIM
    M = np.log2(numerator / denominator)
    return M


def MG_curve(G, DL, DR):

    M = np.piecewise(
        G,
        [G <= DL, G > DL],
        [lambda x: function_before_DL(
            x, DL, DR), lambda x: function_after_DL_V1(x, DL, DR)],
        # [lambda x: function_before_DL(x, DL, DR, LIM), lambda x: function_after_DL_V2(x, DL, DR, LIM)],
    )

    return M


def sigmoid_GM_curve(G, DL, DR, LIM):

    try:
        M = LIM / (1 + np.exp(4*DR / LIM*(DL-G)+2))
    except Warning:
        M = np.array([0, 0, 0, 0, 0])

    return M

# def curve_fitting(
#     index, M_values, generation, transformaed_CS_values, initial_timepoint, useWeightedM=True, cf_metohd="trf"
# ):

    # if useWeightedM:
    #     try:
    #         WM = calculate_weight_averaged_M(
    #             index, M_values, transformaed_CS_values, initial_timepoint
    #             )
    #         Weighted_M = WM["Weighted M"].values[1:]
    #         generations = generation.values[1:]
    #     except KeyError:
    #         statistics = pd.Series()
    #         statistics.loc["Individual insertions"] = np.nan
    #         statistics.loc["Confidence score"] = np.nan
    #         statistics.loc["DL"] = np.nan
    #         statistics.loc["DR"] = np.nan
    #         statistics.loc["DL_pcov"] = np.nan
    #         statistics.loc["DR_pcov"] = np.nan
    #         statistics.loc["linalg_cond_pcov"] = np.nan
    #         statistics.loc["fvec"] = np.nan
    #         return statistics

    # if isinstance(index, tuple):
    #     index = [index]
    # else:
    #     index = index.tolist()

    # sub_Ms = M_values.loc[index]
    # sub_Ms = sub_Ms.stack(level="Sample").stack(level="Timepoint").to_frame("M")

    # generation_index = sub_Ms.index.get_level_values("Timepoint")
    # sub_generation = generation.loc[generation_index].values

    # CS_index = sub_Ms.index.droplevel(level="Timepoint")
    # sub_CS = transformaed_CS_values.loc[CS_index, "Confidence_score"].tolist()

    # sub_Ms["G"] = sub_generation
    # sub_Ms["Confidence_score"] = sub_CS

    # confidence_weights = 1/sub_Ms["Confidence_score"].values

    # init_guess = [0,0]

    # try:
    #     if useWeightedM:
    #         xdata = generations
    #         ydata = Weighted_M
    #     else:
    #         xdata = sub_Ms["G"].values
    #         ydata = sub_Ms["M"].values
    #     popt, pcov, infodict, mesg, ier = curve_fit(MG_curve, xdata, ydata, method=cf_metohd, p0=init_guess, maxfev=100000, full_output=True)
    #     # popt, pcov, infodict, mesg, ier = curve_fit(MG_curve, xdata, ydata, method=cf_metohd, bounds=(0,[10, 3]), maxfev=100000, full_output=True, p0=init_guess)
    #     # popt, pcov, infodict, mesg, ier = curve_fit(sigmoid_GM_curve, xdata, ydata, method=cf_metohd, maxfev=100000, full_output=True, p0=init_guess)

    #     DL, DR = popt
    #     DL_pcov = pcov[0,0]
    #     DR_pcov = pcov[1,1]
    #     linalg_cond_pcov = np.linalg.cond(pcov)
    #     fvec = ",".join(map(str, infodict["fvec"]))
    # except RuntimeError:
    #     DL, DR = np.nan, np.nan
    #     DL_pcov, DR_pcov = np.nan, np.nan
    #     linalg_cond_pcov = np.nan
    #     fvec = np.nan

    # try:
    #     confidence_score = sub_Ms.xs(initial_timepoint, level="Timepoint")[
    #         "Confidence_score"
    #     ]

    #     individual_insertions = confidence_score.shape[0]
    #     sum_confidence_score = confidence_score.sum()
    # except KeyError:
    #     confidence_score = pd.Series()
    #     individual_insertions = 0
    #     sum_confidence_score = 0

    # statistics = pd.Series()
    # statistics.loc["Individual insertions"] = individual_insertions
    # statistics.loc["Confidence score"] = sum_confidence_score
    # statistics.loc["DL"] = DL
    # statistics.loc["DR"] = DR
    # statistics.loc["DL_pcov"] = round(DL_pcov, 3)
    # statistics.loc["DR_pcov"] = round(DR_pcov, 3)
    # statistics.loc["linalg_cond_pcov"] = round(linalg_cond_pcov, 3)
    # statistics.loc["fvec"] = fvec
    # statistics.loc["WM_YES0"] = Weighted_M[0]
    # statistics.loc["WM_YES1"] = Weighted_M[1]
    # statistics.loc["WM_YES2"] = Weighted_M[2]
    # statistics.loc["WM_YES3"] = Weighted_M[3]
    # statistics.loc["WM_YES4"] = Weighted_M[4]

    # return statistics


def curve_fitting(
    index,
    M_values,
    generation,
    transformaed_CS_values,
    initial_timepoint,
    useWeightedM=True,
    cf_method="trf",
    xtol=1e-6,
):

    if useWeightedM:
        try:
            WM = calculate_weight_averaged_M(
                index, M_values, transformaed_CS_values, initial_timepoint
            )
            Weighted_M = WM["Weighted M"].values[1:]
            generations = generation.values[1:]
        except KeyError:
            statistics = pd.Series()
            statistics.loc["Individual insertions"] = np.nan
            statistics.loc["Confidence score"] = np.nan
            statistics.loc["DL"] = np.nan
            statistics.loc["DR"] = np.nan
            statistics.loc["DL_pcov"] = np.nan
            statistics.loc["DR_pcov"] = np.nan
            statistics.loc["linalg_cond_pcov"] = np.nan
            statistics.loc["fvec"] = np.nan
            return statistics

    if isinstance(index, tuple):
        index = [index]
    else:
        index = index.tolist()

    sub_Ms = M_values.loc[index]
    sub_Ms = sub_Ms.stack(level="Sample").stack(
        level="Timepoint").to_frame("M")

    generation_index = sub_Ms.index.get_level_values("Timepoint")
    sub_generation = generation.loc[generation_index].values

    CS_index = sub_Ms.index.droplevel(level="Timepoint")
    sub_CS = transformaed_CS_values.loc[CS_index, "Confidence_score"].tolist()

    sub_Ms["G"] = sub_generation
    sub_Ms["Confidence_score"] = sub_CS

    confidence_weights = 1/sub_Ms["Confidence_score"].values

    init_guess = [0, 0, 0.1]

    try:
        if useWeightedM:
            xdata = generations
            ydata = Weighted_M
        else:
            xdata = sub_Ms["G"].values
            ydata = sub_Ms["M"].values
        # popt, pcov, infodict, mesg, ier = curve_fit(MG_curve, xdata, ydata, method=cf_metohd, p0=init_guess, maxfev=100000, full_output=True)
        # popt, pcov, infodict, mesg, ier = curve_fit(MG_curve, xdata, ydata, method=cf_metohd, bounds=(0,[10, 3]), maxfev=100000, full_output=True, p0=init_guess)
        if cf_method != "lm":
            popt, pcov, infodict, mesg, ier = curve_fit(
                sigmoid_GM_curve, xdata, ydata, method=cf_method, maxfev=100000, full_output=True, bounds=([-1, -0.1, -1], [10, 2.5, 10]), loss="cauchy", p0=init_guess, ftol=xtol)
        else:
            popt, pcov, infodict, mesg, ier = curve_fit(
                sigmoid_GM_curve, xdata, ydata, method=cf_method, maxfev=100000, full_output=True, xtol=xtol)

        DL, DR, LIM = popt
        Final_M = sigmoid_GM_curve(13, DL, DR, LIM)
        # AUC = np.trapz(sigmoid_GM_curve(
        #     np.linspace(0, 13, 100), DL, DR, LIM), np.linspace(0, 13, 100))
        DL_pcov = pcov[0, 0]
        DR_pcov = pcov[1, 1]
        LIM_pcov = pcov[2, 2]
        linalg_cond_pcov = np.linalg.cond(pcov)
        fvec = ",".join(map(str, infodict["fvec"]))
    except RuntimeError:
        DL, DR, LIM, Final_M = np.nan, np.nan, np.nan, np.nan
        DL_pcov, DR_pcov, LIM_pcov = np.nan, np.nan, np.nan
        linalg_cond_pcov = np.nan
        fvec = np.nan

    try:
        confidence_score = sub_Ms.xs(initial_timepoint, level="Timepoint")[
            "Confidence_score"
        ]

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
    statistics.loc["Final_M"] = Final_M
    # statistics.loc["AUC"] = AUC
    statistics.loc["DL_pcov"] = round(DL_pcov, 3)
    statistics.loc["DR_pcov"] = round(DR_pcov, 3)
    statistics.loc["LIM_pcov"] = round(LIM_pcov, 3)
    statistics.loc["linalg_cond_pcov"] = round(linalg_cond_pcov, 3)
    statistics.loc["fvec"] = fvec
    statistics.loc["WM_YES0"] = Weighted_M[0]
    statistics.loc["WM_YES1"] = Weighted_M[1]
    statistics.loc["WM_YES2"] = Weighted_M[2]
    statistics.loc["WM_YES3"] = Weighted_M[3]
    statistics.loc["WM_YES4"] = Weighted_M[4]

    return statistics
