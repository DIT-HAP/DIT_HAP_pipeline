
from util.calculate_weight_averaged_M import calculate_weight_averaged_M
from scipy.optimize import curve_fit, minimize
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import warnings
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

    # Check the condition and apply a penalty if it's met
    penalty = 0
    # avoid division by zero
    Final_M_noZero = np.abs(y_pred[-1])+1e12
    DR_Final_M = DR / Final_M_noZero

    # # if DR_Final_M > 2.5*np.abs(Final_M_noZero): # loss2
    # if (Final_M_noZero - DR_Final_M) < -0.5:  # loss3
    #     # or any other form of penalty you prefer
    #     penalty = np.abs(DR)  # + 1 / (np.abs(DL)+1e-10)
    #     # penalty = np.abs(DR) + 1 / (np.abs(DL)+1e-10)  # loss4

    # loss5
    # if (Final_M_noZero - DR_Final_M) < -0.5:
    #     penalty = np.abs(DR)

    # Add the penalty terms to the loss
    return np.sum(np.log1p(residuals**2)) + penalty

    # loss6


def curve_fitting(
    index,
    M_values,
    generation,
    transformaed_CS_values,
    initial_timepoint,
    xatol,
    useWeightedM=True,
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

    try:
        if useWeightedM:
            xdata = generations
            ydata = Weighted_M
        else:
            xdata = sub_Ms["G"].values
            ydata = sub_Ms["M"].values
        # popt, pcov, infodict, mesg, ier = curve_fit(MG_curve, xdata, ydata, method=cf_metohd, p0=init_guess, maxfev=100000, full_output=True)
        # popt, pcov, infodict, mesg, ier = curve_fit(MG_curve, xdata, ydata, method=cf_metohd, bounds=(0,[10, 3]), maxfev=100000, full_output=True, p0=init_guess)

        init_guess = [0, 0, 0.1]
        xdata = np.array(xdata)
        f = interp1d(xdata, ydata)
        inter_x = (xdata[:-1] + xdata[1:])/2
        xnew = np.append(xdata, inter_x)
        xnew = np.sort(xnew)
        ynew = f(xnew)
        # xnew = xdata
        # ynew = f(xnew)
        # xnew = np.array(list(range(1, 14, 1)))
        # ynew = f(xnew)
        # result = minimize(cauchy_loss_customed, init_guess, args=(xnew, ynew), bounds=[
        #                   (-1, 10), (-0.1, 2.5), (-1, 10)], options={"maxiter": 1000000, 'disp': False, 'fatol': xatol}, method="Nelder-Mead")
        result = minimize(cauchy_loss_customed, init_guess, args=(xnew, ynew), options={
                          "maxiter": 1000000, 'disp': False, 'fatol': xatol}, method="Nelder-Mead")

        sigmoid_DL, sigmoid_DR, sigmoid_LIM = result.x
        ypred = sigmoid_GM_curve(xnew, sigmoid_DL, sigmoid_DR, sigmoid_LIM)

        two_point_paris = zip(xnew[:-1], ynew[:-1], xnew[1:], ynew[1:])
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
    statistics.loc["sigmoid_DL"] = sigmoid_DL
    statistics.loc["sigmoid_DR"] = sigmoid_DR
    statistics.loc["sigmoid_LIM"] = sigmoid_LIM
    statistics.loc["WM_YES0"] = Weighted_M[0]
    statistics.loc["WM_YES1"] = Weighted_M[1]
    statistics.loc["WM_YES2"] = Weighted_M[2]
    statistics.loc["WM_YES3"] = Weighted_M[3]
    statistics.loc["WM_YES4"] = Weighted_M[4]

    return statistics
