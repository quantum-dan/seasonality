"""
This package includes functions to fit three-sine coefficients and generate
timeseries.  The core functionality is written through a class, with top-level
function wrappers to replicate the R interface.
"""

import pandas as pd
import numpy as np


def apply_sin(days, start, domain):
    """
    Apply a truncated sine function, starting from start, within domain
    days: array of days to apply to
    start: integer start day
    domain: array of integer days containing domain
    """
    return np.sin(((days - start) % 365) * 2 * np.pi / len(domain)) * \
        np.array([day in domain for day in days])


def make_trunc_sins(days, spring, summer, fall, winter):
    """
    Build reference truncated sine functions (sin1 = fall-winter, sin2=
                                              spring-summer)
    over days.
    """
    sin1width = round(((winter - fall) % 365)/2)
    sin1dom = np.concatenate((
        np.arange(fall - sin1width, 367, 1),
        np.arange(1, winter + sin1width + 1, 1)))

    sin2width = round((summer - spring)/2)
    sin2dom = np.arange(spring - sin2width,
                        summer + sin2width + 1,
                        1)

    sin1 = -apply_sin(days, fall - sin1width, sin1dom)
    sin2 = -apply_sin(days, spring - sin2width, sin2dom)
    return (sin1, sin2)


def rollmean_circ(x, n):
    """
    Compute a centered, circular rolling mean over x, using n terms.
    We could be a bit smarter here, but over-optimization and all that; this
    will be used for fairly small arrays.
    So the simple solution is: stick x end-to-end three times, run a centered
    rolling mean over the whole thing using a convolution, and then take the
    middle third.
    If n>len(x), that won't work as expected.
    """
    x3 = np.concatenate((x, x, x))
    xconv = np.convolve(x3, np.ones(n)/n, "same")
    return xconv[len(x):(2*len(x))]


def get_anom_date(day, anom, direction, default, dmin=-1, dmax=999):
    """
    Identify the peak date of the anomaly in question, over the dataset
    day and anom (must be the same length).
    Date is bounded between dmin, dmax.
    Direction is either np.min (negative anomaly) or np.max (positive anomaly).
    """
    if len(day) != len(anom):
        raise ValueError("Days and anomalies must have the same length")
    # First, identify eligible days
    el = (day >= dmin) & (day <= dmax)
    eld = day[el]
    # Now get associated anomalies
    ely = anom[el]
    # Now pick the day with the strongest anomaly (or the mean of
    # candidate days)
    return np.mean(eld[ely == direction(ely)]) \
        if len(eld) > 0 else default


class ThreeSine(object):
    def __init__(self,
                 Intercept, Amplitude, FallDay, WinterDay, FallWinter,
                 SpringDay, SummerDay, SpringSummer,
                 R2=None, RMSE=None):
        self.Intercept = Intercept
        self.Amplitude = Amplitude
        self.FallDay = FallDay
        self.WinterDay = WinterDay
        self.FallWinter = FallWinter
        self.SpringDay = SpringDay
        self.SummerDay = SummerDay
        self.SpringSummer = SpringSummer
        self.R2 = R2
        self.RMSE = RMSE

    def from_coefs(coefs):
        # Initialize from a coefficient data frame - everything is just
        # from the first row
        return ThreeSine(
            Intercept=coefs["Intercept"].iloc[0],
            Amplitude=coefs["Amplitude"].iloc[0],
            FallDay=coefs["FallDay"].iloc[0],
            WinterDay=coefs["WinterDay"].iloc[0],
            FallWinter=coefs["FallWinter"].iloc[0],
            SpringDay=coefs["SpringDay"].iloc[0],
            SummerDay=coefs["SummerDay"].iloc[0],
            SpringSummer=coefs["SpringSummer"].iloc[0],
            # Performance *if available*
            R2=coefs["R2"].iloc[0] if "R2" in coefs else None,
            RMSE=coefs["RMSE"].iloc[0] if "RMSE" in coefs else None
        )

    def from_data(data):
        inp = data.dropna().groupby("day").mean("temperature").\
            sort_values("day")
        if inp.shape[0] < 180:
            raise ValueError("Insufficient data coverage for 3-sine fit; >=180 days required")
        day = inp.index.to_numpy()
        Ts = inp["temperature"].to_numpy()
        index = np.cos((day - 210) * 2 * np.pi / 365)
        # Least squares: lstsq(a,b) is equivalent to lm(b~a)
        # index matrix - with a ones column
        imat = np.array([np.ones(len(index)), index]).transpose()
        [meant, amplitude] = np.linalg.lstsq(imat, Ts, None)[0]
        amp_norm = amplitude / meant  # amplitude as fraction of mean
        # predicted cosine values
        cospr = 1 + index * amp_norm
        # anomaly as a fraction of mean temp
        anomaly = Ts/meant - cospr
        # Use a rolling mean (centered) to find zeroes and peaks
        rolled = rollmean_circ(anomaly, 31)
        # Now we start computing seasonal terms; this form applies for each day
        fallt = get_anom_date(day, rolled, np.min, 330, dmin=300)
        wint = get_anom_date(day, rolled, np.max, 80, dmax=110)
        spt = get_anom_date(day, rolled, np.min, 160, dmin=120, dmax=180)
        sumt = get_anom_date(day, rolled, np.max, 220, dmin=200, dmax=240)
        # Next, compute width and domain of sine functions
        (sin1, sin2) = make_trunc_sins(day, spt, sumt, fallt, wint)
        # Now create the full best fit
        y_anom = anomaly * meant  # actual (C-units) anomaly
        xvars = np.array([sin1, sin2]).transpose()
        [fw, ssu] = np.linalg.lstsq(xvars, y_anom, None)[0]
        # Prepare object
        tsfit = ThreeSine(meant, amplitude, fallt, wint, fw, spt, sumt, ssu)
        # Get performance statistics
        pred = tsfit.generate_ts()
        comb = pred.join(inp, on="day")
        r2 = np.corrcoef(comb["actemp"], comb["temperature"])[0,1]**2
        rmse = np.sqrt(np.mean((comb["actemp"] - comb["temperature"])**2))
        tsfit.R2 = r2
        tsfit.RMSE = rmse
        return tsfit
        

    def to_df(self):
        return pd.DataFrame({
            "Intercept": self.Intercept,
            "Amplitude": self.Amplitude,
            "FallDay": self.FallDay,
            "WinterDay": self.WinterDay,
            "SpringDay": self.SpringDay,
            "SummerDay": self.SummerDay,
            "SpringSummer": self.SpringSummer,
            "FallWinter": self.FallWinter,
            "R2": self.R2,
            "RMSE": self.RMSE
            },
            index=pd.Series([0]))

    def generate_ts(self):
        days = np.arange(1, 367, 1, dtype="int")
        index = np.cos((days - 210) * 2 * np.pi / 365)

        # sin1width = round(((self.WinterDay - self.FallDay) % 365)/2)
        # sin1dom = np.concatenate((
        #     np.arange(self.FallDay - sin1width, 367, 1),
        #     np.arange(1, self.WinterDay + sin1width + 1, 1)))

        # sin2width = round((self.SummerDay - self.SpringDay)/2)
        # sin2dom = np.arange(self.SpringDay - sin2width,
        #                     self.SummerDay + sin2width + 1,
        #                     1)

        # sin1 = -apply_sin(days, self.FallDay - sin1width, sin1dom)
        # sin2 = -apply_sin(days, self.SpringDay - sin2width, sin2dom)
        (sin1, sin2) = make_trunc_sins(days, self.SpringDay,
                                       self.SummerDay,
                                       self.FallDay,
                                       self.WinterDay)

        series = self.Intercept + self.Amplitude * index + \
            self.FallWinter * sin1 + self.SpringSummer * sin2
        return pd.DataFrame({"day": days, "actemp": series})


def fit_sins(data):
    return ThreeSine.from_data(data).to_df()


def make_sints(coefs):
    return ThreeSine.from_coefs(coefs).generate_ts()
