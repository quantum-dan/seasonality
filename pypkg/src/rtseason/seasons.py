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
        # TODO: RESUME HERE ...

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
            })

    def generate_ts(self):
        days = np.arange(1, 367, 1, dtype="int")
        index = np.cos((days - 210) * 2 * np.pi / 365)

        sin1width = round(((self.WinterDay - self.FallDay) % 365)/2)
        sin1dom = np.concatenate((
            np.arange(self.FallDay - sin1width, 367, 1),
            np.arange(1, self.WinterDay + sin1width + 1, 1)))

        sin2width = round((self.SummerDay - self.SpringDay)/2)
        sin2dom = np.arange(self.SpringDay - sin2width,
                            self.SummerDay + sin2width + 1,
                            1)

        sin1 = -apply_sin(days, self.FallDay - sin1width, sin1dom)
        sin2 = -apply_sin(days, self.SpringDay - sin2width, sin2dom)

        series = self.Intercept + self.Amplitude * index + \
            self.FallWinter * sin1 + self.SpringSummer * sin2
        return pd.DataFrame({"day": days, "actemp": series})
