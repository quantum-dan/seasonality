from rtseason import *
import pandas as pd


def get_results():
    return pd.read_csv("fit.txt")


def get_data():
    """
    Data were generated in R using dataRetrieval:
        readNWISdv("07086000", "00010", "2010-01-01", "2023-12-31") %>%
        select(Date, temperature=X_00010_00003) %>% drop_na() %>%
        mutate(day = as.integer(format(Date, "%j")))
    """
    return pd.read_csv("data.txt")


def test_fit(data, fit):
    """
    The fit comparison includes performance information from make_sints,
    so it tests both aspects.
    """
    new_fit = fit_sins(data)
    print("Fit comparison (new then reference):")
    print(new_fit)
    print(fit)

