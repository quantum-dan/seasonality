# Stream Seasonal Thermal Regimes

This repository contains reusable code (R script and Python package) for analyzing stream seasonal thermal regimes, in the form of an annual temperature cycle function, using the ``three-sine'' approach.  Notes on research use and reproducibility (for the associated manuscript) and citation information are included at the end of the README, in the [Research Data and Use](#Research-Data-and-Use) section.  This upper portion deals with the use of the R script and the Python package.

## Overview

The stream seasonal thermal regime - i.e., long-term seasonal behavior of stream temperatures - can be described using the annual temperature cycle, or the general pattern of stream temperatures throughout the year.  As a function, mean temperature = f(day of year).  The three-sine approach (Philippus et al., 2024; see [Citation](#Citation)) is a general function that captures the annual temperature cycle with very high accuracy for streams across the United States, including high-elevation/snow-dominated watersheds, arid regions, etc.  Thus, detailed stream seasonal thermal behavior can be summarized using eight coefficients, from which almost any seasonal summary statistic of interest can be extracted (by recreating the annual temperature cycle in full, then computing the summary statistic).  This repository contains Python and R implementations for fitting three-sine coefficients to data from a given stream as well as recreating an annual temperature cycle timeseries from the coefficients.

## Usage

This applies to both Python and R.  There are two key functions for the end user: `fit.sins` (R) / `fit_sins` (Python) and `make.sints` (R) / `make_sints` (Python).

The `fit.sins` function accepts as its main argument a data frame containing columns for `day` (day of year, as an integer) and `temperature` (observed temperature).  This data frame may be day-of-year means already, but it can also be an observed timeseries with multiple entries per day (e.g., by date over several years).  The function will compute day-of-year means.  Observations must be available for at least 180 days.  The R version has several other arguments that were used for generating data for the paper but can be left as default by the end user.

The result of `fit.sins` is a data frame containing the eight three-sine coefficients as well as fit statistics.  The coefficients are labeled Intercept (annual mean), Amplitude (main cosine amplitude), FallWinter, SpringSummer (autumn-winter and spring-summer coefficients), FallDay, WinterDay, SpringDay, and SummerDay (days-of-year of anomaly peaks).  The last two columns are R2 and RMSE for the three-sine fit versus the observed day-of-year means.  In Python, the optional argument `return.object` can be set to `True` to return, instead of a data frame, a three-sine fit object with the above results as properties and a `generate_ts` method for generating a timeseries.

The `make.sints` function accepts a three-sine coefficient data frame (R2/RMSE not required) as an argument and returns a data frame containing `day` (day of year) and `actemp` (annual cycle temperature, or day-of-year mean temperature) columns.

In Python, the actual three-sine functionality is implemented through the `ThreeSine` class with fitting (`ThreeSine.from_data()`) and timeseries-generation (`.generate_ts`) methods, but the above functions are provided so that the interface is (roughly) the same as the R version.

# Research Data and Use

## Reproducibility

This repository contains analysis code  used to develop an analysis of stream seasonal thermal regimes in the United States.  The R Notebook, `analysis.Rmd`, will reproduce all supporting analysis and generate figures when run with the data files downloaded into a `Data` subdirectory of the working directory.  `analysis.Rmd` calls more involved functions from `functions.R`.  The knitted version of the notebook is also included as `analysis.pdf`.

To run the Notebook, in the working directory there must be a `Data` directory containing the data files and a `Figures` directory with a `MovingWindowSeasons` subdirectory, where figures will be stored.  For map generation, `functions.R` will also look for an EPA Level I Ecoregions (https://gaftp.epa.gov/EPADataCommons/ORD/Ecoregions/cec_na/na_cec_eco_l1.zip) shapefile in an `Ecoregions` directory.

The required data can be downloaded from [Hydroshare](http://www.hydroshare.org/resource/7d960b7fdfee480895fd845bade1b75a).

## Citation

The three-sine approach implemented here is introduced and analyzed in the paper ``Improved Annual Temperature Cycle Function for Stream Seasonal Thermal Regimes'' (Philippus, Corona and Hogue 2024), currently in review.  Citation information will be added upon acceptance.
