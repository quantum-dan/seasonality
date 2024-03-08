library(dataRetrieval)
library(tidyverse)
library(R.utils)
library(sf)
source("util.R")

renamer <- c(
  "Intercept" = "Sinusoid Mean (C)",
  "Amplitude" = "Sinusoid Amplitude",
  "temperature" = "Mean Temperature (C)",
  "FallDay" = "Day of Peak Autumn Low",
  "WinterDay" = "Day of Peak Winter High",
  "FallWinter" = "Autumn/Winter\nAnomaly",
  "SpringDay" = "Day of Peak Spring Low",
  "SummerDay" = "Day of Peak Summer High",
  "SpringSummer" = "Spring/Summer\nAnomaly",
  "R2" = "Model R2",
  "N" = "Gage Count"
)
renamer.nn <- c(
  "Intercept" = "Sinusoid Mean (C)",
  "Amplitude" = "Sinusoid Amplitude (C)",
  "temperature" = "Mean Temperature (C)",
  "FallDay" = "Day of Peak Autumn Low",
  "WinterDay" = "Day of Peak Winter High",
  "FallWinter" = "Autumn/Winter Anomaly (C)",
  "SpringDay" = "Day of Peak Spring Low",
  "SummerDay" = "Day of Peak Summer High",
  "SpringSummer" = "Spring/Summer Anomaly (C)",
  "R2" = "Model R2",
  "N" = "Gage Count"
)
renamer.am <- c(
  "amplitude" = "Sinusoid Amplitude",
  "temperature" = "Mean Temperature (C)",
  "FY" = "Autumn Low Magnitude",
  "SpX" = "Day of Maximum Spring Low",
  "SpY" = "Spring Low Magnitude",
  "SuY" = "Summer High Magnitude",
  "WX" = "Day of Maximum Winter High",
  "WY" = "Winter High Magnitude"
)

applysin <- function(days, start, domain) {
  # Apply one of the truncated sine functions, within the domain and offset
  # from start
  sin(((days - start) %% 365) * (2*pi)/length(domain)) * (days %in% domain)
}

freezeit <- function(temperature) {
  temperature[temperature < 0] <- 0
  temperature
}


fit.sins <- function(data, mod=FALSE, aic=FALSE, fixed=FALSE,
                     erv=tibble(
                       Intercept = NA,
                       Amplitude = NA,
                       FallDay = NA,
                       WinterDay = NA,
                       FallWinter = NA,
                       SpringDay = NA,
                       SummerDay = NA,
                       SpringSummer = NA,
                       R2 = NA
                     )) {
  # Fit across all data.  Use three-sinusoid approach.  Data must have day
  # and temperature.
  # mod = return model instead of tibble, unless aic
  # is true, in which case return tibble(sin, anomaly) with AICs
  # Fixed: force summer day to 220, fall day to 330, and spring day to 160
  inp <- data %>%
    group_by(day) %>%
    summarize(temperature = mean(temperature, na.rm=T)) %>%
    drop_na() %>%
    arrange(day)
  
  if (nrow(inp) < 180) {
    warning("Insufficient data coverage for 3-sine fit")
    if (mod) {
      if (aic) {
        tibble(Sin = NA, Anom = NA)
      } else NA
    } else erv
  } else {
    day <- inp$day
    Ts <- inp$temperature
    
    # First, fit the raw cosine.
    index <- cos((day - 210) * 2 * pi / 365)
    cosfit <- lm(Ts ~ index)$coefficients
    meant <- cosfit[[1]]
    amplitude <- cosfit[[2]] / meant
    cospr <- 1 + index * amplitude
    anomaly <- Ts/meant - cospr
    
    # Now identify coordinates for sines.
    # Use rolling mean to find zeros and peaks.
    rolled <- stats::filter(anomaly, rep(1/31, 31), circular=T)
    
    falld <- day[day >= 300]
    fally <- rolled[day %in% falld]
    fallt <- if (length(falld) > 0 & !fixed) mean(falld[fally == min(fally)][1]) else 330
    
    wind <- day[day <= 110]
    winy <- rolled[day %in% wind]
    wint <- if (length(wind) > 0) mean(wind[winy == max(winy)]) else 80
    
    spd <- day[day >= 120 & day <= 180]
    spy <- rolled[day %in% spd]
    spt <- if (length(spd) > 0 & !fixed) mean(spd[spy == min(spy)]) else 160
    
    sumd <- day[day >= 200 & day <= 240]
    sumy <- rolled[day %in% sumd]
    sumt <- if (length(sumd) > 0 & !fixed) mean(sumd[sumy == max(sumy)]) else 220
    
    # Define sine functions
    # Domain = half of the gap between peaks before and after
    sin1width <- round(((wint - fallt) %% 365)/2)
    sin1dom <- c((fallt - sin1width):366, 1:(wint + sin1width))
    
    sin2width <- round((sumt - spt)/2)
    sin2dom <- (spt - sin2width):(sumt + sin2width)
    
    sin1 <- -applysin(day, fallt - sin1width, sin1dom)
    sin2 <- -applysin(day, spt - sin2width, sin2dom)
    
    fullfit <- lm((Ts - meant * cospr) ~ 0 + sin1 + sin2)
    co <- fullfit$coefficients
    
    if (mod) {
      # Return lm for fit comparison (AIC)
      if (aic) {
        tibble(
          Sin = 2*(2 - as.numeric(logLik(lm(Ts ~ index)))),
          Anom = 2*(8 - as.numeric(logLik(fullfit)))
        )
      } else fullfit
    } else {
      tibble(
        Intercept = meant,
        Amplitude = amplitude * meant,
        FallDay = fallt,
        WinterDay = wint,
        FallWinter = if (!is.na(co[[1]])) co[[1]] else 0,
        SpringDay = spt,
        SummerDay = sumt,
        SpringSummer = co[[2]],
        R2 = cor(Ts, (fullfit$fitted.values + meant * cospr))^2
        , RMSE = sqrt(mean(((fullfit$fitted.values + meant * cospr) - Ts)^2))
      )
    }
  }
}

fourier <- function(temps, N) {
  M <- length(temps)
  ft <- fft(temps)
  # Discrete Fourier Transform: the "far end" are the complex conjugates of the
  # "near end" and both ends are needed to reproduce the timeseries.  Only one
  # end would need to be stored (complex conjugates, again), but with a real and
  # complex part it still corresponds to 2N numbers.
  ft <- c(ft[1:N], rep(0, M - 2*N), ft[(M-N+1):M])
  ts <- Re(fft(ft, inverse=TRUE)) / M
  ftr <- ft[1:N]
  R2 <- cor(temps, ts)^2
  RMSE <- sqrt(mean((ts - temps)^2))
  names(ftr) <- paste0("Fourier", 1:N)
  as_tibble_row(c(Mod(ftr), "R2" = R2, "RMSE" = RMSE))
}

normalize.sins <- function(data) {
  # Normalize sin fit data to temperature.  Assume data has sin-fit columns
  # and temperature, and is grouped.
  data %>%
    mutate(
      across(c(Amplitude, FallWinter, SpringSummer),
             ~.x / mean(Intercept))
    )
}

make.sints <- function(sinfit) {
  # Sin fit --> tibble(day, actemp)
  # sinfit: tibble(Intercept, Amplitude, FallDay, WinterDay, FallWinter,
  # SpringDay, SummerDay, SpringSummer, R2)
  # Coefficients fit the actual temperature, not the normalized temperature.
  
  day <- 1:366
  
  index <- cos((day - 210) * 2 * pi / 365)
  
  # Prepare sin functions
  wint <- sinfit$WinterDay
  fallt <- sinfit$FallDay
  sumt <- sinfit$SummerDay
  spt <- sinfit$SpringDay
  
  sin1width <- round(((wint - fallt) %% 365)/2)
  sin1dom <- c((fallt - sin1width):366, 1:(wint + sin1width))
  
  sin2width <- round((sumt - spt)/2)
  sin2dom <- (spt - sin2width):(sumt + sin2width)
  
  sin1 <- -applysin(day, fallt - sin1width, sin1dom)
  sin2 <- -applysin(day, spt - sin2width, sin2dom)
  
  # Compute fit
  tibble(
    day = day,
    actemp = sinfit$Intercept + sinfit$Amplitude * index +
      sinfit$FallWinter * sin1 + sinfit$SpringSummer * sin2
  )
}

add.sinusoid.all <- function(data) {
  # add sinusoid, no ID
  data %>%
    mutate(
      Index = cos((day - 210) * 2 * pi / 365),
      ntemp = temperature / mean(temperature)
    ) %>%
    group_by(day) %>%
    mutate(
      cycle.temp = median(ntemp, na.rm=T)
    ) %>%
    ungroup() %>%
    mutate(
      amplitude = lm(ntemp ~ Index)$coefficients[[2]],
      reftemp = freezeit(1 + Index * amplitude),
      anomaly = cycle.temp - reftemp
    )
}

add.sinusoid <- function(data) {
  data %>%
    group_by(id) %>%
    group_modify(~add.sinusoid.all(.)) %>%
    ungroup()
}

add.resfit <- function(data) {
  # This goes into group_modify and similar.
  # Requires day, reftemp, amplitude, and anomaly (observed residuals).
  # Fits residual prediction to median of anomaly across days.
  #
  # Adds:
  # MedAnom (median anomaly),
  # ResMod (modeled residual),
  # RefMod (modeled reference temperature = reftemp + amplitude * ResMod)
  # And fit characteristics, WX, WY, SpX, SpY, SuY, FY
  inp <- data %>%
    group_by(day) %>%
    summarize(
      across(c(reftemp, amplitude, anomaly),
             ~median(., na.rm=T))
    ) %>%
    ungroup() %>%
    arrange(day)
  fit <- resfit(inp$day, inp$anomaly)
  ResMod <- resinterp(inp$day, fit)
  newdata <- left_join(
    tibble(
      day = inp$day,
      MedAnom = inp$anomaly
    ),
    tibble(
      day = 1:length(ResMod),
      ResMod = ResMod
    ),
    by="day"
  ) %>% 
    mutate(RefMod = inp$reftemp + ResMod * inp$amplitude) %>%
    cbind(fit)
  left_join(data,
            newdata,
            by="day")
}

sinu.resfit <- function(data) {
  # Takes original data; adds sinusoid and residual fit, by gage.
  data %>%
    add.sinusoid %>%
    group_by(id) %>%
    group_modify(~add.resfit(.)) %>%
    mutate(Modeled = RefMod * mean(temperature)) %>%
    ungroup()
}


bucketize <- function(data, bvar, unit, size) {
  # data must have id column.  Grouped by bvar (character).  Labeled
  # with unit.  Size should divide 100.
  # bvar should be constant across id.
  ref <- (data %>% group_by(id) %>% slice_head(n=1))[[bvar]]
  maxes <- c(quantile(ref, (1/size)*(1:(size-1))), max(ref) * 2)
  
  buckets <- tibble(
    Min=lag(maxes, default = 0),
    Max=maxes,
    Bucket=(100/size)*(0:(size-1)),
    BucketName=paste0(
      str_pad((100/size)*(0:(size-1)), width = 2, side = "left", pad = "0"),
      "% (",
      as.integer(Min),
      "+ ",
      unit,
      ")"
    )
  )
  
  data %>%
    group_by(id) %>%
    group_modify(~{
      mutate(.x,
             Bucket = filter(buckets,
                             Min <= first(.x[[bvar]]),
                             Max > first(.x[[bvar]]))$Bucket[1],
             BucketName = filter(buckets,
                                 Min <= first(.x[[bvar]]),
                                 Max > first(.x[[bvar]]))$BucketName[1])
    }) %>%
    drop_na
}


resfit <- function(days, res) {
  # Days and res should be in order.  res --> residuals.
  # Returns a tibble: {WX, WY, SpX, SpY, SuY, FY}
  # Late winter high; can assume it is just the first instance
  lwh <- max(res[days <= 120])
  lwh_x <- days[res == lwh][1]
  # Late spring low; have to limit it to spring days only because winter may have lows
  dsp <- days[days > 120 & days <= 180]
  rsp <- res[days > 120 & days <= 180]
  lsl <- min(rsp)
  lsl_x <- dsp[rsp == lsl]
  # Summer high; fixed X range
  suh <- mean(res[days >= 200 & days <= 240])
  # Late fall low; fixed X, but take a mean just in case
  lfl <- mean(res[days >= 325 & days <= 335])
  
  tibble(
    WX=lwh_x,
    WY=lwh,
    SpX=lsl_x,
    SpY=lsl,
    SuY=suh,
    FY=lfl
  )
}

resinterp <- function(days, fit) {
  # returns the fitted series from fit.
  results <- tryCatch(
    c(
      # Interpolate LFL to LWH
      # 35 --> LWH to Day 0
      fit$FY + (1:fit$WX + 35) / (fit$WX + 35) * (fit$WY - fit$FY),
      # Interpolate LWH to LSL
      fit$WY + ((fit$WX + 1):fit$SpX - fit$WX) / (fit$SpX - fit$WX) * (fit$SpY - fit$WY),
      # Interpolate LSL to SuH (200)
      fit$SpY + ((fit$SpX + 1):200 - fit$SpX) / (200 - fit$SpX) * (fit$SuY - fit$SpY),
      # Summer high
      rep(fit$SuY, 40),
      # SuH to LFL (330)
      fit$SuY + (241:330 - 240) / (330-240) * (fit$FY - fit$SuY),
      # And LFL to LWH - length(days) accounts for leap years etc
      fit$FY + (331:length(days) - 330) / (fit$WX + 35) * (fit$WY - fit$FY)
    ),
    error = function(e) NA)
  results
}


el.splot <- function(data) {
  # seasonality plot for elevation buckets
  plt <- data %>%
    group_by(id) %>%
    mutate(temperature = temperature / mean(temperature)) %>%
    ungroup() %>%
    bucketize("elevation", "m", 10) %>%
    group_by(BucketName, day) %>%
    summarize(Med = median(temperature, na.rm=T)) %>%
    group_by(BucketName) %>%
    arrange(day) %>%
    mutate(
      MedRoll = stats::filter(Med, rep(1/30, 30), circular=T)
    ) %>%
    ggplot(aes(day, MedRoll, color=BucketName)) +
    geom_line(size=2) +
    scale_x_continuous(breaks=30*(0:12)) +
    labs(x="Day of Year",
         y="Median Normalized\nSeasonal Temperature\n(30-Day Mean)",
         color="Elevation\nBucket") +
    scale_color_viridis_d() +
    theme_std() +
    theme(legend.position = "bottom")

  write.img(plt, "Figures/SeasonalBucketsSmooth")
}


##### USGS retrieval functions

get.usgs <- function(id, start, end) {
  # returns data frame with: id, year, time, temperature
  # year and temperature as dbl, id and time as chr
  
  tryCatch(
    withTimeout({
      dv <- 
        readNWISdv(id, "00010", start, end)
      if (nrow(dv) > 0) {
        dv %>%
          select(id=site_no,
                 Date,
                 temperature = X_00010_00003) %>%
          mutate(
            day = as.integer(format(Date, "%j"))
          ) %>%
          group_by(id, Date, day) %>%
          summarize(temperature = mean(temperature, na.rm=T))
      } else {
        tibble(id=id, Date=NA, day=NA, temperature=NA)
      }
    }
    ,
    onTimeout = "error",
    timeout = 600*10^3)
    , error = function(e) {
      print(e)
      # print(paste("Timed out/failed", id))
      write(id, file="failed.txt", append=TRUE)
      tibble(id=id, Date=NA, day=NA, temperature=NA)
    }
  )
}

get.gageinfo <- function(gnums) {
  # Gage numbers -> tibble(id, da, huc, elev) (i.e. drainage area)
  readNWISsite(gnums) %>%
    select(id=site_no, da=drain_area_va,
           huc=huc_cd, elev=alt_va) %>%
    mutate(da = da * 1.6^2,
           elev = elev/3.28)
}

add.gi <- function(data) {
  gi <- get.gageinfo(unique(data$id))
  left_join(gi,
            data,
            by="id")
}

list.gages <- function(state) {
  # List all gages with temperature data.
  whatNWISsites(stateCd=state, parameterCd="00010",
                siteType="ST")  # list sites
}

validate.gages <- function(glist,
                           start="1950-01-01",
                           end="2022-12-31",
                           log=F) {
  if (log)
    cat("\n\n")
  map_dbl(glist,
          ~{
            if (log & (runif(1) > 0.999)) {
              cat(".")
            }
            tryCatch(
              withTimeout({
                nrow(drop_na(readNWISdv(., "00010", start, end)))
              },
              onTimeout = "error",
              timeout = 600*10^3),
              error = function(e) {
                0
              }
            )
          })
}



states <- map_data("state")
eco.raw <- if(!exists("eco.raw")) {
  st_read(".Ecoregions/NA_CEC_Eco_Level1.shp", quiet=TRUE)
} else eco.raw

plot.eco <- function(akhi=F, vir.random=F, legrow=3) {
  ecos <- c("Northern Forests", "Eastern Temperate Forests", "Tropical Wet Forests", "Northwestern Forested Mountains", "Great Plains", "North American Deserts", "Mediterranean California", "Marine West Coast Forest", "Temperate Sierras", "Southern Semiarid Highlands")
  ecoreg <- eco.raw %>% st_transform("WGS84") %>%
    rename(ecoregion = NA_L1NAME) %>%
    mutate(ecoregion = str_to_title(ecoregion))
  ecoreg <- if (akhi) ecoreg else filter(ecoreg, ecoregion %in% ecos)
  plt <- ggplot(ecoreg) +
    geom_sf(aes(fill=if (vir.random) fct_shuffle(as.factor(ecoregion)) else ecoregion), color=NA) +
    scale_y_continuous(limits=c(
      if (akhi) 18 else 25,
      if (akhi) 72 else 50
    )) +
    scale_x_continuous(limits=c(
      if (akhi) -170 else -125,
      -60
    )) +
    (if (vir.random) {
      scale_fill_viridis_d(begin = 0.3)
    } else {
      scale_fill_discrete(
        # order: Arctic Cord,ETF, GP,Hudson Plain, MWCF, MC, NAD, NF, NWFM, SSH,
        # Taiga, TS, Tropical Dry Forests, TWF, Tundra
        type = c(
          "#CCCCFF", "#AA22AA", "#999955", "#6060F0",
          "#BBBBDD", "orange", "#FF6666", "#D5D5BC", "skyblue",
          "lightpink", "#60AA60", "darkgrey", "darkgreen", "lightblue", "#CCEECC",
          "#0000FF"
        ),
        guide = guide_legend(nrow=legrow, title.position = "top")
      )
    })
  
  plt
}
