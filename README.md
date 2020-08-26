# Overview

This repository covers the analysis supporting $R$ estimation used with the these [LMIC projections for COVID19](https://github.com/cmmid/covidm_reports).

# Methods

This analysis uses publicly available data to

 - identify a reasonable point for local interventions to have achieved their effect
 - using that point and the publicly available reported cases, estimate a static $R_t$ for the pre- and post-intervention periods

The $R_t$ estimate is done with the [EpiNow2](https://github.com/epiforecasts/EpiNow2) library.

# HPC Engineering

Since we are doing this calculation for many locations, identically but independently, the computation is an "embarassingly parallel" problem.

As such, we use "array jobs". Each array point is associated with administrative region target for the $R_t$ calculation (e.g., we use the ISO3 country codes). Additionally, `EpiNow2` works with multiple chains; rather than further distribute the chains to array points, we simply use the multiple cores typically available to each node in an HPC setting.