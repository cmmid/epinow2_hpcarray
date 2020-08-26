# Overview

This repository covers the analysis supporting $R$ estimation used with the these [LMIC projections for COVID19](https://github.com/cmmid/covidm_reports).

# Methods

This analysis uses publicly available data to

 - identify a reasonable point for local interventions to have achieved their effect
 - using that point and the publicly available reported cases, estimate a static $R_t$ for the pre- and post-intervention periods

The $R_t$ estimate is done with the [EpiNow2](https://github.com/epiforecasts/EpiNow2) library.