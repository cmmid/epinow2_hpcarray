suppressPackageStartupMessages({
  require(data.table)
  require(wpp2019)
  require(countrycode)
})

.args <- if (interactive()) c(
  "../cases.rds", "ene-ifr.csv", "intros.rds"
) else commandArgs(trailingOnly = TRUE)

#' get all the death data
deaths.dt <- readRDS(.args[1])[, .(deaths), keyby=.(iso3, date) ]

#' identify which countries we can't do infection seeding for
exclude <- deaths.dt[,.(tot = sum(deaths)),by=iso3][tot == 0, iso3]
iso3c <- deaths.dt[,.(tot = sum(deaths)),by=iso3][tot != 0, iso3]
#' need the age pyramids for these countries
include <- countrycode(iso3c, "iso3c", "iso3n")
exclude <- unique(c(exclude, iso3c[is.na(include)]))
include <- include[!is.na(include)]

warning(
  sprintf(
    "excluding the following ISOs with no reported deaths:\n%s",
    paste(exclude, collapse = "\n")
  )
)

#' going to seed from all the deaths from 1st to approximately
#' one generation later (5 days)
death_n_lim <- 5

#' this gives us the initial deaths that we're going to inflate
#' into infections
red.dt <- deaths.dt[!(iso3 %in% exclude)][deaths != 0][,
  .SD[date <= (date[1]+death_n_lim)],
  by = iso3
]

#' IFR from ENE study, as extracted at
#' https://github.com/mbevand/covid19-age-stratified-ifr
ifr.dt <- fread(.args[2])

data(pop)

popreformat <- function(p, gen) as.data.table(p)[
  country_code %in% include,
  .(gender = gen, age = 1:10, pop = {
    res <- rowSums(matrix(`2020`[1:20], byrow = T, nrow = 10))
    res[10] <- res[10] + `2020`[21]
    res*1000
  }),
  keyby=.(iso3 = countrycode(country_code, "iso3n", "iso3c"))
]

pop.dt <- rbind(
  popreformat(popF, "F"),
  popreformat(popM, "M")
)[, .(pop = sum(pop)), keyby=.(iso3, age)]

ifr.exp <- pop.dt[
  ifr.dt[, agen := 1:.N], on=.(age = agen)
][,
  .(ifr.ave = sum(pop*ifr)/sum(pop)/100),
  keyby=iso3
]

excess.dt <- fread(
  "https://github.com/TheEconomist/covid-19-excess-deaths-tracker/raw/master/output-data/excess-deaths/all_weekly_excess_deaths.csv"
)[country %in% c("Brazil", "Mexico", "South Africa", "Turkey")]

ref <- excess.dt[,
  .SD[cumsum(covid_deaths) > 0],
  by=country,
  .SDcols = c("country","week", grep("deaths", names(excess.dt), value = TRUE))
][,
  .SD[
    excess_deaths > 0
  ][1:{
    ny <- excess_deaths_per_100k < covid_deaths_per_100k
    if (any(ny)) which.max(ny) else .N
  }],
  by= country
][,
  .(mul = sum(excess_deaths_per_100k - covid_deaths_per_100k)/sum(covid_deaths_per_100k) + 1),
  by = country
]

mu.inf <- ref[, mean(mul)]

red.dt[ifr.exp, on=.(iso3), intros := round(deaths*mu.inf/ifr.ave)]
red.dt[, intro.date := date - (17+5) ]
red.dt[, intro.day := as.integer(date - date[1]), by=iso3 ]

saveRDS(red.dt, tail(.args, 1))
