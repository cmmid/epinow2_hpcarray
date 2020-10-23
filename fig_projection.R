
suppressPackageStartupMessages({
  require(data.table)
  require(qs)
  require(ggplot2)
})

.debug <- "ZAF"
.args <- if (interactive()) sprintf(c(
  "~/Dropbox/Covid_LMIC/All_Africa_paper/r0_fitting/%s/alt_projection.qs",
#  "~/Downloads/task_20201022_082900_AD_EP_epinow2-hpcarray-preliminary_output_UKO/%s/projection.qs",
  "intros/intros.rds",
  "cases.rds", "altZ_%s.png"
), .debug) else commandArgs(trailingOnly = TRUE)

dt <- qread(.args[1])[
  compartment %in% c("cases", "subclinical","icu_p","nonicu_p","death_o")
]

dt[compartment == "death_o", compartment := "deaths"]

tariso <- gsub(".+([[:upper:]]{3})\\.png$","\\1", tail(.args, 1))

day0 <- readRDS(.args[2])[iso3 == tariso, intro.date[1]]

isoref <- gsub("^.*([[:upper:]]{3}).*$","\\1",.args[1])
cases <- melt(
  readRDS(.args[3])[iso3 == isoref],
  id.vars = c("date"),
  measure.vars = c("cases", "deaths"),
  variable.name = "compartment"
)
cases[order(date), cvalue := cumsum(value), by=compartment ]

dt[, date := day0 + t ]

allage.dt <- dt[,
  .(value = sum(value)),
  keyby=.(scen_id, run, compartment, date = day0 + t)
]
allage.dt[order(date), cvalue := cumsum(value), by=.(scen_id, run, compartment) ]

#' drop a line for intervention turn on
#' area for ignore window
#' area for end of fitting period

p <- ggplot(allage.dt[
  compartment %in% c("cases","subclinical","deaths")
][run == 3]) +
  facet_grid(compartment ~ ., scales = "free_y") +
  aes(date, value, color = factor(scen_id), linetype = c('ll','lo','md','hi','hh')[run], group = interaction(scen_id, run)) +
  geom_line() +
  geom_line(aes(color = "observed", linetype = "md", group = NULL), cases) +
  theme_minimal() +
  scale_x_date(date_breaks = "months", date_labels = "%b") +
  scale_y_log10("Incidence", labels = scales::label_number_si()) +
  scale_color_discrete(guide = "legend") +
  scale_linetype_manual(
    "Outcome\nQuantile (%)",
    labels = c(ll="2.5",lo="25",md="50",hi="75",hh="97.5"),
    values = c(md="solid", lo="dashed", hi="dashed", ll="dotted", hh="dotted")
  )
  
ggsave(tail(.args, 1), p, width = 9, height = 6, units = "in")