
# this simplifies local definitions
# of the '?=' assigned VARIABLES
# n.b.: 1) local.makefile is not required (because of the `-include`)
# 2) local.makefile should *only* define VARIABLES; no targets
-include local.makefile

default: test

# where are we putting analysis inputs?
INDIR ?= ~/Dropbox/input_rt_reports
# where are we putting analysis outputs?
OUTDIR ?= ~/Dropbox/output_rt_reports

# convenience definitions; note that `=` is lazy assignment
# so $^, $*, $@, etc not bound until executing a target
R = Rscript $^ $@
Rpipe = Rscript $^ $| $@
Rstar = Rscript $^ $* $@

# a prerequisite that will always be unsatisfied
FORCE:

${INDIR}/cases.rds: get_ecdc_cases.R | FORCE
	${R}

${INDIR}/iso3.csv: cases_to_csv.R ${INDIR}/cases.rds
	${R}

${INDIR}/interventions.rds: interventions.R
	${R}

${OUTDIR}/%/result.rds: compute.R ${INDIR}/cases.rds
	mkdir -p $(@D)
	${Rstar}

SLURMTEMP ?= slurm.template

${OUTDIR}/jobs.slurm: slurm.R ${SLURMTEMP} ${INDIR}/iso3.csv | ${OUTDIR}
	${Rpipe}

setup: ${INDIR}/cases.rds ${INDIR}/iso3.csv ${INDIR}/interventions.rds ${OUTDIR}/jobs.slurm

test: ${OUTDIR}/KEN/result.rds ${OUTDIR}/GHA/result.rds