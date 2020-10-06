
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
OTHDIR ?= ~/Dropbox/covidm_reports/hpc_inputs
# how many cores available?
NCORES ?= 3

${INDIR} ${OUTDIR}:
	mkdir -p $@

# convenience definitions; note that `=` is lazy assignment
# so $^, $*, $@, etc not bound until executing a target
R = Rscript $(filter-out FORCE,$^) $@
Rpipe = Rscript $(filter-out FORCE,$^) $| $@
Rstar = Rscript $(filter-out FORCE,$^) $* $@

# a prerequisite that will always be unsatisfied
FORCE:

${INDIR}/cases.rds: get_ecdc_cases.R
	${R}

${INDIR}/iso3.csv: cases_to_csv.R ${INDIR}/cases.rds
	${R}

WHOURL := https://www.who.int/docs/default-source/documents/phsm/20200722-phsm-who-int.zip

# requires python tool xlsx2csv
${INDIR}/rawinterventions.csv:
	curl ${WHOURL} --output tmp.zip
	unzip tmp.zip *.xlsx
	mv *.xlsx tmp.xlsx
	xlsx2csv tmp.xlsx > $@
	rm tmp.zip tmp.xlsx

${INDIR}/interventions.rds: interventions.R ${INDIR}/rawinterventions.csv
	${R}

${INDIR}/rt_bounds.rds: rt_bounds.R $(addprefix ${INDIR}/,cases.rds interventions.rds)
	${R}

${OUTDIR}/%/result.rds: compute.R $(addprefix ${INDIR}/,cases.rds rt_bounds.rds)
	mkdir -p $(@D)
	Rscript $(filter-out FORCE,$^) $* ${NCORES} $@

${OUTDIR}/fits/%.rds: param_fit.R ${OUTDIR}/%/result.rds ${OTHDIR}/%/params_set.rds ${OTHDIR}/%/contact_matrices.rds ${OTHDIR}/covidm_fit_yu.qs
	${R}

SLURMTEMP ?= slurm.template

${OUTDIR}/jobs.slurm: slurm.R ${SLURMTEMP} ${INDIR}/iso3.csv | ${OUTDIR}
	${Rpipe}

${OUTDIR}/SI_fig_int_distro.png: SI_fig_int_distro.R $(addprefix ${INDIR}/,interventions.rds cases.rds)
	${R}

sifigs: ${OUTDIR}/SI_fig_int_distro.png

setup: ${INDIR} ${OUTDIR} $(addprefix ${INDIR}/,cases.rds iso3.csv jobs.slurm rawinterventions.csv interventions.rds rt_bounds.rds)

ALLISOS ?= $(shell cat ${INDIR}/iso3.csv)

test: $(patsubst %,${OUTDIR}/%/result.rds,${ALLISOS})