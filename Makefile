
# this simplifies local definitions
# of the '?=' assigned VARIABLES
# n.b.: 1) local.makefile is not required (because of the `-include`)
# 2) local.makefile should *only* define VARIABLES; no targets
-include local.makefile

default: test

# where are we putting analysis inputs?
INDIR ?= ~/Dropbox/input_rt_reports
# where are we putting analysis outputs?
OUTDIR ?= ~/Dropbox/Covid_LMIC/All_Africa_paper/r0_fitting
OTHDIR ?= ~/Dropbox/covidm_reports/hpc_inputs

COVIDMPATH ?= ../covidm

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

${OUTDIR}/%/result.rds: compute.R $(addprefix ${INDIR}/,cases.rds rt_bounds.rds) ${OTHDIR}/%/contact_matrices.rds ${OTHDIR}/%/params_set.rds ${OTHDIR}/covidm_fit_yu.qs
	mkdir -p $(@D)
	Rscript $(filter-out FORCE,$^) $* ${NCORES} $@

${OUTDIR}/eligible.csv: params_filter.R | ${OUTDIR}
	Rscript $< $| $@

${OUTDIR}/fits/%.rds: param_fit.R ${OUTDIR}/%/result.rds ${OTHDIR}/%/params_set.rds ${OTHDIR}/%/contact_matrices.rds ${OTHDIR}/covidm_fit_yu.qs
	${R}

intros/intros.rds intros/urban.rds:
	cd $(@D) && $(MAKE) -C $(@F)

NGM.rda: NGM.R
	${R}

${OUTDIR}/scens/%.rds: gen_scenarios.R ${OUTDIR}/fits/%.rds ${OUTDIR}/%/result.rds intros/intros.rds
	${R}

testscen: ${OUTDIR}/scens/ZAF.rds

${OUTDIR}/%/projection.qs: run_scenario.R ${OUTDIR}/scens/%.rds ${OTHDIR}/%/params_set.rds ${OTHDIR}/covidm_fit_yu.qs ${OUTDIR}/%/result.rds ${OTHDIR}/%/contact_matrices.rds intros/intros.rds intros/urban.rds | ${COVIDMPATH}
	Rscript $^ $| $@

testfit: ${OUTDIR}/fits/ZAF.rds
testproj: ${OUTDIR}/ZAF/projection.qs

.PRECIOUS: ${OUTDIR}/%/result.rds ${OUTDIR}/fits/%.rds ${OUTDIR}/scens/%.rds ${OUTDIR}/%/projection.qs

SLURMTEMP ?= slurm.template

${OUTDIR}/jobs.slurm: slurm.R ${SLURMTEMP} ${INDIR}/iso3.csv | ${OUTDIR}
	${Rpipe}

${OUTDIR}/SI_fig_int_distro.png: SI_fig_int_distro.R $(addprefix ${INDIR}/,interventions.rds cases.rds)
	${R}

sifigs: ${OUTDIR}/SI_fig_int_distro.png

setup: ${INDIR} ${OUTDIR} $(addprefix ${INDIR}/,cases.rds iso3.csv rawinterventions.csv interventions.rds rt_bounds.rds) NGM.rda intros/intros.rds intros/urban.rds

ALLISOS ?= $(shell cat ${INDIR}/iso3.csv)
FIGISOS := $(shell cd ~/Downloads/task_20201022_082900_AD_EP_epinow2-hpcarray-preliminary_output_UKO && ls)

${OUTDIR}/figs:
	mkdir -p $@

${OUTDIR}/figs/%.png: fig_projection.R ~/Downloads/task_20201022_082900_AD_EP_epinow2-hpcarray-preliminary_output_UKO/%/projection.qs intros/intros.rds cases.rds | ${OUTDIR}/figs
	Rscript $^ $| $@

checkiso: $(patsubst %,${OUTDIR}/figs/%.png,${FIGISOS})

test: $(patsubst %,${OUTDIR}/%/result.rds,${ALLISOS})