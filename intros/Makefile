
INDIR ?= ..

R = Rscript $^ $@

intros.rds: intros.R ${INDIR}/cases.rds ene-ifr.csv
	${R}

WBURL := http://api.worldbank.org/v2/en/indicator/SP.URB.TOTL.IN.ZS?downloadformat=csv

urban.rds:
	curl ${WBURL} --output tmp.zip
	unzip tmp.zip API_SP.URB.TOTL.IN.ZS_DS2_*.csv
	mv API_SP.URB.TOTL.IN.ZS_DS2_*.csv tmp.csv
	Rscript urbanization.R tmp.csv $@
	rm tmp.zip tmp.csv