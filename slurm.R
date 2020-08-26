
.args <- if (interactive()) c(
  "slurm.template", "iso3.csv", ".", "jobs.slurm"
) else commandArgs(trailingOnly = TRUE)

# this isn't very flexible
# essentially, hardcoding format as *known*
template <- .args[1]
fmt <- readChar(template, file.info(template)$size)

isosfile <- .args[2]
isocount <- dim(read.csv(isosfile, header = FALSE))[1]

outputdir <- .args[3]
slurmtarget <- tail(.args, 1)

out <- sprintf(fmt, isocount, isosfile, outputdir)
write(out, slurmtarget)