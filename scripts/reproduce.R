#!/usr/bin/env Rscript
# Run from the repository root. Each module runs in a fresh R process so that
# plotting and module order cannot consume another simulation's random stream.
args <- commandArgs(trailingOnly=TRUE)
out <- if (length(args)) args[[1]] else "outputs"
if (!file.exists("VERSION")) stop("Run this script from the ATIG_2026 repository root.")
required <- c("bigsnpr","bigstatsr","jsonlite")
missing <- required[!vapply(required,requireNamespace,FALSE,quietly=TRUE)]
if (length(missing)) stop("Install these R packages first: ",paste(missing,collapse=", "))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
out <- normalizePath(out,mustWork=TRUE)
Sys.setenv(OMP_NUM_THREADS=1,OPENBLAS_NUM_THREADS=1,MKL_NUM_THREADS=1,
           VECLIB_MAXIMUM_THREADS=1)
modules <- c(population_structure="simulate_population_structure.R",
             genomic_control="simulate_genomic_control.R",
             ld_score="simulate_ld_score.R",small_calculations="small_calculations.R")
times <- list()
run <- function(file,arguments) {
  status <- system2(file.path(R.home("bin"),"Rscript"),
                    shQuote(c(file,arguments)))
  if (status != 0) stop("Failed: ",file)
}
for (name in names(modules)) {
  cat("Running",name,"\n")
  elapsed <- system.time(run(file.path("population_structure/R",modules[[name]]),
                             file.path(out,name)))[["elapsed"]]
  times[[name]] <- elapsed
}
times$figures <- system.time(run("population_structure/R/plot_figures.R",out))[["elapsed"]]
run("scripts/validate.R",out)
jsonlite::write_json(list(version=trimws(readLines("VERSION")),
                         UTC=format(Sys.time(),tz="UTC",usetz=TRUE),
                         elapsed_seconds=times, R=R.version.string,
                         packages=setNames(lapply(required,function(p) as.character(packageVersion(p))),required),
                         platform=R.version$platform, system=Sys.info()[c("sysname","release","machine")],
                         BLAS=extSoftVersion()[["BLAS"]], threads=1,
                         data="Synthetic teaching examples and separately labelled published summary values"),
                    file.path(out,"run_manifest.json"),auto_unbox=TRUE,pretty=TRUE)
cat("Complete. Figures, data, checks and provenance are in",out,"\n")
