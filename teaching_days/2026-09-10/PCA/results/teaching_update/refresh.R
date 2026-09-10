# Run from the PCA folder. Re-render teaching additions using the verified fit.
# Refuse changed fitting code or input files; retain the original full-run outputs.
library(bigsnpr)
get_chunk <- function(path, label) {
  x <- readLines(path, warn = FALSE)
  start <- grep(paste0("^```\\{r ", label, "[,}]"), x)
  stopifnot(length(start) == 1L)
  last <- start + which(x[(start + 1L):length(x)] == "```")[1]
  code <- trimws(x[(start + 1L):(last - 1L)])
  code <- trimws(sub("[[:space:]]+#.*$", "", code))
  code[nzchar(code)]
}
fit_chunks <- c("setup", "packages", "data", "relatedness", "initial-pca",
                "outlier-scores", "refit-pca", "projection")
for (label in fit_chunks) {
  original <- get_chunk("results/full_run_sources/PCA_tutorial_answers.Rmd", label)
  current <- get_chunk("PCA_tutorial_answers.Rmd", label)
  if (label == "initial-pca") current <- c(
    get_chunk("PCA_tutorial_answers.Rmd", "reference-samples"), current)
  if (label == "outlier-scores") original <- head(original, length(current))
  stopifnot(identical(original, current))
}
manifest <- jsonlite::read_json("results/source_manifest.json")
for (f in manifest$files) if (startsWith(f$path, "data/")) {
  stopifnot(identical(digest::digest(file = f$path, algo = "sha256"), f$sha256))
}
stopifnot(identical(digest::digest(file = "results/pca_results.rds", algo = "sha256"),
                     manifest$teaching_update$refresh_info$fitted_results_sha256))
started <- Sys.time()
z <- readRDS("results/pca_results.rds")
e <- new.env(parent = globalenv())
e$ATIG_SAVED_RENDER <- TRUE
e$obj.bed <- bed("data/1000G_phase3_common_norel.bed")
e$rel <- read.delim("results/related_pairs.tsv")
e$plink2 <- file.path("data", if (.Platform$OS.type == "windows") "plink2.exe" else "plink2")
e$obj.svd <- z$initial
e$obj.svd2 <- z$final
e$PCs <- z$PCs
e$S <- z$S
for (name in c("ind.rel", "ind.norel", "ind.row")) e[[name]] <- z[[name]]
e$ind.col <- attr(z$initial, "subset")
e$ind.proj <- setdiff(seq_len(nrow(e$obj.bed)), z$ind.row)
# Also fail immediately if the rendering hooks ever try to refit a model.
for (name in c("bed_autoSVD", "snp_plinkKINGQC", "bed_projectSelfPCA"))
  e[[name]] <- function(...) stop("Saved-result refresh must not refit PCA or relatedness.")
skipped <- c("relatedness", "initial-pca", "outlier-scores", "refit-pca", "projection")
knitr::opts_hooks$set(label = function(options) {
  if (options$label %in% skipped) options$eval <- FALSE
  options
})
if (!rmarkdown::pandoc_available())
  Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64")
rmarkdown::render("PCA_tutorial_answers.Rmd", output_format = "html_document",
                  envir = e, quiet = TRUE)
knitr::opts_hooks$restore()
# Reuse the original 13 validation expressions, without rewriting fitted outputs.
runner <- as.list(parse("run_exercise.R"))
with_block <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("with")), runner)
stopifnot(length(with_block) == 1L)
validation <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                       identical(x[[2]], as.name("checks")), as.list(with_block[[1]][[3]]))
stopifnot(length(validation) == 1L)
eval(validation[[1]], envir = e)
stopifnot(all(e$checks))
additional <- c(
  unchanged_PC_scores = identical(e$PCs, z$PCs),
  unchanged_variance = isTRUE(all.equal(e$pve, z$pve, tolerance = 1e-12)),
  unchanged_first8_comparison = isTRUE(all.equal(e$pc_cor, z$pc_cor)),
  unchanged_export_scores = isTRUE(all.equal(as.matrix(e$PCs_df[, 1:20]), z$PCs,
                                             check.attributes = FALSE)),
  cutoff_counts_monotone = all(diff(e$cutoff_counts$excluded) <= 0),
  population_means_complete = nrow(e$pop_means) == length(unique(e$anc_fit$Population)),
  plot_sign_alignment = all(diag(e$score_cor) * e$plot_sign >= 0),
  original_PC19_best_matches_final_PC18 = which.max(abs(e$score_cor[19, ])) == 18L
)
stopifnot(all(additional))
# Run each student starter in the verified analysis context, independently of
# the displayed worksheet's eval=FALSE setting. Keep exercise variables separate.
student_env <- new.env(parent = e)
starter_labels <- c("kinship-starter", "pc-plot-starter", "comparison-starter", "sign-starter")
starter_results <- list()
grDevices::pdf(tempfile(fileext = ".pdf"))
for (label in starter_labels) starter_results[[label]] <- eval(
  parse(text = get_chunk("PCA_tutorial.Rmd", label)), envir = student_env)
grDevices::dev.off()
classroom_checks <- c(
  score_definition = isTRUE(all.equal(e$initial_scores,
                                      sweep(e$obj.svd$u, 2, e$obj.svd$d, "*"))),
  projection_roles = identical(which(e$fit_role != "Reference"), e$ind.proj) &&
    all(e$fit_role[e$ind.rel] == "Projected: related") &&
    all(e$fit_role[e$ind.norel[e$S >= 0.5]] == "Projected: outlier"),
  sign_flip_preserves_distances = isTRUE(e$sign_distances_match),
  student_labels_aligned = identical(student_env$anc_info$sample.ID, e$obj.bed$fam$sample.ID),
  student_plot_ready = inherits(starter_results[["pc-plot-starter"]], "ggplot"),
  student_comparison_matches = isTRUE(all.equal(student_env$C, e$score_cor)),
  student_sign_starter = isTRUE(starter_results[["sign-starter"]])
)
stopifnot(all(classroom_checks))
checks <- c(e$checks, additional, classroom_checks)
write.table(data.frame(check = names(checks), passed = unname(checks)),
            "results/teaching_update/validation.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(e$cutoff_counts, "results/teaching_update/cutoff_sensitivity.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(e$pop_means, "results/teaching_update/population_means.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(e$score_cor, "results/teaching_update/PC_correlations.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
knitr::opts_chunk$set(eval = FALSE)
rmarkdown::render("PCA_tutorial.Rmd", output_format = "html_document", quiet = TRUE)
knitr::opts_chunk$restore()
jsonlite::write_json(list(
  date = as.character(Sys.Date()), elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
  mode = "Teaching update rendered from saved, previously validated full-data fit; no PCA refit.",
  unchanged_fit_chunks = fit_chunks,
  fitting_code_note = "The original initial-pca code is now split into reference-samples and initial-pca; their concatenated R statements match the full-run source.",
  student_starters_executed = starter_labels, checks_passed = length(checks),
  fitted_results_sha256 = digest::digest(file = "results/pca_results.rds", algo = "sha256")
), "results/teaching_update/refresh_info.json", auto_unbox = TRUE, pretty = TRUE)
cat("Teaching refresh passed", length(checks), "checks.\n")
