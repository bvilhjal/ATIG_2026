# Run from the PCA folder. RStudio supplies Pandoc for HTML rendering.
if (!requireNamespace("rmarkdown", quietly = TRUE) ||
    !requireNamespace("knitr", quietly = TRUE)) {
  stop("Install rmarkdown and knitr first; see README.md.")
}
if (!rmarkdown::pandoc_available()) {
  candidates <- c(
    "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64",
    "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64"
  )
  candidates <- candidates[file.exists(file.path(candidates, "pandoc"))]
  if (length(candidates)) Sys.setenv(RSTUDIO_PANDOC = candidates[1])
}
if (!rmarkdown::pandoc_available()) stop("Pandoc is required; run from RStudio or install Pandoc.")
dir.create("results", showWarnings = FALSE)
started <- Sys.time()
run_env <- new.env(parent = globalenv())
rmarkdown::render("PCA_tutorial_answers.Rmd", output_format = "html_document",
                  envir = run_env, quiet = FALSE)
with(run_env, {
  checks <- c(
    sample_accounting = length(ind.norel) + length(ind.rel) == nrow(obj.bed),
    reference_accounting = length(ind.row) + sum(S >= 0.5) == length(ind.norel),
    projection_accounting = length(ind.row) + length(ind.proj) == nrow(obj.bed),
    sample_ids_aligned = identical(anc_info$sample.ID, obj.bed$fam$sample.ID),
    initial_dimensions = identical(dim(obj.svd$u), c(length(ind.norel), 20L)),
    final_dimensions = identical(dim(PCs), c(nrow(obj.bed), 20L)),
    finite_scores = all(is.finite(PCs)),
    orthogonal_initial = max(abs(crossprod(obj.svd$u) - diag(20))) < 1e-6,
    orthogonal_final = max(abs(crossprod(obj.svd2$u) - diag(20))) < 1e-6,
    centered_initial = max(abs(colMeans(predict(obj.svd)))) < 1e-6,
    positive_variance = all(pve > 0) && sum(pve) <= 1 + 1e-8,
    ordered_variance = all(diff(pve) <= 1e-8),
    reference_scores_retained = isTRUE(all.equal(PCs[ind.row, ], predict(obj.svd2)))
  )
  stopifnot(all(checks))
  write.table(data.frame(check = names(checks), passed = unname(checks)),
              "results/validation.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cbind(anc_info, PCs_df[paste0("PC", 1:20)]),
              "results/PCs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(PC = 1:20, pve = pve, cumulative_pve = cumulative_pve),
              "results/variance_explained.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(rel, "results/related_pairs.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(sample.ID = obj.bed$fam$sample.ID[ind.norel], S = S,
                         reference = S < 0.5),
              "results/outlier_scores.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  # Save fitted results, not the memory-mapped genotype object or temporary files.
  saveRDS(list(initial = obj.svd, final = obj.svd2, PCs = PCs,
               ind.rel = ind.rel, ind.norel = ind.norel, ind.row = ind.row,
               pve = pve, pc_cor = pc_cor, S = S, checks = checks),
          "results/pca_results.rds")
  summary <- c(samples = nrow(obj.bed), input_variants = ncol(obj.bed),
               related_pairs = nrow(rel), flagged_individuals = length(ind.rel),
               higher_threshold_pairs = nrow(rel2), initial_samples = length(ind.norel),
               initial_variants = length(attr(obj.svd, "subset")),
               outliers = sum(S >= 0.5), final_reference_samples = length(ind.row),
               final_variants = length(attr(obj.svd2, "subset")),
               projected_samples = length(ind.proj), PCs = ncol(PCs),
               PC1_PVE = pve[1], first20_PVE = sum(pve), cores = ncores)
  write.table(data.frame(metric = names(summary), value = unname(summary)),
              "results/summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  capture.output(sessionInfo(), file = "results/sessionInfo.txt")
  capture.output(system2(plink2, "--version", stdout = TRUE), file = "results/plink_version.txt")
})
writeLines(c(paste("Started:", format(started, tz = "UTC", usetz = TRUE)),
             paste("Finished:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
             paste("Elapsed seconds:", round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2)),
             "RNG seeds: 2026 (initial PCA), 2027 (refitted PCA)",
             paste("Pandoc:", rmarkdown::pandoc_version()),
             paste("Thread limits:", paste(Sys.getenv(c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")), collapse = ","))),
           "results/run_info.txt")
print(read.delim("results/summary.tsv"))
cat("All", length(run_env$checks), "numerical and output checks passed.\n")
