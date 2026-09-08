#!/usr/bin/env Rscript
# Figure 28a. A small population LD model makes every contribution explicit.
# Within each block, a fair binary ancestral allele is flipped independently
# at SNP m with probability (1 - u_m)/2. Thus cor(G_i, G_m) = u_i*u_m
# off the diagonal. Summing two independent haplotypes preserves correlation.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1]] else "."
dir.create(out, recursive = TRUE, showWarnings = FALSE)
library(grid)

R <- diag(12)
R[1:7, 1:7] <- tcrossprod(c(.2, .6, .8, 1, .9, .7, .1))
R[8:11, 8:11] <- tcrossprod(c(.8, 1, .9, .6))
diag(R) <- 1
R2 <- R^2
score <- rowSums(R2)
focal <- 4L
cutoff <- .20
friends <- which(R2[focal, ] >= cutoff & seq_len(12) != focal)
stopifnot(isSymmetric(R), min(eigen(R, symmetric=TRUE)$values) > 0,
          all(diag(R2) == 1), all(R2 >= 0 & R2 <= 1),
          abs(score[4] - 3.35) < 1e-12, score[12] == 1,
          identical(friends, c(2L, 3L, 5L, 6L)))
write.csv(R, file.path(out, "population_genotype_correlation.csv"), row.names=FALSE)
write.csv(R2, file.path(out, "population_squared_correlation.csv"), row.names=FALSE)
write.csv(data.frame(SNP=1:12, r_squared_to_SNP4=R2[4, ], LD_score=score,
                     friend_of_SNP4_at_example_cutoff=(1:12 %in% friends)),
          file.path(out, "ld_scores_and_contributions.csv"), row.names=FALSE)
jsonlite::write_json(list(
  figure="28a", version="1", data="Exact illustrative population LD; not observed data",
  model="Independent binary haplotype blocks; allele flip probabilities (1-u)/2; two independent haplotypes per diploid genotype",
  block_1_loadings=c(.2,.6,.8,1,.9,.7,.1), block_2_loadings=c(.8,1,.9,.6),
  independent_SNP=12L, focal_SNP=focal, focal_score=unname(score[focal]),
  self_contribution=1, other_contribution=unname(score[focal]-1),
  illustrative_r_squared_cutoff=cutoff, other_friends=as.integer(friends),
  convention="Scores sum all 12 variants, including the focal variant; no threshold is used in the LD score",
  caveat="The r-squared cutoff illustrates a friend count. GCTA --ld uses a significance threshold; this is not a reproduction of that implementation. Finite-sample bias-corrected LD-score estimates need not be >=1.",
  references=c("https://doi.org/10.1038/ng.3211", "https://yanglab.westlake.edu.cn/software/gcta/"),
  minimum_eigenvalue=min(eigen(R, symmetric=TRUE)$values),
  R_version=R.version.string
), file.path(out, "model.json"), auto_unbox=TRUE, pretty=TRUE, digits=15)
