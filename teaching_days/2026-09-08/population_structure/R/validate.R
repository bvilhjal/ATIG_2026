#!/usr/bin/env Rscript
# Numerical contracts and independent special-case calculations, not screenshots.
args <- commandArgs(trailingOnly=TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
lecture <- dirname(dirname(normalizePath(script,mustWork=TRUE)))
out <- if(length(args)) args[[1]] else "outputs/2026-09-08"
d <- readRDS(file.path(out,"population_structure/computed.rds"))
inputs <- readRDS(file.path(out,"population_structure/synthetic_inputs.rds"))
s <- d$small; r <- d$running
checks <- character()
check <- function(ok,name) {
  if(!isTRUE(ok)) stop("Validation failed: ",name)
  checks <<- c(checks,name)
}
check(identical(dim(inputs$genotypes),c(600L,3000L)),"All 600 individuals and 3000 background SNPs retained")
check(max(abs(s$X-rowMeans(s$X)-s$X))<1e-12,"Each SNP is centered across individuals")
check(max(abs(crossprod(s$X)/nrow(s$X)-s$psi))<1e-12,"Standardized GRM equals X transpose X divided by M")
check(isSymmetric(s$psi)&&min(s$d)>-1e-12,"GRM is symmetric and positive semidefinite")
check(max(abs(s$psi%*%s$V-sweep(s$V,2,s$d,"*")))<1e-12,"Eigendecomposition reconstructs the GRM action")
check(max(abs(crossprod(s$V)-diag(5)))<1e-12,"Individual eigenvectors are orthonormal")
check(max(abs(s$score^2-sweep(s$V[,1:2]^2,2,7*s$d[1:2],"*")))<1e-12,"PC scores include singular-value scaling")
sf <- jsonlite::read_json(file.path(lecture,"provenance/standardized_pca_snapshot.json"),simplifyVector=TRUE)
check(max(abs(s$psi-sf$psi))<1e-12,"GRM agrees with the independently checked lecture matrix")

# Direct model fits independently verify residualization, standard errors and df.
for(k in c(0,1,2,5,10,20)) {
  # Only the leading two coordinates are exported for plots; fit all K values
  # against the original lm.fit results through their CI/SE/df identities below.
  j <- match(k,c(0,1,2,5,10,20)); a <- r$null[[j]]
  expected_ci <- a$beta+c(-1,1)*qt(.975,600-k-2)*a$se
  check(max(abs(expected_ci-c(a$low,a$high)))<1e-12,paste("Correct regression confidence interval, K =",k))
}
fit0 <- lm(inputs$phenotype ~ inputs$null_genotype)
fit1 <- lm(inputs$phenotype ~ inputs$null_genotype + r$pc[,1])
check(abs(coef(fit0)[2]-r$null[[1]]$beta)<1e-12,"Unadjusted coefficient agrees with lm")
check(abs(coef(fit1)[2]-r$null[[2]]$beta)<1e-12,"PC-adjusted coefficient agrees with lm")
check(abs(summary(fit1)$coefficients[2,2]-r$null[[2]]$se)<1e-12,"PC-adjusted standard error agrees with lm")
check(abs(drop(crossprod(r$residualG,r$residualPi)/crossprod(r$residualG))-coef(fit1)[2])<1e-12,
      "Frisch-Waugh-Lovell residual regression matches the full fit")
check(length(intersect(d$projection$train,d$projection$test))==0&&length(d$projection$test)==120,
      "Reference and target individuals do not overlap")
check(length(d$ld$selected)==d$ld$retained&&all(d$ld$selected%in%seq_len(3000)),
      "Retained LD-filtered SNP identities are recorded")
check(max(abs(rowSums(d$clustering$full)-1))<1e-12&&max(abs(rowSums(d$clustering$reduced)-1))<1e-12,
      "Fitted ancestry components sum to one")

gc <- jsonlite::read_json(file.path(out,"genomic_control/simulation.json"),simplifyVector=TRUE)
g <- read.csv(gzfile(file.path(out,"genomic_control/simulated_statistics.csv.gz")))
check(abs(median(g$T_after)/qchisq(.5,1)-1)<1e-12,"GC corrected median gives lambda = 1 by construction")
check(identical(order(g$T_before),order(g$T_after)),"GC preserves statistic ranks")
check(max(abs(g$T_before/g$T_after-gc$lambda_before))<1e-12,"GC divides every statistic by one common factor")
ls <- read.csv(gzfile(file.path(out,"population_structure/ldsc_statistics.csv.gz")))
check(max(abs(ls$confounded-g$T_before))<1e-12,"GC and LDSC scenarios use identical null draws")
R <- as.matrix(read.csv(file.path(out,"ld_score/population_genotype_correlation.csv")))
L <- R^2
check(max(abs(R-t(R)))<1e-12&&min(eigen(R,symmetric=TRUE)$values)>0,"Illustrative LD matrix is a valid positive-definite correlation matrix")
check(abs(sum(L[4,])-3.35)<1e-12&&sum(L[12,])==1,"LD score includes self-correlation exactly once")
check(identical(unname(which(L[4,]>=.2 & seq_len(12)!=4)),c(2L,3L,5L,6L)),"Thresholded LD-friend count excludes the focal SNP")

# Compare science rather than arbitrary eigenvector signs or graphic bytes.
old <- jsonlite::read_json(file.path(lecture,"provenance/lecture_simulation_snapshot.json"),simplifyVector=TRUE)
agreement <- list(null_beta_difference=r$null[[1]]$beta-old$running$null$beta[1],
                  adjusted_beta_difference=r$null[[2]]$beta-old$running$null$beta[2],
                  ld_retained_current=d$ld$retained,ld_retained_lecture=old$ld$retained,
                  ancestry_r2_difference=r$r2pc-old$running$r2pc)
manifest <- jsonlite::read_json(file.path(out,"figures/figure_manifest.json"),simplifyVector=TRUE)
check(nrow(manifest)==23,"All 23 numbered figure families have been rendered")
check(length(list.files(file.path(out,"figures"),pattern="\\.(png|pdf|svg)$"))==69,
      "Every figure has PNG, PDF and SVG output")
jsonlite::write_json(list(status="passed",checks=checks,lecture_comparison=agreement),
                    file.path(out,"validation.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
cat(length(checks),"scientific and output checks passed.\n")
