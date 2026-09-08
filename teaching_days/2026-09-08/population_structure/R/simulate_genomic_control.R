# Figure 26. Genomic control under uniform null inflation.
# Reuse the null draws in the lecture's subsequent LDSC example.
library(jsonlite)
args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args)) args[1] else 'genomic_control_simulation'
dir.create(B, recursive = TRUE, showWarnings = FALSE)
set.seed(904)
n <- 100000L
ld_score <- runif(n, 5, 95) # Consumed to match the subsequent lecture example.
z <- rnorm(n)
T_before <- 1.4 * z^2
null_median <- qchisq(0.5, df = 1)
lambda <- median(T_before) / null_median
T_after <- T_before / lambda
lambda_after <- median(T_after) / null_median

# Q-Q coordinates use every statistic. Thin only the displayed interior;
# retain all 200 largest values and identical ranks in both panels.
rank <- seq_len(n)
expected <- qchisq((rank - 0.5) / n, df = 1)
ord <- order(T_before)
ix <- unique(as.integer(round(c(seq(1, n - 1000, length.out = 180),
       seq(n - 999, n - 200, length.out = 160), (n - 199):n))))
display <- data.frame(rank = ix, expected = expected[ix],
                     before = T_before[ord][ix], after = T_after[ord][ix])
stopifnot(abs(null_median - qnorm(0.75)^2) < 1e-12,
          abs(lambda_after - 1) < 1e-12,
          identical(order(T_before), order(T_after)),
          max(abs(T_before / T_after - lambda)) < 1e-12)
p_before <- pchisq(T_before, 1, lower.tail = FALSE)
p_after <- pchisq(T_after, 1, lower.tail = FALSE)
stopifnot(all(p_after >= p_before))

write.csv(data.frame(snp = rank, z = z, T_before = T_before,
                    T_after = T_after), gzfile(file.path(B, 'simulated_statistics.csv.gz')),
          row.names = FALSE)
write.csv(display, file.path(B, 'qq_coordinates.csv'), row.names = FALSE)
result <- list(n = n, seed = 904, R = R.version.string,
  model = 'z_i iid N(0,1); T_i = 1.4 z_i^2; all SNPs are null',
  lambda_true = 1.4, lambda_before = lambda, lambda_after = lambda_after,
  null_median = null_median, median_before = median(T_before),
  median_after = median(T_after),
  false_positive_rate_0_05_before = mean(p_before < 0.05),
  false_positive_rate_0_05_after = mean(p_after < 0.05),
  max_expected = max(expected), max_before = max(T_before),
  max_after = max(T_after), displayed_points = nrow(display),
  coordinates = as.list(display),
  notes = paste('Original teaching simulation, not participant data or a genotype simulation.',
     'Identical null draws to the following lecture LDSC example.',
     'The corrected median equals the null median by construction.',
     'Uniform null inflation is an illustrative assumption; polygenicity can also inflate statistics.',
     'Display thinning does not affect estimation; all 200 largest statistics are shown.'))
write_json(result, file.path(B, 'simulation.json'), auto_unbox = TRUE,
           pretty = TRUE, digits = 16)
cat(sprintf('n=%d; lambda=%.8f -> %.8f; median=%.8f -> %.8f; FPR(.05)=%.5f -> %.5f; maxima=%.3f,%.3f,%.3f; plotted=%d\n',
  n, lambda, lambda_after, median(T_before), median(T_after),
  mean(p_before < .05), mean(p_after < .05), max(expected),
  max(T_before), max(T_after), nrow(display)))
