# Exact classroom examples; no random draws and no participant data.
args <- commandArgs(trailingOnly=TRUE)
out <- if (length(args)) args[[1]] else "outputs/2026-09-08/small_calculations"
dir.create(out, recursive=TRUE, showWarnings=FALSE)

G <- rbind(c(0,1,0,2,1), c(1,0,1,2,1), c(2,1,1,0,2), c(0,2,1,0,1))
raw_cor <- cor(G)
X <- G - rowMeans(G)
Psi <- crossprod(X)/nrow(X)
G_flip <- G
G_flip[1,] <- 2-G_flip[1,]
X_flip <- G_flip-rowMeans(G_flip)
stopifnot(max(abs(crossprod(X_flip)/nrow(X_flip)-Psi)) < 1e-12)

# A two-individual centered matrix with an analytic eigendecomposition.
toyX <- cbind(2: -2, -(2: -2))
toyPsi <- crossprod(toyX)/5
toyV <- matrix(c(1,-1,1,1),2,2)/sqrt(2)
toyD <- diag(c(4,0))
stopifnot(max(abs(toyPsi-toyV%*%toyD%*%t(toyV))) < 1e-12)

# Full binomial likelihoods include the coefficient for each heterozygote.
pa <- c(.25,.57,.29,.38)
pb <- c(.40,.32,.84,.22)
g <- c(2,0,1,1)
loglik <- function(alpha) sum(dbinom(g,2,alpha*pa+(1-alpha)*pb,log=TRUE))
opt <- optimize(loglik,c(0,1),maximum=TRUE,tol=1e-12)
alpha <- seq(0,1,length.out=1001)
L1 <- exp(loglik(1)); L2 <- exp(loglik(0))
stopifnot(abs(L1-.002242376599)<1e-12, abs(L2-.00682518380544)<1e-12,
          abs(opt$maximum-.2176550693)<1e-7)
result <- list(raw_genotypes=G,raw_correlations=raw_cor,centered_gram=Psi,
               toyX=toyX,toyPsi=toyPsi,toyV=toyV,toyD=toyD,
               pa=pa,pb=pb,g=g,alpha=alpha,loglik=vapply(alpha,loglik,0.0),
               alpha_hat=opt$maximum,L1=L1,L2=L2,ratio=L2/L1,
               provenance="Deterministic examples from the supplied lecture; calculations checked independently.")
saveRDS(result,file.path(out,"calculations.rds"))
jsonlite::write_json(result,file.path(out,"calculations.json"),auto_unbox=TRUE,
                     digits=16,pretty=TRUE)
cat(sprintf("Four-SNP ancestry likelihood: L1=%.12f, L2=%.12f, alpha=%.9f\n",L1,L2,opt$maximum))
