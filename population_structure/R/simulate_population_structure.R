# Synthetic population-structure teaching examples. See docs/METHODS.md.
library(bigsnpr)
library(jsonlite)
options(bigstatsr.ncores.max = 2L, bigstatsr.ncores = 1L, bigstatsr.check.parallel.blas = FALSE)
set.seed(20260907)
args <- commandArgs(trailingOnly = TRUE)
B <- if (length(args)) args[[1]] else 'outputs/population_structure'
dir.create(B, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(B, 'fbm'), showWarnings = FALSE)
mk <- function(x, name) {
  unlink(file.path(B, 'fbm', paste0(name, c('.bk', '.rds'))))
  z <- FBM.code256(nrow(x), ncol(x), code = CODE_012,
                  backingfile = file.path(B, 'fbm', name))
  z[,] <- x
  z
}
align <- function(o, a) {
  for (k in seq_len(ncol(o$u))) {
    sg <- if (k == 1) sign(cor(o$u[, k], a)) else sign(o$u[which.max(abs(o$u[, k])), k])
    if (sg == 0) sg <- 1
    o$u[, k] <- o$u[, k] * sg
    o$v[, k] <- o$v[, k] * sg
  }
  o
}
N <- 600; M <- 3000
a <- rep(0:1, each = N / 2)
base <- runif(M, .2, .8); delta <- sample(c(-1, 1), M, TRUE) * .065
freq <- outer(a - .5, 2 * delta, '*') + matrix(base, N, M, byrow = TRUE)
G <- matrix(rbinom(N * M, 2, freq), N, M)
fbm <- mk(G, 'background')
svd <- align(big_randomSVD(fbm, fun.scaling = snp_scaleBinom(), k = 20, ncores = 1), a)
v <- svd$u
g0 <- rbinom(N, 2, ifelse(a == 0, .2, .6))
g1 <- rbinom(N, 2, .4)
pi <- 1.2 * (a - .5) + .25 * g1 + rnorm(N)
fit <- function(g, K) {
  C <- if (K == 0) matrix(1, N, 1) else cbind(1, v[, seq_len(K), drop = FALSE])
  z <- lm.fit(cbind(C, g), pi)
  k <- ncol(C) + 1
  se <- sqrt(sum(z$residuals^2) / z$df.residual * chol2inv(z$qr$qr[seq_len(k), seq_len(k), drop = FALSE])[k, k])
  b <- unname(z$coefficients[k]); ci <- b + c(-1, 1) * qt(.975, z$df.residual) * se
  list(K = K, beta = b, se = se, low = ci[1], high = ci[2], p = 2 * pt(-abs(b / se), z$df.residual))
}
fit0 <- lapply(c(0, 1, 2, 5, 10, 20), function(k) fit(g0, k))
fit1 <- lapply(c(0, 1, 2, 5, 10, 20), function(k) fit(g1, k))
C <- cbind(1, v[, 1]); resg <- lm.fit(C, g0)$residuals; respi <- lm.fit(C, pi)$residuals
jitter <- runif(N, -.13, .13)
running <- list(N=N, M=M, ancestry=a, genotype=g0, jitter=jitter, phenotype=pi,
  pc=v[,1:2], score=predict(svd)[,1:2], residualG=resg, residualPi=respi,
  null=fit0, causal=fit1, truth=c(0,.25), r2pc=cor(v[,1], a)^2)

# A large correlated region overwhelms the unfiltered PCA.
set.seed(901)
block <- rbinom(N, 2, .4)
GLD <- G
for (j in 1201:1700) GLD[,j] <- ifelse(runif(N) < .985, block, G[,j])
FLD <- mk(GLD, 'ld_example')
raw <- big_randomSVD(FLD, fun.scaling=snp_scaleBinom(), k=4, ncores=1)
clean <- snp_autoSVD(FLD, infos.chr=rep(1:3, each=1000),
   infos.pos=rep(seq_len(1000)*10000, 3), k=4, ncores=1,
   roll.size=20, int.min.size=10, verbose=FALSE)
clean <- align(clean, a)
ld <- list(raw=raw$u[,1:2], clean=clean$u[,1:2], ancestry=a, block=block,
   rawLoad=raw$v[,1], selected=attr(clean,'subset'),
   rawAncestryR2=cor(raw$u[,1],a)^2, cleanAncestryR2=cor(clean$u[,1],a)^2,
   rawBlockR2=cor(raw$u[,1],block)^2,
   retained=length(attr(clean,'subset')), regions=attr(clean,'lrldr'))

# A second simulation creates a technical shift independent of ancestry.
set.seed(902)
batch <- sample(rep(0:1, length.out=N))
batchfreq <- pmax(.02, pmin(.98, freq + outer(batch-.5, rep(.15,M), '*')))
GB <- matrix(rbinom(N*M, 2, batchfreq), N, M)
FB <- mk(GB, 'batch_example')
sb <- align(big_randomSVD(FB, fun.scaling=snp_scaleBinom(), k=4, ncores=1), a)
batchout <- list(pc=sb$u[,1:2], ancestry=a, batch=batch,
   r2ancestry=apply(sb$u[,1:2],2,function(x)cor(x,a)^2),
   r2batch=apply(sb$u[,1:2],2,function(x)cor(x,batch)^2))

# Use an independent fitting subset and the actual OADP projection implementation.
set.seed(903)
train <- sort(c(sample(which(a==0),240),sample(which(a==1),240)))
test <- setdiff(seq_len(N),train)
sp <- snp_autoSVD(fbm, infos.chr=rep(1:3,each=1000),
   infos.pos=rep(seq_len(1000)*10000,3), ind.row=train,
   k=4, ncores=1, verbose=FALSE)
proj <- snp_projectSelfPCA(sp,fbm,ind.row=test,ncores=1)
projection <- list(train=train,test=test,reference=predict(sp)[,1:2],
  target=proj$OADP_proj[,1:2], ancestry=a, methods=names(proj),
  selected=length(attr(sp,'subset')))

# Same five individuals for every step of the small PCA example.
smallG <- rbind(c(1,1,1,0,0),c(0,1,2,1,2),c(2,1,1,0,1),
 c(0,0,1,2,2),c(2,1,1,0,0),c(0,0,1,1,1),c(2,2,1,1,0))
# The final lecture uses binomial standardization, not centering alone.
smallP <- rowMeans(smallG) / 2
smallScale <- sqrt(2 * smallP * (1 - smallP))
smallX <- (smallG - 2 * smallP) / smallScale
psi <- crossprod(smallX)/nrow(smallG)
es <- eigen(psi,symmetric=TRUE); V <- es$vectors
for(k in 1:5) if(V[which.max(abs(V[,k])),k]<0) V[,k] <- -V[,k]
score <- sweep(V[,1:2],2,sqrt(7*es$values[1:2]),'*')
theta <- 0.52
rotation <- matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2,2)
small <- list(G=smallG,X=smallX,psi=psi,V=V,d=es$values,
  score=score,geometry=score %*% t(rotation),axes=rotation,
  means=rowMeans(smallG),frequencies=smallP,scales=smallScale,explained=es$values/sum(es$values))
stopifnot(max(abs(psi%*%V-sweep(V,2,es$values,'*')))<1e-10,
 max(abs(crossprod(V)-diag(5)))<1e-10)

# Matched inflation under the LDSC working model; not simulated genotypes.
set.seed(904)
ns <- 100000; ell <- runif(ns,5,95); z <- rnorm(ns)
tc <- 1.4*z^2; tp <- (1+.008*ell)*z^2
qq <- function(t) {
  o <- sort(-log10(pmax(pchisq(t,1,lower.tail=FALSE),1e-300)))
  e <- sort(-log10((seq_len(ns)-.5)/ns))
  ix <- unique(round(c(seq(1,ns-1000,length.out=120),seq(ns-999,ns,length.out=120))))
  list(x=e[ix],y=o[ix])
}
bins <- cut(ell,breaks=seq(5,95,length.out=21),include.lowest=TRUE)
means <- function(t) as.numeric(tapply(t,bins,mean))
ldsc <- list(confQQ=qq(tc), polyQQ=qq(tp), ell=means(ell), conf=means(tc), poly=means(tp),
  confFit=unname(coef(lm(tc~ell))),polyFit=unname(coef(lm(tp~ell))),
  confLambda=median(tc)/qchisq(.5,1),polyLambda=median(tp)/qchisq(.5,1),
  confMean=mean(tc),polyMean=mean(tp),Nstats=ns)

# Supervised ancestry likelihood: the same 30 genomes under two reference panels.
set.seed(905)
nm <- 1200; pa <- runif(nm,.15,.85); pb <- runif(nm,.15,.85)
pc <- pmax(.03,pmin(.97,(pa+pb)/2+rnorm(nm,0,.09)))
pf <- cbind(pa,pb,pc); truegroup <- rep(1:3,each=10)
gc <- sapply(truegroup,function(q)rbinom(nm,2,pf[,q]))
grid3 <- as.matrix(expand.grid(a=seq(0,1,.05),b=seq(0,1,.05)))
grid3 <- grid3[rowSums(grid3)<=1+1e-10,,drop=FALSE]
grid3 <- cbind(grid3,1-rowSums(grid3))
grid2 <- cbind(seq(0,1,.01),seq(1,0,-.01),0)
infer <- function(grid) {
  pred <- pf %*% t(grid)
  ll <- t(gc)%*%log(pred)+t(2-gc)%*%log(1-pred)
  grid[max.col(ll,ties.method='first'),,drop=FALSE]
}
clustering <- list(full=infer(grid3),reduced=infer(grid2),group=truegroup,M=nm,N=30,
                  reference_frequencies=pf, genotypes=gc)
out <- list(running=running,ld=ld,batch=batchout,projection=projection,
 small=small,ldsc=ldsc,clustering=clustering,
 provenance=list(seed=20260907, module_seeds=c(main=20260907, LD=901, batch=902, projection=903, LDSC=904, ancestry=905), R=R.version.string,
 bigsnpr=as.character(packageVersion('bigsnpr')),bigstatsr=as.character(packageVersion('bigstatsr')),
 source='population_structure/R/simulate_population_structure.R', notes='All teaching simulations are synthetic. No participant data were downloaded.'))
write_json(out,file.path(B,'computed.json'),auto_unbox=TRUE,digits=10,pretty=FALSE)
saveRDS(out,file.path(B,'computed.rds'))
cat('Running null before/after:',fit0[[1]]$beta,fit0[[2]]$beta,'; PC1 R2',running$r2pc,'\n')
cat('LD before/after ancestry R2:',ld$rawAncestryR2,ld$cleanAncestryR2,'retained',ld$retained,'\n')
cat('Batch PC R2:',batchout$r2batch,'; projection',names(proj),'\n')
cat('LDSC intercepts:',ldsc$confFit[1],ldsc$polyFit[1],'\n')

# Preserve synthetic inputs for practical exercises and independent checking.
saveRDS(list(genotypes=G, ancestry=a, phenotype=pi, null_genotype=g0,
             causal_genotype=g1, background_frequencies=freq),
        file.path(B, 'synthetic_inputs.rds'), compress='xz')
write.csv(data.frame(ld_score=ell, z=z, confounded=tc, polygenic=tp),
          gzfile(file.path(B, 'ldsc_statistics.csv.gz')), row.names=FALSE)
write.csv(data.frame(id=seq_len(N), group=a, null_dosage=g0, causal_dosage=g1, phenotype=pi),
          file.path(B, 'individuals.csv'), row.names=FALSE)
writeLines(capture.output(sessionInfo()), file.path(B, 'sessionInfo.txt'))
