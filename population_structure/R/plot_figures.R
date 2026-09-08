#!/usr/bin/env Rscript
# Portable publication/lecture figures. Data are calculated in separate modules.
# R's plotmath renders mathematical labels; the optional LaTeX workflow is only
# needed to reproduce the original LD heatmap's exact typography and layout.
args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[[1]] else "outputs"
figdir <- file.path(out,"figures")
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
d <- readRDS(file.path(out,"population_structure/computed.rds"))
sm <- readRDS(file.path(out,"small_calculations/calculations.rds"))
gc <- jsonlite::read_json(file.path(out,"genomic_control/simulation.json"),simplifyVector=TRUE)
ld <- as.matrix(read.csv(file.path(out,"ld_score/population_squared_correlation.csv")))
colnames(ld) <- rownames(ld) <- as.character(1:12)
ink <- "#16313F"; gray <- "#53656F"; pale <- "#E5EBEE"
blue <- "#0072B2"; orange <- "#E69F00"; green <- "#009E73"; pink <- "#CC79A7"
pal <- c(blue,orange,green,pink,"#56B4E9")
manifest <- list()
use_quartz <- identical(unname(Sys.info()["sysname"]),"Darwin") && capabilities("aqua")
pdftocairo <- Sys.getenv("PDFTOCAIRO",Sys.which("pdftocairo"))
if(use_quartz && !nzchar(pdftocairo)) stop("On macOS, install Poppler (pdftocairo) for SVG export, or set PDFTOCAIRO.")

figure <- function(id,name,title,caption,draw,width=14,height=8,kind="original synthetic illustration") {
  stem <- paste0("F",id,"_",name)
  for(ext in if(use_quartz) c("png","pdf") else c("png","pdf","svg")) {
    file <- file.path(figdir,paste0(stem,".",ext))
    if(ext=="png") png(file,width=width,height=height,units="in",res=160,type=if(use_quartz) "quartz" else "cairo",bg="white")
    if(ext=="pdf") {
      if(use_quartz) quartz(type="pdf",file=file,width=width,height=height,family="Arial",bg="white")
      else cairo_pdf(file,width=width,height=height,family="sans",bg="white")
    }
    if(ext=="svg") svg(file,width=width,height=height,family="sans",bg="white")
    par(family="sans",fg=ink,col.axis=gray,col.lab=ink,col.main=ink,
        cex=1.2,las=1,mar=c(4.3,4.6,1.5,1.2),oma=c(3.4,0.7,4.2,0.7),mgp=c(2.8,.7,0))
    draw()
    mtext(title,side=3,outer=TRUE,line=1.3,adj=0,font=2,cex=1.5,col=ink)
    mtext(paste(strwrap(paste0("Figure ",as.integer(id),". ",caption),width=145),collapse="\n"),side=1,outer=TRUE,line=1.6,adj=0,cex=.88,col=gray)
    dev.off()
  }
  if(use_quartz) {
    status <- system2(pdftocairo,shQuote(c("-svg",file.path(figdir,paste0(stem,".pdf")),file.path(figdir,paste0(stem,".svg")))))
    if(status!=0) stop("SVG conversion failed for ",stem)
  }
  manifest[[length(manifest)+1L]] <<- list(id=paste0("F",id),file=stem,title=title,caption=caption,kind=kind)
}
scatter <- function(x,y,group,xlab,ylab,...) {
  plot(x,y,pch=16,cex=.62,col=adjustcolor(pal[group+1L],alpha.f=.65),
       xlab=xlab,ylab=ylab,bty="l",...)
}
groupkey <- function(labels=c("Group A","Group B"),cols=pal[1:2],where="topright") {
  legend(where,labels,col=cols,pch=16,bty="n",cex=.9)
}
heat <- function(z,labels=rownames(z),range=base::range(z),numbers=FALSE) {
  n <- nrow(z); m <- ncol(z)
  colors <- colorRampPalette(c("#F3F7F9","#B9D9E9",blue,"#003E68"))(256)
  plot.new(); plot.window(c(.5,m+.5),c(n+.5,.5),asp=1)
  for(i in seq_len(n)) for(j in seq_len(m)) {
    cc <- 1L+round(255*(z[i,j]-range[1])/diff(range))
    rect(j-.5,i-.5,j+.5,i+.5,col=colors[max(1,min(256,cc))],border="white")
    if(numbers) text(j,i,sprintf("%.2f",z[i,j]),cex=.85,
                     col=if(z[i,j]>mean(range)) "white" else ink)
  }
  axis(1,seq_len(m),labels=if(length(labels)==m) labels else seq_len(m),tick=FALSE)
  axis(2,seq_len(n),labels=labels,tick=FALSE)
  box(col=pale)
}
blank <- function() {plot.new();plot.window(c(0,1),c(0,1));}
edge <- function(x,y,x2,y2,col=gray,lty=1) arrows(x,y,x2,y2,length=.12,lwd=2,col=col,lty=lty)
chromosomes <- function(x=.17,y=.65,w=.7,h=.075,rows=list(c(1,1,1,1,1,2,2,2),c(1,1,1,2,2,2,2,2))) {
  for(k in 1:2) {
    yy <- y-(k-1)*.18
    for(b in seq_along(rows[[k]])) rect(x+(b-1)*w/8,yy,x+b*w/8,yy+h,
                                       col=pal[rows[[k]][b]],border="white")
    text(x-.035,yy+h/2,paste("Copy",k),adj=1)
  }
}

figure("01","standardized_pca","A single genotype example links the PCA calculations",
       "Seven SNPs, five individuals; binomial-standardized genotypes. Colors identify individuals.",function(){
  par(mfrow=c(1,2))
  heat(d$small$psi,as.character(1:5),numbers=TRUE)
  title(expression(Psi == X^T*X/M),line=1)
  a <- d$small$geometry; axes <- d$small$axes; lim <- max(abs(a))*1.35
  plot(a,pch=16,col=pal,cex=1.2,asp=1,xlim=c(-lim,lim),ylim=c(-lim,lim),
       xlab="Rotated leading-PC plane",ylab="",axes=FALSE)
  abline(h=0,v=0,col=pale)
  for(k in 1:2) {
    edge(-lim*.8*axes[1,k],-lim*.8*axes[2,k],lim*.8*axes[1,k],lim*.8*axes[2,k],ink)
    text(lim*.9*axes[1,k],lim*.9*axes[2,k],paste0("PC",k),font=2)
  }
  for(j in 1:5) {
    xy <- d$small$score[j,1]*axes[,1]
    segments(a[j,1],a[j,2],xy[1],xy[2],lty=2,col=pal[j])
    text(a[j,1],a[j,2],j,pos=3,col=pal[j],font=2)
  }
  title(sprintf("PC1: %.1f%%; PC2: %.1f%% of variation",100*d$small$explained[1],100*d$small$explained[2]))
},kind="deterministic worked example")

figure("02","ancestry_likelihood","A four-SNP likelihood for an ancestry proportion",
       "Known reference allele frequencies and four genotypes; independent SNPs and Hardy–Weinberg equilibrium assumed.",function(){
  par(mfrow=c(1,2))
  blank()
  text(.02,.92,"Reference allele frequencies",adj=0,font=2)
  text(.18,.77,"SNP");text(.43,.77,"Reference 1");text(.68,.77,"Reference 2");text(.9,.77,"Dosage")
  for(i in 1:4) { yy <- .67-(i-1)*.13; text(.18,yy,i);text(.43,yy,sm$pa[i]);text(.68,yy,sm$pb[i]);text(.9,yy,sm$g[i]) }
  text(.05,.08,sprintf("L₁ = %.6f; L₂ = %.6f; ratio = %.2f",sm$L1,sm$L2,sm$ratio),adj=0,cex=.9)
  plot(sm$alpha,exp(sm$loglik-max(sm$loglik)),type="l",lwd=3,col=blue,bty="l",
       xlab=expression(alpha[j]),ylab="Likelihood relative to its maximum",ylim=c(0,1.15))
  abline(v=sm$alpha_hat,lty=2,col=orange,lwd=2)
  text(sm$alpha_hat,.1,sprintf("Maximum at %.3f",sm$alpha_hat),pos=4,col=ink)
},kind="deterministic worked example")

figure("03","reference_panels","The reference panel changes the inferred components",
       "The same 30 synthetic genomes are fitted against three reference profiles (top) or two (bottom).",function(){
  par(mfrow=c(2,1),mar=c(2.4,4.6,2,1))
  for(k in 1:2) {
    q <- if(k==1) d$clustering$full else d$clustering$reduced
    barplot(t(q),col=pal[1:3],border=NA,space=.07,ylim=c(0,1),ylab="Proportion",
            main=if(k==1) "References A + B + C" else "References A + B",cex.main=1)
    axis(1,at=c(5.35,16.05,26.75),labels=c("Profile A","Profile B","Profile C"),tick=FALSE)
  }
})

figure("04","local_global_ancestry","Local ancestry labels sum to global ancestry",
       "A deterministic chromosome schematic; lengths, not the number of inferred tracts, determine global ancestry.",function(){
  blank(); chromosomes()
  legend("topright",c("Reference 1","Reference 2"),fill=pal[1:2],bty="n")
  text(.5,.28,expression(alpha[j] == frac(1,2*L)*sum(sum(L[b]*I(A[jhb]==1),b==1,B),h==1,2)),cex=1.8)
  text(.5,.10,"8 of 16 equally long copy-segments: α = 0.50",font=2,cex=1.35)
},kind="original deterministic schematic")

figure("05","stratification_diamond","Population structure can produce a spurious SNP association",
       "Original schematic for the running simulation. The tested SNP has zero direct effect; a separate causal SNP affects the phenotype.",function(){
  blank()
  text(.5,.92,"Group / ancestry",font=2,cex=1.4)
  text(.13,.51,"Null SNP dosage",font=2,cex=1.3)
  text(.85,.51,"Environmental shift",font=2,cex=1.3)
  text(.5,.13,expression(Phenotype~pi),font=2,cex=1.4)
  edge(.41,.84,.15,.60,blue);edge(.59,.84,.84,.60,orange);edge(.84,.41,.59,.19,orange)
  segments(.17,.4,.42,.2,lty=2,lwd=2,col=gray)
  text(.23,.26,"No direct effect\nβᵢ = 0",cex=.95,col=gray)
  text(.13,.70,"Allele frequency:\n0.20 in A; 0.60 in B",cex=.95)
  text(.82,.27,"Mean difference: 1.2",cex=.95)
},kind="original deterministic schematic")

r <- d$running
figure("06","pooled_association","A null SNP appears associated in the pooled sample",
       "600 synthetic individuals; horizontal jitter is for display only. The true tested-SNP effect is zero.",function(){
  scatter(r$genotype+r$jitter,r$phenotype,r$ancestry,"Null SNP dosage (jittered)",expression(Phenotype~pi))
  groupkey(); f <- r$null[[1]]
  legend("topleft",sprintf("β̂ = %.3f\n95%% CI [%.3f, %.3f]\nP = %.2g",f$beta,f$low,f$high,f$p),bty="n")
})
figure("07","background_pca","PC1 captures the simulated ancestry contrast",
       "PCA uses 3,000 background SNPs; both tested SNPs are excluded. All 600 individuals are shown.",function(){
  scatter(r$pc[,1],r$pc[,2],r$ancestry,expression(v[1]~coordinate),expression(v[2]~coordinate))
  groupkey();mtext(sprintf("R²(PC1, group) = %.3f",r$r2pc),side=3,line=.3)
})
figure("08","residual_association","After PC adjustment, compare residual with residual",
       "Both dosage and phenotype are residualized on an intercept and PC1, using the same 600 individuals.",function(){
  scatter(r$residualG,r$residualPi,r$ancestry,"Residual SNP dosage","Residual phenotype")
  groupkey();f <- r$null[[2]]
  legend("topleft",sprintf("β̂ = %.3f\n95%% CI [%.3f, %.3f]\nP = %.2g",f$beta,f$low,f$high,f$p),bty="n")
})
figure("09","known_effects","Compare estimates with the effects used to generate the data",
       "A single simulation: estimates and 95% t intervals before and after PC1 adjustment; vertical lines show the true effects.",function(){
  par(mfrow=c(1,2),mar=c(4.3,7,2,1))
  for(k in 1:2) {
    a <- if(k==1) r$null else r$causal
    plot(NA,xlim=c(-.3,.65),ylim=c(.5,2.5),yaxt="n",xlab=expression(hat(beta)[i]),ylab="",bty="l",
         main=if(k==1) "Null SNP: βᵢ = 0" else "Causal SNP: βᵢ = 0.25")
    abline(v=if(k==1) 0 else .25,col=gray,lty=2)
    axis(2,1:2,c("PC1 adjusted","Unadjusted"),tick=FALSE)
    for(j in 1:2) { y <- 3-j; segments(a[[j]]$low,y,a[[j]]$high,y,col=pal[c(2,1)[j]],lwd=3);points(a[[j]]$beta,y,pch=16,col=pal[c(2,1)[j]],cex=1.4) }
  }
})
figure("10","pc_sensitivity","How many PCs? Inspect sensitivity of the effect estimate",
       "Null-SNP point estimates and 95% t intervals. Stability in this example does not determine K for another study.",function(){
  a <- do.call(rbind,lapply(r$null,as.data.frame))
  plot(a$K,a$beta,type="b",pch=16,col=blue,lwd=2,ylim=range(a$low,a$high),bty="l",
       xlab="Number of PCs K",ylab=expression(hat(beta)[i]))
  segments(a$K,a$low,a$K,a$high,col=blue,lwd=2);abline(h=0,col=gray,lty=2)
})

figure("11","ld_pca","A large LD region can dominate the first PC",
       "Synthetic correlated block at SNPs 1,201–1,700; big_randomSVD before filtering, snp_autoSVD after filtering.",function(){
  par(mfrow=c(1,2),mar=c(4.3,4.6,3.4,1.2))
  scatter(d$ld$raw[,1],d$ld$raw[,2],d$ld$ancestry,expression(v[1]),expression(v[2]))
  title("Before LD filtering",line=1.7)
  groupkey(where="center");mtext(sprintf("R²(PC1, group) = %.3f",d$ld$rawAncestryR2),side=3,line=.3,cex=.85)
  scatter(d$ld$clean[,1],d$ld$clean[,2],d$ld$ancestry,expression(v[1]),expression(v[2]))
  title("After LD filtering",line=1.7)
  groupkey(where="center");mtext(sprintf("R² = %.3f; %d SNPs retained",d$ld$cleanAncestryR2,d$ld$retained),side=3,line=.3,cex=.85)
})
figure("12","ld_loadings","Inspect SNP loadings before interpreting a PC",
       "All 3,000 PC1 loadings from the unfiltered synthetic example; the highlighted region was deliberately made correlated.",function(){
  lim <- range(d$ld$rawLoad)
  plot(seq_along(d$ld$rawLoad),d$ld$rawLoad,type="n",xlab="SNP in genomic order",ylab="PC1 SNP loading",bty="l")
  rect(1200,lim[1],1700,lim[2],col=adjustcolor(orange,.15),border=NA)
  points(d$ld$rawLoad,pch=16,cex=.4,col=blue)
})
figure("13","batch_pca","The same PCs can track ancestry and technical batch",
       "Identical coordinates in both panels; only the coloring changes. Batch was assigned independently of ancestry.",function(){
  par(mfrow=c(1,2))
  scatter(d$batch$pc[,1],d$batch$pc[,2],d$batch$ancestry,expression(v[1]),expression(v[2]),main="Colored by group");groupkey()
  scatter(d$batch$pc[,1],d$batch$pc[,2],d$batch$batch,expression(v[1]),expression(v[2]),main="Colored by batch");groupkey(c("Batch 1","Batch 2"))
})
figure("14","oadp_projection","Fit the axes once, then project held-out individuals",
       "snp_autoSVD fits 480 reference individuals; snp_projectSelfPCA supplies OADP scores for 120 held-out individuals.",function(){
  a <- d$projection
  plot(a$reference,pch=16,cex=.6,col="#B3BEC5",xlab="PC1 score",ylab="PC2 score",bty="l")
  points(a$target,pch=18,cex=.9,col=green)
  groupkey(c("Reference (480)","Projected (120)"),c("#B3BEC5",green))
})

figure("15","genomic_control","Genomic control rescales the same association statistics",
       "100,000 synthetic null statistics. The corrected median equals the null median by construction; the same ranks are plotted.",function(){
  par(mfrow=c(1,2));a <- gc$coordinates
  lim <- c(0,ceiling(max(a$before,a$expected)))
  for(k in 1:2) {
    plot(a$expected,if(k==1) a$before else a$after,pch=16,cex=.5,col=if(k==1) orange else blue,
         xlim=lim,ylim=lim,asp=1,bty="l",xlab=expression(Expected~chi[1]^2),ylab=expression(Observed~chi[1]^2),
         main=if(k==1) sprintf("Before: λGC = %.3f",gc$lambda_before) else "After: λGC = 1.000")
    abline(0,1,col=gray,lty=2)
  }
})
figure("16","matched_qq","Similar Q–Q inflation can have different causes",
       "100,000 matched normal draws under two LDSC working-model scenarios; both have expected mean χ² = 1.4.",function(){
  par(mfrow=c(1,2))
  for(k in 1:2) {
    a <- if(k==1) d$ldsc$confQQ else d$ldsc$polyQQ
    plot(a$x,a$y,pch=16,cex=.6,col=pal[k],xlim=c(0,5.4),ylim=c(0,8),bty="l",
         xlab=expression(Expected~-log[10](P)),ylab=expression(Observed~-log[10](P)),main=paste("Scenario",c("A","B")[k]))
    abline(0,1,col=gray,lty=2)
  }
})
figure("17","ldsc_scenarios","LD-score relationships separate these model scenarios",
       "The same statistics as Figure 16. Points are means in 20 bins; lines show generating means, not a fitted real GWAS.",function(){
  par(mfrow=c(1,2))
  for(k in 1:2) {
    plot(d$ldsc$ell,if(k==1) d$ldsc$conf else d$ldsc$poly,pch=16,col=pal[k],ylim=c(.8,2),xlim=c(0,100),bty="l",
         xlab="LD score",ylab=expression(Mean~chi[1]^2),main=c("A: LD-independent inflation","B: LD-dependent component")[k])
    if(k==1) abline(h=1.4,col=gray,lwd=2) else abline(1,.008,col=gray,lwd=2)
  }
})
figure("18","ld_score_heatmap","LD score sums squared correlations with neighboring SNPs",
       "Exact illustrative population LD. SNP 4: score 3.35, including self-LD 1; four other SNPs exceed the example cutoff 0.20.",function(){
  par(mfrow=c(1,2))
  heat(ld,1:12,range=c(0,1));rect(.5,3.5,12.5,4.5,border=orange,lwd=3)
  title(expression(A:~r[im]^2~heatmap),line=1)
  bp <- barplot(ld[4,],col=replace(rep(blue,12),4,orange),border=NA,ylim=c(0,1.16),
                xlab="Neighbor SNP m",ylab=expression(r[4*m]^2),names.arg=1:12,main="B: Contributions from the highlighted row")
  abline(h=.2,col=gray,lty=2);text(bp,ld[4,]+.04,sprintf("%.2f",ld[4,]),cex=.75)
  mtext("LD score = 1 (self) + 2.35 (other SNPs) = 3.35",side=3,line=-1.3,cex=.9)
},kind="deterministic valid population-LD model")

figure("19","pedigree_relationships","Relatives create correlated observations",
       "An expected additive relationship matrix for two unrelated, non-inbred parents and two full siblings; realized genomic values vary.",function(){
  par(mfrow=c(1,2));blank()
  points(c(.25,.75),c(.78,.78),pch=c(15,16),cex=4,col=pal[1:2])
  segments(.25,.78,.75,.78);segments(.5,.78,.5,.5);segments(.25,.5,.75,.5)
  segments(c(.25,.75),.5,c(.25,.75),.25);points(c(.25,.75),c(.25,.25),pch=c(15,16),cex=4,col=pal[3:4])
  text(c(.25,.75),.93,c("Father","Mother"));text(c(.25,.75),.1,c("Child 1","Child 2"))
  K <- matrix(.5,4,4);diag(K)<-1;K[1,2]<-K[2,1]<-0
  heat(K,c("F","M","C1","C2"),range=c(0,1),numbers=TRUE)
},kind="original deterministic schematic")
figure("20","parental_transmission","Parental transmission and family environment are different pathways",
       "Original teaching DAG. Conditioning on parental genotypes uses within-family segregation; causal interpretation still requires assumptions.",function(){
  blank();text(.5,.92,"Parental genotypes",font=2,cex=1.4)
  text(.14,.5,"Transmitted alleles\nin the child",font=2,cex=1.2)
  text(.86,.5,"Family environment",font=2,cex=1.2)
  text(.5,.10,expression(Child~phenotype~pi),font=2,cex=1.4)
  edge(.40,.83,.14,.62,blue);edge(.60,.83,.86,.62,orange)
  edge(.14,.36,.40,.17,blue);edge(.86,.36,.60,.17,orange)
  text(.19,.22,"Direct pathway",col=blue);text(.83,.22,"Indirect pathway",col=orange)
},kind="original deterministic schematic")
figure("21","local_ancestry_model","Local ancestry supplies information beyond a global proportion",
       "Original two-copy schematic for an ancestry-specific association model; it is not output from Tractor-Mix or observed genotype data.",function(){
  blank();chromosomes(rows=list(c(2,2,1,1,1,1,1,2),c(2,2,2,1,1,1,2,2)))
  xx <- .17+6.5*.7/8
  segments(xx,.4,xx,.85,lty=2,lwd=2);text(xx,.9,"Tested SNP",font=2)
  text(.25,.20,"Ancestry-specific dosage",font=2);edge(.48,.20,.60,.20,ink)
  text(.8,.20,"Ancestry-specific effects\n+ a relatedness model",font=2)
},kind="original deterministic schematic")

rare <- read.csv("population_structure/data/rare_variant_inflation.csv")
moba <- read.csv("population_structure/data/moba_attenuation.csv")
figure("22","rare_variant_summary","Rare-variant structure can survive common-PC adjustment",
       "Published summary values, not simulations: Hanson et al. (2026), Nature Communications, doi:10.1038/s41467-026-73776-9.",function(){
  par(mar=c(4.3,13,2,1));yy <- 5:1
  plot(rare$lambda,yy,pch=16,col=blue,yaxt="n",ylim=c(.5,5.5),xlim=c(0,7.5),bty="l",xlab=expression(lambda[GC]),ylab="")
  axis(2,yy,rare$adjustment,tick=FALSE);abline(v=1,lty=2,col=gray)
  text(rare$lambda,yy,sprintf("%.2f",rare$lambda),pos=4)
},kind="redrawn published numerical summary")
figure("23","moba_summary","Parental adjustment attenuates some SNP associations",
       "Published estimates and 95% CIs, not simulations: Corfield et al. (2026), Nature, doi:10.1038/s41586-026-10926-5.",function(){
  par(mar=c(4.3,13,2,1));yy <- 4:1
  plot(moba$attenuation_percent,yy,pch=16,col=pal[1:4],yaxt="n",ylim=c(.5,4.5),xlim=c(0,55),bty="l",
       xlab="Attenuation of SNP effect estimates (%)",ylab="")
  axis(2,yy,paste0(moba$trait,"\nN = ",format(moba$n,big.mark=",")),tick=FALSE)
  segments(moba$lower_95,yy,moba$upper_95,yy,col=pal[1:4],lwd=3)
  points(moba$attenuation_percent,yy,pch=16,col=pal[1:4],cex=1.3)
},kind="redrawn published numerical summary")

jsonlite::write_json(manifest,file.path(figdir,"figure_manifest.json"),auto_unbox=TRUE,pretty=TRUE)
cat(length(manifest),"figures written in PNG, PDF and SVG formats.\n")
