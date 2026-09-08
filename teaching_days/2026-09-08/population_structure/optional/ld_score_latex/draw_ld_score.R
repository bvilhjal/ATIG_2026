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

# Coordinates are in points on a 1600 x 1000 canvas, measured from the top.
W <- 1600; H <- 1000
ink <- "#16313F"; muted <- "#53656F"; pale <- "#E5EBEE"
blue <- "#0072B2"; orange <- "#E69F00"; soft <- "#F3F7F9"
txt <- function(label, x, y, size=22, col=ink, face="plain", just="left", rot=0) {
  grid.text(label, x=unit(x,"pt"), y=unit(H-y,"pt"), just=just, rot=rot,
            gp=gpar(fontfamily="Arial", fontsize=size, col=col, fontface=face))
}
rect <- function(x,y,w,h,fill=NA,col=NA,lwd=1,lty=1) {
  grid.rect(x=unit(x,"pt"),y=unit(H-y,"pt"),width=unit(w,"pt"),height=unit(h,"pt"),
            just=c("left","top"),gp=gpar(fill=fill,col=col,lwd=lwd,lty=lty))
}
line <- function(x1,y1,x2,y2,col=pale,lwd=1,lty=1) {
  grid.lines(x=unit(c(x1,x2),"pt"),y=unit(H-c(y1,y2),"pt"),gp=gpar(col=col,lwd=lwd,lty=lty))
}
palette <- colorRampPalette(c("#F3F7F9", "#B9D9E9", "#5FA9CF", blue, "#003E68"))(256)
ld_col <- function(x) palette[1L+round(x*255)]

if (identical(unname(Sys.info()["sysname"]),"Darwin") && capabilities("aqua")) {
  quartz(type="pdf",file=file.path(out,"base.pdf"),width=W/72,height=H/72,pointsize=22,family="Arial",bg="white")
} else {
  cairo_pdf(filename=file.path(out,"base.pdf"),width=W/72,height=H/72,pointsize=22,family="sans",bg="white")
}
grid.newpage()
txt("LD score sums a SNP’s LD with its neighbours",90,63,44,face="bold")
txt("Strong LD with many nearby SNPs gives a larger score.",90,116,25,col=muted)
line(90,155,1510,155,lwd=1.5)
txt("A   LD heatmap",90,192,27,face="bold")
txt("12 SNPs in genomic order; cells show r²",90,229,22,col=muted)
txt("B   Contributions to SNP 4’s LD score",795,192,27,face="bold")

# Heatmap. The orange row and square preserve the same r-squared colour scale.
hx <- 128; hy <- 283; cell <- 42; side <- 12*cell
for (i in 1:12) for (m in 1:12) {
  rect(hx+(m-1)*cell, hy+(i-1)*cell, cell, cell, fill=ld_col(R2[i,m]), col="white", lwd=.6)
}
for (i in 1:12) {
  txt(as.character(i),hx-20,hy+(i-.5)*cell,20,
      col=if(i==4) "#B77B00" else muted, face=if(i==4) "bold" else "plain",just="right")
  txt(as.character(i),hx+(i-.5)*cell,hy-18,20,col=muted,just="centre")
}
txt("SNP",hx-21,hy-18,18,col=muted,just="right")
rect(hx,hy+3*cell,side,cell,col=orange,lwd=3)
rect(hx+3*cell,hy+3*cell,cell,cell,col=orange,lwd=4)
txt("1",hx+3.5*cell,hy+3.5*cell,22,col="white",face="bold",just="centre")
rect(hx+11*cell,hy+11*cell,cell,cell,col=muted,lwd=2.5)
txt("1",hx+11.5*cell,hy+11.5*cell,22,col="white",face="bold",just="centre")
grid.lines(x=unit(c(hx+side+7,hx+side+46),"pt"),
           y=unit(rep(H-(hy+3.5*cell),2),"pt"),
           arrow=arrow(length=unit(8,"pt"),type="closed"), gp=gpar(col=orange,lwd=2.2))
txt("SNP 4",hx+side+20,hy+3.5*cell-29,19,col="#B77B00",face="bold",just="centre")

# Compact colour key aligned under the heatmap.
cx <- hx; cy <- 818; cw <- 310
for (k in 1:256) rect(cx+(k-1)*cw/256,cy,cw/256+.25,13,fill=palette[k])
for (a in c(0,.5,1)) txt(format(a,trim=TRUE),cx+a*cw,850,18,col=muted,just="centre")
txt(expression(r^2),cx+cw+26,825,23)
txt("Diagonal: each SNP is perfectly correlated with itself.",90,892,22,col=muted)

# Bar heights exactly match the highlighted heatmap row. Self-LD is orange.
bx <- 835; by <- 335; bw <- 660; bh <- 275; base <- by+bh
step <- bw/12
for (v in seq(0,1,by=.2)) {
  yy <- base-v*bh
  line(bx,yy,bx+bw,yy,col=pale,lwd=1)
  txt(sprintf("%.1f",v),bx-15,yy,18,col=muted,just="right")
}
for (m in 1:12) {
  value <- R2[4,m]; xc <- bx+(m-.5)*step
  if (value>0) rect(xc-17,base-value*bh,34,value*bh,
                    fill=if(m==4) orange else blue)
  if (value==0) line(xc-12,base,xc+12,base,col="#B8C5CD",lwd=3)
  txt(sprintf("%.2f",value),xc,base-value*bh-15,18,
      col=if(m==4) "#A46E00" else ink, face=if(m==4) "bold" else "plain",just="centre")
  txt(as.character(m),xc,base+26,20,
      col=if(m==4) "#A46E00" else muted,face=if(m==4) "bold" else "plain",just="centre")
}
line(bx,base-.2*bh,bx+bw,base-.2*bh,col=muted,lwd=1.5,lty="dashed")
rect(bx+7.0*step,base-.2*bh-39,270,27,fill="white")
txt("Example cutoff: 0.20",bx+7.2*step,base-.2*bh-26,20,col=muted)
txt("SNP",bx+bw/2,base+62,21,col=muted,just="centre")
txt(expression(r[4*m]^2),bx-57,by+bh/2,23,just="centre",rot=90)

# LaTeX equations are added to the blank regions by build_figure.py.
rect(785,693,735,88,fill="#F5F8FA")
txt("Illustrative LD-friend count: 4 other SNPs",805,820,24,face="bold")
txt("SNPs 2, 3, 5 and 6 have r² ≥ 0.20.",805,854,23,col=muted)
txt("Weak LD also contributes to the score.",805,892,22,col=muted)
line(90,922,1510,922,lwd=1)
txt("Figure 28a. Illustrative population LD in one window; the score includes the focal SNP.",90,947,21,col=muted)
txt("Bulik-Sullivan et al. (2015a), Nature Genetics  •  LD friends: GCTA documentation",90,980,19,col=muted)
dev.off()
cat(sprintf("SNP 4 LD score = %.2f; other friends at r² >= %.2f: %s\n",score[4],cutoff,paste(friends,collapse=", ")))
