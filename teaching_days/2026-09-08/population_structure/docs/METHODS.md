# Models and notation

Let there be $M$ SNPs and $N$ individuals. SNP $i$ is tested, $m$ indexes another
SNP, $j$ indexes an individual, and $k$ indexes a principal component. This is
the lecture's notation. The lecture matrices have SNPs in rows and individuals
in columns. The `bigsnpr` API stores individuals in rows; that transpose is
intentional.

## Standardized genotype PCA

With allele dosage $g_{ij}$ and estimated counted-allele frequency $p_i$, write

$$x_{ij}=\frac{g_{ij}-2p_i}{\sqrt{2p_i(1-p_i)}},\qquad
\boldsymbol\Psi=\frac{1}{M}\mathbf X^{\mathsf T}\mathbf X.\tag{1}$$

$\boldsymbol\Psi$ (capital Greek **Psi**) is the genetic relationship Gram
matrix. The spectral decomposition, also called eigendecomposition, is

$$\boldsymbol\Psi=\mathbf V\mathbf D\mathbf V^{\mathsf T},\qquad
\mathbf V^{\mathsf T}\mathbf V=\mathbf I.\tag{2}$$

The columns $\mathbf v_k$ are unit eigenvectors; individual PC scores are
$\sqrt{M d_k}\mathbf v_k$. Their signs are arbitrary. The seven-SNP,
five-individual example is standardized exactly as in equation (1). Its GRM
diagonal need not be exactly one: binomial scaling is not per-individual sample
variance normalization, and the sample is small. The original centered-only
version remains only as an explicitly historical snapshot.

The additional four-SNP raw-vector correlation table illustrates why correlating
each individual's raw dosage vector across SNPs is different from genotype PCA.
An allele-coding reversal leaves the centered Gram matrix unchanged. The
two-individual appendix example also has an analytic eigendecomposition, checked
against the numerical matrix.

Background: [Patterson et al. (2006), PLOS Genetics](https://doi.org/10.1371/journal.pgen.0020190).

## Running stratification example

Seed **20260907** generates 600 individuals, equally split between two groups,
and 3,000 background SNPs. Baseline frequencies are uniform on $[0.2,0.8]$;
each background SNP differs by $+0.13$ or $-0.13$ between groups. Conditional on
group, dosages are independent binomial draws. PCA uses binomial scaling and
`bigstatsr::big_randomSVD` with 20 components.

The tested null SNP has frequencies 0.20 and 0.60 in the two groups. A separate
causal SNP has frequency 0.40 in both groups and effect 0.25. Neither tested SNP
enters PCA. Phenotype is generated as

$$\pi_j=1.2(a_j-0.5)+0.25g_{\mathrm{causal},j}+\varepsilon_j,
\qquad \varepsilon_j\stackrel{\mathrm{iid}}\sim N(0,1).\tag{3}$$

Here $a_j\in\{0,1\}$ is a simulation group indicator; it is not the ancestry
mixture proportion $\alpha_j$. Each association fit includes an intercept, the
tested dosage, and $K\in\{0,1,2,5,10,20\}$ PCs. Confidence intervals use the
residual degrees of freedom $N-K-2$. Residual plots adjust both the dosage and
phenotype for an intercept and PC1. Horizontal dosage jitter changes display
positions only.

This is an original teaching simulation inspired by the adjustment principle in
[Price et al. (2006), Nature Genetics](https://doi.org/10.1038/ng1847), not a
reproduction of that paper's data or simulation study.

## PCA diagnostics and projection

1. **LD-driven PC, seed 901.** SNPs 1,201–1,700 mostly copy one genotype vector
   generated independently of group. Each copied entry retains its background
   genotype with probability 0.015. Compare `big_randomSVD` with
   `bigsnpr::snp_autoSVD`; the latter performs LD clumping and long-range-LD
   detection. The retained SNP indices are saved.
2. **Technical batch, seed 902.** A balanced random batch indicator is generated
   independently of group. Conditional binomial frequencies shift by 0.15
   between batches, clipped to $[0.02,0.98]$. The same PCs are colored by group
   and batch. Association with a real batch label would be diagnostic, not by
   itself proof of a technical artifact.
3. **Held-out projection, seed 903.** A balanced reference subset of 480
   individuals fits `snp_autoSVD`; 120 distinct individuals are projected with
   `snp_projectSelfPCA` using `OADP_proj`. Scores rather than unit eigenvectors
   are plotted. All genotypes come from the same synthetic matrix, so allele
   matching is exact. An external cohort would additionally require allele and
   variant harmonization.

Software reference: [bigsnpr PCA documentation](https://privefl.github.io/bigsnpr/articles/bedpca.html).

## Ancestry likelihood and chromosomes

The four-SNP calculation evaluates a binomial likelihood under two known
reference-frequency profiles. For a mixture proportion $\alpha_j$, the modeled
allele frequency is $\alpha_jp_{i1}+(1-\alpha_j)p_{i2}$. Full binomial likelihoods
include the combinatorial coefficients; those coefficients cancel in ratios
and optimization over ancestry. The optimum is approximately 0.217655.

Seed **905** generates 30 diploid genomes at 1,200 independent SNPs from three
frequency profiles, ten genomes per profile. Profile C has intermediate A/B
frequencies plus independent perturbations. The same genomes are fitted to a
three-profile simplex grid (spacing 0.05) or an A/B grid (spacing 0.01). This is
a supervised likelihood demonstration, **not STRUCTURE, ADMIXTURE or FRAPPE
output**. A mixed fitted bar need not establish a recent admixture event.

The chromosome illustrations are deterministic schematics. Their colors refer
to chosen reference populations. Global ancestry averages assigned local
ancestry over both chromosome copies, weighted by analyzed segment length.
They are not empirical ancestry reconstructions. Saved image-generation prompts
describe the separate raster illustrations; those outputs are not deterministic
and are not evidence of genotype measurements.

## Genomic control and LDSC working model

Seed **904** first generates $\ell_i\sim U(5,95)$, then $z_i\sim N(0,1)$ for
100,000 independent draws. This order deliberately matches the lecture's
original random stream. The two scenarios are

$$T_i^{(A)}=1.4z_i^2,\qquad
T_i^{(B)}=(1+0.008\ell_i)z_i^2.\tag{4}$$

Both have expected mean statistic 1.4. Figure 16 plots selected empirical
quantiles; Figure 17 plots means in 20 equal-width LD-score bins. Lines show
the generating conditional means. Scenario A illustrates LD-independent
inflation; scenario B illustrates an LD-dependent component under the LDSC
working model. There are no underlying genotype data for these two scenarios.

Genomic control uses the exact same scenario-A draws:

$$\lambda_{\mathrm{GC}}=\frac{\operatorname{median}_i(T_i)}
{\operatorname{median}(\chi_1^2)},\qquad
T_i^{\mathrm{GC}}=T_i/\lambda_{\mathrm{GC}}.\tag{5}$$

The corrected median gives $\lambda_{\mathrm{GC}}=1$ by construction. Rank order
is unchanged. Q–Q display thinning retains all 200 largest statistics; estimation
uses all 100,000 draws. Uniform null inflation is the illustration's assumption,
not a guarantee that genomic control resolves realistic structure. Polygenicity
can inflate statistics too. This simulation does not estimate SNP effects.

Background: [Bulik-Sullivan et al. (2015a), Nature Genetics](https://doi.org/10.1038/ng.3211).

## LD score and an illustrative friend count

The 12-SNP LD example is an exact, valid population-correlation model. Within a
block, let a fair binary ancestral allele be independently flipped at SNP $m$
with probability $(1-u_m)/2$. The off-diagonal correlation is $u_i u_m$; summing
two independent haplotypes preserves the correlation. The two independent
blocks have loadings $(0.2,0.6,0.8,1,0.9,0.7,0.1)$ and $(0.8,1,0.9,0.6)$;
SNP 12 is independent. The diagonal is one.

$$\ell_i=\sum_{m=1}^{M}r_{im}^{\,2}
=1+\sum_{m\ne i}r_{im}^{\,2}.\tag{6}$$

SNP 4 has score $3.35=1+2.35$. Four other SNPs exceed the illustrative
$r^2\geq0.20$ cutoff. The cutoff is **not used in equation (6)**. GCTA's
LD-friends routine instead tests significance of LD; the pedagogical threshold
count is not a reimplementation of that routine. Finite-sample bias-corrected
LD-score estimates need not be at least one.

Source conventions: [GCTA documentation](https://yanglab.westlake.edu.cn/software/gcta/)
and the [LDSC FAQ](https://github.com/bulik/ldsc/wiki/FAQ).

## Published numerical summaries and family diagrams

`rare_variant_inflation.csv` transcribes five reported genomic-inflation values
for a birthplace analysis of 279,390 UK-born participants of inferred European
ancestry. These are alternative covariate adjustments. Near-one global inflation
did not eliminate all spatially confounded loci. Values and sample definition
were checked against the primary paper on 8 September 2026:
[Hanson et al. (2026), Nature Communications](https://doi.org/10.1038/s41467-026-73776-9).

`moba_attenuation.csv` transcribes reported proportional attenuation of SNP
effect estimates, sample sizes and 95% confidence intervals. These are not
changes in SNP heritability, and the intervals are not reconstructed from rounded
standard errors. Values were checked against the primary paper on 8 September
2026: [Corfield et al. (2026), Nature](https://doi.org/10.1038/s41586-026-10926-5).

The pedigree relationship matrix and parental-transmission diagrams are original
schematics. The former represents expected additive relationships, not estimated
genomic relatedness. A covariance model for related observations does not, by
itself, identify direct genetic effects. The lecture's additional mixed-model
slides are explanatory mathematics, not results from an extra simulation here.

## Reproducibility and scale

Main reference versions are R 4.6.0, bigsnpr 1.12.21, bigstatsr 1.6.2 and jsonlite
2.0.0. The main seed feeds randomized SVD as well as subsequent draws, so package
changes can alter particular realizations. Each module uses a fresh R process;
BLAS/OpenMP threading is limited to one. File-backed genotype arrays are created
only under the output folder. No code downloads individual data.

The simulations are modest teaching examples, not biobank benchmarks. A single
selected realization supports the displayed arithmetic, not broad conclusions
about calibration or power. Such conclusions would need prespecified repeated
simulations across sample structure, LD, architectures and model violations.
