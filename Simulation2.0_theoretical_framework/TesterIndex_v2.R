# HybridTools

## ------------------------------------------------------------
## Toy "good tester" framework (QTL truth) for 2 pools
## UPDATED:
##  (1) Directional dominance: favorable allele has mostly + additive effects,
##      and dominance deviations are mostly positive (toward the favorable allele).
##  (2) Adds theoretical expected tester alignment curve based on
##      sigma_g^2 and sigma_S^2 and number of testers:
##        Rel(t) = sigma_g^2 / (sigma_g^2 + sigma_S^2 / t)
##        Acc(t) = sqrt(Rel(t))
##      where sigma_g^2 is variance of true GCA among A lines and sigma_S^2 is
##      variance of SCA interaction from the full A×B table decomposition.
##
## - Two pools A,B: 20 inbreds each
## - 100 biallelic sites, with LD along genome
## - 20 QTL among those sites, with additive (a) + dominance (d)
## - True hybrid value H_centre(i,j) (add + dom, no noise)
## - True target for A: gca_true[i] = mean_j∈B H_centre(i,j)
## - Tester set T ⊂ B of size t: gca_testcross[i] = mean_j∈T H_centre(i,j)
## - Alignment: rho(T)=cor(gca_true, gca_testcross)
## - Diagnostics: representativeness, redundancy, discrimination, etc.
## - Optimisation:
##    (a) GA subset selection for each t (alignment curve)
##    (b) Greedy conditional selection (tester2 conditional on tester1, etc.)
## - Visualisation: PCA (+ optional UMAP) with chosen tester(s) highlighted
## ------------------------------------------------------------

## Packages (install if needed)
library(GA)
library(ggplot2)
library(Matrix)
library(FieldSimR)
library(AlphaSimR)

set.seed(1)

## --------------------------
## 1) Simulate LD genotypes
## --------------------------
nA <- 20; nB <- 20
m  <- 100                 # sites
nQTL <- 20
rho_ld <- 0.8            # LD "strength" along genome (AR1-ish)
pA0 <- 0.3               # baseline allele freq in pool A
pB0 <- 0.6               # baseline allele freq in pool B
p_sd <- 0.10              # freq heterogeneity along sites

# AR(1) latent for site-specific allele frequency variation + mild drift between pools
ar1_latent <- function(m, rho) {
  z <- numeric(m); z[1] <- rnorm(1)
  for(k in 2:m) z[k] <- rho*z[k-1] + sqrt(1-rho^2)*rnorm(1)
  z
}
z <- ar1_latent(m, rho_ld)

# site-specific frequencies for each pool (clipped to [0.02,0.98])
pA <- pmin(pmax(pA0 + p_sd*scale(z)[,1] + 2*rnorm(m,0,0.02), 0.02), 0.98)
pB <- pmin(pmax(pB0 + p_sd*scale(z)[,1] + 2*rnorm(m,0,0.02), 0.02), 0.98)
plot(pA, pB); abline(a=0, b=1)

# simulate inbred haplotypes with LD via a simple Markov chain of alleles along sites
# Each inbred genotype is 0/2 coded (fully homozygous at each site)
sim_inbred_pool <- function(n, p_vec, rho_switch=0.92) {
  # rho_switch ~ probability to "stay" with previous allele state; induces LD
  X <- matrix(0, n, length(p_vec))
  for(i in 1:n){
    a <- rbinom(1, 1, p_vec[1])
    X[i,1] <- 2*a
    for(k in 2:length(p_vec)){
      if(runif(1) < rho_switch){
        a <- X[i,k-1]/2
      } else {
        a <- rbinom(1, 1, p_vec[k])
      }
      X[i,k] <- 2*a
    }
  }
  X
}

XA <- sim_inbred_pool(nA, pA, rho_switch=0.93)
XA[1:10, 1:20]
XB <- sim_inbred_pool(nB, pB, rho_switch=0.93)

rownames(XA) <- paste0("A",1:nA)
rownames(XB) <- paste0("B",1:nB)
colnames(XA) <- colnames(XB) <- paste0("S",1:m)

## LD proxy among sites (correlation of 0/2 genotypes across all inbreds)
Xall <- rbind(XA, XB)
LD <- cor(Xall)  # 100x100
plot_matrix(LD, order = F)
plot(LD[,2])

## --------------------------
## 2) Choose QTL + effects
##    Directional dominance
## --------------------------
qtl_idx <- sort(sample(1:m, nQTL))
a <- rep(0, m); d <- rep(0, m)

## Directional additive effects:
##   - define allele coded as "2" (reference allele count=2) as favorable at QTL
##   - set additive effects mostly positive
a_q <- abs(rnorm(nQTL, mean=0.8, sd=0.4))  # mostly positive
a[qtl_idx] <- a_q
mean(a)

## Directional dominance:
##   - dominance deviations mostly positive and correlated with additive magnitude
##   - mild noise, clipped at >= 0
d_q <- 0.6*a_q + rnorm(nQTL, mean=0.2, sd=0.15)
d_q <- pmax(d_q, 0)             # enforce directional (positive) dominance
d_q <- pmin(d_q, 1.8)           # cap
d[qtl_idx] <- d_q
mean(d)

## weights for representativeness: total "importance" at QTL
w_qtl <- (a^2 + d^2)  # 0 on non-QTL 
w_a <- a^2 # additive only importance
w_d <- d^2 # dominance only importance

## --------------------------
## 3) True hybrid value H_centre(i,j)
## --------------------------
# Hybrid genotypic value across loci (add + dom, no noise):
# additive: a_l*(m_ijl - 1), where m_ijl = (x_il + x_jl)/2 ∈ {0,1,2}
# dominance: d_l * I(heterozygous), hetero if x_il != x_jl (parents are homozygous)
hybrid_value <- function(x_i, x_j, a, d) {
  m_ij <- (x_i + x_j)/2               # 0,1,2
  add  <- sum(a * (m_ij - 1))
  het  <- as.numeric(x_i != x_j)      # 0/1
  dom  <- sum(d * het)
  add + dom
}

# Build H matrix for all A x B (nA x nB)
H <- matrix(0, nA, nB, dimnames=list(rownames(XA), rownames(XB)))
for(i in 1:nA){
  for(j in 1:nB){
    H[i,j] <- hybrid_value(XA[i,], XB[j,], a, d)
  }
}
H

## --------------------------
## 4) Variance-component decomposition from full A×B table
##    Model: H_ij = mu + gA_i + gB_j + s_ij
##    (truth, no residual)
## --------------------------
(mu_hat <- mean(H))
H_centre <- H - mu_hat
gcaA_true <- rowMeans(H_centre)
gcaB_true <- colMeans(H_centre)
sca_true  <- H_centre - outer(gcaA_true, rep(1,nB)) - outer(rep(1,nA), gcaB_true)

(sigma_gA2 <- popVar(gcaA_true))          # Var of A GCA effects across A
(sigma_gB2 <- popVar(gcaB_true))          # Var of B GCA effects across B
(sigma_S2  <- popVar(as.vector(sca_true))) # Var of SCA interaction across all A×B

# Expected alignment (accuracy) of estimating A GCA using t common testers (no residual):
# Rel(t) = sigma_gA^2 / (sigma_gA^2 + sigma_S^2 / t)
# Acc(t) = sqrt(Rel(t))
gca_acc_expected <- function(t) sqrt(sigma_gA2 / (sigma_gA2 + sigma_S2 / t))
plot(unlist(lapply(1:20, function(x) gca_acc_expected(x))), type = "l")
var_gca_expected <- function(t) sigma_gA2 + sigma_S2 / t
plot(unlist(lapply(1:20, function(x) var_gca_expected(x))), type = "l", ylim=c(sigma_gA2, 15)); abline(h = sigma_gA2, col="red")

sca_acc_expected <- function(t) sqrt(1 - 1 / t)
plot(unlist(lapply(1:20, function(x) sca_acc_expected(x))), type = "l")
var_sca_expected <- function(t) sigma_S2 * (1 - 1/t)
plot(unlist(lapply(1:20, function(x) var_sca_expected(x))), type = "l", ylim=c(0, 15)); abline(h = sigma_S2, col="red")

# helper: proxy g for a tester set T (indices into B)
gca_testcross <- function(T_idx) rowMeans(H_centre[, T_idx, drop=FALSE])
sca_testcross <- function(T_idx){
  Hsub <- H_centre[, T_idx, drop = FALSE]
  gca_i <- rowMeans(Hsub)
  gca_j <- colMeans(Hsub)
  sca <- sweep(Hsub, 1, gca_i, "-")
  sca <- sweep(sca, 2, gca_j, "-")
  return(sca)
}

# alignment rho(T): Corr over A candidates between true and proxy
gca_align <- function(T_idx) cor(gcaA_true, gca_testcross(T_idx))
summary(gca_alignment <- unlist(lapply(1:20, function(x) gca_align(x))))
plot(gca_alignment, type = "l"); abline(h=gca_acc_expected(1))

gca_var <- function(T_idx) popVar(gca_testcross(T_idx))
summary(gca_variance <- unlist(lapply(1:20, function(x) gca_var(x))))
plot(gca_variance, type = "l"); abline(h = sigma_gA2, col="blue")

plot(gca_variance, gca_alignment); abline(v = sigma_gA2, col="blue"); abline(h=gca_acc_expected(1), col="red");text(gca_variance, gca_alignment, labels=rownames(XB), pos=3, cex=0.5)

# need to consider pairs for testcross SCA alignment and variance
sca_var <- function(T_idx) mean(diag(popVar(sca_testcross(T_idx))))

summary(sca_variance2 <- apply(combn(1:20,2), 2, function(x) sca_var(x)))
plot(sca_variance2, type = "l"); abline(h = sigma_S2, col="blue")
summary(gca_alignment2 <- apply(combn(1:20,2), 2, function(x) gca_align(x)))
plot(sca_variance2, gca_alignment2); abline(v = sigma_S2, col="blue"); abline(h=gca_acc_expected(2), col="red")
summary(gca_variance2 <- apply(combn(1:20,2), 2, function(x) gca_var(x)))
plot(sca_variance2, gca_variance2); abline(v = sigma_S2,col="blue"); abline(h = sigma_gA2,col="blue")
plot(gca_variance2, gca_alignment2); abline(v = sigma_gA2,col="blue"); abline(h=gca_acc_expected(2),col="red")

# trios
summary(sca_variance3 <- apply(combn(1:20,3), 2, function(x) sca_var(x)))
plot(sca_variance3, type = "l"); abline(h = sigma_S2, col="blue")
summary(gca_alignment3 <- apply(combn(1:20,3), 2, function(x) gca_align(x)))
plot(sca_variance3, gca_alignment3); abline(v = sigma_S2,col="blue"); abline(h=gca_acc_expected(3), col="red")
summary(gca_variance3 <- apply(combn(1:20,3), 2, function(x) gca_var(x)))
plot(sca_variance3, gca_variance3); abline(v = sigma_S2,col="blue"); abline(h = sigma_gA2,col="blue")
plot(gca_variance3, gca_alignment3); abline(v = sigma_gA2,col="blue"); abline(h=gca_acc_expected(3),col="red")

# DT up to here

## --------------------------
## 5) Diagnostics for a tester set T
## --------------------------
# Representativeness: weighted squared distance between allele freqs in T and pool B
# (lower is better)
rep_score <- function(weights, T_idx) {
  pT <- colMeans(XB[T_idx,,drop=FALSE] / 2)
  pB_hat <- colMeans(XB / 2)
  -sum(weights * (pT - pB_hat)^2)
}

# Discrimination among A using selected testers: Var of proxy g (bigger = more separation)
disc_power <- function(T_idx) var(gca_testcross(T_idx))

# Redundancy: average pairwise relatedness among chosen testers 
# (lower is better)
# Here: mean pairwise correlation across sites
redundancy <- function(T_idx){
  if(length(T_idx) < 2) return(0)
  G <- cor(t(XB[T_idx,,drop=FALSE]))
  # weighted correlation?
  mean(G[upper.tri(G)])
}

# Combine into an index (tune lambdas as you like)
tester_index <- function(T_idx, lam=c(alignment=1, rep=0.1, disc=0.1, red=0.1)){
    lam["alignment"]*gca_align(T_idx) + # not available in practice
    lam["rep"]*rep_score(weights = 1, T_idx) +
    lam["disc"]*log1p(disc_power(T_idx)) -
    lam["red"]*redundancy(T_idx)
}

## --------------------------
## 6) Optimisation methods for choosing testers from B
## --------------------------

## (A) GA subset selection for fixed t: maximise alignment or index
ga_select_testers <- function(t, objective=c("alignment","index"),
                              lam=c(alignment=1, rep=0.1, disc=0.1, red=0.1),
                              popSize=60, maxiter=80, run=30, seed=NULL){
  objective <- match.arg(objective)
  if(!is.null(seed)) set.seed(seed)
  n <- nB
  
  fitness <- function(bitstring){
    if(sum(bitstring) != t) return(-1e9)
    T_idx <- which(bitstring == 1)
    if(objective == "alignment") return(gca_align(T_idx))
    tester_index(T_idx, lam=lam)
  }
  
  GA::ga(type="binary", nBits=n, fitness=fitness,
         popSize=popSize, maxiter=maxiter, run=run, pmutation=0.15,
         elitism=max(1, round(0.05*popSize)), keepBest=TRUE, monitor=FALSE)
}

## (B) Greedy conditional selection: pick tester1, then tester2 given tester1, etc.
greedy_select_testers <- function(t, objective=c("alignment","index"),
                                  lam=c(alignment=1, rep=0.1, disc=0.1, red=0.1)){
  objective <- match.arg(objective)
  chosen <- integer(0)
  remaining <- 1:nB
  for(k in 1:t){
    scores <- sapply(remaining, function(j){
      T <- c(chosen, j)
      if(objective == "alignment") gca_align(T) else tester_index(T, lam)
    })
    best_j <- remaining[which.max(scores)]
    chosen <- c(chosen, best_j)
    remaining <- setdiff(remaining, best_j)
  }
  chosen
}

## --------------------------
## 7) Alignment curve vs #testers
##    - GA empirical best rho(T)
##    - Greedy conditional best rho(T)
##    - Theoretical expected accuracy sqrt(Rel(t)) from sigma_gA2, sigma_S2
## --------------------------
t_max <- 10
curve_ga <- data.frame(t=1:t_max, rho=NA, method="GA")
curve_gr <- data.frame(t=1:t_max, rho=NA, method="Greedy")
curve_th <- data.frame(t=1:t_max, rho=sapply(1:t_max, gca_acc_expected), method="Theory")

best_sets_ga <- vector("list", t_max)
best_sets_gr <- vector("list", t_max)

for(t in 1:t_max){
  ga_fit <- ga_select_testers(t, objective="alignment", popSize=70, maxiter=1000, run=40, seed=100+t)
  best_bits <- GA::summary(ga_fit)$solution[1,]
  T_ga <- which(best_bits == 1)
  best_sets_ga[[t]] <- T_ga
  curve_ga$rho[t] <- gca_align(T_ga)
  
  T_gr <- greedy_select_testers(t, objective="alignment")
  best_sets_gr[[t]] <- T_gr
  curve_gr$rho[t] <- gca_align(T_gr)
}

curve <- rbind(curve_ga, curve_gr, curve_th)

p_curve <- ggplot(curve, aes(t, rho, color=method, linetype=method)) +
  geom_line(linewidth=1) + geom_point(size=2) + scale_x_continuous(breaks = 1:t_max) +
  ylim(0,1) +
  labs(
    title="Alignment curve: Corr(true GCA vs tester-set proxy) + theoretical expectation",
    subtitle=paste0("Directional dominance; sigma_gA^2=", round(sigma_gA2,3),
                    ", sigma_S^2=", round(sigma_S2,3)),
    x="# testers in set", y="Alignment (rho)"
  ) +
  theme_minimal()
print(p_curve)

## Choose a t to visualise (e.g. t=4)
t_show <- 5
T_star <- best_sets_ga[[t_show]]  # GA best set of size t_show
cat("Chosen testers (GA, t=", t_show, "): ", paste(rownames(XB)[T_star], collapse=", "), "\n", sep="")

## --------------------------
## 8) PCA / (optional) UMAP plots for genotype space
## --------------------------
pca <- prcomp(Xall, center=TRUE, scale.=TRUE)
scores <- as.data.frame(pca$x[,1:2])
scores$id <- rownames(scores)
scores$pool <- c(rep("A", nA), rep("B", nB))
scores$role <- "Other"
scores$role[scores$pool=="B" & scores$id %in% rownames(XB)[T_star]] <- "ChosenTester"

p_pca <- ggplot(scores, aes(PC1, PC2, shape=pool, color=role)) + 
  geom_point(size=3, alpha=0.9) +
  geom_text(data = subset(scores, role == "ChosenTester"),
            aes(label = id),
            vjust = -0.8, size = 3)+ 
  labs(title=paste0("PCA of genotypes (chosen testers highlighted), t=", t_show),
       x="PC1", y="PC2") +
  theme_minimal()
print(p_pca)

# Optional UMAP
  emb <- uwot::umap(scale(Xall), n_neighbors=10, min_dist=0.2, metric="euclidean") # , random_state=1
  um <- data.frame(UMAP1=emb[,1], UMAP2=emb[,2], id=rownames(Xall),
                   pool=c(rep("A", nA), rep("B", nB)))
  um$role <- "Other"
  um$role[um$pool=="B" & um$id %in% rownames(XB)[T_star]] <- "ChosenTester"
  
  p_umap <- ggplot(um, aes(UMAP1, UMAP2, shape=pool, color=role)) +
    geom_point(size=3, alpha=0.9) +
    labs(title=paste0("UMAP of genotypes (chosen testers highlighted), t=", t_show)) +
    theme_minimal()
print(p_umap)

#####################################################################
#####################################################################
# up to here 4 March...


## --------------------------
## 9) Additional diagnostics for the chosen set + variance components
## --------------------------
cat("\n--- Variance components from full A×B true table ---\n")
cat("sigma_gA^2 (Var GCA_A): ", round(sigma_gA2, 4), "\n")
cat("sigma_gB^2 (Var GCA_B): ", round(sigma_gB2, 4), "\n")
cat("sigma_S^2  (Var SCA):   ", round(sigma_S2, 4), "\n")

cat("\n--- Chosen set diagnostics (GA) ---\n")
cat("Alignment rho(T):           ", round(gca_align(T_star), 3), "\n")
cat("Expected rho from theory:   ", round(gca_acc_expected(length(T_star)), 3), "\n")
cat("Representativeness R(T):    ", round(rep_score(T_star), 3), "\n")
cat("Redundancy(T):              ", round(redundancy(T_star), 3), "\n")
cat("Discrimination Var(gca_testcross):", round(disc_power(T_star), 3), "\n")

## --------------------------
## 10) Notes / next refinements
## --------------------------
# 1) This is true: QTL (a,d) known and used for weighting.
#    Later: replace (a,d) with marker-effect BLUPs or Bayesian posteriors.
# 2) The theoretical curve sqrt(sigma_gA^2/(sigma_gA^2 + sigma_S^2/t)) is a
#    benchmark under the classic model with common testers and SCA as “noise”.
#    GA/Greedy can exceed it if the tester subset is not “random” but chosen to
#    reduce effective SCA noise for the specific A panel (true advantage).
# 3) LD-aware representativeness: replace sum w*(pT-pB)^2 by Mahalanobis using LD
#    among QTL sites:
#      (pT-pB)' * (LD_QTL + eps*I)^(-1) * (pT-pB)
# ------------------------------------------------------------


## ------------------------------------------------------------
## Some add-ons
## ------------------------------------------------------------

## 1) Replace (a,d) true weights with marker-effect BLUP weights (minimal placeholder)
##    Here we simulate "estimated" additive & dominance effects by adding noise to truth.
##    Later you would replace a_hat/d_hat with effects from your fitted model.
set.seed(123)
a_hat <- a + rnorm(m, 0, 0.4) * (a != 0)   # noisy estimates on QTL only (minimal)
d_hat <- d + rnorm(m, 0, 0.3) * (d != 0)
w_hat <- (a_hat^2 + d_hat^2)

rep_score_hat <- function(T_idx) {
  pT <- colMeans(XB[T_idx,,drop=FALSE] / 2)
  pB_hatfreq <- colMeans(XB / 2)
  -sum(w_hat * (pT - pB_hatfreq)^2)
}

## If you want the index to use these weights:
tester_index_hatW <- function(T_idx, lam=c(alignment=1, rep=0.1, disc=0.1, red=0.1)){
  lam["alignment"]*gca_align(T_idx) +
    lam["rep"]*rep_score_hat(T_idx) +
    lam["disc"]*log1p(disc_power(T_idx)) -
    lam["red"]*redundancy(T_idx)
}

## 2) “Chosen testers can exceed the theory” check: compare GA rho vs theory, and
##    compute RANDOM tester baseline to show theory aligns with random sampling.
random_alignment <- function(t, nrep=2000){
  rhos <- replicate(nrep, {
    T <- sample.int(nB, t)
    gca_align(T)
  })
  c(mean=mean(rhos), sd=sd(rhos), q05=quantile(rhos,0.05), q50=quantile(rhos,0.5), q95=quantile(rhos,0.95))
}

rand_tbl <- do.call(rbind, lapply(1:t_max, function(t){
  s <- random_alignment(t, nrep=1000)
  data.frame(t=t, rho_mean=s["mean"], rho_sd=s["sd"], rho_q05=s["q05"], rho_q50=s["q50"], rho_q95=s["q95"],
             rho_theory=gca_acc_expected(t))
}))
print(rand_tbl)

## 3) LD-aware representativeness via Mahalanobis distance on QTL sites
##    Uses LD (correlation) matrix computed earlier as a proxy for marker covariance.
qtl <- qtl_idx
LD_qtl <- LD[qtl, qtl, drop=FALSE]

# regularize to ensure invertible
eps <- 1e-4
LD_qtl_reg <- LD_qtl + diag(eps, nrow(LD_qtl))

inv_LD_qtl <- solve(LD_qtl_reg)

rep_score_mahal <- function(T_idx, weights=c("truth","hat","none")){
  weights <- match.arg(weights)
  pT <- colMeans(XB[T_idx, qtl, drop=FALSE] / 2)
  pB_hatfreq <- colMeans(XB[, qtl, drop=FALSE] / 2)
  diff <- pT - pB_hatfreq
  
  md2 <- as.numeric(t(diff) %*% inv_LD_qtl %*% diff) # Mahalanobis distance^2
  
  if(weights == "truth"){
    return(-md2)
  } else if(weights == "truth"){
    w <- w_qtl[qtl]
  } else {
    w <- w_hat[qtl]
  }
  # weight by QTL importance (simple scaling)
  return( -sum(w) * md2 )
}

## Example: swap representativeness term in the index to the LD-aware one
tester_index_LDrep <- function(T_idx, lam=c(alignment=1, rep=0.1, disc=0.1, red=0.1), weights="truth"){
  lam["alignment"]*gca_align(T_idx) +
    lam["rep"]*rep_score_mahal(T_idx, weights=weights) +
    lam["disc"]*log1p(disc_power(T_idx)) -
    lam["red"]*redundancy(T_idx)
}

## 4) Conditional GA (tester 2+ chosen conditional on fixed tester 1, etc.)
##    Minimal: fix chosen indices and GA over remaining using an offset bitstring.
ga_select_conditional <- function(t_total, fixed_idx,
                                  objective=c("alignment","index","index_LD"),
                                  lam=c(alignment=1, rep=0.1, disc=0.1, red=0.1),
                                  popSize=60, maxiter=80, run=30, seed=NULL){
  objective <- match.arg(objective)
  if(!is.null(seed)) set.seed(seed)
  fixed_idx <- sort(unique(fixed_idx))
  stopifnot(all(fixed_idx %in% 1:nB))
  stopifnot(length(fixed_idx) <= t_total)
  
  remaining <- setdiff(1:nB, fixed_idx)
  nBits <- length(remaining)
  t_need <- t_total - length(fixed_idx)
  
  fitness <- function(bitstring){
    if(sum(bitstring) != t_need) return(-1e9)
    add_idx <- remaining[which(bitstring == 1)]
    T_idx <- c(fixed_idx, add_idx)
    
    if(objective == "alignment") return(gca_align(T_idx))
    if(objective == "index")     return(tester_index(T_idx, lam=lam))
    tester_index_LDrep(T_idx, lam=lam, weights="truth")
  }
  
  GA::ga(type="binary", nBits=nBits, fitness=fitness,
         popSize=popSize, maxiter=maxiter, run=run, pmutation=0.15,
         elitism=max(1, round(0.05*popSize)), keepBest=TRUE, monitor=FALSE)
}

## Example conditional run:
##  - pick best single tester by brute force
best1 <- which.max(sapply(1:nB, function(j) gca_align(j)))
cat("Best single tester:", rownames(XB)[best1], "rho=", round(gca_align(best1),3), "\n")

##  - then choose 4-testers total conditional on that best1
ga_cond <- ga_select_conditional(t_total=4, fixed_idx=best1, objective="alignment",
                                 popSize=70, maxiter=100, run=40, seed=777)
best_bits <- GA::summary(ga_cond)$solution[1,]
rem <- setdiff(1:nB, best1)
T_cond <- c(best1, rem[which(best_bits == 1)])
cat("Conditional GA testers (t=4):", paste(rownames(XB)[T_cond], collapse=", "),
    "rho=", round(gca_align(T_cond),3), "\n")