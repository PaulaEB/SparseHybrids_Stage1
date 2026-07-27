# Below we will consider two examples,
# The first will have -1 genotype favoured in pool 1 and +1 favoured in pool 2 reflected in the GCA.
# The second will have +1 genotype favoured in both pools in terms of GCA, but the -1 genotype in pool 1
# is actually desired for long term improvement.

# Both examples will have the below in common, but different allele frequencies
a <- 1
d <- 2
np <- 10

######################################
# Example 1
(p1 <- 0.1)
q1 <- 1-p1
p2 <- 0.9
(q2 <- 1-p2)

set.seed(2)
table(X1 <- cbind(sample(c(-1,1), np, prob = c(q1, p1), replace = T)))
table(X2 <- cbind(sample(c(-1,1), np, prob = c(q2, p2), replace = T)))

(X1[rep(1:np, each = np),])
(X2[rep(1:np, times = np),])
X12 <- matrix((X1[rep(1:np, each = np),] + X2[rep(1:np, times = np),])/2, ncol = np)
X12
X1[1]
X2[1]
(W12 <- -1*abs(X12) + 1)

(H <- a*X12 + d*W12)
colnames(H) <- rownames(H) <- 1:np

# centring
Hcent <- H - mean(H)
mean(Hcent)
X12cent <- X12 - mean(X12)
mean(X12cent)
table(H)
table(Hcent)
# plot
plot(X12cent + mean(X12), Hcent + rnorm(np^2, 0, 0.01))
plot(X12, Hcent + rnorm(np^2, 0, 0.01))
# Now we consider the breeding values and dominance deviations in the hybrids (ignoring GCA/SCA for now)
(alpha <- cov(c(X12cent),c(Hcent))/var(c(X12cent)))
# 0.408

# overall allele frequencies in the hybrids
(p <- mean(X12+1)/2) # 0.55
(q <- 1-p) # 0.45
a + (q-p)*d #alpha not considering F, does not match 
# 0.8
(F <- 1 - mean(X12==0) / (2*p*q)) # accounting for non-random mating (i.e., structured mating between pools and not within) Falconer eq. 3.15
# -0.4949495
a + (q-p)*d*(1-F)/(1+F) # check via theory alpha considering F
a + d*cov(c(X12),c(W12))/var(c(X12)) # another way which respects the across-pool mating
# i.e., the fundamental issue causing a distinction with just using a + (q-p)*d is assortative or structured mating between pools, rather than inbreeding itself.
# 0.408
(mu <- mean(Hcent) - alpha*mean(c(X12cent)))
# 0
# insert allele substitution effect, alpha via regression line
plot(X12cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1)); abline(a=mu, b=alpha)
# plot(X12cent + mean(X12), Hcent + rnorm(np^2, 0, 0.01)); abline(a=mu-mean(X12), b=alpha)
(bv <- X12cent*alpha)
points(unique(c(X12cent)), unique(c(bv)), pch = 1, cex = 2) # insert breeding values of the three types of hybrids
(dd <- Hcent - X12cent*alpha) # insert dominance deviations
lines(unique(c(X12cent))[c(1,1)], c(unique(c(bv))[1], unique(c(bv))[1] + unique(c(dd))[1]), col = "darkblue")
lines(unique(c(X12cent))[c(2,2)], c(unique(c(bv))[2], unique(c(bv))[2] + unique(c(dd))[2]), col = "darkred")
lines(unique(c(X12cent))[c(3,3)], c(unique(c(bv))[3], unique(c(bv))[3] + unique(c(dd))[3]), col = "darkorange")
plot(Hcent, bv + dd); abline(a=0,b=1) # check

# GCA's and sca's
# pop 1
q1 # 0.9 # so we want to move towards the -1 genotype in pool 1, and we are starting very close to that
(X1cent <- X1[rep(1:np, each = np),] - mean(X1))
plot(X1cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1))
# plot(X1cent + mean(X1), Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1))
(alpha1 <- cov(c(X1cent),c(Hcent))/var(c(X1cent)))
# -0.3
(p1_obs <- mean(X1 == 1)) # observed allele frequencies
(q1_obs <- mean(X1 == -1))
(p2_obs <- mean(X2 == 1))
(q2_obs <- mean(X2 == -1))
(a + (q2_obs-p2_obs)*d)/2 # check with theory
# -0.3
(F1 <- 1 - mean(W12)/(p1_obs*q2_obs + q1_obs*p2_obs)) # 1 - observed/expected hets = 0 in this case since we make all crosses
(a + (q2_obs-p2_obs)*d*(1-F1)/(1+F1))/2
a/2 + d*cov(c(X1cent), c(W12))/var(c(X1cent)) # another way to account for non-random mating (when it exists)
# -0.3
(mu1 <- mean(Hcent) - alpha1*mean(c(X1cent)))
#-1.822986e-16
plot(X1cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1)); abline(a=mu1, b=alpha1)
# gcas of pop1
(X1cent2 <- X1  - mean(X1))
(gca1 <- X1cent2*alpha1)
cbind(colMeans(Hcent), gca1)# check with hybrid table
points(unique(c(X1cent)), unique(c(round(X1cent*alpha1,3))), pch = 1, cex = 2)
# so here we can see that the -1 genotype in pool 1 is favoured, as we would expect.

# pop 2
p2 # 0.9 # so we want to move towards the +1 genotype in pool 2, and we are starting very close to that
(X2cent <- X2[rep(1:np, times = np),] - mean(X2))
plot(X2cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1))
(alpha2 <- cov(c(X2cent),c(Hcent))/var(c(X2cent)))
# 1.1
(a + (q1_obs-p1_obs)*d)/2 # check with theory
# 1.1
(F2 <- F1) # same as before
(a + (q1_obs-p1_obs)*d*(1-F2)/(1+F2))/2
a/2 + d*cov(c(X2cent), c(W12))/var(c(X2cent)) # another way to account for non-random mating (when it exists)
# 1.1
(mu2 <- mean(Hcent) - alpha2*mean(c(X2cent)))
plot(X2cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1), xlim=c(-2,2)); abline(a=mu2, b=alpha2);abline(a=mu1,b=alpha1, col="red")
# gcas of pop1
(X2cent2 <- X2  - mean(X2))
(gca2 <- X2cent2*alpha2)
cbind(rowMeans(Hcent), gca2)# check with hybrid table
points(unique(c(X2cent)), unique(c(round(X2cent*alpha2,3))), pch = 1, cex = 2)
# points(unique(c(X1cent)), unique(c(round(X1cent*alpha1,3))), pch = 1, cex = 2, col="red")
# so here we can see that the +1 genotype in pool 2 is favoured, as we would expect.

# SCA
# deviation from GCA and SCA
sca <- Hcent - (gca1[rep(1:np, each = np)] + gca2[rep(1:np, times = np)])
plot(Hcent, gca1[rep(1:np, each = np)] + gca2[rep(1:np, times = np)] + sca); abline(a=0,b=1) # check
# nice.

# lets plot this...
plot(X12cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1))
# Firstly, we need to conert our alphas to the scale of the hybrid genotypes
unique(c(X12)) # need to swap 1 and 2
# GCA1
points(unique(c(X12cent))[c(2,1,3)], (c(-1,-1,1) - mean(X1))*alpha1, pch = 1, cex = 1, col = "darkorange") # GCA1
lines(unique(c(X12cent))[c(2,1,3)], (c(-1,-1,1) - mean(X1))*alpha1, pch = 1, col = "darkorange") # GCA1
lines(unique(c(X12cent))[c(2,1,3)], (c(-1,0,1) - mean(X1))*alpha1, lty = 3, lwd = 0.5, col = "darkorange") # GCA1
# GCA2
points(unique(c(X12cent))[c(2,1,3)], (c(-1,1,1) - mean(X2))*alpha2, pch = 1, cex = 1, col = "darkred") # GCA2
lines(unique(c(X12cent))[c(2,1,3)], (c(-1,1,1) - mean(X2))*alpha2, pch = 1, col = "darkred") # GCA2
lines(unique(c(X12cent))[c(2,1,3)], (c(-1,0,1) - mean(X2))*alpha2, lty = 3, lwd = 0.5, col = "darkred") # GCA2
# combined
(gca <- (c(-1,-1,1) - mean(X1))*alpha1 + (c(-1,1,1) - mean(X2))*alpha2)
points(unique(c(X12cent))[c(2,1,3)], gca, pch = 1, cex = 1, col = "darkblue") # combined GCA line (GCA1 + GCA2)
lines(unique(c(X12cent))[c(2,1,3)], gca, pch = 1, col = "darkblue") # combined GCA line (GCA1 + GCA2)
# sca deviations (from combined GCA line)
lines(unique(c(X12cent))[c(2,2)], c(gca[1], gca[1] + unique(c(sca))[2]))
lines(unique(c(X12cent))[c(3,3)], c(gca[3], gca[3] + unique(c(sca))[3]))
lines(unique(c(X12cent))[c(1,1)], c(gca[2], gca[2] + unique(c(sca))[1]))
# The reason we dont have straight lines here is because at the heterozygote, the gca of the parents from
# pools 1 and 2 are equal to their values at -1 and 1, respectively (because those give the hybrid at 0).

# Note that the other way to get the hybrid is if the parents are +1 and -1 (reciprocals of above)
# in this case the plot is given by:
plot(X12cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1))
# GCA1
points(unique(c(X12cent))[c(2,1,3)], (c(-1,1,1) - mean(X1))*alpha1, pch = 1, cex = 1, col = "darkorange") # GCA1
lines(unique(c(X12cent))[c(2,1,3)], (c(-1,1,1) - mean(X1))*alpha1, pch = 1, col = "darkorange") # GCA1
lines(unique(c(X12cent))[c(2,1,3)], (c(-1,0,1) - mean(X1))*alpha1, lty = 3, lwd = 0.5, col = "darkorange") # GCA1
# GCA2
points(unique(c(X12cent))[c(2,1,3)], (c(-1,-1,1) - mean(X2))*alpha2, pch = 1, cex = 1, col = "darkred") # GCA2
lines(unique(c(X12cent))[c(2,1,3)], (c(-1,-1,1) - mean(X2))*alpha2, pch = 1, col = "darkred") # GCA2
 lines(unique(c(X12cent))[c(2,1,3)], (c(-1,0,1) - mean(X2))*alpha2, lty = 3, lwd = 0.5, col = "darkred") # GCA2
# combined
(gca <- (c(-1,1,1) - mean(X1))*alpha1 + (c(-1,-1,1) - mean(X2))*alpha2)
points(unique(c(X12cent))[c(2,1,3)], gca, pch = 1, cex = 1, col = "darkblue") # combined GCA line (GCA1 + GCA2)
lines(unique(c(X12cent))[c(2,1,3)], gca, pch = 1, col = "darkblue") # combined GCA line (GCA1 + GCA2)
lines(unique(c(X12cent))[c(2,2)], c(gca[1], gca[1] + unique(c(sca))[2])) # exactly the same as above
lines(unique(c(X12cent))[c(3,3)], c(gca[3], gca[3] + unique(c(sca))[3])) # exactly the same as above
lines(unique(c(X12cent))[c(1,1)], c(gca[2], gca[2] + unique(c(sca))[4])) # new line for the reciprocal hybrid

# also, we can see the above favours moving pool 1 to -1 genotype and pool 2 to +1 genotype
# below, we will see what happens when GCAs do not favour this, although the locus still has d = 2 > a = 1 like here

# assume pop2 (rows) is testers
Hcent_T <- t(scale(t(H), scale = FALSE))
rowMeans(Hcent_T)
p_T <- as.numeric(X2 == 1)
q_T <- 1- p_T
(alpha_T <- a+(q_T - p_T)*d)
apply(rep(alpha_T, 10)*X12, 1, var)
# 0.1777778 0.1777778 0.1777778 0.1777778 0.1777778 0.1777778 1.6000000 0.1777778 0.1777778 0.1777778
# variance of testcrosses, showing that tester with -1 has highest variance
# lets show this below logically...
# a=1, d= 2, the favourable testers have genotype = 1, given pool 1 has allele frequencies p = 0.2, q = 0.8

X1
# aa  Aa
# 0.8 0.2
gca1 # so 5 & 6 which are +1 genotype have lower gca's than the -1 genotyps, 
# so we want to show gca(g = -1) > gca(g = +1)
# This means that if we have a tester which is AA, it produces hybrids
# Aa     AA
# at frequencies
# 0.8 and 0.2
# d      a
# meaning that a pool 1 parent with aa has hybrids Aa with d =2, 
# while a parent with AA has hybrids AA with a = 1.
# so the gap is 1.

# Alternatively, if we have a tester which is aa, it produces hybrids
# Aa     aa
# at frequencies
# 0.2 and 0.8
# d      -a
# meaning that a pool 1 parent with aa has hybrids aa with -a =-1, 
# while a parent with AA has hybrids Aa with d = 2.
# so the gap is 3.

# in both cases, the rankings of parent GCA's are correct, but a tester
# with -1 will give greater variance as shown by below
apply(rep(alpha_T, 10)*X12, 1, var)*9/10 # note we need to multiply by 9/10 since var uses 1/(n-1) instead of 1/n
0.5*(1+F11)*p1_obs*q1_obs*alpha_T^2

# ok, so the above works for a single tester. If you are interested in the full factorial, i'll leave it to you.
# but the above shows the tester case perfectly.
# I would suggest doing the reciprocal crosses (testers are pool 1), since it will show a different tester state which is favorable...


# additive variance in the hybrids
(F11 <- 1 - mean(X1==0) / (2*p1_obs*q1_obs))
(sigma_tc = 0.5*(1+F11)*p1_obs*q1_obs*((alpha_T)^2))

(mu_tc <- a * (p1_obs * p2_obs - q1_obs * q2_obs) + d * (p1_obs * q2_obs + p2_obs * q1_obs))
(mean(H))

######################################
# Example 2
a # same as above
d # same as above
(p1 <- 0.3) # I changed here now
q1 <- 1-p1
p2 <- 0.7 # I changed here now
(q2 <- 1-p2)

# The code below is exactly the same as above, now using different allele frequencies
set.seed(2)
table(X1 <- cbind(sample(c(-1,1), np, prob = c(q1, p1), replace = T)))
table(X2 <- cbind(sample(c(-1,1), np, prob = c(q2, p2), replace = T)))

(X1[rep(1:np, each = np),])
(X2[rep(1:np, times = np),])
X12 <- matrix((X1[rep(1:np, each = np),] + X2[rep(1:np, times = np),])/2, ncol = np)
X12
X1[1]
X2[1]
(W12 <- -1*abs(X12) + 1)

(H <- a*X12 + d*W12)
colnames(H) <- rownames(H) <- 1:np

# centring
Hcent <- H - mean(H)
mean(Hcent)
X12cent <- X12 - mean(X12)
mean(X12cent)
# plot
plot(X12cent + mean(X12), Hcent + rnorm(np^2, 0, 0.01))
# Now we consider the breeding values and dominance deviations in the hybrids (ignoring GCA/SCA for now)
(alpha <- cov(c(X12cent),c(Hcent))/var(c(X12cent)))
# 0.76 (0.408 before)
# overall allele frequencies in the hybrids
(p <- mean(X12+1)/2) # 0.55 (same as before)
(q <- 1-p) # 0.45
(F <- 1 - mean(X12==0) / (2*p*q)) # accounting for non-random mating (i.e., structured mating between pools and not within)
# -0.09090909 (-0.4949495 before)
a + (q-p)*d*(1-F)/(1+F) # check via theory
a + d*cov(c(X12),c(W12))/var(c(X12)) # another way which respects the across-pool mating
# i.e., the fundamental issue causing a distinction with just using a + (q-p)*d is assortative or structured mating between pools, rather than inbreeding itself.
# 0.76 (0.408 before)
(mu <- mean(Hcent) - alpha*mean(c(X12cent)))
# 0
# insert allele substitution effect, alpha via regression line
plot(X12cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1)); abline(a=mu, b=alpha)
# plot(X12cent + mean(X12), Hcent + rnorm(np^2, 0, 0.01)); abline(a=mu-mean(X12), b=alpha)
(bv <- X12cent*alpha)
points(unique(c(X12cent)), unique(c(bv)), pch = 1, cex = 2) # insert breeding values of the three types of hybrids
(dd <- Hcent - X12cent*alpha) # insert dominance deviations
lines(unique(c(X12cent))[c(1,1)], c(unique(c(bv))[1], unique(c(bv))[1] + unique(c(dd))[1]), col = "darkblue")
lines(unique(c(X12cent))[c(2,2)], c(unique(c(bv))[2], unique(c(bv))[2] + unique(c(dd))[2]), col = "darkred")
lines(unique(c(X12cent))[c(3,3)], c(unique(c(bv))[3], unique(c(bv))[3] + unique(c(dd))[3]), col = "darkorange")
plot(Hcent, bv + dd); abline(a=0,b=1) # check

# GCA's and sca's
# pop 1
q1 # 0.7 # so we still want to move towards the -1 genotype in pool 1, but we are now starting further away.
(X1cent <- X1[rep(1:np, each = np),] - mean(X1))
plot(X1cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1))
(alpha1 <- cov(c(X1cent),c(Hcent))/var(c(X1cent)))
# 0.1 (-0.3 before)
(p1_obs <- mean(X1 == 1)) # observed allele frequencies
(q1_obs <- mean(X1 == -1))
(p2_obs <- mean(X2 == 1))
(q2_obs <- mean(X2 == -1))
(a + (q2_obs-p2_obs)*d)/2 # check with theory
# 0.1 (-0.3 before)
(F1 <- 1 - mean(W12)/(p1_obs*q2_obs + q1_obs*p2_obs)) # 1 - observed/expected hets = 0 in this case since we make all crosses
(a + (q2_obs-p2_obs)*d*(1-F1)/(1+F1))/2 #0.1
a/2 + d*cov(c(X1cent), c(W12))/var(c(X1cent)) # 0.1 another way to account for non-random mating (when it exists)
# 0.1 (-0.3 before)
(mu1 <- mean(Hcent) - alpha1*mean(c(X1cent)))
plot(X1cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1)); abline(a=mu1, b=alpha1)
# gcas of pop1
(X1cent2 <- X1  - mean(X1))
(gca1 <- X1cent2*alpha1)
cbind(colMeans(Hcent), gca1)# check with hybrid table
points(unique(c(X1cent)), unique(c(round(X1cent*alpha1,3))), pch = 1, cex = 2)
# so here we can see that the -1 genotype in pool 1 is no longer favoured via the gca effects
# Despite long-term we do want that -1 genotype in pool 1, the key issue here as that the
# gca effects reflect the average combining ability of the parents, so with a change in allele frequency
# in pool 2 (from 0.9 for +1 to 0.7), the average looks worse for the -1 genotype, and hence the gca effects
# want to push the pool towards +1, favouring a local maximum at a = 1, ignoring the long-term potential to
# get towards d (via favouring the -1 genotype).

# pop 2
q2 # 0.3 # so we still want to move towards the 1 genotype in pool 3, but we are now starting further away.
(X2cent <- X2[rep(1:np, times = np),] - mean(X2))
plot(X2cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1))
(alpha2 <- cov(c(X2cent),c(Hcent))/var(c(X2cent)))
# 0.7 (1.1 before)
(a + (q1_obs-p1_obs)*d)/2 # check with theory
# 0.1 (-0.3 before)
(F2 <- F1) # 1 - observed/expected hets = 0 in this case since we make all crosses
(a + (q1_obs-p1_obs)*d*(1-F2)/(1+F2))/2 #0.7
a/2 + d*cov(c(X2cent), c(W12))/var(c(X2cent)) # 0.7 another way to account for non-random mating (when it exists)
# 0.1 (-0.3 before)
(mu2 <- mean(Hcent) - alpha2*mean(c(X2cent)))
plot(X2cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1)); abline(a=mu2, b=alpha2)

# gcas of pop2
(X2cent2 <- X2  - mean(X2))
(gca2 <- X2cent2*alpha2)
cbind(rowMeans(Hcent), gca2)# check with hybrid table
points(unique(c(X2cent)), unique(c(round(X2cent*alpha2,3))), pch = 1, cex = 2)

# lets take the same a and d values, but change the frequency of the pool 2 alleles, and check the alpha1
# and the gca for parents in pool1
p22 <- seq(0,1,0.01)
q22 <- 1 - p22
(alpha11 <- a + (q22-p22)*d)
plot(p22, alpha11, ylim = c(-3,3)); abline(a=0,b=0)
lines(p22, -1*alpha11, col = "red") # -1 genotype's gca in pool 1
lines(p22, 1*alpha11, col = "blue") # +1 genotype's gca in pool 1
# so what allele frequency do we need to get a positive gca for the -1 genotype?
# assume -1*alpha11 >= 0
# ie.. alpha11 <= 0
# a + (q22-p22)*d <= 0
# (q22-p22) <= -a/d
# so when a=1,d=2
# we want (q22-p22) <= -1/2
# given q = 1-p
# 1-2*p22 <= -1/2
gca1_x <- cbind(alpha11, p22)
gca1_x[which.min(abs(gca1_x[, "alpha11"])), "p22"]
lines(c(0.75,0.75), c(-3,3))
# p22 >= 0.75 (which we had before when p22 = 0.9, but not now as p22 = 0.7)

# pop 2
p2 # 0.7 # so we still want to move towards the +1 genotype in pool 2, but we are now starting further away.
(X2cent <- X2[rep(1:np, times = np),] - mean(X2))
plot(X2cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1))
(alpha2 <- cov(c(X2cent),c(Hcent))/var(c(X2cent)))
# 0.7 (1.1 before)
(a + (q1_obs-p1_obs)*d)/2 # check with theory
# 0.7 (1.1 before)
(F2 <- F1) # same as before
(a + (q1_obs-p1_obs)*d*(1-F2)/(1+F2))/2
a/2 + d*cov(c(X2cent), c(W12))/var(c(X2cent)) # another way to account for non-random mating (when it exists)
# 0.7 (1.1 before)
(mu2 <- mean(Hcent) - alpha2*mean(c(X2cent)))
plot(X2cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1)); abline(a=mu2, b=alpha2)
# gcas of pop1
(X2cent2 <- X2  - mean(X2))
(gca2 <- X2cent2*alpha2)
cbind(rowMeans(Hcent), gca2)# check with hybrid table
points(unique(c(X2cent)), unique(c(round(X2cent*alpha2,3))), pch = 1, cex = 2)
# so here we can see that the +1 genotype in pool 2 is still favoured

# SCA
# deviation from GCA and SCA
sca <- Hcent - (gca1[rep(1:np, each = np)] + gca2[rep(1:np, times = np)])
plot(Hcent, gca1[rep(1:np, each = np)] + gca2[rep(1:np, times = np)] + sca); abline(a=0,b=1) # check
# nice.

# lets plot this...
plot(X12cent, Hcent + rnorm(np^2, 0, 0.01), ylim = c(-3,1))
# Firstly, we need to conert our alphas to the scale of the hybrid genotypes
unique(c(X12)) # need to swap 1 and 2
# GCA1
points(unique(c(X12cent))[c(2,1,3)], (c(-1,-1,1) - mean(X1))*alpha1, pch = 1, cex = 1, col = "darkorange") # GCA1
lines(unique(c(X12cent))[c(2,1,3)], (c(-1,-1,1) - mean(X1))*alpha1, pch = 1, col = "darkorange") # GCA1
# lines(unique(c(X12cent))[c(2,1,3)], (c(-1,0,1) - mean(X1))*alpha1, lty = 3, lwd = 0.5, col = "darkorange") # GCA1
# GCA2
points(unique(c(X12cent))[c(2,1,3)], (c(-1,1,1) - mean(X2))*alpha2, pch = 1, cex = 1, col = "darkred") # GCA2
lines(unique(c(X12cent))[c(2,1,3)], (c(-1,1,1) - mean(X2))*alpha2, pch = 1, col = "darkred") # GCA2
# lines(unique(c(X12cent))[c(2,1,3)], (c(-1,0,1) - mean(X2))*alpha2, lty = 3, lwd = 0.5, col = "darkred") # GCA2
# combined
(gca <- (c(-1,-1,1) - mean(X1))*alpha1 + (c(-1,1,1) - mean(X2))*alpha2)
points(unique(c(X12cent))[c(2,1,3)], gca, pch = 1, cex = 1, col = "darkblue") # combined GCA line (GCA1 + GCA2)
lines(unique(c(X12cent))[c(2,1,3)], gca, pch = 1, col = "darkblue") # combined GCA line (GCA1 + GCA2)
# sca deviations (from combined GCA line)
lines(unique(c(X12cent))[c(2,2)], c(gca[1], gca[1] + unique(c(sca))[2]))
lines(unique(c(X12cent))[c(1,1)], c(gca[2], gca[2] + unique(c(sca))[1]))
lines(unique(c(X12cent))[c(3,3)], c(gca[3], gca[3] + unique(c(sca))[3]))
# The reason we dont have straight lines here is because at the heterozygote, the gca of the parents from
# pools 1 and 2 are equal to their values at -1 and 1, respectively (because those give the hybrid at 0).

# i've left out the reciprocal hybrids for brevity.

# so this all shows that the short-term target via GCA at this locus will actually be a, not d (which is >a).
# provided that p22 < 0.75, given a=1, d=2.

# so what can we do? well we can calculate the expected gain in the hybrids after one round of parent selection

# firstly, assume low selection intensity over two cycles,
# lets say ~80% each time (this will show we will get a gain in hybrids due to less hybrids at -1 genotype but then as pool 2 becomes fixed,
# the selection focus via gca will shift back to the -1 genotype)
# round 1 ------------------------
cbind(X1, gca1)
cbind(X2, gca2) # gca favours +1 genotype in both pools
(X11 <- cbind(X1[rank(gca1, ties.method = "first") > 2, ]))
(X22 <- cbind(X2[rank(gca2, ties.method = "first") > 2, ]))
# create hybrids
(X122 <- matrix((X11[rep(1:8, each = 8)] + X22[rep(1:8, times = 8)])/2, ncol = 8))
(W122 <- -1*abs(X122) + 1)
(H12 <- a*X122 + d*W122)
mean(H12)
# 1.375
# new observed frequencies
(p1_obs2 <- mean(X11 == 1)) # 0.5 (0.3 before)
(q1_obs2 <- mean(X11 == -1)) # so we've actually reduced from 0.7 to 0.5
(p2_obs2 <- mean(X22 == 1)) # now 0.875, so we've gotten past that 0.75 mark.
(q2_obs2 <- mean(X22 == -1))
(gca11 <- colMeans(H12 - mean(H12)))
cbind(X11, gca11) # so now gca has switched to favouring the -1 genotype in pool 1
(gca22 <- rowMeans(H12 - mean(H12)))
cbind(X22, gca22) # same as above, still favouring +1 genotype in pool 2
# round 2 ------------------------
(X111 <- cbind(X11[rank(gca11, ties.method = "first") > 2, ]))
(X222 <- cbind(X22[rank(gca22, ties.method = "first") > 2, ]))
# create hybrids
(X1222 <- matrix((X111[rep(1:6, each = 6)] + X222[rep(1:6, times = 6)])/2, ncol = 6))
(W1222 <- -1*abs(X1222) + 1)
(H122 <- a*X1222 + d*W1222)
mean(H122)
# 1.375 for H12 and 1.666667 for H122...
# new observed frequencies
mean(X111 == 1) # 0.3333333 back towards the favourable direction
mean(X111 == -1)
mean(X222 == 1) # now 1, fixed
mean(X222 == -1)
(gca111 <- colMeans(H122 - mean(H122)))
cbind(X111, gca111) # still favouring -1 genotype, as expected
(gca222 <- rowMeans(H122 - mean(H122))) # all zero, fixed

# you could go and do a third round, to see the hybrid mean increasing again towards 2.

# round 3 ------------------------
(X1111 <- cbind(X111[rank(gca111, ties.method = "first") > 2, ]))
(X2222 <- cbind(X222[rank(gca222, ties.method = "first") > 2, ]))
# create hybrids
(X12222 <- matrix((X1111[rep(1:4, each = 4)] + X2222[rep(1:4, times = 4)])/2, ncol = 4))
(W12222 <- -1*abs(X12222) + 1)
(H1222 <- a*X12222 + d*W12222)
mean(H1222)
# 2 mean...
# new observed frequencies
mean(X1111 == 1) # 0.3333333 back towards the favourable direction
mean(X1111 == -1)
mean(X2222 == 1) # now 1, fixed
mean(X2222 == -1)
(gca1111 <- colMeans(H1222 - mean(H1222)))
cbind(X1111, gca1111) # still favouring -1 genotype, as expected
(gca2222 <- rowMeans(H1222 - mean(H1222))) # all zero, fixed
cbind(X2222, gca2222)

plot(1:4, c(mean(H), mean(H12),  mean(H122), mean(H1222))) # gain over cycles
lines(1:4, c(mean(H), mean(H12),  mean(H122), mean(H1222)))

# now reverse the gca value by *-1 in the first cycle for pool 1
# round 1 ------------------------
cbind(X1, -gca1) # force gca to favour -1 genotype in pool 1
cbind(X2, gca2) # gca favours +1 genotype in pool 2
(X11 <- cbind(X1[rank(-gca1, ties.method = "first") > 2, ]))
(X22 <- cbind(X2[rank(gca2, ties.method = "first") > 2, ]))
# create hybrids
(X122 <- matrix((X11[rep(1:8, each = 8)] + X22[rep(1:8, times = 8)])/2, ncol = 8))
(W122 <- -1*abs(X122) + 1)
(H12 <- a*X122 + d*W122)
mean(H12)
# 1.5 # higher than before
# new observed frequencies
mean(X11 == 1) # 0.25 (0.3 before)
mean(X11 == -1) # so we've now increased from 0.7 to 0.75, which is desirable
mean(X22 == 1) # now 0.875, so we've gotten past that 0.75 mark (just like in the above)
mean(X22 == -1)
(gca11 <- colMeans(H12 - mean(H12)))
cbind(X11, gca11) # same as above, still favoruing -1 genotype in pool 1
(gca22 <- rowMeans(H12 - mean(H12)))
cbind(X22, gca22) # same as above, still favoruing +1 genotype in pool 2
# round 2 ------------------------
(X111 <- cbind(X11[rank(gca11, ties.method = "first") > 2, ]))
(X222 <- cbind(X22[rank(gca22, ties.method = "first") > 2, ]))
# create hybrids
(X1222 <- matrix((X111[rep(1:6, each = 6)] + X222[rep(1:6, times = 6)])/2, ncol = 6))
(W1222 <- -1*abs(X1222) + 1)
(H122 <- a*X1222 + d*W1222)
mean(H122)
# 2 yes, we have created favourable hybrids quickly
# new observed frequencies
mean(X111 == 1)
mean(X111 == -1) # now 1, fixed
mean(X222 == 1) # still 1, fixed
mean(X222 == -1)
(gca111 <- colMeans(H122 - mean(H122)))
cbind(X111, gca111) # all zero, fixed
(gca222 <- rowMeans(H122 - mean(H122))) # all zero, fixed
cbind(X222, gca222)
plot(1:3, c(mean(H), 1.375, 1.375), ylim = c(1,2.5)) # gain over cycles
lines(1:3, c(mean(H), 1.375, 1.375))
points(1:3, c(mean(H), 1.5, 2), col = "green")
lines(1:3, c(mean(H), 1.5, 2), col = "green")

# we can also do the above via manipulation of selection index theory (for later).

# In other words, selecting on GCA alone might find local maxima (a), which does not exploit the overdominance.
# I.e. it attempts to boost short-term gain which sacrifices time getting long-term complementarity and genetic gain.
# In this example, it was actually much better for short and long term improvement by forcing -gca1.


