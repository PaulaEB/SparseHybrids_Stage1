# some code for Sven
n <- 100 # genotypes
m <- 1000 # markers
set.seed(123)
a <- scale(rnorm(m))
hist(a)
d <- 1 + scale(rnorm(m))*0.5 # directional
hist(d)
mean(d)
table(cut(d, breaks = c(-Inf, 0, 1, Inf), 
          labels = c("under", "partial/complete", "over")))

# pool 1
(X1 <- matrix(sample(c(-1,1), n*m, prob = c(0.7,0.3), replace = T), ncol = m))
hist(colMeans(X1+1)/2)
(mu1 <- colMeans(X1)) # means
plot(sigma1 <- apply(X1, 2, sd), mu1) # sds

# pool 2
(X2 <- matrix(sample(c(-1,1), n*m, prob = c(0.3,0.7), replace = T), ncol = m)) # pool 2
hist(colMeans(X2+1)/2)
(mu2 <- colMeans(X2)) # means
plot(sigma2 <- apply(X2, 2, sd), mu2) # sds

# compare parameters
plot(colMeans(X1+1)/2, colMeans(X2+1)/2, col = c("blue","red"))
cor(colMeans(X1+1)/2, colMeans(X2+1)/2) # -0.02564634
hist(colMeans(X1+1)/2, xlim = c(0, 1), col="blue"); hist(colMeans(X2+1)/2, add = T, col = "red")
plot(mu1, mu2, col = c("blue","red"))
plot(sigma1, sigma2, col = c("blue","red")); abline(a=0, b=1)

# create hybrids and calculate GCA and SCA effects
# 10,000 hybrids...
X12 <- (X1[rep(1:n, each = n), ] + X2[rep(1:n, times = n), ])/2 # P1 cycles slowly, P2 cycles quickly
W12 <- -1*abs(X12) + 1
hybrid_mat <- matrix(X12 %*% a + W12 %*% d, ncol = n)
colnames(hybrid_mat) <- 1:n
rownames(hybrid_mat) <- (n + 1):(2*n)
hybrid_mat[1:10, 1:10]
hist(gca1 <- colMeans(hybrid_mat) - mean(hybrid_mat))
hist(gca2 <- rowMeans(hybrid_mat) - mean(hybrid_mat))
hist(sca12 <- hybrid_mat  - mean(hybrid_mat) - 
                           matrix(rep(gca1, each = n), ncol = n) - 
                           matrix(rep(gca2, times = n), ncol = n))
mean(sca12) # 3.979126e-14

# GCA variance
var(gca1) # 339.871
var(gca2) # 235.8704

# SCA variance
diag(var(sca12)) # heterogeneity in the variances across P2 WHEN CHOOSING TESTERS p1
var(sca12[,1])
var(sca12[1,]) # heterogeneity in the variances across P1 WHEN CHOOSING TESTERS p2

hist(diag(var(sca12))) # heterogeneity in the variances across P2 WHEN CHOOSING TESTERS p1
hist(diag(var(t(sca12)))) # heterogeneity in the variances across P1 WHEN CHOOSING TESTERS p2
mean(diag(var(sca12))); mean(diag(var(t(sca12)))) # 221.2509

# regress onto gca effects to identify good parents 
hist(a<-cor(sca12, gca2))
rev(order(a))[1:10] # 69 80 31 76 39 83 63 45 36 89
# plot(gca2, sca12[,which(b1 == max(b1))]); abline(a=0, b=max(b1)) # good tester i98
hist(b1 <- cov(sca12, gca2)/var(gca2)) # fair bit of heterogeneity
cov(sca12[,1], gca2)/var(gca2)# -0.028426189
cov(sca12[,69], gca2)/var(gca2) #0.2623899
cov(sca12[,69], gca2)
# plot(gca2, sca12[,which(b1 == min(b1))]); abline(a=0, b=min(b1)) # bad tester i98

cov(sca12[,98], gca2)/var(gca2) #-0.2243718
cov(sca12[,98], gca2) #-52.92266
rev(order(gca1)) # best testers GCA from P1, are not related to best b1
plot(gca1,b1)
hist(hybrid_mat[,order(gca1)[1:10]])
hist(hybrid_mat[,order(b1)[1:10]])
hist(hybrid_mat[,rev(order(gca1))[1:10]])
rev(order(b1))[1:10] #69 80 31 76 39 83 63 45 36 89
hist(b1)
hist(hybrid_mat);abline(v=mean(hybrid_mat), col="blue",lwd=3) #complete
hist(hybrid_mat[,rev(order(b1))[1:10]]) #highest b1
hist(hybrid_mat[,order(b1)[1:10]]) #lowest b1

hist(cor(t(sca12), gca1))
hist(b2 <- cov(t(sca12), gca1)/var(gca1)) # fair bit of heterogeneity
plot(gca1, sca12[which(b2 == max(b2)),]); abline(a=0, b=max(b2)) # so i45 from P2 produces high hybrids when GCA from P1 is high
plot(gca1, sca12[which(b2 == min(b2)),]); abline(a=0, b=min(b2)) # bad tester i61
hist(cor(t(sca12), gca2))
hist(hybrid_mat[,order(gca2)[1:10]])
hist(hybrid_mat[,order(b2)[1:10]])
hist(hybrid_mat[,rev(order(gca2))[1:10]])
hist(hybrid_mat[,rev(order(b2))[1:10]])

#Careful as the ids of b1 and b2 are different but which extracts positions
plot(b2, b1) # top or right
which(b1 > 0.15) 
# 31 39 63 69 76 80 83
rev(order(b1))[1:10]
# 69 80 31 76 39 83 63 45 36 89
which(b2 > 0.15)
# 42 45 48 59 69
rev(order(b2))[1:10]
# 45 42 48 59 69 81  8 78 90 24
hist(hybrid_mat) # looking at > 650
hist(hybrid_mat,col = adjustcolor("blue",0.5))
hist(hybrid_mat[b2 > 0.15, b1 > 0.15],col = adjustcolor("green",0.5),xlim = c(min(hybrid_mat),max(hybrid_mat)))
abline(v=mean(hybrid_mat), col="blue",lwd=3); abline(v=mean(hybrid_mat[b2 > 0.15, b1 > 0.15]), col="darkgreen",lwd=3);abline(v=mean(hybrid_mat[b2 < -0.2, b1 < -0.1]), col="darkred",lwd=3)# nice one. gets the highest
hist(hybrid_mat[b2 < -0.2, b1 < -0.1],add=T,col=adjustcolor("red",0.5)) # nice one. gets the lowest

hist(hybrid_mat) # looking at > 650
hist(hybrid_mat,col = adjustcolor("blue",0.5))
hist(hybrid_mat[ b1 > 0.15],col = adjustcolor("green",0.5),xlim = c(min(hybrid_mat),max(hybrid_mat)), ylim=c(0, 250))
abline(v=mean(hybrid_mat), col="blue",lwd=3); abline(v=mean(hybrid_mat[b1 > 0.15]), col="darkgreen",lwd=3);abline(v=mean(hybrid_mat[b2 < -0.2, b1 < -0.1]), col="darkred",lwd=3)# nice one. gets the highest
hist(hybrid_mat[, b1 < -0.1],add=T,col=adjustcolor("red",0.5),xlim=c(0,280)) # nice one. gets the lowest

hist(gca1)
which(hybrid_mat == max(hybrid_mat), arr.ind = TRUE)
which(colSums(hybrid_mat == max(hybrid_mat)) == 1) # 97
gca1[97] # 35.49576 
b1[97]  # -0.04517902
hist(gca2)
which(rowSums(hybrid_mat == max(hybrid_mat)) == 1) # 33
gca2[33] # 38.67975
b2[33] #-0.02644912
sca12[33,97] #37.22081
mean(hybrid_mat) + gca1[97] + gca2[33] + sca12[33,97] #693.3316
max(hybrid_mat) #693.3316

which(gca1 == max(gca1), arr.ind = TRUE) #48
which(gca2 == max(gca2), arr.ind = TRUE) #33
which(sca12 == max(sca12), arr.ind = TRUE) #row 25 col 59
which(b1 == max(b1), arr.ind = TRUE) # 69
which(b2 == max(b2), arr.ind = TRUE) # 45

# particularly good combination, despite low b1&b2
# think about the above in terms of tester choice for hybrids (i.e., high b1/b2, or parent selection)

# Take single testers
# GCA alignment
# consider testers from P1
plot(gca1_var <- diag(var(gca2 + sca12)), gca1_cor <- cor(gca2 + sca12, gca2)) # top right, high alignment and discriminating ability
points(gca1_var[69], gca1_cor[69], col = "red")
points(gca1_var[89], gca1_cor[89], col = "green")
which(gca1_cor > 0.8) # 69 89    69 was as above!! In max b1 but 89 was not in b1 >0.15
rev(order(gca1_cor))[1:10] #89 69 31 26 17 82 80 39 14 49

# gca1_var_hmat <- diag(var(hybrid_mat))
# gca1_cor_hmat<-cor(hybrid_mat, gca2)
# all.equal(as.numeric(gca1_var), as.numeric(gca1_var_hmat)) #TRUE
# all.equal(as.numeric(gca1_cor), as.numeric(gca1_cor_hmat)) #TRUE

plot(b1, sqrt(gca1_cor^2*gca1_var))
points(b1[69], sqrt(gca1_cor^2*gca1_var)[69], col = "red")
points(b1[89], sqrt(gca1_cor^2*gca1_var)[89], col = "green")

plot((1-gca1_cor^2)*gca1_var, gca1_cor^2*gca1_var) # top, and left
points(((1-gca1_cor^2)*gca1_var)[69], (gca1_cor^2*gca1_var)[69], col = "red")
points(((1-gca1_cor^2)*gca1_var)[89], (gca1_cor^2*gca1_var)[89], col = "green")

# consider testers from P2
plot(gca2_var <- diag(var(t(hybrid_mat))), gca2_cor <- cor(t(hybrid_mat), gca1))
which(gca2_cor > 0.84) # 42 45 69 89   45,42,69 was as above!! In max b2 but 89 was not in b2 >0.15
plot((1-gca2_cor^2)*gca2_var, gca2_cor^2*gca2_var) # top, and left
plot(b2, sqrt(gca2_cor^2*gca2_var))

# good parents vs good testers
which(b1 == max(b1)) # 69
which(gca1 == max(gca1)) # 48
plot(b1, gca1) # top right would be a good tester and good parent
points(b1[69], gca1[69], col = "red") # you can see our good tester, i69, is not a good parent :(
points(b1[89], gca1[89], col = "green") # you can see our good tester, i48, is not a good parent :(
points(b1[48], gca1[48], col = "blue") # you can see our good tester, i48, is not a good parent :(

######################################################################################################################################
# Lets compare the true alignment to that expected from a heritability type expression
# firstly, assuming the classical form
hist(gca1_cor_exp1a <- sqrt(var(gca2)/(var(gca2) + diag(var(sca12))))) 
# so a range here, driven by the variation in the sca variances
# compare these to the true alignments...
plot(gca1_cor, gca1_cor_exp1a) # correlated, but not perfect...
# so we can do better...
(b1 <- cov(sca12, gca2)/var(gca2)) # recall
plot(gca1_cor, gca1_cor_exp1b <- sqrt((1+b1)^2*var(gca2)/((1+b1)^2*var(gca2) + (diag(var(sca12)) - b1^2*var(gca2))))); abline(a=0, b=1)
plot(gca1_cor, sqrt((1+b1)^2*var(gca2)/((1+2*b1)*var(gca2) + diag(var(sca12))))); abline(a=0, b=1) # just another way of writing it
# perfect!
# below is just a toy example demonstrating the above for any two random variables
xx <- rnorm(10000)
yy <- rnorm(10000) + xx*1
(r <- cor(xx, yy))
zz <- xx + yy
cor(xx, zz)
(b <- cov(xx, yy) / var(xx))   # signal in yy aligned with xx
popVar(e <- residuals(lm(yy ~ xx)))  # noise part of yy
popVar(yy) - (b)^2 *popVar(xx)
sqrt((1 + b)^2 * popVar(xx) / ((1 + 2*b) * popVar(xx) + popVar(e)))
# # but what about when averaging over testers/environments? TO-DO

################################################################################
#GCA regression
# this is the alpha_T that depends on the opposite pool allele frequencies
# Tester 14 of pool 1 across pool 2: here 
plot(X1[,2], hybrid_mat[,14], col="red");abline(lm(hybrid_mat[,14]~X1[,2]), col="red");
points(X2[,2], hybrid_mat[,14], col="blue",pch=4);abline(lm(hybrid_mat[,14] ~ X2[,2]), col="blue")
################################################################################
##graphs report 9 months

##BEST B1 AND WORST B1
par(mar = c(3.3, 3.3, 1.2, 0.8),
    mgp = c(1.7, 0.5, 0),
    tck = -0.015,
    bty = "l",
    cex.axis = 0.75,
    cex.lab  = 0.95)

plot(gca2, sca12[, which.min(b1)],
     pch = 1, col = "darkred", cex = 0.6,
     xlab = expression("Candidate lines GCA (" * g[1] * ")"),
     ylab = expression("Tester-specific SCA (" * s[1~italic(j)] * ")"))

abline(a = 0, b = min(b1))

legend("topright",
       legend = bquote(italic(j) : T[.(which.min(b1))]),
       bty = "n", text.col = "grey20", cex = 0.8)

par(mar = c(3.3, 3.3, 1.2, 0.8),
    mgp = c(1.7, 0.5, 0),
    tck = -0.015,
    bty = "l",
    cex.axis = 0.75,
    cex.lab  = 0.95)

plot(gca2, sca12[, which.max(b1)],
     pch = 1, col = "royalblue", cex = 0.6,
     xlab = expression("Candidate lines GCA (" * g[1] * ")"),
     ylab = expression("Tester-specific SCA (" * s[1~italic(j)] * ")"))

abline(a = 0, b = max(b1))

legend("topright",
       legend = bquote(italic(j) : T[.(which.max(b1))]),
       bty = "n", text.col = "grey20", cex = 0.8)

## HYBRID VALUES
par(cex.axis = 0.75,
    cex.lab  = 0.95)
hist(as.numeric(hybrid_mat),
     col = adjustcolor("royalblue", 0.4), border = NA,
     xlab = "Hybrid genetic value", ylab = "Frequency",ylim=c(0, 2800))

abline(v = mean(as.numeric(hybrid_mat)),col = "blue", lwd = 1)
abline(v = mean(as.numeric(hybrid_mat[, rev(order(b1))[1:10]])), col = "darkgreen", lwd = 1)
abline(v = mean(as.numeric(hybrid_mat[, order(b1)[1:10]])), col = "darkred",   lwd = 1)

legend(par("usr")[2], par("usr")[4], xjust = 1, yjust = 1,  legend = expression("Full factorial", "Top 10 " * beta[j], "Bottom 10 " * beta[j]), 
       col = c("blue", "darkgreen", "darkred"), lwd = 1, bty = "n", cex = 0.6,seg.len = 1, x.intersp = 0.6)

par(cex.axis = 0.75,
    cex.lab  = 0.95)
hist(as.numeric(hybrid_mat),
     col = adjustcolor("royalblue", 0.4), border = NA,
     xlab = "Hybrid genetic value", ylab = "Frequency",ylim=c(0, 2800))

abline(v = mean(as.numeric(hybrid_mat)),col = "blue", lwd = 1)
abline(v = mean(as.numeric(hybrid_mat[, rev(order(gca1))[1:10]])), col = "darkgreen", lwd = 1)
abline(v = mean(as.numeric(hybrid_mat[, order(gca1)[1:10]])), col = "darkred",   lwd = 1)

legend(par("usr")[2], par("usr")[4], xjust = 1, yjust = 1,  legend = expression("Full factorial", "Top 10 " * g[2], "Bottom 10 " * g[2]), 
       col = c("blue", "darkgreen", "darkred"), lwd = 1, bty = "n", cex = 0.6,seg.len = 1, x.intersp = 0.6)

## testers for P1 gca1_cor is line 498 report
rwhich(b1 == max(b1)) # 69
which(gca1 == max(gca1)) # Higher g2 (púrpura) 48
rev(order(gca1_cor))[1:10] # Higher gca to sca1j correlation  89 69 31 26 17 82 80 39 14 49
rev(order(b1))[1:10] # High beta1 (azul) 69 80 31 76 39 83 63 45 36 89
rev(order(gca1_var))[1:10] #76 69 45 80 42 83 95 94 63 65
rev(order(gca1))[1:10] # Higher g2 (púrpura) 48 58 97 29 38 22 94 17 92 51
which(gca1_cor > 0.8) # 69 89    69 was as above!! In max b1 but 89 was not in b1 >0.15

par(mar = c(3.3, 3.3, 1.2, 0.8),
    mgp = c(1.7, 0.5, 0),
    tck = -0.015,
    bty = "l",
    cex.axis = 0.8,
    cex.lab  = 0.95)

# 1. Variables
gca1_var <- diag(var(gca2 + sca12))
gca1_cor <- cor(gca2 + sca12, gca2)

plot(gca1_var, gca1_cor, cex = 0.6, pch = 16,
     col = adjustcolor("grey30", 0.5),
     xlab = expression("Testcross variance (" * sigma[tc]^2 * ")"),
     ylab = expression("Ranking accuracy (" * r[italic(g[1]) * "," * italic(g[1~(j)])] * ")"))

# Puntos destacados
points(gca1_var[69], gca1_cor[69], pch = 19, cex = 0.8, col = "#0C66C1")    # 1
points(gca1_var[89], gca1_cor[89], pch = 19, cex = 0.8, col = "#C71D0B")    # 2
points(gca1_var[48], gca1_cor[48], pch = 19, cex = 0.8, col = "purple4")   # 3
points(gca1_var[38], gca1_cor[38], pch = 19, cex = 0.8, col = "darkgreen") # 4

# Etiquetas de texto sobre los puntos (desfase en el eje Y para legibilidad)
text(gca1_var[69], gca1_cor[69], labels = "1", cex = 0.75, font = 2, pos = 3, offset = 0.3, col = "#0C66C1")
text(gca1_var[89], gca1_cor[89], labels = "2", cex = 0.75, font = 2, pos = 3, offset = 0.3, col = "#C71D0B")
text(gca1_var[48], gca1_cor[48], labels = "3", cex = 0.75, font = 2, pos = 3, offset = 0.3, col = "purple4")
text(gca1_var[38], gca1_cor[38], labels = "4", cex = 0.75, font = 2, pos = 3, offset = 0.3, col = "darkgreen")

usr <- par("usr")
x_coor <- usr[2] - (usr[2] - usr[1]) * 0.02
y_coor <- usr[3] + (usr[4] - usr[3]) * 0.02

legend(x = x_coor, y = y_coor,
       legend = expression("1: Higher "*beta[1], 
                           "2: Higher "*r^2,
                           "3: Higher "*g[2], 
                           "4: High "*beta[1]*" + "*g[2]),
       col    = c("#0C66C1", "#C71D0B", "purple4", "darkgreen"),
       pch    = c(19, 19, 19, 19),
       pt.cex = c(0.8, 0.8, 0.8, 0.8),
       bty = "n", cex = 0.8, 
       x.intersp = 0.6, y.intersp = 0.85,
       xjust = 1, yjust = 0)
# good parents vs good testers

# testers with different profiles

which(b1 == max(b1)) # 69
which(gca1 == max(gca1)) # Higher g2 (púrpura) 48
rev(order(gca1_cor))[1:10] # Higher gca to sca1j correlation  89 69 31 26 17 82 80 39 14 49
rev(order(b1))[1:10] # High beta1 (azul) 69 80 31 76 39 83 63 45 36 89
rev(order(gca1_var))[1:10] #76 69 45 80 42 83 95 94 63 65
rev(order(gca1))[1:10] # Higher g2 (púrpura) 48 58 97 29 38 22 94 17 92 51

par(mar = c(3.3, 3.3, 1.2, 0.8),
    mgp = c(1.7, 0.5, 0),
    tck = -0.015,
    bty = "l",
    cex.axis = 0.8,
    cex.lab  = 0.95)

# Gráfico base
plot(b1, gca1, cex = 0.6, pch = 16,
     col = adjustcolor("grey30", 0.5),
     xlab = expression("Regression slope (" * beta[1] * ")"),
     ylab = expression("Tester GCA (" * g[2] * ")"))

# Puntos destacados con sus colores asignados
points(b1[69], gca1[69], pch = 19, cex = 0.8, col = "#0C66C1")    # 1: High beta1 (Azul)
points(b1[89], gca1[89], pch = 19, cex = 0.8, col = "#C71D0B")    # 2: Higher r2 (Rojo)
points(b1[48], gca1[48], pch = 19, cex = 0.8, col = "purple4")   # 3: Higher g2 (Púrpura)
points(b1[38], gca1[38], pch = 19, cex = 0.8, col = "darkgreen") # 4: High beta1 + g2 (Verde)

# Etiquetas numéricas idénticas sobre los puntos (con pos = 3 para que queden arriba)
text(b1[69], gca1[69], labels = "1", cex = 0.75, font = 2, pos = 3, offset = 0.3, col = "#0C66C1")
text(b1[89], gca1[89], labels = "2", cex = 0.75, font = 2, pos = 3, offset = 0.3, col = "#C71D0B")
text(b1[48], gca1[48], labels = "3", cex = 0.75, font = 2, pos = 3, offset = 0.3, col = "purple4")
text(b1[38], gca1[38], labels = "4", cex = 0.75, font = 2, pos = 3, offset = 0.3, col = "darkgreen")

# Leyenda unificada incluyendo la numeración de los perfiles
legend("bottomright",
       legend = expression("1: Higher "*beta[1], 
                           "2: Higher "*r^2,
                           "3: Higher "*g[2], 
                           "4: High "*beta[1]*" + "*g[2]),
       col    = c("#0C66C1", "#C71D0B", "purple4", "darkgreen"),
       pch    = c(19, 19, 19, 19),
       pt.cex = c(0.8, 0.8, 0.8, 0.8),
       bty = "n", cex = 0.8, 
       x.intersp = 0.7, y.intersp = 0.9,
       inset = c(0.01, 0.01))

which(gca1 > 30 & b1 > 0.1) # 38
##############################################################################
# using allele frequencies
hist(mu1_cor <- cor(t(X1), mu1 <- colMeans(X1)))
which(mu1_cor > 0.1)
# 10 20 45
plot(mu1_cor, gca1_cor)

w <- rep(1/m, m) # assuming equal weights across sites
(Ex1 <- sqrt(((X1 - rep(mu1, each = n))^2) %*% w))
plot(Ex1, mu1_cor) # top left
which(Ex1 < 0.89)
# 57 63 97
plot(Ex1, gca1_cor)
plot(Ex1, gca1_cor^2*gca1_var)

# using weighted allele frequencies
w1 <- a + (1 - 2*colMeans(X2+1)/2)*d # using marker GCA as weights
plot(gca1, (X1 %*% w1 - 164.19)/2); abline(a=0, b=1)
(Ex1w <- sqrt((((X1 - rep(mu1, each = n)) %*% w1)^2)))
plot(Ex1w, mu1_cor) # top left
plot(Ex1w, gca1_cor)
plot(Ex1w, gca1_cor^2*gca1_var)

# Its difficult for 1 tester, since the best GCA cor will occur when the SCA var is small...
# for example, take 69 89 from P1
plot(apply(sca12, 2, var), gca1_cor)
apply(sca12, 2, var)[c(69, 89)]
# i do wonder why 69 has such a high correlation - this is interesting and worth finding out...
(gca1_cor^2*gca1_var)[c(69, 89)]
plot(apply(sca12, 2, var), gca1_cor^2*gca1_var)
# Lets leave this for now, but i want to know how to exploit high GCA alignment and discriminating ability,
# so we need to know what drives these two features...
plot(b1, gca1_cor)
plot(b1, gca1_var)
plot(gca1_var, gca1_cor)
plot(b1, gca1_cor^2*gca1_var)
# So if we know b1, then job done,
# but what if we dont??
# what if we only have access to allele frequency data?
# Can we get marker effects for b1?
wb1 <- t(X1) %*% solve(X1 %*% t(X1)) %*% b1
plot(w1, wb1) # good check for outliers... work back to gxe and getting outlier genotypes/environments
# also, GWAS for what makes a good tester...
plot(b1, X1 %*% wb1); abline(a=0, b=1)
Ex1wb <- X1 %*% wb1
plot(Ex1wb, gca1_cor)
plot(Ex1wb, sqrt(gca1_cor^2*gca1_var))
# OK, so assuming we can get b estimated well, then we can identify the 
# best testers. recall, the b will be population dependent though...
# Expectations to update slopes...
# and what about GxE (later).

# so the above captures mean and variance, nice...
# we can also tailor the search for b based on the type of tester.
# e.g., high b1 will also give good hybrids!!
# (although it may not identify the very best...)

# lastly, lets write out our hybrid effects in terms of GCA and SCA
range(hybrid_mat - (mean(hybrid_mat) +
                      matrix(rep(gca1, each = n), ncol = n) + 
                      matrix(rep(gca2, times = n), ncol = n) +
                      sca12))
# and now we partition the "non-crossover" from the SCA...
hybrid_mat[1,1] - mean(hybrid_mat) - (1 + b2[1])*gca1[1] - (1 + b1[1])*gca2[1] - (sca12[1,1] - b2[1]*gca1[1] - b1[1]*gca2[1])
  
sca12_adj <- sca12 - gca2 %*% t(b1)
cor(sca12_adj, gca2)
sca12_adj <- sca12 - t(gca1 %*% t(b2))
cor(t(sca12_adj), gca1)

sca12_adj <- sca12 - gca2 %*% t(b1) - t(gca1 %*% t(b2))
range(cor(t(sca12_adj), gca1))
range(cor(sca12_adj, gca2))
cor(rep(gca1, each = n), rep(gca2, times = n)) # zero by definition

range(hybrid_mat - (mean(hybrid_mat) +
                      # 
                      t(gca1 %*% t(1+ b2)) + # b2 here is the part perfectly correlated with gca1 when crossed to P2. 
                      gca2 %*% t(1+b1) +     # b1 here is the part perfectly correlated with gca2 when crossed to P1.
                      # 
                      sca12_adj))
hist(1 + b1)
which(b1 == max(b1)) # 69 in P1 is best at discriminating other
hist(1 + b2)
which(b2 == max(b2)) # 45 in P2 is best at discriminating other
# with a pretty picture for Sven

# First, we show three inbreds from P1, and highlight that tester 45 from P2 is the best at discriminating these
p1_ids <- c(which(gca1 == min(gca1)), which(gca1 == max(gca1)), which(round(gca1,0) == 0)[1])
plot_reaction(cov.mat = 1+b2, slopes = gca1, ids = p1_ids, envs = "45") + labs(y = "GCA for pool 1", x = "Tester from Pool 2")

# Next, we show three inbreds from P2, and highlight that tester 69 from P1 is the best at discriminating these
p2_ids <- c(which(gca2 == min(gca2)), which(gca2 == max(gca2)), which(round(gca2,0) == 0)[1])
plot_reaction(cov.mat = 1+b1, slopes = gca2, ids = p2_ids, envs = "69") + labs(y = "GCA for pool 1", x = "Tester from Pool 2")

# We should also include the deviations due to SCA into the above...
# pool 1
plot_reaction(cov.mat = 1+b2, slopes = gca1, deviations = sca12_adj,
              ids = p1_ids, envs = "45") + labs(y = "GCA for pool 1", x = "Tester from Pool 2")
# cool
# GCA vs SCA plot
plot(colMeans(sca12_adj^2),gca1) # would hybrids coming from top left be more stable across envs?
# i.e, they originate from a half-sib family which performs similarly...
# Ultimately, the best performing single hybrid will come from top right
# application to swards in forages.

# pool 2
plot_reaction(cov.mat = 1+b1, slopes = gca2, deviations = t(sca12_adj),
              ids = p2_ids, envs = "69") + labs(y = "GCA for pool 2", x = "Tester from Pool 1")
# plot(diag(var(t(sca12_adj)))[order(b1)])
# op vs rmsd type plot
plot(rowMeans((sca12_adj)^2),gca2) 



###################################
# take some tester sets...
t <- 3    # testers per set
s <- 5000 # sets

# Pool 1
set.seed(123)
X1sample <- t(matrix(unlist(lapply(1:s, function(x) sample(1:n, t))), ncol = s))
X1sample[1,] <- c(69, 80, 31)
hist(table(X1sample))

# Pool 2
X2sample <- t(matrix(unlist(lapply(1:s, function(x) sample(1:n, t))), ncol = s))
hist(table(X2sample))

# so the question is, which set/s give the best GCA alignment and discriminating ability.
# is it the testers with the highest b's?, plus sca parts sum to zero... (dcross)
# so individually they can give high hibrid values and collectively they give good gca alignment.

# start with which sets give what gca / sca
hybrid_mat_k <- apply(X1sample[,, drop = FALSE], 1, function(i) hybrid_mat[, i, drop = FALSE], simplify = F)
gca2_k <- lapply(hybrid_mat_k, function(x) rowMeans(x) - mean(x))
hist(unlist(lapply(gca2_k, function(x) cor(x, gca2))))
hist(unlist(lapply(gca2_k, function(x) cor(x, gca2)^2*var(x))))

b1_k <- apply(X1sample[,, drop = FALSE], 1, function(i) b1[i], simplify = F)
plot(b1_mean_k <- unlist(lapply(b1_k, function(x) mean(x))), unlist(lapply(gca2_k, function(x) cor(x, gca2)*sd(x))))
plot(b1_sd_k <- unlist(lapply(b1_k, function(x) sd(x))), unlist(lapply(gca2_k, function(x) cor(x, gca2)*sd(x))))
# ok, so if we get the mean across testers high, then we also get the correct GCA, but what about SCA?

plot(b1_sd_k, b1_mean_k)
# whats the benefit of getting variance in the bk?
# more diversity captured?

sca2_k <- lapply(hybrid_mat_k, function(x) x - rep(colMeans(x), each = n) - rep(rowMeans(x), times = t) + mean(x))
which(b1_mean_k == max(b1_mean_k))
sca2_k[[1]]
# 69 45 80
rev(order(b1))
# 69  80  31

plot_reaction(cov.mat = 1+b1, slopes = gca2, deviations = t(sca12_adj),
              ids = 1:10, envs = c("69","80","31")) + labs(y = "GCA for pool 1", x = "Tester from Pool 2")
ss1 <- svd(scale(X1, scale = F))
plot(cumsum(ss1$d^2)/sum(ss1$d^2))
plot(ss1$u[,2], ss1$u[,1])
points(ss1$u[69,2], ss1$u[69,1], col = "purple", cex = 2)
points(ss1$u[80,2], ss1$u[80,1], col = "darkorange", cex = 2)
points(ss1$u[31,2], ss1$u[31,1], col = "navy", cex = 2)

ss1 <- svd(scale(rbind(X1, X2), scale = F))
plot(cumsum(ss1$d^2)/sum(ss1$d^2))
plot(ss1$u[,2], ss1$u[,1])
points(ss1$u[69,2], ss1$u[69,1], col = "purple", cex = 2)
points(ss1$u[80,2], ss1$u[80,1], col = "darkorange", cex = 2)
points(ss1$u[31,2], ss1$u[31,1], col = "navy", cex = 2)

# what is the trade-off, i.e., optimisation
plot(gca1, b1) # top right...
plot(gca2, b2) # top right...

# OP vs rmsd 
# op vs rmsd type plot
plot(rowMeans((sca12_adj)^2),gca2) 

# Now lets create the b's for each set. 

plot(apply(hybrid_mat - mean(hybrid_mat), 2, max), gca1)
plot(apply(hybrid_mat - mean(hybrid_mat), 2, max), b1)

