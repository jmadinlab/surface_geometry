# Code from paper

### Scope and grain parameters ----
L <- 2     # Extent, 2 by 2 m reef patches
L0 <- 2/32 # Grain, reoslution of processing ~ 6 cm

### Data preparation
source("R/functions.R")
dat <- read.csv("output/master_20200709.csv", as.is=TRUE)

### Surface descriptor assocations ----

pdf("figs/fig2.pdf", width = 8, height = 8)

layout(matrix(c(1, 4, 4, 2, 4, 4, 3, 4, 4), 3))
par(mar=c(4, 6, 2, 1))

plot(R2_log10 ~ D_theory, dat, axes=FALSE, ylab=expression(paste("Rugosity (", italic(R)^2 - 1, ")")), xlab=expression(paste("Fractal dimension (", italic(D), ")")), ylim=c(-2, 0.7), xlim=c(2, 2.6), col=rgb(0,0,0,0.3), cex=1, pch=1)
# points(R2_log10 ~ D, dat[dat$site=="Megaplot",], pch=20, col=rgb(0,0,1,0.5))
axis(2, at=log10(c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5)), labels=c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5), las=2)
axis(1, at=c(2, 2.2, 2.4, 2.6))
pts <- chull(dat$R2_log10, dat$D_theory)
pts <- c(pts, pts[1])
# lines(dat$D_theory[pts], dat$R2_log10[pts], lty=2)
mtext(expression(bold(a)), 3, 0, adj=0, cex=1.2)
text(2.1, 0.6, substitute(expression(italic(r)^2 == MYVALUE), list(MYVALUE=format(cor(dat$D_theory, dat$R2_log10)^2, digits = 3)))[[2]])

plot(R2_log10 ~ HL0_log10, dat, axes=FALSE, xlab=expression(paste("Height range (", frac(Delta * italic(H), sqrt(2) * italic(L)[0]), ")")), ylab=expression(paste("Rugosity (", italic(R)^2 - 1, ")")), xlim=log10(c(2, 50)), ylim=c(-2, 0.7), col=rgb(0,0,0,0.3), cex=1, pch=1)
# points(R2_log10 ~ HL0_log10, dat[dat$site=="Megaplot",], pch=20, col=rgb(0,0,1,0.5))
axis(1, at=log10(c(2, 5, 10, 20, 50)), labels=c(2, 5, 10, 20, 50))
axis(2, at=log10(c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5)), labels=c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5), las=2)
pts <- chull(dat$HL0_log10, dat$R2_log10)
pts <- c(pts, pts[1])
# lines(dat$HL0_log10[pts], dat$R2_log10[pts], lty=2)
mtext(expression(bold(b)), 3, 0, adj=0, cex=1.2)
text(0.5, 0.6, substitute(expression(italic(r)^2 == MYVALUE), list(MYVALUE=format(cor(dat$HL0_log10, dat$R2_log10)^2, digits = 3)))[[2]])

plot(HL0_log10 ~ D_theory, dat, axes=FALSE, ylab=expression(paste("Height range (", frac(Delta * italic(H), sqrt(2) * italic(L)[0]), ")")), xlab=expression(paste("Fractal dimension (", italic(D), ")")), ylim=log10(c(2, 50)), xlim=c(2, 2.6), col=rgb(0,0,0,0.3), cex=1, pch=1)
# points(HL0_log10 ~ D, dat[dat$site=="Megaplot",], pch=20, col=rgb(0,0,1,0.5))
axis(2, at=log10(c(2, 5, 10, 20, 50)), labels=c(2, 5, 10, 20, 50))
axis(1, at=c(2, 2.2, 2.4, 2.6))
pts <- chull(dat$D_theory, dat$HL0_log10)
pts <- c(pts, pts[1])
# lines(dat$D_theory[pts], dat$HL0_log10[pts], lty=2)
mtext(expression(bold(c)), 3, 0, adj=0, cex=1.2)
text(2.4, 1.6, substitute(expression(italic(r)^2 == MYVALUE), list(MYVALUE=format(cor(dat$HL0_log10, dat$D_theory)^2, digits = 3)))[[2]])

par(mar=c(1, 1, 1, 1))

grid.lines <- 350 #350
R <- dat$R2_log10 # Rugosity
H <- dat$HL0_log10  # Height range
D <- dat$D  # Height range
M <- dat$mid  # Height range

R.pred <- seq(min(R)-0.25, max(R)+0.25, length=grid.lines)
H.pred <- seq(min(H)-0.25, max(H)+0.25, length=grid.lines)
mat <- expand.grid(R=R.pred, H=H.pred)
D.mat <- matrix(D_func((sqrt(2) * L0)*(10^mat$H), sqrt(10^mat$R+1), L, L0), nrow = grid.lines, ncol = grid.lines)
dat$R2_log10 <- log10(dat$R_theory^2 - 1)
dat$HL0_log10 <- log10(dat$H / (L0 * sqrt(2)))

scatter3D(R, H, D, pch = 20, col=rev(hcl.colors(100, "reds")), cex = 1.2, xlim=c(min(R)-0.25, max(R)+0.25), ylim=c(min(H)-0.25, max(H)+0.25), zlim=c(2, max(D)+0.05), ylab="Height range", xlab="Rugosity", zlab="Fractal dimension", surf=list(x=R.pred, y=H.pred, z=D.mat, facets=NA, col=rgb(0,0,0,0.02), fitpoints=dat$D_theory), theta=215, phi=0, colkey=list(side = 4, length = 0.5, width = 1, line.clab = 1, dist=-0.1), clab=expression(italic(D)))

ss_tot <- sum((dat$D - mean(dat$D_theory))^2)
ss_f <- sum((dat$D_theory - mean(dat$D_theory))^2)
ss_res <- sum((dat$D - dat$D_theory)^2)
r_2 <- 1 - (ss_res / ss_tot)
text3D(1.8, 1.2, 2.6, labels=expression(italic(r)^2 == 0.983), add = TRUE, cex=1)

mtext(expression(bold(d)), 3, -2, adj=0.2, cex=1.2)

dev.off()

### Geometric-biodiversity coupling ----

megaplot <- dat[dat$site=="megaplot",]
megaplot$abd_sqrt <- sqrt(megaplot$abd)
megaplot$spp_sqrt <- sqrt(megaplot$spp)
megaplot$pie_asin <- asin(megaplot$pie)

megaplotA <- megaplot[megaplot$abd > 0,] # Test results removing zero abundance patches
megaplotR <- megaplot[megaplot$R2_log10 < 0.2,] # Test results removing high rugosity patches

### PIE Statistics ----

pieR_gam <- gam(pie_asin ~ s(R2_log10), data=megaplot, method = "REML")
summary(pieR_gam)
# par(mfrow = c(2,2))
# gam.check(pieR_gam)
pieR_gam2 <- gam(pie_asin ~ s(R2_log10), data=megaplotR, method = "REML")
summary(pieR_gam2)
pieR_lm <- lm(pie_asin ~ R2_log10 + R2_log10_sq, megaplot)
drop1(pieR_lm, test="Chisq")
# summary(pieR_lm)$r.squared

pieD_gam <- gam(pie_asin ~ s(D_theory), data=megaplot, method = "REML")
summary(pieD_gam)
# par(mfrow = c(2,2))
# gam.check(pieD_gam)
pieD_lm <- lm(pie_asin ~ D_theory + D_theory_sq, megaplot)
drop1(pieD_lm, test="Chisq")
# summary(pieD_lm)$r.squared

pieH_gam <- gam(pie_asin ~ s(HL0_log10), data=megaplot, method = "REML")
summary(pieH_gam)
# par(mfrow = c(2,2))
# gam.check(pieH_gam)
pieH_lm <- lm(pie_asin ~ HL0_log10 + HL0_log10_sq, megaplot)
drop1(pieH_lm, test="Chisq")
# summary(pieH_lm)$r.squared

pie_gam <- gam(pie_asin ~ s(R2_log10) + s(HL0_log10) + s(D_theory), data=megaplot, method = "REML")
summary(pie_gam)
# par(mfrow = c(2,2))
# gam.check(pie_gam)

pie_gam2 <- gam(pie_asin ~ s(R2_log10) + s(HL0_log10) + D_theory, data=megaplot, method = "REML")
summary(pie_gam2)
# par(mfrow = c(2,2))
# gam.check(pie_gam2)

### Richness statistics ----

sppR_gam <- gam(spp_sqrt ~ s(R2_log10), data=megaplot, method = "REML")
summary(sppR_gam)
# par(mfrow = c(2,2))
# gam.check(sppR_gam)
sppR_gam2 <- gam(spp_sqrt ~ s(R2_log10), data=megaplotR, method = "REML")
summary(sppR_gam2)
sppR_lm <- lm(spp_sqrt ~ R2_log10 + R2_log10_sq, megaplot)
drop1(sppR_lm, test="Chisq")
# summary(sppR_lm)$r.squared

sppD_gam <- gam(spp_sqrt ~ s(D_theory), data=megaplot, method = "REML")
summary(sppD_gam)
# par(mfrow = c(2,2))
# gam.check(sppD_gam)
sppD_lm <- lm(spp_sqrt ~ D_theory + D_theory_sq, megaplot)
drop1(sppD_lm, test="Chisq")
# summary(sppD_lm)$r.squared

sppH_gam <- gam(spp_sqrt ~ s(HL0_log10), data=megaplot, method = "REML")
summary(sppH_gam)
# par(mfrow = c(2,2))
# gam.check(sppH_gam)
sppH_lm <- lm(spp_sqrt ~ HL0_log10 + HL0_log10_sq, megaplot)
drop1(sppH_lm, test="Chisq")
# summary(sppH_lm)$r.squared

spp_gam <- gam(spp_sqrt ~ s(R2_log10) + s(HL0_log10) + s(D_theory), data=megaplot, method = "REML")
summary(spp_gam)
# par(mfrow = c(2,2))
# gam.check(spp_gam)

spp_gam2 <- gam(spp_sqrt ~ s(R2_log10) + s(HL0_log10) + D_theory, data=megaplot, method = "REML")
summary(spp_gam2)
# par(mfrow = c(2,2))
# gam.check(spp_gam2)

### Abundance statistics ----

abdR_gam <- gam(abd_sqrt ~ s(R2_log10), data=megaplot, method = "REML")
summary(abdR_gam)
# par(mfrow = c(2,2))
# gam.check(abdR_gam)
abdR_gam2 <- gam(abd_sqrt ~ s(R2_log10), data=megaplotR, method = "REML")
summary(abdR_gam2)
abdR_lm <- lm(abd_sqrt ~ R2_log10 + R2_log10_sq, megaplot)
drop1(abdR_lm, test="Chisq")
# summary(abdR_lm)$r.squared

abdD_gam <- gam(abd_sqrt ~ s(D_theory), data=megaplot, method = "REML")
summary(abdD_gam)
# par(mfrow = c(2,2))
# gam.check(abdD_gam)
abdD_lm <- lm(abd_sqrt ~ D_theory + D_theory_sq, megaplot)
drop1(abdD_lm, test="Chisq")
# summary(abdD_lm)$r.squared

abdH_gam <- gam(abd_sqrt ~ s(HL0_log10), data=megaplot, method = "REML")
summary(abdH_gam)
# par(mfrow = c(2,2))
# gam.check(abdH_gam)
abdH_lm <- lm(abd_sqrt ~ HL0_log10 + HL0_log10_sq, megaplot)
drop1(abdH_lm, test="Chisq")
# summary(abdH_lm)$r.squared

abd_gam <- gam(abd_sqrt ~ s(R2_log10) + s(HL0_log10) + s(D_theory), data=megaplot, method = "REML")
summary(abd_gam)
# par(mfrow = c(2,2))
# gam.check(abd_gam)

abd_gam2 <- gam(abd_sqrt ~ s(R2_log10) + s(HL0_log10) + D_theory, data=megaplot, method = "REML")
summary(abd_gam2)
# par(mfrow = c(2,2))
# gam.check(abd_gam2)

### Plots geometry-biodiversity ---- 

x_R <- seq(min(megaplot$R2_log10), max(megaplot$R2_log10), length.out = 100)
x_R2 <- seq(min(megaplotR$R2_log10), max(megaplotR$R2_log10), length.out = 100)
x_D <- seq(min(megaplot$D_theory), max(megaplot$D_theory), length.out = 100)
x_H <- seq(min(megaplot$HL0_log10), max(megaplot$HL0_log10), length.out = 100)

png("figs/figS7.png", width = 8, height = 8, units = 'in', res = 300)
par(mfrow=c(3, 3))

# PIE

plot(pie_asin ~ R2_log10, megaplot, col="grey", axes=FALSE, xlab=expression(paste(italic(R)^2 - 1, " (log scale)")), ylab="Probability of interspecific encounter")
axis(1); axis(2, las=2); mtext("a", 3, 0, adj=0)
txt <- substitute(expression(italic(r)^2 == MYVALUE), list(MYVALUE=format(summary(pieR_gam)$r.sq, digits = 3)))[[2]]
mtext(txt, 3, 0, adj=1, cex=0.8)

y_pred <- predict(pieR_gam, data.frame(R2_log10 = x_R), se.fit=TRUE, type="response")
polygon(c(x_R, rev(x_R)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,0,0.2))
lines(x_R, y_pred$fit, lty=2, col="black")

y_pred <- predict(pieR_lm, data.frame(R2_log10 = x_R, R2_log10_sq = x_R^2), se.fit=TRUE, type="response")
polygon(c(x_R, rev(x_R)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,1,0.2))
lines(x_R, y_pred$fit, lty=1, col="blue")

y_pred <- predict(pieR_gam2, data.frame(R2_log10 = x_R2), se.fit=TRUE, type="response")
polygon(c(x_R2, rev(x_R2)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(1,0,0,0.2))
lines(x_R2, y_pred$fit, lty=1, col="red")

legend("bottomright", legend=c("GAM","LM (poly)", "GAM (subset R)"), lty=c(2, 1, 1), col=c("black", "blue", "red"), bty="n")

plot(pie_asin ~ D_theory, megaplot, col="grey", axes=FALSE, xlab="D", ylab="Probability of interspecific encounter")
axis(1); axis(2, las=2); mtext("b", 3, 0, adj=0)

y_pred <- predict(pieD_gam, data.frame(D_theory = x_D), se.fit=TRUE, type="response")
polygon(c(x_D, rev(x_D)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,0,0.2))
lines(x_D, y_pred$fit, lty=2)

y_pred <- predict(pieD_lm, data.frame(D_theory = x_D, D_theory_sq = x_D^2), se.fit=TRUE, type="response")
polygon(c(x_D, rev(x_D)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,1,0.2))
lines(x_D, y_pred$fit, lty=1, col="blue")

txt <- substitute(expression(italic(r)^2 == MYVALUE), list(MYVALUE=format(summary(pieD_gam)$r.sq, digits = 3)))[[2]]
mtext(txt, 3, 0, adj=1, cex=0.8)

plot(pie_asin ~ HL0_log10, megaplot, col="grey", axes=FALSE, xlab=expression(paste(frac(Delta * italic(H), sqrt(2) * italic(L)[0]), " (log scale)")), ylab="Probability of interspecific encounter")
axis(1); axis(2, las=2); mtext("c", 3, 0, adj=0)

y_pred <- predict(pieH_gam, data.frame(HL0_log10 = x_H), se.fit=TRUE, type="response")
polygon(c(x_H, rev(x_H)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,0,0.2))
lines(x_H, y_pred$fit, lty=2)

y_pred <- predict(pieH, data.frame(HL0_log10 = x_H, HL0_log10_sq = x_H^2), se.fit=TRUE, type="response")
polygon(c(x_H, rev(x_H)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,1,0.2))
lines(x_H, y_pred$fit, lty=1, col="blue")

txt <- substitute(expression(italic(r)^2 == MYVALUE), list(MYVALUE=format(summary(pieH_gam)$r.sq, digits = 3)))[[2]]
mtext(txt, 3, 0, adj=1, cex=0.8)

# RICHNESS

plot(spp_sqrt ~ R2_log10, megaplot, col="grey", axes=FALSE, xlab=expression(paste(italic(R)^2 - 1, " (log scale)")), ylab="Richness")
axis(1); axis(2, las=2); mtext("d", 3, 0, adj=0)

y_pred <- predict(sppR_gam, data.frame(R2_log10 = x_R), se.fit=TRUE)
polygon(c(x_R, rev(x_R)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,0,0.2))
lines(x_R, y_pred$fit, lty=2, col="black")

y_pred <- predict(sppR_lm, data.frame(R2_log10 = x_R, R2_log10_sq = x_R^2), se.fit=TRUE)
polygon(c(x_R, rev(x_R)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,1,0.2))
lines(x_R, y_pred$fit, lty=1, col="blue")

y_pred <- predict(sppR_gam2, data.frame(R2_log10 = x_R2), se.fit=TRUE, type="response")
polygon(c(x_R2, rev(x_R2)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(1,0,0,0.2))
lines(x_R2, y_pred$fit, lty=1, col="red")

txt <- substitute(expression(italic(r)^2 == MYVALUE), list(MYVALUE=format(summary(sppR_gam)$r.sq, digits = 3)))[[2]]
mtext(txt, 3, 0, adj=1, cex=0.8)


plot(spp_sqrt ~ D_theory, megaplot, col="grey", axes=FALSE, xlab="D", ylab="Richness")
axis(1); axis(2, las=2); mtext("e", 3, 0, adj=0)

y_pred <- predict(sppD_gam, data.frame(D_theory = x_D), se.fit=TRUE)
polygon(c(x_D, rev(x_D)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,0,0.2))
lines(x_D, y_pred$fit, lty=2)

y_pred <- predict(sppD_lm, data.frame(D_theory = x_D, D_theory_sq = x_D^2), se.fit=TRUE)
polygon(c(x_D, rev(x_D)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,1,0.2))
lines(x_D, y_pred$fit, lty=1, col="blue")

txt <- substitute(expression(italic(r)^2 == MYVALUE), list(MYVALUE=format(summary(sppD_gam)$r.sq, digits = 3)))[[2]]
mtext(txt, 3, 0, adj=1, cex=0.8)


plot(spp_sqrt ~ HL0_log10, megaplot, col="grey", axes=FALSE, xlab=expression(paste(frac(Delta * italic(H), sqrt(2) * italic(L)[0]), " (log scale)")), ylab="Richness")
axis(1); axis(2, las=2); mtext("f", 3, 0, adj=0)

y_pred <- predict(sppH_gam, data.frame(HL0_log10 = x_H), se.fit=TRUE)
polygon(c(x_H, rev(x_H)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,0,0.2))
lines(x_H, y_pred$fit, lty=2)

y_pred <- predict(sppH_lm, data.frame(HL0_log10 = x_H, HL0_log10_sq = x_H^2), se.fit=TRUE)
polygon(c(x_H, rev(x_H)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,1,0.2))
lines(x_H, y_pred$fit, lty=1, col="blue")

txt <- substitute(expression(italic(r)^2 == MYVALUE), list(MYVALUE=format(summary(sppH_gam)$r.sq, digits = 3)))[[2]]
mtext(txt, 3, 0, adj=1, cex=0.8)

# ABUNDANCE

plot(abd_sqrt ~ R2_log10, megaplot, col="grey", axes=FALSE, xlab=expression(paste(italic(R)^2 - 1, " (log scale)")), ylab="Abundance")
axis(1); axis(2, las=2); mtext("g", 3, 0, adj=0)

y_pred <- predict(abdR_gam, data.frame(R2_log10 = x_R), se.fit=TRUE)
polygon(c(x_R, rev(x_R)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,0,0.2))
lines(x_R, y_pred$fit, lty=2)

y_pred <- predict(abdR_lm, data.frame(R2_log10 = x_R, R2_log10_sq = x_R^2), se.fit=TRUE)
polygon(c(x_R, rev(x_R)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,1,0.2))
lines(x_R, y_pred$fit, lty=1, col="blue")

y_pred <- predict(abdR_gam2, data.frame(R2_log10 = x_R2), se.fit=TRUE, type="response")
polygon(c(x_R2, rev(x_R2)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(1,0,0,0.2))
lines(x_R2, y_pred$fit, lty=1, col="red")

txt <- substitute(expression(italic(r)^2 == MYVALUE), list(MYVALUE=format(summary(abdR_gam)$r.sq, digits = 3)))[[2]]
mtext(txt, 3, 0, adj=1, cex=0.8)


plot(abd_sqrt ~ D_theory, megaplot, col="grey", axes=FALSE, xlab="D", ylab="Abundance")
axis(1); axis(2, las=2); mtext("h", 3, 0, adj=0)

y_pred <- predict(abdD_gam, data.frame(D_theory = x_D), se.fit=TRUE)
polygon(c(x_D, rev(x_D)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,0,0.2))
lines(x_D, y_pred$fit, lty=2)

y_pred <- predict(abdD_lm, data.frame(D_theory = x_D, D_theory_sq = x_D^2), se.fit=TRUE)
polygon(c(x_D, rev(x_D)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,1,0.2))
lines(x_D, y_pred$fit, lty=1, col="blue")

txt <- substitute(expression(italic(r)^2 == MYVALUE), list(MYVALUE=format(summary(abdD_gam)$r.sq, digits = 3)))[[2]]
mtext(txt, 3, 0, adj=1, cex=0.8)


plot(abd_sqrt ~ HL0_log10, megaplot, col="grey", axes=FALSE, xlab=expression(paste(frac(Delta * italic(H), sqrt(2) * italic(L)[0]), " (log scale)")), ylab="Abundance")
axis(1); axis(2, las=2); mtext("i", 3, 0, adj=0)

y_pred <- predict(abdH_gam, data.frame(HL0_log10 = x_H), se.fit=TRUE)
polygon(c(x_H, rev(x_H)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,0,0.2))
lines(x_H, y_pred$fit, lty=2)

y_pred <- predict(abdH_lm, data.frame(HL0_log10 = x_H, HL0_log10_sq = x_H^2), se.fit=TRUE)
polygon(c(x_H, rev(x_H)), c(y_pred$fit + 2*y_pred$se.fit, rev(y_pred$fit - 2*y_pred$se.fit)), border=NA, col=rgb(0,0,1,0.2))
lines(x_H, y_pred$fit, lty=1, col="blue")

txt <- substitute(expression(italic(r)^2 == MYVALUE), list(MYVALUE=format(summary(abdH_gam)$r.sq, digits = 3)))[[2]]
mtext(txt, 3, 0, adj=1, cex=0.8)

dev.off()

# Diagnostics

png("figs/figS8_diag.png", width = 8, height = 12, units = 'in', res = 300)
par(mfcol = c(4, 3))
gam.check(sppR_gam)
gam.check(sppD_gam)
gam.check(sppH_gam)
dev.off()

### Site examples plot ----


D <- seq(2, 2.6, length=50)
R2_log10 <- seq(min(dat$R2_log10)-0.1, max(dat$R2_log10)+0.1, length=50)

mat_HL0 <- matrix(NA, length(D), length(R2_log10))
mat_pie <- matrix(NA, length(D), length(R2_log10))

for (i in 1:length(D)) {
  for (j in 1:length(R2_log10)) {
    HL0_log10 <- HL0_func(D[i], sqrt((10^R2_log10[j]) + 1), L, L0)
    mat_HL0[i, j] <- HL0_log10
    mat_pie[i, j] <- predict(pie_gam2, data.frame(R2_log10=R2_log10[j], D_theory=D[i], HL0_log10=HL0_log10))
  }
}

mat <- matrix(1, length(D), length(R2_log10))
matx <- matrix(D, length(D), length(R2_log10))
maty <- matrix(R2_log10, length(D), length(R2_log10), byrow=TRUE)

mat_pie[!in.chull(matx, maty, dat$D_theory[chull(dat$R2_log10, dat$D_theory)], dat$R2_log10[chull(dat$R2_log10, dat$D_theory)])] <- NA
mat_HL0[!in.chull(matx, maty, dat$D_theory[chull(dat$R2_log10, dat$D_theory)], dat$R2_log10[chull(dat$R2_log10, dat$D_theory)])] <- NA


# png("figs/fig3.png", width = 8.8, height = 4, units = "in", res = 300)
pdf("figs/fig3t.pdf", width = 8.8, height = 4)

layout(matrix(c(c(1, 2, rep(3, 3), 8, 9), c(4, 5, rep(3, 3), 10, 11), c(6, 7, rep(3, 3), 12, 13)), nrow=3, byrow=TRUE))

pp <- 0.5

par(mar=c(pp, pp, pp, 0))
plot(1,1, axes=FALSE, xlab="", ylab="", type="n", xlim=c(0, 1), ylim=c(0, 1), asp=1, yaxs="i", xaxs="i")
rasterImage(readJPEG("data/images_for_figures/01_northreef_dem.jpg"), 0, 0, 1, 1)
mtext("a", 3, -1.5, adj=0.1, cex=1.2, col="black")
par(mar=c(pp, 0, pp, pp))
plot(1,1, axes=FALSE, xlab="", ylab="", type="n", xlim=c(0, 1), ylim=c(0, 1), asp=1, yaxs="i", xaxs="i")
rasterImage(readJPEG("data/images_for_figures/01_northreef_mos.jpg"), 0, 0, 1, 1)

par(mar=c(5, 4.5, 1, 1))

image(D, R2_log10, mat_HL0, col=NA, axes=FALSE, ylab=expression(paste("Rugosity (", italic(R)^2 - 1, ")")), xlab=expression(paste("Fractal dimension (", italic(D), ")")), xlim=c(2, 2.6), ylim=c(-1.5, 0.7))
axis(2, at=log10(c(0.05, 0.1, 0.2, 0.5, 1, 2, 5)), labels=c(0.05, 0.1, 0.2, 0.5, 1, 2, 5))
axis(1, at=seq(2, 2.6, 0.2))
# mtext(expression(bold(a)), 3, 0, adj=0, cex=1.2)

temp <- dat[dat$site=="Mermaid Cove",]
pts <- chull(temp$D_theory, temp$R2_log10)
pts <- c(pts, pts[1])
polygon(temp$D_theory[pts], temp$R2_log10[pts], col=hcl.colors(6, alpha=0.3)[2], border=NA)
points(R2_log10 ~ D_theory, temp, col=hcl.colors(6, alpha=1)[2], pch=20)
text(2.07, log10(0.065), "a", cex=1.2)

text(2.09, 0.2, "Greater height\nrange", cex=0.7)
text(2.53, -1.3, "Smaller height\nrange", cex=0.7)

temp <- dat[dat$site=="Osprey",]
pts <- chull(temp$D_theory, temp$R2_log10)
pts <- c(pts, pts[1])
polygon(temp$D_theory[pts], temp$R2_log10[pts], col=hcl.colors(6, alpha=0.3)[1], border=NA)
points(R2_log10 ~ D_theory, temp, col=hcl.colors(6, alpha=1)[1], pch=20)
text(2.44, log10(0.22), "b", cex=1.2)

temp <- dat[dat$site=="Lagoon 2",]
pts <- chull(temp$D_theory, temp$R2_log10)
pts <- c(pts, pts[1])
polygon(temp$D_theory[pts], temp$R2_log10[pts], col=hcl.colors(6, alpha=0.3)[3], border=NA)
points(R2_log10 ~ D_theory, temp, col=hcl.colors(6, alpha=1)[3], pch=20)
text(2.5, log10(2.1), "c", cex=1.2)

temp <- dat[dat$site=="Resort",]
pts <- chull(temp$D_theory, temp$R2_log10)
pts <- c(pts, pts[1])
polygon(temp$D_theory[pts], temp$R2_log10[pts], col=hcl.colors(6, alpha=0.3)[4], border=NA)
points(R2_log10 ~ D_theory, temp, col=hcl.colors(6, alpha=1)[4], pch=20)
text(2.27, log10(3.5), "d", cex=1.2)

temp <- dat[dat$site=="South Island",]
pts <- chull(temp$D_theory, temp$R2_log10)
pts <- c(pts, pts[1])
polygon(temp$D_theory[pts], temp$R2_log10[pts], col=hcl.colors(6, alpha=0.3)[5], border=NA)
points(R2_log10 ~ D_theory, temp, col=hcl.colors(6, alpha=1)[5], pch=20)
text(2.1, log10(0.8), "e", cex=1.2)

temp <- dat[dat$site=="Horseshoe",]
pts <- chull(temp$D_theory, temp$R2_log10)
pts <- c(pts, pts[1])
polygon(temp$D_theory[pts], temp$R2_log10[pts], col=hcl.colors(6, alpha=0.3)[6], border=NA)
points(R2_log10 ~ D_theory, temp, col=hcl.colors(6, alpha=1)[6], pch=20)
text(2.325, log10(1), "f", cex=1.2)


par(mar=c(pp, pp, pp, 0))
plot(1,1, axes=FALSE, xlab="", ylab="", type="n", xlim=c(0, 1), ylim=c(0, 1), asp=1, yaxs="i", xaxs="i")
rasterImage(readJPEG("data/images_for_figures/02_osprey_dem.jpg"), 0, 0, 1, 1)
mtext("b", 3, -1.5, adj=0.1, cex=1.2, col="black")
par(mar=c(pp, 0, pp, pp))
plot(1,1, axes=FALSE, xlab="", ylab="", type="n", xlim=c(0, 1), ylim=c(0, 1), asp=1, yaxs="i", xaxs="i")
rasterImage(readJPEG("data/images_for_figures/02_osprey_mos.jpg"), 0, 0, 1, 1)

par(mar=c(pp, pp, pp, 0))
plot(1,1, axes=FALSE, xlab="", ylab="", type="n", xlim=c(0, 1), ylim=c(0, 1), asp=1, yaxs="i", xaxs="i")
rasterImage(readJPEG("data/images_for_figures/03_lagoon02_dem.jpg"), 0, 0, 1, 1)
mtext("c", 3, -1.5, adj=0.1, cex=1.2, col="black")
par(mar=c(pp, 0, pp, pp))
plot(1,1, axes=FALSE, xlab="", ylab="", type="n", xlim=c(0, 1), ylim=c(0, 1), asp=1, yaxs="i", xaxs="i")
rasterImage(readJPEG("data/images_for_figures/03_lagoon02_mos.jpg"), 0, 0, 1, 1)

par(mar=c(pp, pp, pp, 0))
plot(1,1, axes=FALSE, xlab="", ylab="", type="n", xlim=c(0, 1), ylim=c(0, 1), asp=1, yaxs="i", xaxs="i")
rasterImage(readJPEG("data/images_for_figures/04_resort_dem.jpg"), 0, 0, 1, 1)
mtext("d", 3, -1.5, adj=0.1, cex=1.2, col="black")
par(mar=c(pp, 0, pp, pp))
plot(1,1, axes=FALSE, xlab="", ylab="", type="n", xlim=c(0, 1), ylim=c(0, 1), asp=1, yaxs="i", xaxs="i")
rasterImage(readJPEG("data/images_for_figures/04_resort_mos04.jpg"), 0, 0, 1, 1)

par(mar=c(pp, pp, pp, 0))
plot(1,1, axes=FALSE, xlab="", ylab="", type="n", xlim=c(0, 1), ylim=c(0, 1), asp=1, yaxs="i", xaxs="i")
rasterImage(readJPEG("data/images_for_figures/05_southis_dem.jpg"), 0, 0, 1, 1)
mtext("e", 3, -1.5, adj=0.1, cex=1.2, col="black")
par(mar=c(pp, 0, pp, pp))
plot(1,1, axes=FALSE, xlab="", ylab="", type="n", xlim=c(0, 1), ylim=c(0, 1), asp=1, yaxs="i", xaxs="i")
rasterImage(readJPEG("data/images_for_figures/05_southis_mos.jpg"), 0, 0, 1, 1)

par(mar=c(pp, pp, pp, 0))
plot(1,1, axes=FALSE, xlab="", ylab="", type="n", xlim=c(0, 1), ylim=c(0, 1), asp=1, yaxs="i", xaxs="i")
rasterImage(readJPEG("data/images_for_figures/06_horseshoe_dem.jpg"), 0, 0, 1, 1)
mtext("f", 3, -1.5, adj=0.1, cex=1.2, col="black")
par(mar=c(pp, 0, pp, pp))
plot(1,1, axes=FALSE, xlab="", ylab="", type="n", xlim=c(0, 1), ylim=c(0, 1), asp=1, yaxs="i", xaxs="i")
rasterImage(readJPEG("data/images_for_figures/06_horseshoe_mos.jpg"), 0, 0, 1, 1)

dev.off()

### Geometric/biodiversity coupling ----

pdf("figs/fig4.pdf", width = 6, height = 8)

# layout(matrix(c(rep(1, 12), rep(2, 12)), nrow=5, byrow=TRUE))

par(mfrow=c(2, 1), mar=c(2, 4, 2, 2))

image(D, R2_log10, sin(mat_pie), col=rev(hcl.colors(100, "reds",  alpha=0.7)), axes=FALSE, ylab="", xlab="", xlim=c(2, 2.6), ylim=c(-1.7, 0.7))

axis(2, at=log10(c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5)), labels=c(0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5))
axis(1, at=seq(2, 2.6, 0.2))

# points(dat$D[dat$site=="Megaplot"], dat$R2_log10[dat$site=="Megaplot"], col=rgb(0,0,0,0.3), pch=20)
contour(D, R2_log10, sin(mat_pie), levels=seq(0, 1, 0.1), add=TRUE)
mtext(expression(bold(a)), 3, 0, adj=0, cex=1.2)
mtext("Diversity (probability of interspecific encounter)", 3, 0, adj=0.1)
txt <- substitute(expression(italic(r)^2 == MYVALUE), list(MYVALUE=format(summary(pie_gam)$r.sq, digits = 3)))[[2]]
mtext(txt, 3, 0, adj=1, cex=0.8)

text(2.09, 0.2, "Greater height\nrange", cex=0.7)
text(2.53, -1.3, "Smaller height\nrange", cex=0.7)
mtext(expression(paste("Rugosity (", italic(R)^2 - 1, ")")), 2, 2)
mtext(expression(italic(D)), 1, 2, cex=0.9)

par(mar=c(2, 2, 2, 2))
plot(1,1, axes=FALSE, xlab="", ylab="", type="n", xlim=c(0, 1.67802), ylim=c(0, 1), asp=1, yaxs="i", xaxs="i")
rasterImage(readJPEG("data/images_for_figures/Trimodal_dem4a.jpg"), 0, 0, 1.67802, 1)
mtext(expression(bold(b)), 3, 0, adj=0.1, cex=1.2, col="black")

dev.off()


### Reef record maps ----

png("figs/figS2.png", width = 6.5, height = 14, units = 'in', res = 300)

par(mfrow=c(7, 3), mar=c(2, 2, 1, 0), oma=c(1.2, 1.2, 0, 0))
# keep <- data.frame()
breakpoints <- seq(-5, -1.3, length.out=50)

records <- dat[dat$site != "megaplot",]

for (i in 1:21) {
  map <- raster(paste0("data/records/", unique(records$rec)[i], ".tif"))
  plot(map, main="", breaks=breakpoints, col=hcl.colors(50), asp=1, axes=FALSE, box=FALSE, legend=FALSE)
  x0 <- mean(map@extent[1:2]) - 4
  y0 <- mean(map@extent[3:4]) - 4
  rect(rep(seq(x0, x0+6, 2), 4), rep(seq(y0, y0+6, 2), each=4), rep(seq(x0, x0+6, 2), 4)+2, rep(seq(y0, y0+6, 2), each=4)+2, border="white")
  axis(1, at=seq(x0, x0+8, 2), labels=seq(0, 8, 2))
  axis(2, at=seq(y0, y0+8, 2), labels=seq(0, 8, 2))
  mtext(unique(temp$site)[i], 3, -1, cex=0.8, adj=0)
}
mtext("Easting (m)", 1, 0, outer=TRUE)
mtext("Northing (m)", 2, 0, outer=TRUE)

dev.off()

### Abundance versus richness ----

png("figs/figS6.png", width = 5, height = 5, units = 'in', res = 300)

plot(spp ~ abd, megaplot, col="grey", axes=FALSE, xlab="Abundance", ylab="Richness")
axis(1)
axis(2, las=2)
taxabd <- lm(spp_sqrt ~ abd_sqrt, data=megaplot)
summary(taxabd)
x_new <- seq(min(megaplot$abd_sqrt), max(megaplot$abd_sqrt), length.out = 100)
y_pred <- predict(taxabd, data.frame(abd_sqrt = x_new), interval="confidence")
polygon(c(x_new, rev(x_new))^2, c(y_pred[,2], rev(y_pred[,3]))^2, border=NA, col=rgb(0,0,0,0.2))
lines((x_new)^2, y_pred[,1]^2, lty=2, col="black")

dev.off()

### Fractal dimension residuals ----

png("figs/figS5.png", width = 4, height = 4, units = 'in', res = 300)

hist(dat$D_theory - dat$D, breaks=50, xlab="D theory - D empirical", main=NA, prob=TRUE, axes=FALSE)
axis(1, las=2)
axis(2, las=2)
lines(density(dat$D_theory - dat$D))
abline(v=0, lwd=2, col="red")
# mtext(expression(bold(b)), 3, 0, adj=0, cex=1.2)

dev.off()
